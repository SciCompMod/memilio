/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Elisabeth Kluth, Khoa Nguyen, David Kerkmann, Rene Schmieding
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#include "abm/model_functions.h"
#include "abm/location.h"
#include "abm/person.h"
#include "abm/random_events.h"
#include "abm/infection.h"
#include "abm/virus_variant.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/logging.h"

namespace mio
{
namespace abm
{

ScalarType total_exposure_by_contacts(const ContactExposureRates& rates, const CellIndex cell_index,
                                      const VirusVariant virus, const AgeGroup age_receiver,
                                      const LocalInfectionParameters& params)
{
    assert(age_receiver < rates.size<AgeGroup>());
    ScalarType total_exposure = 0;
    for (AgeGroup age_transmitter(0); age_transmitter < rates.size<AgeGroup>(); ++age_transmitter) {
        total_exposure +=
            rates[{cell_index, virus, age_transmitter}] * params.get<ContactRates>()[{age_receiver, age_transmitter}];
    }
    return total_exposure;
}

ScalarType total_exposure_by_air(const AirExposureRates& rates, const CellIndex cell_index, const VirusVariant virus,
                                 const Parameters& global_params)
{
    return rates[{cell_index, virus}] * global_params.get<AerosolTransmissionRates>()[{virus}];
}

void interact(PersonalRandomNumberGenerator& personal_rng, Person& person, const Location& location,
              const AirExposureRates& local_air_exposure, const ContactExposureRates& local_contact_exposure,
              const TimePoint t, const TimeSpan dt, const Parameters& global_parameters)
{
    // make sure all dimensions are set correctly and all indices are valid
    assert(location.get_cells().size() == local_air_exposure.size<CellIndex>().get());
    assert(location.get_cells().size() == local_contact_exposure.size<CellIndex>().get());
    assert(local_contact_exposure.size<VirusVariant>() == local_air_exposure.size<VirusVariant>());
    assert(local_contact_exposure.size<VirusVariant>() == VirusVariant::Count);
    assert(local_contact_exposure.size<AgeGroup>().get() == global_parameters.get_num_groups());
    assert(person.get_age() < local_contact_exposure.size<AgeGroup>());
    assert(std::all_of(person.get_cells().begin(), person.get_cells().end(), [&](const auto& cell) {
        return cell < location.get_cells().size();
    }));

    if (person.get_infection_state(t) == InfectionState::Susceptible) {
        auto& local_parameters = location.get_infection_parameters();
        // TODO: we need to define what a cell is used for, as the loop may lead to incorrect results for multiple cells
        auto age_receiver          = person.get_age();
        ScalarType mask_protection = person.get_mask_protective_factor(global_parameters);
        assert(person.get_cells().size() && "Person is in multiple cells. Interact logic is incorrect at the moment.");
        for (auto cell_index : person.get_cells()) {
            std::pair<VirusVariant, ScalarType> local_indiv_expected_trans[static_cast<uint32_t>(VirusVariant::Count)];
            for (uint32_t v = 0; v != static_cast<uint32_t>(VirusVariant::Count); ++v) {
                VirusVariant virus = static_cast<VirusVariant>(v);
                ScalarType exposed_viral_shed =
                    (total_exposure_by_contacts(local_contact_exposure, cell_index, virus, age_receiver,
                                                local_parameters) +
                     total_exposure_by_air(local_air_exposure, cell_index, virus, global_parameters)) *
                    (1 - mask_protection) * (1 - person.get_protection_factor(t, virus, global_parameters));
                ScalarType infection_rate =
                    global_parameters.get<InfectionRateFromViralShed>()[{virus}] * exposed_viral_shed;
                local_indiv_expected_trans[v] = std::make_pair(virus, infection_rate);
            }
            VirusVariant virus =
                random_transition(personal_rng, VirusVariant::Count, dt,
                                  local_indiv_expected_trans); // use VirusVariant::Count for no virus submission
            if (virus != VirusVariant::Count) {
                person.add_new_infection(Infection(personal_rng, virus, age_receiver, global_parameters, t + dt / 2,
                                                   mio::abm::InfectionState::Susceptible,
                                                   person.get_latest_protection(t + dt / 2),
                                                   false)); // Starting time in second order approximation
            }
        }
    }
    person.add_time_at_location(dt);
}

void add_exposure_contribution(AirExposureRates& local_air_exposure, ContactExposureRates& local_contact_exposure,
                               const Person& person, const Location& location, const Parameters& params,
                               const TimePoint t, const TimeSpan dt)
{
    if (person.get_location() != location.get_id()) {
        mio::log_debug("In add_exposure_contribution: Person {} is not at Location {}", person.get_id().get(),
                       location.get_id().get());
    }

    if (person.is_infected(t)) {
        auto& infection = person.get_infection();
        auto virus      = infection.get_virus_variant();
        auto age        = person.get_age();
        // average viral shed over the time step to second order accuracy using midpoint rule
        const auto viral_shed = infection.get_viral_shed(t + dt / 2);
        const auto quarantine_factor =
            person.is_in_quarantine(t, params) ? (1.0 - params.get<QuarantineEffectiveness>()) : 1.0;

        for (CellIndex cell : person.get_cells()) {
            auto air_contribution     = viral_shed * quarantine_factor;
            auto contact_contribution = viral_shed * quarantine_factor;

            if (location.get_infection_parameters().get<UseLocationCapacityForTransmissions>()) {
                air_contribution *= location.get_cells()[cell.get()].compute_space_per_person_relative();
            }

            local_air_exposure[{cell, virus}] += air_contribution;
            local_contact_exposure[{cell, virus, age}] += contact_contribution;
        }
    }
}

bool change_location(Person& person, const Location& destination, const TransportMode mode,
                     const std::vector<uint32_t>& cells)
{
    assert(std::all_of(cells.begin(), cells.end(), [&](const auto& cell) {
        return cell < destination.get_cells().size();
    })); // make sure cell indices are valid

    if (person.get_location() != destination.get_id()) {
        person.set_location(destination.get_type(), destination.get_id(), destination.get_model_id());
        person.get_cells() = cells;
        person.set_last_transport_mode(mode);

        return true;
    }
    else {
        mio::log_debug("In change_location: Person {} already is at Location {}", person.get_id().get(),
                       destination.get_id().get());
        return false;
    }
}

void adjust_contact_rates(Location& location, size_t num_agegroups)
{
    if (location.get_infection_parameters().get<MaximumContacts>() == std::numeric_limits<ScalarType>::max()) {
        return;
    }
    for (auto contact_from = AgeGroup(0); contact_from < AgeGroup(num_agegroups); contact_from++) {
        ScalarType total_contacts = 0.;
        // slizing would be preferred but is problematic since both Tags of ContactRates are AgeGroup
        for (auto contact_to = AgeGroup(0); contact_to < AgeGroup(num_agegroups); contact_to++) {
            total_contacts += location.get_infection_parameters().get<ContactRates>()[{contact_from, contact_to}];
        }
        if (total_contacts > location.get_infection_parameters().get<MaximumContacts>()) {
            for (auto contact_to = AgeGroup(0); contact_to < AgeGroup(num_agegroups); contact_to++) {
                location.get_infection_parameters().get<ContactRates>()[{contact_from, contact_to}] =
                    location.get_infection_parameters().get<ContactRates>()[{contact_from, contact_to}] *
                    location.get_infection_parameters().get<MaximumContacts>() / total_contacts;
            }
        }
    }
}

} // namespace abm
} // namespace mio
