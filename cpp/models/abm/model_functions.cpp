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

ScalarType daily_transmissions_by_contacts(const ContactExposureRates& rates, const CellIndex cell_index,
                                           const VirusVariant virus, const AgeGroup age_receiver,
                                           size_t age_receiver_group_size, const LocalInfectionParameters& params)
{
    assert(age_receiver < rates.size<AgeGroup>());
    ScalarType prob = 0;
    for (AgeGroup age_transmitter(0); age_transmitter < rates.size<AgeGroup>(); ++age_transmitter) {
        if (age_receiver == age_transmitter &&
            age_receiver_group_size > 1) // adjust for the person not meeting themself
        {
            std::cout << "Adjusting for self-contact: age_receiver " << age_receiver.get()
                      << ", age_receiver_group_size " << age_receiver_group_size << " with rate "
                        << rates[{cell_index, virus, age_transmitter}] << "and contact rate " << params.get<ContactRates>()[{age_receiver, age_transmitter}] <<" and further factors " <<
                         (double)(age_receiver_group_size) << "\n";

            prob += rates[{cell_index, virus, age_transmitter}] *params.get<ContactRates>()[{age_receiver, age_transmitter}] * age_receiver_group_size /
                    (age_receiver_group_size - 1);
               
        }
        else {
            prob += rates[{cell_index, virus, age_transmitter}] *
                    params.get<ContactRates>()[{age_receiver, age_transmitter}];
        }
    }
    return prob;
}

ScalarType daily_transmissions_by_air(const AirExposureRates& rates, const CellIndex cell_index,
                                      const VirusVariant virus, const Parameters& global_params)
{
    return rates[{cell_index, virus}] * global_params.get<AerosolTransmissionRates>()[{virus}];
}

void interact(PersonalRandomNumberGenerator& personal_rng, Person& person, const Location& location,
              const PopulationByAge& local_population_by_age, const AirExposureRates& local_air_exposure,
              const ContactExposureRates& local_contact_exposure, const TimePoint t, const TimeSpan dt,
              const Parameters& global_parameters)
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
        for (CellIndex cell_index : person.get_cells()) {
            std::pair<VirusVariant, ScalarType> local_indiv_trans_prob[static_cast<uint32_t>(VirusVariant::Count)];
            for (uint32_t v = 0; v != static_cast<uint32_t>(VirusVariant::Count); ++v) {
                VirusVariant virus                      = static_cast<VirusVariant>(v);
                size_t local_population_by_age_receiver = local_population_by_age[{cell_index, age_receiver}];
                auto daily_contacts_exposure = daily_transmissions_by_contacts(
                    local_contact_exposure, cell_index, virus, age_receiver, local_population_by_age_receiver,
                    local_parameters);
                auto daily_air_exposure = daily_transmissions_by_air(local_air_exposure, cell_index, virus, global_parameters);
                auto mask=  (1 - mask_protection);
                auto personal_protection =
                    (1 - person.get_protection_factor(t, virus, global_parameters));
                ScalarType local_indiv_trans_prob_v =
                    (daily_contacts_exposure + daily_air_exposure) *
                    mask * personal_protection
                    *global_parameters.get<InfectionRateFromViralShed>()[{virus}];

                local_indiv_trans_prob[v] = std::make_pair(virus, local_indiv_trans_prob_v);
                std::cout << "Person " << person.get_id().get() <<"has local_indiv_trans_prob_v: "
                          << local_indiv_trans_prob_v << " daily_contacts_exposure = exposed_viral_shed/n_s = viral_shed_normiert: " << daily_contacts_exposure
                          << " daily_air_exposure: " << daily_air_exposure << " mask factor: " << mask
                          << " personal protection: " << personal_protection << "\n";
            }
            VirusVariant virus =
                random_transition(personal_rng, VirusVariant::Count, dt,
                                  local_indiv_trans_prob); // use VirusVariant::Count for no virus submission
            if (virus != VirusVariant::Count) {
                person.add_new_infection(Infection(personal_rng, virus, age_receiver, global_parameters, t + dt / 2,
                                                   mio::abm::InfectionState::Exposed, person.get_latest_protection(),
                                                   false)); // Starting time in first approximation
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
        // average infectivity over the time step to second order accuracy using midpoint rule
        const auto infectivity = infection.get_infectivity(t + dt / 2);
        const auto quarantine_factor =
            person.is_in_quarantine(t, params) ? (1.0 - params.get<QuarantineEffectiveness>()) : 1.0;
        std::cout << "Person " << person.get_id().get() << " adds exposure contribution for virus " << static_cast<uint32_t>(virus)
                  << ": infectivity " << infectivity << ", quarantine factor " << quarantine_factor << "\n";

        for (CellIndex cell : person.get_cells()) {
            auto air_contribution     = infectivity * quarantine_factor;
            auto contact_contribution = infectivity * quarantine_factor;

            if (location.get_infection_parameters().get<UseLocationCapacityForTransmissions>()) {
                air_contribution *= location.get_cells()[cell.get()].compute_space_per_person_relative();
            }

            local_air_exposure[{cell, virus}] += air_contribution;
            local_contact_exposure[{cell, virus, age}] += contact_contribution;
        }
    }
}

void normalize_exposure_contribution(ContactExposureRates& local_contact_exposure,
                                     const PopulationByAge& local_population_by_age)
{
    // make sure all dimensions are set correctly and all indices are valid
    assert(local_population_by_age.size<AgeGroup>() == local_contact_exposure.size<AgeGroup>());
    assert(local_population_by_age.size<CellIndex>() == local_contact_exposure.size<CellIndex>());
    assert(local_contact_exposure.size<VirusVariant>() == VirusVariant::Count);

    for (auto index : make_index_range(local_contact_exposure.size())) {
        auto age_index = reduce_index<Index<CellIndex, AgeGroup>>(index);
        if (local_population_by_age[age_index] > 0) {
            // this instruction is not and does not need to be atomic
            local_contact_exposure[index] = local_contact_exposure[index] / local_population_by_age[age_index];
            std::cout << "Normalized contact exposure " << local_contact_exposure[index] << " with population "
                      << local_population_by_age[age_index] << " to "
                      << local_contact_exposure[index] << "\n";
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
