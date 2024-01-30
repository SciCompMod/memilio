#include "functions.h"

#include "abm/location.h"
#include "abm/person.h"
#include "abm/random_events.h"
#include "abm/infection.h"
#include "abm/trip_list.h"
#include "abm/virus_variant.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include <algorithm>
#include <cassert>
#include <cstddef>

namespace mio
{

namespace abm
{

// TODO: on argument order: maybe personal_rng first, as it always(?) is a non-const reference

ScalarType daily_transmissions_by_contacts(const Location::ContactExposureRates& rates, CellIndex cell_index,
                                           VirusVariant virus, AgeGroup age_receiver,
                                           const LocalInfectionParameters& params)
{
    assert(age_receiver < rates.size<AgeGroup>());
    ScalarType prob = 0;
    for (AgeGroup age_transmitter(0); age_transmitter < rates.size<AgeGroup>(); ++age_transmitter) {
        prob +=
            rates[{cell_index, virus, age_transmitter}] * params.get<ContactRates>()[{age_receiver, age_transmitter}];
    }
    return prob;
}

ScalarType daily_transmissions_by_air(const Location::AirExposureRates& rates, CellIndex cell_index, VirusVariant virus,
                                      const Parameters& global_params)
{
    return rates[{cell_index, virus}] * global_params.get<AerosolTransmissionRates>()[{virus}];
}

void interact(Person& person, const Location& location, const Location::AirExposureRates& local_air_exposure,
              const Location::ContactExposureRates& local_contact_exposure, const TimePoint t, const TimeSpan dt,
              const Parameters& global_parameters, Person::RandomNumberGenerator& personal_rng)
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
        for (auto cell_index :
             person.get_cells()) { // TODO: the logic here is incorrect in case a person is in multiple cells
            std::pair<VirusVariant, ScalarType> local_indiv_trans_prob[static_cast<uint32_t>(VirusVariant::Count)];
            for (uint32_t v = 0; v != static_cast<uint32_t>(VirusVariant::Count); ++v) {
                VirusVariant virus = static_cast<VirusVariant>(v);
                ScalarType local_indiv_trans_prob_v =
                    (std::min(local_parameters.get<MaximumContacts>(),
                              daily_transmissions_by_contacts(local_contact_exposure, cell_index, virus, age_receiver,
                                                              local_parameters)) +
                     daily_transmissions_by_air(local_air_exposure, cell_index, virus, global_parameters)) *
                    dt.days() * (1 - mask_protection) * (1 - person.get_protection_factor(t, virus, global_parameters));

                local_indiv_trans_prob[v] = std::make_pair(virus, local_indiv_trans_prob_v);
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

void add_exposure_contribution(Location::AirExposureRates& local_air_exposure,
                               Location::ContactExposureRates& local_contact_exposure, const Person& person,
                               Location& location, TimePoint t, TimeSpan dt)
{ // TODO: location should be const, but there is no const cell accessor
    assert([&]() {
        if (person.get_location() != location.get_id()) {
            mio::log_warning("Person with id {} is not at Location with id {}", person.get_person_id(),
                             location.get_index());
        }
        return true;
    }());
    if (person.is_infected(t)) {
        auto& infection = person.get_infection();
        auto virus      = infection.get_virus_variant();
        auto age        = person.get_age();
        // average infectivity over the time step to second order accuracy using midpoint rule
        for (CellIndex cell : person.get_cells()) {
            if (location.get_infection_parameters().get<UseLocationCapacityForTransmissions>()) {
                local_air_exposure[{cell, virus}] +=
                    infection.get_infectivity(t + dt / 2) *
                    location.get_cells()[cell.get()].compute_space_per_person_relative();
            }
            else {
                local_air_exposure[{cell, virus}] += infection.get_infectivity(t + dt / 2);
            }
            local_contact_exposure[{cell, virus, age}] += infection.get_infectivity(t + dt / 2);
        }
    }
}

void interact(Person& person, Location& location, const std::vector<Person>& local_population, const TimePoint t,
              const TimeSpan dt, const Parameters& global_parameters, Person::RandomNumberGenerator& personal_rng)
{ // TODO: make location const&, need to fix add_exposure_contribution
    Location::AirExposureRates local_air_exposure{{CellIndex(location.get_cells().size()), VirusVariant::Count}, 0.};
    Location::ContactExposureRates local_contact_exposure{
        {CellIndex(location.get_cells().size()), VirusVariant::Count, AgeGroup(global_parameters.get_num_groups())},
        0.};
    for (const Person& p : local_population) {
        add_exposure_contribution(local_air_exposure, local_contact_exposure, p, location, t, dt);
    }
    interact(person, location, local_air_exposure, local_contact_exposure, t, dt, global_parameters, personal_rng);
}

bool migrate(Person& person, const Location& destination, const std::vector<uint32_t>& cells, const TransportMode mode)
{
    assert(std::all_of(cells.begin(), cells.end(), [&](const auto& cell) {
        return cell < destination.get_cells().size();
    })); // make sure cell indices are valid

    if (person.get_location() != destination.get_id()) {
        person.set_location(destination);
        person.get_cells() = cells;
        person.set_last_transport_mode(mode);

        return true;
    }
    else {
        return false;
    }
}

} // namespace abm

} // namespace mio