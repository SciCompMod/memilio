/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Elisabeth Kluth, Carlotta Gerstein, Martin J. Kuehn, Khoa Nguyen, David Kerkmann
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
#include "abm/mask_type.h"
#include "abm/mask.h"
#include "abm/location.h"
#include "abm/random_events.h"
#include "abm/infection.h"
#include "memilio/utils/random_number_generator.h"
#include <mutex>
#include <numeric>

namespace mio
{
namespace abm
{

Location::Location(LocationId loc_id, size_t num_agegroups, uint32_t num_cells)
    : m_id(loc_id)
    , m_capacity_adapted_transmission_risk(false)
    , m_parameters(num_agegroups)
    , m_cells(num_cells, num_agegroups)
    , m_required_mask(MaskType::Community)
    , m_npi_active(false)
{
    assert(num_cells > 0 && "Number of cells has to be larger than 0.");
}

Location Location::copy_location_without_persons(size_t num_agegroups)
{
    Location copy_loc  = Location(*this);
    copy_loc.m_persons = std::vector<observer_ptr<Person>>();
    copy_loc.m_cells   = std::vector<Cell>{num_agegroups};
    for (uint32_t idx = 0; idx < m_cells.size(); idx++) {
        copy_loc.set_capacity(get_capacity(idx).persons, get_capacity(idx).volume, idx);
        copy_loc.get_cached_exposure_rate_contacts(idx) = get_cached_exposure_rate_contacts(idx);
        copy_loc.get_cached_exposure_rate_air(idx)      = get_cached_exposure_rate_air(idx);
    }
    return copy_loc;
}

ScalarType Location::transmission_contacts_per_day(uint32_t cell_index, VirusVariant virus, AgeGroup age_receiver,
                                                   size_t num_agegroups) const
{
    assert(age_receiver.get() < num_agegroups);
    ScalarType prob = 0;
    for (uint32_t age_transmitter = 0; age_transmitter != num_agegroups; ++age_transmitter) {
        prob += m_cells[cell_index].m_cached_exposure_rate_contacts[{virus, static_cast<AgeGroup>(age_transmitter)}] *
                m_parameters.get<ContactRates>()[{age_receiver, static_cast<AgeGroup>(age_transmitter)}];
    }
    return prob;
}

ScalarType Location::transmission_air_per_day(uint32_t cell_index, VirusVariant virus,
                                              const Parameters& global_params) const
{
    return m_cells[cell_index].m_cached_exposure_rate_air[{virus}] *
           global_params.get<AerosolTransmissionRates>()[{virus}];
}

void Location::interact(Person::RandomNumberGenerator& rng, Person& person, TimePoint t, TimeSpan dt,
                        const Parameters& global_params) const
{
    // TODO: we need to define what a cell is used for, as the loop may lead to incorrect results for multiple cells
    auto age_receiver          = person.get_age();
    ScalarType mask_protection = person.get_mask_protective_factor(global_params);
    assert(person.get_cells().size() && "Person is in multiple cells. Interact logic is incorrect at the moment.");
    for (auto cell_index :
         person.get_cells()) { // TODO: the logic here is incorrect in case a person is in multiple cells
        std::pair<VirusVariant, ScalarType> local_indiv_trans_prob[static_cast<uint32_t>(VirusVariant::Count)];
        for (uint32_t v = 0; v != static_cast<uint32_t>(VirusVariant::Count); ++v) {
            VirusVariant virus = static_cast<VirusVariant>(v);
            ScalarType local_indiv_trans_prob_v =
                (std::min(
                     m_parameters.get<MaximumContacts>(),
                     transmission_contacts_per_day(cell_index, virus, age_receiver, global_params.get_num_groups())) +
                 transmission_air_per_day(cell_index, virus, global_params)) *
                (1 - mask_protection) * dt.days() * (1 - person.get_protection_factor(t, virus, global_params));

            local_indiv_trans_prob[v] = std::make_pair(virus, local_indiv_trans_prob_v);
        }
        VirusVariant virus =
            random_transition(rng, VirusVariant::Count, dt,
                              local_indiv_trans_prob); // use VirusVariant::Count for no virus submission
        if (virus != VirusVariant::Count) {
            person.add_new_infection(Infection(rng, virus, age_receiver, global_params, t + dt / 2,
                                               mio::abm::InfectionState::Exposed, person.get_latest_protection(),
                                               false)); // Starting time in first approximation
        }
    }
}

void Location::cache_exposure_rates(TimePoint t, TimeSpan dt, size_t num_agegroups)
{
    //cache for next step so it stays constant during the step while subpopulations change
    //otherwise we would have to cache all state changes during a step which uses more memory
    for (auto& cell : m_cells) {
        cell.m_cached_exposure_rate_contacts = {{VirusVariant::Count, AgeGroup(num_agegroups)}, 0.};
        cell.m_cached_exposure_rate_air      = {{VirusVariant::Count}, 0.};
        for (auto&& p : cell.m_persons) {
            if (p->is_infected(t)) {
                auto& inf  = p->get_infection();
                auto virus = inf.get_virus_variant();
                auto age   = p->get_age();
                /* average infectivity over the time step 
                 *  to second order accuracy using midpoint rule
                */
                cell.m_cached_exposure_rate_contacts[{virus, age}] += inf.get_infectivity(t + dt / 2);
                cell.m_cached_exposure_rate_air[{virus}] += inf.get_infectivity(t + dt / 2);
            }
        }
        if (m_capacity_adapted_transmission_risk) {
            cell.m_cached_exposure_rate_air.array() *= cell.compute_space_per_person_relative();
        }
    }
}

void Location::add_person(Person& p, std::vector<uint32_t> cells)
{
    std::lock_guard<std::mutex> lk(m_mut);
    m_persons.push_back(&p);
    for (uint32_t cell_idx : cells)
        m_cells[cell_idx].m_persons.push_back(&p);
}

void Location::remove_person(Person& p)
{
    std::lock_guard<std::mutex> lk(m_mut);
    m_persons.erase(std::remove(m_persons.begin(), m_persons.end(), &p), m_persons.end());
    for (auto&& cell : m_cells) {
        cell.m_persons.erase(std::remove(cell.m_persons.begin(), cell.m_persons.end(), &p), cell.m_persons.end());
    }
}

size_t Location::get_number_persons() const
{
    return m_persons.size();
}

/*
For every cell in a location we have a transmission factor that is nomalized to m_capacity.volume / m_capacity.persons of 
the location "Home", which is 66. We multiply this rate with the individual size of each cell to obtain a "space per person" factor.
*/
ScalarType Cell::compute_space_per_person_relative()
{
    if (m_capacity.volume != 0) {
        return 66.0 / m_capacity.volume;
    }
    else {
        return 1.0;
    }
}

size_t Cell::get_subpopulation(TimePoint t, InfectionState state) const
{
    return count_if(m_persons.begin(), m_persons.end(), [&](observer_ptr<Person> p) {
        return p->get_infection_state(t) == state;
    });
}

size_t Location::get_subpopulation(TimePoint t, InfectionState state) const
{
    return count_if(m_persons.begin(), m_persons.end(), [&](observer_ptr<Person> p) {
        return p->get_infection_state(t) == state;
    });
}

} // namespace abm
} // namespace mio
