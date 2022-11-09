/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Elisabeth Kluth
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
#include "memilio/utils/random_number_generator.h"
#include "abm/location.h"
#include "abm/person.h"
#include "abm/random_events.h"

#include <numeric>

namespace mio
{
namespace abm
{

Location::Location(LocationType type, uint32_t index, uint32_t num_cells)
    : m_type(type)
    , m_index(index)
    , m_subpopulations{}
    , m_cached_exposure_rate({AgeGroup::Count, VaccinationState::Count})
    , m_testing_scheme()
    , m_cells(std::vector<Cell>(num_cells))
{
}

boost::optional<Infection> Location::interact(const Person& person, const TimePoint& t, const TimeSpan dt,
                              const GlobalInfectionParameters& /*global_params*/) const
{
    auto infection         = person.get_infection(); // always susceptible
    auto vaccination_state = person.get_vaccination_state();
    auto age               = person.get_age();

    std::shared_ptr<Virus> virus; // subject to change
    *virus = Virus(1,1,1);   // need to determine virus from other infected
                             // persons at location/cell
    
    Infection* temp;
    if (!person.get_cells().empty()) {
        for (auto cell_index : person.get_cells()) {
            return random_transition({}, dt,
            {*new Infection(virus, t), m_cells[cell_index].cached_exposure_rate[{age, vaccination_state}]});
        }
        
    }
    else {
        return random_transition(infection, dt,
            {{*new Infection(virus, t), m_cached_exposure_rate[{age, vaccination_state}]}});
    }
}

void Location::begin_step(TimeSpan /*dt*/, const GlobalInfectionParameters& global_params)
{
    //cache for next step so it stays constant during the step while subpopulations change
    //otherwise we would have to cache all state changes during a step which uses more memory
    if (m_cells.empty() && m_num_persons == 0) {
        m_cached_exposure_rate = {{AgeGroup::Count, VaccinationState::Count}, 0.};
    }
    else if (m_cells.empty()) {
        auto num_carriers              = get_subpopulation(InfectionState::Carrier);
        auto num_infected              = get_subpopulation(InfectionState::Infected);
        m_cached_exposure_rate.array() = std::min(m_parameters.get<MaximumContacts>(), double(m_num_persons)) /
                                         m_num_persons *
                                         (global_params.get<SusceptibleToExposedByCarrier>().array() * num_carriers +
                                          global_params.get<SusceptibleToExposedByInfected>().array() * num_infected);
    }
    else {
        for (auto& cell : m_cells) {
            if (cell.num_people == 0) {
                cell.cached_exposure_rate = {{AgeGroup::Count, VaccinationState::Count}, 0.};
            }
            else {
                cell.cached_exposure_rate.array() =
                    std::min(m_parameters.get<MaximumContacts>(), double(cell.num_people)) / cell.num_people *
                    (global_params.get<SusceptibleToExposedByCarrier>().array() * cell.num_carriers +
                     global_params.get<SusceptibleToExposedByInfected>().array() * cell.num_infected);
            }
        }
    }
}

void Location::add_person(const Person& p, const TimePoint& t)
{
    ++m_num_persons;
    InfectionState s = p.get_infection_state(t);
    change_subpopulation(s, +1);
    if (!m_cells.empty()) {
        for (auto i : p.get_cells()) {
            ++m_cells[i].num_people;
            if (s == InfectionState::Carrier) {
                ++m_cells[i].num_carriers;
            }
            else if (s == InfectionState::Infected || s == InfectionState::Infected_Severe ||
                     s == InfectionState::Infected_Critical) {
                ++m_cells[i].num_infected;
            }
        }
    }
}

void Location::remove_person(const Person& p, const TimePoint& t)
{
    --m_num_persons;
    InfectionState s = p.get_infection_state(t);
    change_subpopulation(s, -1);
    for (auto i : p.get_cells()) {
        --m_cells[i].num_people;
        if (s == InfectionState::Carrier) {
            --m_cells[i].num_carriers;
        }
        else if (s == InfectionState::Infected || s == InfectionState::Infected_Severe ||
                 s == InfectionState::Infected_Critical) {
            --m_cells[i].num_infected;
        }
    }
}

void Location::changed_state(const Person& p, InfectionState old_infection_state, const TimePoint& t)
{
    change_subpopulation(old_infection_state, -1);
    change_subpopulation(p.get_infection_state(t), +1);
    auto current_infection_state = p.get_infection_state(t);
    for (auto i : p.get_cells()) {
        if (old_infection_state == InfectionState::Carrier) {
            --m_cells[i].num_carriers;
        }
        else if (old_infection_state == InfectionState::Infected ||
                 old_infection_state == InfectionState::Infected_Severe ||
                 old_infection_state == InfectionState::Infected_Critical) {
            --m_cells[i].num_infected;
        }
        if (current_infection_state == InfectionState::Carrier) {
            ++m_cells[i].num_carriers;
        }
        else if (current_infection_state == InfectionState::Infected ||
                 current_infection_state == InfectionState::Infected_Severe ||
                 current_infection_state == InfectionState::Infected_Critical) {
            ++m_cells[i].num_infected;
        }
    }
}

void Location::change_subpopulation(InfectionState s, int delta)
{
    m_subpopulations[size_t(s)] += delta;
    assert(m_subpopulations[size_t(s)] >= 0 && "subpopulations must be non-negative");
}

int Location::get_subpopulation(InfectionState s) const
{
    return m_subpopulations[size_t(s)];
}

Eigen::Ref<const Eigen::VectorXi> Location::get_subpopulations() const
{
    return Eigen::Map<const Eigen::VectorXi>(m_subpopulations.data(), m_subpopulations.size());
}

} // namespace abm
} // namespace mio
