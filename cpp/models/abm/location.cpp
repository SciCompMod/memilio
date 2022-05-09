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

Location::Location(LocationType type, uint32_t index, uint32_t num_cells)
    : m_type(type)
    , m_index(index)
    , m_subpopulations{}
    , m_cached_exposure_rate({AbmAgeGroup::Count, mio::VaccinationState::Count})
    , m_testing_scheme()
    , m_cells(std::vector<Cell>(num_cells))
{
    set_default_capacity(type);
}

InfectionState Location::interact(const Person& person, TimeSpan dt,
                                  const GlobalInfectionParameters& global_params) const
{
    auto infection_state   = person.get_infection_state();
    auto vaccination_state = person.get_vaccination_state();
    auto age               = person.get_age();
    switch (infection_state) {
    case InfectionState::Susceptible:
        if (!person.get_cells().empty()) {
            for (auto cell_index : person.get_cells()) {
                InfectionState new_state = random_transition(
                    infection_state, dt,
                    {{InfectionState::Exposed, m_cells[cell_index].cached_exposure_rate[{age, vaccination_state}]}});
                if (new_state != infection_state) {
                    return new_state;
                }
            }
            return infection_state;
        }
        else {
            return random_transition(infection_state, dt,
                                     {{InfectionState::Exposed, m_cached_exposure_rate[{age, vaccination_state}]}});
        }
    case InfectionState::Carrier:
        return random_transition(
            infection_state, dt,
            {{InfectionState::Infected, global_params.get<CarrierToInfected>()[{age, vaccination_state}]},
             {InfectionState::Recovered_Carrier, global_params.get<CarrierToRecovered>()[{age, vaccination_state}]}});
    case InfectionState::Infected:
        return random_transition(
            infection_state, dt,
            {{InfectionState::Recovered_Infected, global_params.get<InfectedToRecovered>()[{age, vaccination_state}]},
             {InfectionState::Infected_Severe, global_params.get<InfectedToSevere>()[{age, vaccination_state}]}});
    case InfectionState::Infected_Severe:
        return random_transition(
            infection_state, dt,
            {{InfectionState::Recovered_Infected, global_params.get<SevereToRecovered>()[{age, vaccination_state}]},
             {InfectionState::Infected_Critical, global_params.get<SevereToCritical>()[{age, vaccination_state}]}});
    case InfectionState::Infected_Critical:
        return random_transition(
            infection_state, dt,
            {{InfectionState::Recovered_Infected, global_params.get<CriticalToRecovered>()[{age, vaccination_state}]},
             {InfectionState::Dead, global_params.get<CriticalToDead>()[{age, vaccination_state}]}});
    case InfectionState::Recovered_Carrier: //fallthrough!
    case InfectionState::Recovered_Infected:
        return random_transition(
            infection_state, dt,
            {{InfectionState::Susceptible, global_params.get<RecoveredToSusceptible>()[{age, vaccination_state}]}});
    default:
        return infection_state; //some states don't transition
    }
}

void Location::begin_step(TimeSpan /*dt*/, const GlobalInfectionParameters& global_params, bool generalized)
{
    //cache for next step so it stays constant during the step while subpopulations change
    //otherwise we would have to cache all state changes during a step which uses more memory
    if (m_cells.empty() && m_num_persons == 0) {
        m_cached_exposure_rate = {{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.};
    }
    else if (m_cells.empty()) {
        auto num_carriers              = get_subpopulation(InfectionState::Carrier);
        auto num_infected              = get_subpopulation(InfectionState::Infected);
        auto rel_risk_loc              = compute_rel_risk(generalized);
        m_cached_exposure_rate.array() = std::min(m_parameters.get<MaximumContacts>(), double(m_num_persons)) /
                                         m_num_persons *
                                         (global_params.get<SusceptibleToExposedByCarrier>().array() * num_carriers +
                                          global_params.get<SusceptibleToExposedByInfected>().array() * num_infected) *
                                         rel_risk_loc;
    }
    else {
        for (auto& cell : m_cells) {
            if (cell.num_people == 0) {
                cell.cached_exposure_rate = {{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.};
            }
            else {
                auto rel_risk_loc = compute_rel_risk(generalized);
                cell.cached_exposure_rate.array() =
                    std::min(m_parameters.get<MaximumContacts>(), double(cell.num_people)) / cell.num_people *
                    (global_params.get<SusceptibleToExposedByCarrier>().array() * cell.num_carriers +
                     global_params.get<SusceptibleToExposedByInfected>().array() * cell.num_infected) *
                    rel_risk_loc;
            }
        }
    }
}

void Location::add_person(const Person& p)
{
    ++m_num_persons;
    InfectionState s = p.get_infection_state();
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

void Location::remove_person(const Person& p)
{
    --m_num_persons;
    InfectionState s = p.get_infection_state();
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

void Location::changed_state(const Person& p, InfectionState old_infection_state)
{
    change_subpopulation(old_infection_state, -1);
    change_subpopulation(p.get_infection_state(), +1);
    for (auto i : p.get_cells()) {
        if (old_infection_state == InfectionState::Carrier) {
            --m_cells[i].num_carriers;
        }
        else if (old_infection_state == InfectionState::Infected ||
                 old_infection_state == InfectionState::Infected_Severe ||
                 old_infection_state == InfectionState::Infected_Critical) {
            --m_cells[i].num_infected;
        }
        if (p.get_infection_state() == InfectionState::Carrier) {
            ++m_cells[i].num_carriers;
        }
        else if (p.get_infection_state() == InfectionState::Infected ||
                 p.get_infection_state() == InfectionState::Infected_Severe ||
                 p.get_infection_state() == InfectionState::Infected_Critical) {
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

void Location::set_default_capacity(LocationType type)
{
    if (type == LocationType::Home) {
        m_capacity.persons = 5;
        m_capacity.volume  = 104;
    }
    else if (type == LocationType::School) {
        m_capacity.persons = 600;
        m_capacity.volume  = 150;
    }
    else if (type == LocationType::Work) {
        m_capacity.persons = 10;
        m_capacity.volume  = 518;
    }
    else if (type == LocationType::SocialEvent) {
        m_capacity.persons = 100;
        m_capacity.volume  = 150;
    }
    else if (type == LocationType::BasicsShop) {
        m_capacity.persons = 50;
        m_capacity.volume  = 2000;
    }
    else if (type == LocationType::Hospital) {
        m_capacity.persons = 4;
        m_capacity.volume  = 100;
    }
    else if (type == LocationType::ICU) {
        m_capacity.persons = 3;
        m_capacity.volume  = 100;
    }
    else if (type == LocationType::Car) {
        m_capacity.persons = 5;
        m_capacity.volume  = 4;
    }
    else if (type == LocationType::PublicTransport) {
        m_capacity.persons = 150;
        m_capacity.volume  = 75;
    }
    else {
        m_capacity.persons = 0;
        m_capacity.volume  = 0;
    }
}

int Location::get_capacity(bool use_volume)
{
    if (use_volume) {
        if (m_capacity.volume > 0) {
            return m_capacity.volume;
        }
        else {
            mio::log_warning("Volume capacity not set but called. Using capacity in person number.");
            return m_capacity.persons;
        }
    }
    else {
        return m_capacity.persons;
    }
}

double Location::compute_rel_risk(bool generalized)
{
    double rel_risk_loc = 1.0;
    if (generalized) {
        rel_risk_loc = (m_capacity.volume * m_num_persons) / (8.0 * m_capacity.persons * m_capacity.persons);
    }
    return rel_risk_loc;
}

} // namespace mio
