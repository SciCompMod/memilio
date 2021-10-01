/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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

namespace epi
{

Location::Location(LocationType type, uint32_t index)
    : m_type(type)
    , m_index(index)
    , m_subpopulations{}
    , m_cached_exposure_rate({AbmAgeGroup::Count, epi::VaccinationState::Count})
    , m_testing_scheme()
{
}

InfectionState Location::interact(const Person& person, TimeSpan dt, const GlobalInfectionParameters& global_params) const
{
    auto infection_state = person.get_infection_state();
    auto vaccination_state = person.get_vaccination_state();
    auto age = person.get_age();
        switch (infection_state) {
        case InfectionState::Susceptible:
                return random_transition(infection_state, dt, {{InfectionState::Exposed, m_cached_exposure_rate[{age,vaccination_state}]}});
        case InfectionState::Carrier:
            return random_transition(
                infection_state, dt,
                {{InfectionState::Infected, global_params.get<CarrierToInfected>()[{age,vaccination_state}]},
                 {InfectionState::Recovered_Carrier, global_params.get<CarrierToRecovered>()[{age,vaccination_state}]}});
        case InfectionState::Infected:
            return random_transition(
                infection_state, dt,
                {{InfectionState::Recovered_Infected, global_params.get<InfectedToRecovered>()[{age,vaccination_state}]},
                 {InfectionState::Infected_Severe, global_params.get<InfectedToSevere>()[{age,vaccination_state}]}});
        case InfectionState::Infected_Severe:
            return random_transition(
                infection_state, dt,
                {{InfectionState::Recovered_Infected, global_params.get<SevereToRecovered>()[{age,vaccination_state}]},
                 {InfectionState::Infected_Critical, global_params.get<SevereToCritical>()[{age,vaccination_state}]}});
        case InfectionState::Infected_Critical:
            return random_transition(
                infection_state, dt,
                {{InfectionState::Recovered_Infected, global_params.get<CriticalToRecovered>()[{age,vaccination_state}]},
                 {InfectionState::Dead, global_params.get<CriticalToDead>()[{age,vaccination_state}]}});
        case InfectionState::Recovered_Carrier: //fallthrough!
        case InfectionState::Recovered_Infected:
            return random_transition(
                infection_state, dt, {{InfectionState::Susceptible, global_params.get<RecoveredToSusceptible>()[{age,vaccination_state}]}});
        default:
            return infection_state; //some states don't transition
        }
    
}

void Location::begin_step(TimeSpan /*dt*/, const GlobalInfectionParameters& global_params)
{
    //cache for next step so it stays constant during the step while subpopulations change
    //otherwise we would have to cache all state changes during a step which uses more memory
    auto num_carriers = get_subpopulation(InfectionState::Carrier);
    auto num_infected = get_subpopulation(InfectionState::Infected);
    if (m_num_persons == 0){
        m_cached_exposure_rate = {{epi::AbmAgeGroup::Count, epi::VaccinationState::Count}, 0.};
    } 
    else{
        m_cached_exposure_rate.array()
            = std::min(m_parameters.get<EffectiveContacts>(), double(m_num_persons)) / m_num_persons *
                             (global_params.get<SusceptibleToExposedByCarrier>().array() * num_carriers +
                              global_params.get<SusceptibleToExposedByInfected>().array() * num_infected);
    }
}


void Location::add_person(const Person& p)
{
    ++m_num_persons;
    InfectionState s = p.get_infection_state();
    change_subpopulation(s, +1);
}

void Location::remove_person(const Person& p)
{
    --m_num_persons;
    InfectionState s = p.get_infection_state();
    change_subpopulation(s, -1);
}

void Location::changed_state(const Person& p, InfectionState old_infection_state)
{
    change_subpopulation(old_infection_state, -1);
    change_subpopulation(p.get_infection_state(), +1);
}

void Location::change_subpopulation(InfectionState s, int delta)
{
    m_subpopulations[size_t(s)] += delta;
    assert(m_subpopulations[size_t(s)]>=0 && "subpopulations must be non-negative");
}

int Location::get_subpopulation(InfectionState s) const
{
    return m_subpopulations[size_t(s)];
}

Eigen::Ref<const Eigen::VectorXi> Location::get_subpopulations() const
{
    return Eigen::Map<const Eigen::VectorXi>(m_subpopulations.data(), m_subpopulations.size());
}


} // namespace epi
