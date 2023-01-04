/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Elisabeth Kluth, Carlotta Gerstein, Martin J. Kuehn
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
    , m_capacity(LocationCapacity())
    , m_capacity_adapted_transmission_risk(false)
    , m_subpopulations{}
    , m_cached_exposure_rate({AgeGroup::Count, VaccinationState::Count})
    , m_cells(std::vector<Cell>(num_cells))
    , m_required_mask(MaskType::Community)
    , m_npi_active(false)
{
}

InfectionState Location::interact(const Person& person, TimeSpan dt,
                                  const GlobalInfectionParameters& global_params) const
{
    auto infection_state   = person.get_infection_state();
    auto vaccination_state = person.get_vaccination_state();
    auto age               = person.get_age();
    double mask_protection = person.get_protective_factor(global_params);
    switch (infection_state) {
    case InfectionState::Susceptible:
        if (!person.get_cells().empty()) {
            for (auto cell_index : person.get_cells()) {
                InfectionState new_state = random_transition(
                    infection_state, dt,
                    {{InfectionState::Exposed,
                      (1 - mask_protection) * m_cells[cell_index].cached_exposure_rate[{age, vaccination_state}]}});
                if (new_state != infection_state) {
                    return new_state;
                }
            }
            return infection_state;
        }
        else {
            return random_transition(
                infection_state, dt,
                {{InfectionState::Exposed, (1 - mask_protection) * m_cached_exposure_rate[{age, vaccination_state}]}});
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

void Location::begin_step(TimeSpan /*dt*/, const GlobalInfectionParameters& global_params)
{
    //cache for next step so it stays constant during the step while subpopulations change
    //otherwise we would have to cache all state changes during a step which uses more memory
    if (m_cells.empty() && m_num_persons == 0) {
        m_cached_exposure_rate = {{AgeGroup::Count, VaccinationState::Count}, 0.};
    }
    else if (m_cells.empty()) {
        auto num_carriers               = get_subpopulation(InfectionState::Carrier);
        auto num_infected               = get_subpopulation(InfectionState::Infected);
        auto relative_transmission_risk = compute_relative_transmission_risk();
        m_cached_exposure_rate.array()  = relative_transmission_risk *
                                         std::min(m_parameters.get<MaximumContacts>(), double(m_num_persons)) /
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
                auto relative_transmission_risk = compute_relative_transmission_risk();
                cell.cached_exposure_rate.array() =
                    std::min(m_parameters.get<MaximumContacts>(), double(cell.num_people)) / cell.num_people *
                    (global_params.get<SusceptibleToExposedByCarrier>().array() * cell.num_carriers +
                     global_params.get<SusceptibleToExposedByInfected>().array() * cell.num_infected) *
                    relative_transmission_risk;
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

/*
For every location we have a transmission factor that is nomalized to m_capacity.volume / m_capacity.persons of 
the location "Home", which is 66. We multiply this rate with the factor m_capacity.persons / m_capacity.volume of the 
considered location. For the dependency of the transmission risk on the occupancy we use a linear approximation which 
leads to: 66 * (m_capacity.persons / m_capacity.volume) * (m_num_persons / m_capacity.persons).
*/
double Location::compute_relative_transmission_risk()
{
    if (m_capacity.volume != 0 && m_capacity_adapted_transmission_risk) {
        return 66.0 * m_num_persons / m_capacity.volume;
    }
    else {
        return 1.0;
    }
}

} // namespace abm
} // namespace mio
