/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#include "memilio/utils/random_number_generator.h"
#include "abm/random_events.h"

#include <numeric>

namespace mio
{
namespace abm
{

Location::Location(LocationType type, uint32_t index, uint32_t num_cells)
    : m_type(type)
    , m_index(index)
    , m_capacity_adapted_transmission_risk(false)
    , m_subpopulations(Eigen::Index(InfectionState::Count))
    , m_cached_exposure_rate({AgeGroup::Count, VaccinationState::Count})
    , m_cells(std::vector<Cell>(num_cells))
    , m_required_mask(MaskType::Community)
    , m_npi_active(false)
{
    m_subpopulations.add_time_point(0);
    m_subpopulations.get_last_value().setZero();
}

VirusVariant Location::interact(const Person& person, const TimePoint& t, const TimeSpan& dt,
                                const GlobalInfectionParameters& global_params) const
{
    auto infection_state       = person.get_infection_state();
    auto vaccination_state     = person.get_vaccination_state();
    auto age                   = person.get_age();
    ScalarType mask_protection = person.get_protective_factor(global_params);
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
            local_indiv_trans_prob[v] = std::make_pair(virus, local_indiv_trans_prob_v);
        }

        return random_transition(VirusVariant::Count, dt,
                                 local_indiv_trans_prob); // use VirusVariant::Count for no virus submission
    }
    // we need to define what a cell is used for, as the loop may lead to incorrect results for multiple cells
    return VirusVariant::Count;
}

void Location::begin_step(TimePoint t, TimeSpan dt)
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
        m_cached_exposure_rate.array() =
            relative_transmission_risk * std::min(m_parameters.get<MaximumContacts>(), ScalarType(m_num_persons)) /
            m_num_persons *
            (global_params.get<SusceptibleToExposedByCarrier>().array().cast<ScalarType>() * num_carriers +
             global_params.get<SusceptibleToExposedByInfected>().array().cast<ScalarType>() * num_infected);
    }
    else {
        for (auto& cell : m_cells) {
            if (cell.num_people == 0) {
                cell.cached_exposure_rate = {{AgeGroup::Count, VaccinationState::Count}, 0.};
            }
            else {
                auto relative_transmission_risk = compute_relative_transmission_risk();
                cell.cached_exposure_rate.array() =
                    std::min(m_parameters.get<MaximumContacts>(), ScalarType(cell.num_people)) / cell.num_people *
                    (global_params.get<SusceptibleToExposedByCarrier>().array().cast<ScalarType>() * cell.num_carriers +
                     global_params.get<SusceptibleToExposedByInfected>().array().cast<ScalarType>() *
                         cell.num_infected) *
                    relative_transmission_risk;
            }
        }
        if (m_capacity_adapted_transmission_risk) {
            cell.m_cached_exposure_rate_air.array() *= cell.compute_relative_cell_size();
        }
    }
}

void Location::add_person(const Person& p, const uint32_t cell_idx)
{
    m_cells[cell_idx].m_persons.push_back(p);
}

void Location::remove_person(const Person& p)
{
    for (auto&& cell : m_cells) {
        auto it = std::remove(cell.m_persons.begin(), cell.m_persons.end(), p);
        cell.m_persons.erase(it, cell.m_persons.end());
    }
}

int Cell::get_subpopulation(const TimePoint& t, const InfectionState& state) const
{
    return count_if(m_persons.begin(), m_persons.end(), [&](Person p) {
        return p.get_infection_state(t) == state;
    });
}

int Location::get_subpopulation(const TimePoint& t, const InfectionState& state, const uint32_t cell_idx) const
{
    return m_cells[cell_idx].get_subpopulation(t, state);
}

int Location::get_subpopulation(const TimePoint& t, const InfectionState& state) const
{
    int n_persons = 0;
    for (auto&& cell : m_cells) {
        n_persons += cell.get_subpopulation(t, state);
    }
    return n_persons;
}

int Location::get_number_infected_total(TimePoint t) const
{
    m_subpopulations.get_last_value()[size_t(s)] += delta;
    assert(m_subpopulations.get_last_value()[size_t(s)] >= 0 && "subpopulations must be non-negative");
}

Eigen::Ref<const Eigen::VectorXi> Location::get_subpopulations(TimePoint t) const
{
    return (int)m_subpopulations.get_last_value()[size_t(s)];
}
/*
For every cell in a location we have a transmission factor that is nomalized to m_capacity.volume / m_capacity.persons of 
the location "Home", which is 66. We multiply this rate with the individul size of each cell to obtain a "space per person" factor.
*/
ScalarType Location::compute_relative_transmission_risk()
{
    if (m_capacity.volume != 0) {
        return 66.0 / m_capacity.volume;
    }
    else {
        return 1.0;
    }
}

void Location::add_subpopulations_timepoint(const TimePoint& t)
{
    // Get the previous time point index.
    // Since index starts from 0, we need to -1 from the current time point.
    auto last_idx = m_subpopulations.get_num_time_points() - 1;
    m_subpopulations.add_time_point(t.days());
    m_subpopulations.get_last_value() = m_subpopulations.get_value(last_idx);
}

void Location::initialize_subpopulation(const TimePoint& t)
{
    // We make a copy of subpoulation to the new TimePoint t and delete the first TimePoint.
    add_subpopulations_timepoint(t);
    m_subpopulations.remove_time_point(0);
}

} // namespace abm
} // namespace mio
