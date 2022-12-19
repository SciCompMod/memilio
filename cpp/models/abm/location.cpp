/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Elisabeth Kluth, David Kerkmann
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
    , m_cells(std::vector<Cell>(num_cells))
{
}

VirusVariant Location::interact(const Person& person, const TimePoint& t, const TimeSpan& dt) const
{
    auto age = person.get_age();
    for (auto cell_index : person.get_cells()) { // should be changed so that persons can only be in one cell
        std::pair<VirusVariant, double> local_indiv_trans_prob[static_cast<uint32_t>(VirusVariant::Count)];

        for (uint32_t v = 0; v != static_cast<uint32_t>(VirusVariant::Count); ++v) {
            VirusVariant virus              = static_cast<VirusVariant>(v);
            double local_indiv_trans_prob_v = 0;
            for (uint32_t a = 0; a != static_cast<uint32_t>(AgeGroup::Count); ++a) {
                local_indiv_trans_prob_v += (m_cells[cell_index].m_cached_exposure_rate_contacts[{virus, age}] *
                                                 m_parameters.get<ContactRates>()[{age, static_cast<AgeGroup>(a)}] +
                                             m_cells[cell_index].m_cached_exposure_rate_air[{virus}] *
                                                 m_parameters.get<AerosolTransmissionRates>()[{virus}]) *
                                            dt.days() / days(1).days() * person.get_protection_factor(virus, t);
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
    for (auto& cell : m_cells) {
        cell.m_cached_exposure_rate_contacts = {{VirusVariant::Count, AgeGroup::Count}, 0.};
        cell.m_cached_exposure_rate_air      = {{VirusVariant::Count}, 0.};
        for (auto&& p = cell.m_persons.begin(); p != cell.m_persons.end(); ++p) {
            if (p->is_infected(t)) {
                auto inf   = p->get_infection();
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
    return get_subpopulation(t, InfectionState::Carrier) + get_subpopulation(t, InfectionState::Infected) +
           get_subpopulation(t, InfectionState::Infected_Critical) +
           get_subpopulation(t, InfectionState::Infected_Severe);
}

Eigen::Ref<const Eigen::VectorXi> Location::get_subpopulations(TimePoint t) const
{
    std::array<int, size_t(InfectionState::Count)> subpopulations;
    for (uint32_t i = 0; i < subpopulations.size(); ++i) {
        for (auto&& cell : m_cells) {
            subpopulations[i] = cell.get_subpopulation(t, static_cast<InfectionState>(i));
        }
    }
    return Eigen::Map<const Eigen::VectorXi>(subpopulations.data(), subpopulations.size());
}
/*
For every cell in a location we have a transmission factor that is nomalized to m_capacity.volume / m_capacity.persons of 
the location "Home", which is 66. We multiply this rate with the individul size of each cell to obtain a "space per person" factor.
*/
double Cell::compute_relative_cell_size()
{
    if (m_capacity.volume != 0) {
        return 66.0 / m_capacity.volume;
    }
    else {
        return 1.0;
    }
}
} // namespace abm
} // namespace mio
