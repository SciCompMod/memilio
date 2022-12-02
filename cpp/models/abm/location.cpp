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
    , m_cached_exposure_rate({AgeGroup::Count, VirusVariant::Count})
    , m_testing_scheme()
    , m_cells(std::vector<Cell>(num_cells))
{
}

VirusVariant Location::interact(const Person& person, const TimePoint& t, const TimeSpan& dt) const
{
    //auto vaccination_state = person.get_vaccination_state(); // change to get immunity level and use!
    auto age = person.get_age();

    // remove duplicated code!
    // add individual infection probability through immunity level
    if (!person.get_cells().empty()) {
        for (auto cell_index : person.get_cells()) {
            std::pair<VirusVariant, double> local_indiv_trans_prob[static_cast<uint32_t>(VirusVariant::Count)];

            for (uint32_t v = 0; v != static_cast<uint32_t>(VirusVariant::Count); ++v) {
                VirusVariant virus              = static_cast<VirusVariant>(v);
                double local_indiv_trans_prob_v = 0;
                for (uint32_t a = 0; a != static_cast<uint32_t>(AgeGroup::Count); ++a) {
                    local_indiv_trans_prob_v += m_cells[cell_index].m_cached_exposure_rate[{virus, age}] *
                                                m_parameters.get<ContactRates>()[{age, static_cast<AgeGroup>(a)}] *
                                                dt.days() / days(1).days() *
                                                person.get_immunity_level().get_protection_factor(virus, t);
                }
                local_indiv_trans_prob[v] = std::make_pair(virus, local_indiv_trans_prob_v);
            }

            return random_transition(VirusVariant::Count, dt,
                                     local_indiv_trans_prob); // use VirusVariant::Count for no virus submission
        }
        // we need to define what a cell is used for, as the loop may lead to incorrect results for multiple cells
        return VirusVariant::Count;
    }
    else {
        std::pair<VirusVariant, double> local_indiv_trans_prob[static_cast<uint32_t>(VirusVariant::Count)];

        for (uint32_t v = 0; v != static_cast<uint32_t>(VirusVariant::Count); ++v) {
            VirusVariant virus              = static_cast<VirusVariant>(v);
            double local_indiv_trans_prob_v = 0;
            for (uint32_t a = 0; a != static_cast<uint32_t>(AgeGroup::Count); ++a) {
                local_indiv_trans_prob_v += m_cached_exposure_rate[{virus, age}] *
                                            m_parameters.get<ContactRates>()[{age, static_cast<AgeGroup>(a)}] *
                                            dt.days() / days(1).days() *
                                            person.get_immunity_level().get_protection_factor(virus, t);
            }
            local_indiv_trans_prob[v] = std::make_pair(virus, local_indiv_trans_prob_v);
        }

        return random_transition(VirusVariant::Count, dt,
                                 local_indiv_trans_prob); // use VirusVariant::Count for no virus submission
    }
}

void Location::begin_step(TimePoint t, TimeSpan dt)
{
    //cache for next step so it stays constant during the step while subpopulations change
    //otherwise we would have to cache all state changes during a step which uses more memory

    if (m_cells.empty()) {
        m_cached_exposure_rate = {{VirusVariant::Count, AgeGroup::Count}, 0.};
        for (auto&& p = m_persons.begin(); p != m_persons.end(); ++p) {
            if (p->is_infected(t)) {
                auto inf   = p->get_infection();
                auto virus = inf.get_virus_variant();
                auto age   = p->get_age();
                m_cached_exposure_rate[{virus, age}] +=
                    inf.get_infectivity(t + dt / 2); // average infectivity over the time step to second order accuracy
            }
        }
    }
    else {
        for (auto& cell : m_cells) {
            cell.m_cached_exposure_rate = {{VirusVariant::Count, AgeGroup::Count}, 0.};
            for (auto&& p = cell.m_persons.begin(); p != m_persons.end(); ++p) {
                if (p->is_infected(t)) {
                    auto inf   = p->get_infection();
                    auto virus = inf.get_virus_variant();
                    auto age   = p->get_age();
                    cell.m_cached_exposure_rate[{virus, age}] += inf.get_infectivity(
                        t + dt / 2); // average infectivity over the time step to second order accuracy
                }
            }
        }
    }
}

void Location::add_person(const Person& p)
{
    m_persons.push_back(p);
    if (!m_cells.empty()) {
        for (auto i : p.get_cells()) {
            m_cells[i].m_persons.push_back(p);
        }
    }
}

void Location::remove_person(const Person& p)
{
    auto it = std::remove(m_persons.begin(), m_persons.end(), p);
    m_persons.erase(it, m_persons.end());
}

int Location::get_subpopulation(TimePoint t, InfectionState state) const
{
    // cells are not accounted for
    return count_if(m_persons.begin(), m_persons.end(), [&](Person p) {
        return p.get_infection_state(t) == state;
    });
}

int Location::get_number_infected_total(TimePoint t) const
{
    // cells are not accounted for
    return get_subpopulation(t, InfectionState::Carrier) + get_subpopulation(t, InfectionState::Infected) +
           get_subpopulation(t, InfectionState::Infected_Critical) +
           get_subpopulation(t, InfectionState::Infected_Severe);
}

Eigen::Ref<const Eigen::VectorXi> Location::get_subpopulations(TimePoint t) const
{
    std::array<int, size_t(InfectionState::Count)> subpopulations;
    for (uint32_t i = 0; i < subpopulations.size(); ++i) {
        subpopulations[i] = get_subpopulation(t, static_cast<InfectionState>(i));
    }
    return Eigen::Map<const Eigen::VectorXi>(subpopulations.data(), subpopulations.size());
}

} // namespace abm
} // namespace mio
