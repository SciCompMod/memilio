/*
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: Elisabeth Kluth, David Kerkmann, Sascha Korf, Martin J. Kuehn
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

#include "abm/testing_strategy.h"


namespace mio
{
namespace abm
{

TestingCriteria::TestingCriteria(const std::vector<AgeGroup>& ages, const std::vector<LocationType>& location_types,
                                 const std::vector<InfectionState>& infection_states)
    : m_ages(ages)
    , m_location_types(location_types)
    , m_infection_states(infection_states)
{
}

bool TestingCriteria::operator==(TestingCriteria other) const
{
    auto to_compare_ages             = this->m_ages;
    auto to_compare_infection_states = this->m_infection_states;
    auto to_compare_location_types   = this->m_location_types;

    std::sort(to_compare_ages.begin(), to_compare_ages.end());
    std::sort(other.m_ages.begin(), other.m_ages.end());
    std::sort(to_compare_infection_states.begin(), to_compare_infection_states.end());
    std::sort(other.m_infection_states.begin(), other.m_infection_states.end());
    std::sort(to_compare_location_types.begin(), to_compare_location_types.end());
    std::sort(other.m_location_types.begin(), other.m_location_types.end());

    return to_compare_ages == other.m_ages && to_compare_location_types == other.m_location_types &&
           to_compare_infection_states == other.m_infection_states;
}

void TestingCriteria::add_age_group(const AgeGroup age_group)
{
    if (std::find(m_ages.begin(), m_ages.end(), age_group) == m_ages.end()) {
        m_ages.push_back(age_group);
    }
}

void TestingCriteria::remove_age_group(const AgeGroup age_group)
{
    auto last = std::remove(m_ages.begin(), m_ages.end(), age_group);
    m_ages.erase(last, m_ages.end());
}

void TestingCriteria::add_location_type(const LocationType location_type)
{
    if (std::find(m_location_types.begin(), m_location_types.end(), location_type) == m_location_types.end()) {
        m_location_types.push_back(location_type);
    }
}
void TestingCriteria::remove_location_type(const LocationType location_type)
{
    auto last = std::remove(m_location_types.begin(), m_location_types.end(), location_type);
    m_location_types.erase(last, m_location_types.end());
}

void TestingCriteria::add_infection_state(const InfectionState infection_state)
{
    if (std::find(m_infection_states.begin(), m_infection_states.end(), infection_state) == m_infection_states.end()) {
        m_infection_states.push_back(infection_state);
    }
}

void TestingCriteria::remove_infection_state(const InfectionState infection_state)
{
    auto last = std::remove(m_infection_states.begin(), m_infection_states.end(), infection_state);
    m_infection_states.erase(last, m_infection_states.end());
}

} // namespace abm
} // namespace mio
