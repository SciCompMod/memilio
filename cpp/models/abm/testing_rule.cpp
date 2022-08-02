/*
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: David Kerkmann
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

#include "abm/testing_rule.h"

namespace mio
{
namespace abm
{

TestingRule::TestingRule(const std::vector<AgeGroup> ages, const std::vector<LocationType> location_types,
                         const std::vector<InfectionState> infection_states)
    : m_ages(ages)
    , m_location_types(location_types)
    , m_infection_states(infection_states)
{
}

void TestingRule::add_age_group(const AgeGroup ageGroup)
{
    m_ages.push_back(ageGroup);
    std::unique(m_ages.begin(), m_ages.end());
}

void TestingRule::remove_age_group(const AgeGroup ageGroup)
{
    std::remove(m_ages.begin(), m_ages.end(), ageGroup);
}

void TestingRule::add_location_type(const LocationType locationType)
{
    m_location_types.push_back(locationType);
    std::unique(m_location_types.begin(), m_location_types.end());
}
void TestingRule::remove_location_type(const LocationType locationType)
{
    std::remove(m_location_types.begin(), m_location_types.end(), locationType);
}

void TestingRule::add_infection_state(const InfectionState infection_state)
{
    m_infection_states.push_back(infection_state);
    std::unique(m_infection_states.begin(), m_infection_states.end());
}

void TestingRule::remove_infection_state(const InfectionState infection_state)
{
    std::remove(m_infection_states.begin(), m_infection_states.end(), infection_state);
}

bool TestingRule::evaluate(const Person& p, const Location& l) const
{
    return has_requested_age(p) && is_requested_location_type(l) && has_requested_infection_state(p);
}

bool TestingRule::has_requested_age(const Person& p) const
{
    if (m_ages.empty()) {
        return true; // no condition on the age
    }
    return std::find(m_ages.begin(), m_ages.end(), p.get_age()) != m_ages.end();
}

bool TestingRule::is_requested_location_type(const Location& l) const
{
    if (m_location_types.empty()) {
        return true; // no condition on the location
    }
    return std::find(m_location_types.begin(), m_location_types.end(), l.get_type()) != m_location_types.end();
}

bool TestingRule::has_requested_infection_state(const Person& p) const
{
    if (m_infection_states.empty()) {
        return true; // no condition on infection state
    }
    return std::find(m_infection_states.begin(), m_infection_states.end(), p.get_infection_state()) !=
           m_infection_states.end();
}

} // namespace abm
} // namespace mio
