/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Elisabeth Kluth, David Kerkmann, Sascha Korf, Martin J. Kuehn, Khoa Nguyen
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
#include "memilio/utils/random_number_generator.h"

namespace mio
{
namespace abm
{

TestingCriteria::TestingCriteria(const std::vector<AgeGroup>& ages, const std::vector<InfectionState>& infection_states)
{
    for (auto age : ages) {
        m_ages.set(static_cast<size_t>(age), true);
    }
    for (auto infection_state : infection_states) {
        m_infection_states.set(static_cast<size_t>(infection_state), true);
    }
}

bool TestingCriteria::operator==(const TestingCriteria& other) const
{
    return m_ages == other.m_ages && m_infection_states == other.m_infection_states;
}

void TestingCriteria::add_age_group(const AgeGroup age_group)
{
    m_ages.set(static_cast<size_t>(age_group), true);
}

void TestingCriteria::remove_age_group(const AgeGroup age_group)
{
    m_ages.set(static_cast<size_t>(age_group), false);
}

void TestingCriteria::add_infection_state(const InfectionState infection_state)
{
    m_infection_states.set(static_cast<size_t>(infection_state), true);
}

void TestingCriteria::remove_infection_state(const InfectionState infection_state)
{
    m_infection_states.set(static_cast<size_t>(infection_state), false);
}


} // namespace abm
} // namespace mio
