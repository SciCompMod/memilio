/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: Sascha Korf
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
#include "abm/testing_scheme.h"

namespace mio
{
namespace abm
{

TestingScheme::TestingScheme(const std::vector<TestRule>& rules)
    : m_test_rules(rules)
    , : m_active(true)
{
}

TestingScheme::add_rule(const TestRule& rule)
{
    m_test_rules.push_back(rule);
    if (!m_active) { //Inactive Testing scheme becomes active again with a new rule
        m_active = true;
    }
}

TestingScheme::check_rules_for_activeness() const
{
    m_active = false;
    for (const TestRule& rule : m_test_rules) {
        if (rule.get_active_status()) {
            m_active = true;
        }
    }
}

TestingScheme::run_scheme(const Person& person, const Time& time, const LocationType& location_type)
{
}

} // namespace abm
} // namespace mio
