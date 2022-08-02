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

#include "abm/testing_scheme.h"
#include "memilio/utils/random_number_generator.h"

namespace mio
{
namespace abm
{

TestingScheme::TestingScheme(const std::vector<TestingRule> testing_rules, const TimeSpan interval,
                             const double probability)
    : m_testing_rules(testing_rules)
    , m_testing_frequency(interval)
    , m_probability(probability)
{
}

const TimeSpan& TestingScheme::get_interval() const
{
    return m_testing_frequency;
}

double TestingScheme::get_probability() const
{
    return m_probability;
}

void TestingScheme::set_interval(TimeSpan t)
{
    m_testing_frequency = t;
}

void TestingScheme::set_probability(double p)
{
    m_probability = p;
}

void TestingScheme::add_testing_rule(const TestingRule rule)
{
    m_testing_rules.push_back(rule);
    std::unique(m_testing_rules.begin(), m_testing_rules.end());
}

void TestingScheme::remove_testing_rule(const TestingRule rule)
{
    std::remove(m_testing_rules.begin(), m_testing_rules.end(), rule);
}

// const TimePoint& TestingScheme::get_start_date() const
// {
//     return m_start_date;
// }

// const TimePoint& TestingScheme::get_end_date() const
// {
//     return m_end_date;
// }

// const TimeSpan TestingScheme::get_duration() const
// {
//     return TimeSpan(m_end_date.seconds() - m_start_date.seconds());
// }

bool TestingScheme::is_active() const
{
    return m_is_active;
}
void TestingScheme::update_activity_status(const TimePoint t)
{
    m_is_active = (m_start_date <= t && t <= m_end_date);
}

bool TestingScheme::run_scheme(Person& person, const Location& location) const
{
    if (person.get_time_since_negative_test() > m_testing_frequency) {
        double random = UniformDistribution<double>::get_instance()();
        if (random < m_probability) {
            if (std::any_of(m_testing_rules.begin(), m_testing_rules.end(), [person, location](TestingRule tr) {
                    return tr.evaluate(person, location);
                })) {
                return !person.get_tested(m_test_type.get_default());
            }
        }
    }
    return true;
}

} // namespace abm
} // namespace mio
