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

TestingRule::TestingRule(const std::vector<AgeGroup> ages, const std::vector<LocationType> location_types,
                         const std::vector<InfectionState> infection_states)
    : m_ages(ages)
    , m_location_types(location_types)
    , m_infection_states(infection_states)
{
}

void TestingRule::add_age_group(const AgeGroup age_group)
{
    if (std::find(m_ages.begin(), m_ages.end(), age_group) == m_ages.end()) {
        m_ages.push_back(age_group);
    }
}

void TestingRule::remove_age_group(const AgeGroup age_group)
{
    auto last = std::remove(m_ages.begin(), m_ages.end(), age_group);
    m_ages.erase(last, m_ages.end());
}

void TestingRule::add_location_type(const LocationType location_type)
{
    if (std::find(m_location_types.begin(), m_location_types.end(), location_type) == m_location_types.end()) {
        m_location_types.push_back(location_type);
    }
}
void TestingRule::remove_location_type(const LocationType location_type)
{
    auto last = std::remove(m_location_types.begin(), m_location_types.end(), location_type);
    m_location_types.erase(last, m_location_types.end());
}

void TestingRule::add_infection_state(const InfectionState infection_state)
{
    if (std::find(m_infection_states.begin(), m_infection_states.end(), infection_state) == m_infection_states.end()) {
        m_infection_states.push_back(infection_state);
    }
}

void TestingRule::remove_infection_state(const InfectionState infection_state)
{
    auto last = std::remove(m_infection_states.begin(), m_infection_states.end(), infection_state);
    m_infection_states.erase(last, m_infection_states.end());
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

TestingScheme::TestingScheme(const std::vector<TestingRule> testing_rules, const TimeSpan testing_frequency,
                             TimePoint start_date, TimePoint end_date, const double probability,
                             const GenericTest& test_type)
    : m_testing_rules(testing_rules)
    , m_testing_frequency(testing_frequency)
    , m_probability(probability)
    , m_start_date(start_date)
    , m_end_date(end_date)
    , m_test_type(test_type)
{
}

void TestingScheme::add_testing_rule(const TestingRule rule)
{
    if (std::find(m_testing_rules.begin(), m_testing_rules.end(), rule) == m_testing_rules.end()) {
        m_testing_rules.push_back(rule);
    }
}

void TestingScheme::remove_testing_rule(const TestingRule rule)
{
    auto last = std::remove(m_testing_rules.begin(), m_testing_rules.end(), rule);
    m_testing_rules.erase(last, m_testing_rules.end());
}

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
