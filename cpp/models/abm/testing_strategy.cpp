/*
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: Elisabeth Kluth, David Kerkmann, Sascha Korf
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

TestingCriteria::TestingCriteria(const std::vector<AgeGroup>& ages, const std::vector<LocationType>& location_types,
                                 const std::vector<InfectionState>& infection_states)
    : m_ages(ages)
    , m_location_types(location_types)
    , m_infection_states(infection_states)
{
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

bool TestingCriteria::evaluate(const Person& p, const Location& l) const
{
    return has_requested_age(p) && is_requested_location_type(l) && has_requested_infection_state(p);
}

bool TestingCriteria::has_requested_age(const Person& p) const
{
    if (m_ages.empty()) {
        return true; // no condition on the age
    }
    return std::find(m_ages.begin(), m_ages.end(), p.get_age()) != m_ages.end();
}

bool TestingCriteria::is_requested_location_type(const Location& l) const
{
    if (m_location_types.empty()) {
        return true; // no condition on the location
    }
    return std::find(m_location_types.begin(), m_location_types.end(), l.get_type()) != m_location_types.end();
}

bool TestingCriteria::has_requested_infection_state(const Person& p) const
{
    if (m_infection_states.empty()) {
        return true; // no condition on infection state
    }
    return std::find(m_infection_states.begin(), m_infection_states.end(), p.get_infection_state()) !=
           m_infection_states.end();
}

TestingScheme::TestingScheme(const std::vector<TestingCriteria>& testing_criteria,
                             TimeSpan minimal_time_since_last_test, TimePoint start_date, TimePoint end_date,
                             double probability, const GenericTest& test_type)
    : m_testing_criteria(testing_criteria)
    , m_minimal_time_since_last_test(minimal_time_since_last_test)
    , m_probability(probability)
    , m_start_date(start_date)
    , m_end_date(end_date)
    , m_test_type(test_type)
{
}

void TestingScheme::add_testing_criteria(const TestingCriteria criteria)
{
    if (std::find(m_testing_criteria.begin(), m_testing_criteria.end(), criteria) == m_testing_criteria.end()) {
        m_testing_criteria.push_back(criteria);
    }
}

void TestingScheme::remove_testing_criteria(const TestingCriteria criteria)
{
    auto last = std::remove(m_testing_criteria.begin(), m_testing_criteria.end(), criteria);
    m_testing_criteria.erase(last, m_testing_criteria.end());
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
    if (person.get_time_since_negative_test() > m_minimal_time_since_last_test) {
        double random = UniformDistribution<double>::get_instance()();
        if (random < m_probability) {
            if (std::any_of(m_testing_criteria.begin(), m_testing_criteria.end(),
                            [person, location](TestingCriteria tr) {
                                return tr.evaluate(person, location);
                            })) {
                return !person.get_tested(m_test_type.get_default());
            }
        }
    }
    return true;
}

TestingStrategy::TestingStrategy(const std::vector<TestingScheme>& testing_schemes)
    : m_testing_schemes(testing_schemes)
{
}

void TestingStrategy::add_testing_scheme(const TestingScheme& scheme)
{
    if (std::find(m_testing_schemes.begin(), m_testing_schemes.end(), scheme) == m_testing_schemes.end()) {
        m_testing_schemes.push_back(scheme);
    }
}

void TestingStrategy::remove_testing_scheme(const TestingScheme& scheme)
{
    auto last = std::remove(m_testing_schemes.begin(), m_testing_schemes.end(), scheme);
    m_testing_schemes.erase(last, m_testing_schemes.end());
}

bool TestingStrategy::run_strategy(Person& person, const Location& location) const
{
    return std::all_of(m_testing_schemes.begin(), m_testing_schemes.end(), [&person, location](TestingScheme ts) {
        if (ts.is_active()) {
            return ts.run_scheme(person, location);
        }
        return true;
    });
}

void TestingStrategy::update_testing_scheme_activity_status(const TimePoint t)
{
    for (auto& ts : m_testing_schemes) {
        ts.update_activity_status(t);
    }
}

} // namespace abm
} // namespace mio
