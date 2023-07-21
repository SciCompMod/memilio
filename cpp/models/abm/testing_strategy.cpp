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

bool TestingCriteria::evaluate(const Person& p, const Location& l, TimePoint t) const
{
    return has_requested_age(p) && is_requested_location_type(l) && has_requested_infection_state(p, t);
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

bool TestingCriteria::has_requested_infection_state(const Person& p, TimePoint t) const
{
    if (m_infection_states.empty()) {
        return true; // no condition on infection state
    }
    return std::find(m_infection_states.begin(), m_infection_states.end(), p.get_infection_state(t)) !=
           m_infection_states.end();
}

TestingScheme::TestingScheme(const std::vector<TestingCriteria>& testing_criteria,
                             TimeSpan minimal_time_since_last_test, TimePoint start_date, TimePoint end_date,
                             const GenericTest& test_type, double probability)
    : m_testing_criteria(testing_criteria)
    , m_minimal_time_since_last_test(minimal_time_since_last_test)
    , m_start_date(start_date)
    , m_end_date(end_date)
    , m_test_type(test_type)
    , m_probability(probability)
{
}

bool TestingScheme::operator==(const TestingScheme& other) const
{
    return this->m_testing_criteria == other.m_testing_criteria &&
           this->m_minimal_time_since_last_test == other.m_minimal_time_since_last_test &&
           this->m_start_date == other.m_start_date && this->m_end_date == other.m_end_date &&
           this->m_test_type.get_default().sensitivity == other.m_test_type.get_default().sensitivity &&
           this->m_test_type.get_default().specificity == other.m_test_type.get_default().specificity &&
           this->m_probability == other.m_probability;
    //To be adjusted and also TestType should be static.
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
void TestingScheme::update_activity_status(TimePoint t)
{
    m_is_active = (m_start_date <= t && t <= m_end_date);
}

bool TestingScheme::run_scheme(Person& person, const Location& location, TimePoint t) const
{
    if (person.get_time_since_negative_test() > m_minimal_time_since_last_test) {
        double random = UniformDistribution<double>::get_instance()();
        if (random < m_probability) {
            if (std::any_of(m_testing_criteria.begin(), m_testing_criteria.end(),
                            [person, location, t](TestingCriteria tr) {
                                return tr.evaluate(person, location, t);
                            })) {
                return !person.get_tested(t, m_test_type.get_default());
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

void TestingStrategy::update_activity_status(TimePoint t)
{
    for (auto& ts : m_testing_schemes) {
        ts.update_activity_status(t);
    }
}

bool TestingStrategy::run_strategy(Person& person, const Location& location, TimePoint t) const
{
    // Person who is in quarantine but not yet home should go home. Otherwise they can't because they test positive.
    if (location.get_type() == mio::abm::LocationType::Home && person.is_in_quarantine()) {
        return true;
    }
    return std::all_of(m_testing_schemes.begin(), m_testing_schemes.end(), [&person, location, t](TestingScheme ts) {
        if (ts.is_active()) {
            return ts.run_scheme(person, location, t);
        }
        return true;
    });
}

} // namespace abm
} // namespace mio
