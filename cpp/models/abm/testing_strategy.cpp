/*
* Copyright (C) 2020-2025 MEmilio
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
#include "abm/location_id.h"
#include "memilio/utils/random_number_generator.h"
#include <utility>

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

bool TestingCriteria::evaluate(const Person& p, TimePoint t) const
{
    // An empty vector of ages or none bitset of #InfectionStates% means that no condition on the corresponding property is set.
    return (m_ages.none() || m_ages[static_cast<size_t>(p.get_age())]) &&
           (m_infection_states.none() || m_infection_states[static_cast<size_t>(p.get_infection_state(t))]);
}

TestingScheme::TestingScheme(const TestingCriteria& testing_criteria, TimeSpan validity_period, TimePoint start_date,
                             TimePoint end_date, TestParameters test_parameters, ScalarType probability)
    : m_testing_criteria(testing_criteria)
    , m_validity_period(validity_period)
    , m_start_date(start_date)
    , m_end_date(end_date)
    , m_test_parameters(test_parameters)
    , m_probability(probability)
{
    assert(start_date <= end_date && "Start date must be before or equal to end date");
}

bool TestingScheme::operator==(const TestingScheme& other) const
{
    return this->m_testing_criteria == other.m_testing_criteria && this->m_validity_period == other.m_validity_period &&
           this->m_start_date == other.m_start_date && this->m_end_date == other.m_end_date &&
           this->m_test_parameters.sensitivity == other.m_test_parameters.sensitivity &&
           this->m_test_parameters.specificity == other.m_test_parameters.specificity &&
           this->m_probability == other.m_probability;
    //To be adjusted and also TestType should be static.
}

bool TestingScheme::is_active(TimePoint t) const
{
    return (m_start_date <= t && t < m_end_date);
}

bool TestingScheme::run_and_test(PersonalRandomNumberGenerator& rng, Person& person, TimePoint t) const
{
    if (!is_active(t)) { // If the scheme is not active, do nothing; early return
        return false;
    }

    auto test_result = person.get_test_result(m_test_parameters.type);
    // If the agent has a test result valid until now, use the result directly
    if ((test_result.time_of_testing > TimePoint(std::numeric_limits<int>::min())) &&
        (test_result.time_of_testing + m_validity_period >= t)) {
        return test_result.result; // If the test is positive, the entry is not allowed, and vice versa
    }
    if (person.get_compliance(InterventionType::Testing) <
            1.0 && // Dont need to draw a random number if the person is compliant either way
        !person.is_compliant(
            rng, InterventionType::Testing)) { // If the person is not compliant with the testing intervention
        return true;
    }
    // Otherwise, the time_of_testing in the past (i.e. the agent has already performed it).
    if (m_testing_criteria.evaluate(person, t - m_test_parameters.required_time)) {
        ScalarType random = UniformDistribution<ScalarType>::get_instance()(rng);
        if (random < m_probability) {
            bool result = person.get_tested(rng, t - m_test_parameters.required_time, m_test_parameters);
            person.add_test_result(t, m_test_parameters.type, result);
            return result; // If the test is positive, the entry is not allowed, and vice versa
        }
    }
    // If the test is not performed, the entry is allowed
    return false;
}

TestingStrategy::TestingStrategy(const std::vector<LocalStrategy>& location_to_schemes_id,
                                 const std::vector<LocalStrategy>& location_to_schemes_type)
    : m_testing_schemes_at_location_id(location_to_schemes_id.begin(), location_to_schemes_id.end())
    , m_testing_schemes_at_location_type(location_to_schemes_type.begin(), location_to_schemes_type.end())
{
}

void TestingStrategy::add_scheme(const LocationId& loc_id, const TestingScheme& scheme)
{
    if (loc_id.get() >= m_testing_schemes_at_location_id.size()) {
        m_testing_schemes_at_location_id.resize(loc_id.get() + 1);
    }
    m_testing_schemes_at_location_id[loc_id.get()].schemes.push_back(scheme);
}

void TestingStrategy::add_scheme(const LocationType& loc_type, const TestingScheme& scheme)
{
    if ((size_t)loc_type >= m_testing_schemes_at_location_type.size()) {
        m_testing_schemes_at_location_type.resize((size_t)loc_type + 1);
    }
    m_testing_schemes_at_location_type[(size_t)loc_type].schemes.push_back(scheme);
}

bool TestingStrategy::run_and_check(PersonalRandomNumberGenerator& rng, Person& person, const Location& location,
                                    TimePoint t)
{
    // Early return if no scheme defined for this location or type
    auto loc_id   = location.get_id().get();
    auto loc_type = static_cast<size_t>(location.get_type());

    bool has_id_schemes =
        loc_id < m_testing_schemes_at_location_id.size() && !m_testing_schemes_at_location_id[loc_id].schemes.empty();

    bool has_type_schemes = loc_type < m_testing_schemes_at_location_type.size() &&
                            !m_testing_schemes_at_location_type[loc_type].schemes.empty();

    if (!has_id_schemes && !has_type_schemes) {
        return true; // No applicable schemes
    }

    bool entry_allowed = true; // Assume entry is allowed unless a scheme denies it
    // Check schemes for specific location id
    if (has_id_schemes) {
        for (const auto& scheme : m_testing_schemes_at_location_id[loc_id].schemes) {
            if (scheme.run_and_test(rng, person, t)) {
                entry_allowed = false; // Deny entry
            }
        }
    }

    // Check schemes for location type
    if (has_type_schemes) {
        for (const auto& scheme : m_testing_schemes_at_location_type[loc_type].schemes) {
            if (scheme.run_and_test(rng, person, t)) {
                entry_allowed = false; // Deny entry
            }
        }
    }

    // If the location is a home, entry is always allowed regardless of testing, no early return here because we still need to test

    return entry_allowed;
}

} // namespace abm
} // namespace mio
