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
                             TimePoint end_date, TestParameters test_parameters, double probability)
    : m_testing_criteria(testing_criteria)
    , m_validity_period(validity_period)
    , m_start_date(start_date)
    , m_end_date(end_date)
    , m_test_parameters(test_parameters)
    , m_probability(probability)
{
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

void TestingScheme::run_scheme(PersonalRandomNumberGenerator& rng, Person& person, TimePoint t) const
{
    auto test_result = person.get_test_result(m_test_parameters.type);
    // If the agent has a test result valid until now, use the result directly
    if ((test_result.time_of_testing > TimePoint(std::numeric_limits<int>::min())) &&
        (test_result.time_of_testing + m_validity_period >= t)) {
            return
    }
    // Otherwise, the time_of_testing in the past (i.e. the agent has already performed it).
    if (m_testing_criteria.evaluate(person, t - m_test_parameters.required_time)) {
        double random = UniformDistribution<double>::get_instance()(rng);
        if (random < m_probability) {
            bool result = person.get_tested(rng, t - m_test_parameters.required_time, m_test_parameters);
            person.add_test_result(t, m_test_parameters.type, result);
        }
    }
}

TestingStrategy::TestingStrategy(const std::vector<LocalStrategy>& location_to_schemes_map)
    : m_location_to_schemes_map(location_to_schemes_map.begin(), location_to_schemes_map.end())
{
}

void TestingStrategy::add_testing_scheme(const LocationType& loc_type, const LocationId& loc_id,
                                         const TestingScheme& scheme)
{
    auto iter_schemes =
        std::find_if(m_location_to_schemes_map.begin(), m_location_to_schemes_map.end(), [&](const auto& p) {
            return p.type == loc_type && p.id == loc_id;
        });
    if (iter_schemes == m_location_to_schemes_map.end()) {
        //no schemes for this location yet, add a new list with one scheme
        m_location_to_schemes_map.push_back({loc_type, loc_id, std::vector<TestingScheme>(1, scheme)});
    }
    else {
        //add scheme to existing vector if the scheme doesn't exist yet
        auto& schemes = iter_schemes->schemes;
        if (std::find(schemes.begin(), schemes.end(), scheme) == schemes.end()) {
            schemes.push_back(scheme);
        }
    }
}

void TestingStrategy::run_strategy(PersonalRandomNumberGenerator& rng, Person& person, const Location& location,
                                   TimePoint t)
{
    if (m_testing_schemes_per_location.at(location.get_id().get()).empty()) {
        return;
    }

    for (auto&& scheme : m_testing_schemes_per_location.at(location.get_id().get())) {
        scheme.run_scheme(rng, person, t)
    }
}

} // namespace abm
} // namespace mio
