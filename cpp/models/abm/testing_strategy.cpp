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

bool TestingCriteria::evaluate(const Person& p, TimePoint t) const
{
    // An empty vector of ages or none bitset of #InfectionStates% means that no condition on the corresponding property is set.
    return (m_ages.none() || m_ages[static_cast<size_t>(p.get_age())]) &&
           (m_infection_states.none() || m_infection_states[static_cast<size_t>(p.get_infection_state(t))]);
}

TestingScheme::TestingScheme(const TestingCriteria& testing_criteria, TimeSpan minimal_time_since_last_test,
                             TimePoint start_date, TimePoint end_date, const GenericTest& test_type, double probability)
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

bool TestingScheme::is_active() const
{
    return m_is_active;
}

void TestingScheme::update_activity_status(TimePoint t)
{
    m_is_active = (m_start_date <= t && t <= m_end_date);
}

bool TestingScheme::run_scheme(PersonalRandomNumberGenerator& rng, Person& person, TimePoint t) const
{
    if (t - person.get_time_of_last_test() > m_minimal_time_since_last_test) {
        if (m_testing_criteria.evaluate(person, t)) {
            double random = UniformDistribution<double>::get_instance()(rng);
            if (random < m_probability) {
                return !person.get_tested(rng, t, m_test_type.get_default());
            }
        }
    }
    return true;
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

void TestingStrategy::remove_testing_scheme(const LocationType& loc_type, const LocationId& loc_id,
                                            const TestingScheme& scheme)
{
    auto iter_schemes =
        std::find_if(m_location_to_schemes_map.begin(), m_location_to_schemes_map.end(), [&](const auto& p) {
            return p.type == loc_type && p.id == loc_id;
        });
    if (iter_schemes != m_location_to_schemes_map.end()) {
        //remove the scheme from the list
        auto& schemes_vector = iter_schemes->schemes;
        auto last            = std::remove(schemes_vector.begin(), schemes_vector.end(), scheme);
        schemes_vector.erase(last, schemes_vector.end());
        //delete the list of schemes for this location if no schemes left
        if (schemes_vector.empty()) {
            m_location_to_schemes_map.erase(iter_schemes);
        }
    }
}

void TestingStrategy::update_activity_status(TimePoint t)
{
    for (auto& [_type, _id, testing_schemes] : m_location_to_schemes_map) {
        for (auto& scheme : testing_schemes) {
            scheme.update_activity_status(t);
        }
    }
}

bool TestingStrategy::run_strategy(PersonalRandomNumberGenerator& rng, Person& person, const Location& location,
                                   TimePoint t)
{
    // A Person is always allowed to go home and this is never called if a person is not discharged from a hospital or ICU.
    if (location.get_type() == mio::abm::LocationType::Home) {
        return true;
    }

    //lookup schemes for this specific location as well as the location type
    //lookup in std::vector instead of std::map should be much faster unless for large numbers of schemes
    for (auto key : {std::make_pair(location.get_type(), location.get_id()),
                     std::make_pair(location.get_type(), LocationId::invalid_id())}) {
        auto iter_schemes =
            std::find_if(m_location_to_schemes_map.begin(), m_location_to_schemes_map.end(), [&](const auto& p) {
                return p.type == key.first && p.id == key.second;
            });
        if (iter_schemes != m_location_to_schemes_map.end()) {
            //apply all testing schemes that are found
            auto& schemes = iter_schemes->schemes;
            if (!std::all_of(schemes.begin(), schemes.end(), [&rng, &person, t](TestingScheme& ts) {
                    return !ts.is_active() || ts.run_scheme(rng, person, t);
                })) {
                return false;
            }
        }
    }
    return true;
}

} // namespace abm
} // namespace mio
