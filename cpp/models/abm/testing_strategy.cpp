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
#include "abm/location.h"

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

void TestingCriteria::add_infection_state(const InfectionState infection_state)
{
    m_infection_states.set(static_cast<size_t>(infection_state), true);
}

bool TestingCriteria::evaluate(const Person& p, TimePoint t) const
{
    // An empty vector of ages or none bitset of #InfectionStates% means that no condition on the corresponding property is set.
    return (m_ages.none() || m_ages[static_cast<size_t>(p.get_age())]) &&
           (m_infection_states.none() || m_infection_states[static_cast<size_t>(p.get_infection_state(t))]);
}

TestingScheme::TestingScheme(const TestingCriteria& testing_criteria, TimeSpan minimal_time_since_last_test,
                             TimePoint start_date, TimePoint end_date, TestParameters test_parameters,
                             double probability)
    : m_testing_criteria(testing_criteria)
    , m_minimal_time_since_last_test(minimal_time_since_last_test)
    , m_start_date(start_date)
    , m_end_date(end_date)
    , m_test_parameters(test_parameters)
    , m_probability(probability)
{
}

bool TestingScheme::operator==(const TestingScheme& other) const
{
    return this->m_testing_criteria == other.m_testing_criteria &&
           this->m_minimal_time_since_last_test == other.m_minimal_time_since_last_test &&
           this->m_start_date == other.m_start_date && this->m_end_date == other.m_end_date &&
           this->m_test_parameters.sensitivity == other.m_test_parameters.sensitivity &&
           this->m_test_parameters.specificity == other.m_test_parameters.specificity &&
           this->m_probability == other.m_probability;
    //To be adjusted and also TestType should be static.
}

bool TestingScheme::is_active(const TimePoint t) const
{
    return (m_start_date <= t && t <= m_end_date);
}

bool TestingScheme::run_scheme(Person::RandomNumberGenerator& rng, Person& person, TimePoint t) const
{
    if (t - person.get_time_of_last_test() > m_minimal_time_since_last_test) {
        if (m_testing_criteria.evaluate(person, t)) {
            double random = UniformDistribution<double>::get_instance()(rng);
            if (random < m_probability) {
                return person.get_tested(rng, t, m_test_parameters);
            }
        }
    }
    return false;
}

TestingStrategy::TestingStrategy(
    const std::unordered_map<LocationId, std::vector<TestingScheme>>& location_to_schemes_map)
    : m_location_to_schemes_map(location_to_schemes_map.begin(), location_to_schemes_map.end())
{
}

void TestingStrategy::add_testing_scheme(const LocationId& loc_id, const TestingScheme& scheme)
{
    auto iter_schemes =
        std::find_if(m_location_to_schemes_map.begin(), m_location_to_schemes_map.end(), [loc_id](auto& p) {
            return p.first == loc_id;
        });
    if (iter_schemes == m_location_to_schemes_map.end()) {
        //no schemes for this location yet, add a new list with one scheme
        m_location_to_schemes_map.emplace_back(loc_id, std::vector<TestingScheme>(1, scheme));
    }
    else {
        //add scheme to existing vector if the scheme doesn't exist yet
        auto& schemes = iter_schemes->second;
        if (std::find(schemes.begin(), schemes.end(), scheme) == schemes.end()) {
            schemes.push_back(scheme);
        }
    }
}

bool TestingStrategy::entry_allowed_testing_schemes(Person::RandomNumberGenerator& rng, Person& person, LocationId id,
                                                    const mio::abm::TimePoint t)
{
    if (m_testing_schemes_per_location.at(id.index).empty()) {
        return true;
    }

    for (auto&& scheme : m_testing_schemes_per_location.at(id.index)) {
        if (scheme.run_scheme(rng, person, t)) {
            return false;
        }
    }
    return true;
}

void TestingStrategy::update_location_testing_schemes(
    TimePoint t, mio::Range<std::pair<ConstLocationIterator, ConstLocationIterator>> locations)
{
    if (m_testing_schemes_per_location.empty()) {
        for (auto i = size_t(0); i < locations.size(); ++i) {
            m_testing_schemes_per_location.push_back({});
        }
        for (auto&& schemes : m_location_to_schemes_map) {
            for (auto&& scheme : schemes.second) {
                auto start_date = scheme.get_start_date();
                auto end_date   = scheme.get_end_date();
                if (std::find_if(m_update_ts_scheduler.begin(), m_update_ts_scheduler.end(), [start_date](auto& p) {
                        return p == start_date;
                    }) == m_update_ts_scheduler.end()) {
                    m_update_ts_scheduler.push_back(start_date);
                }
                if (std::find_if(m_update_ts_scheduler.begin(), m_update_ts_scheduler.end(), [end_date](auto& p) {
                        return p == end_date;
                    }) == m_update_ts_scheduler.end()) {
                    m_update_ts_scheduler.push_back(end_date);
                }
            }
        }
        std::sort(m_update_ts_scheduler.begin(), m_update_ts_scheduler.end());
    }

    if (std::find(m_update_ts_scheduler.begin(), m_update_ts_scheduler.end(), t) != m_update_ts_scheduler.end()) {
        std::cout << "Updating testing schemes at time " << t.days() << std::endl;
        for (auto i = size_t(0); i < locations.size(); ++i) {
            auto& location = locations[i];
            auto loc_id    = location.get_index();
            auto& schemes  = m_location_to_schemes_map;

            auto iter_schemes = std::find_if(schemes.begin(), schemes.end(), [loc_id, &location](auto& p) {
                return p.first.index == loc_id ||
                       (p.first.index == INVALID_LOCATION_INDEX && p.first.type == location.get_type());
            });
            if (iter_schemes != schemes.end()) {
                m_testing_schemes_per_location[i].clear();
                for (auto&& scheme : iter_schemes->second) {
                    if (scheme.is_active(t)) {
                        m_testing_schemes_per_location[i].push_back(scheme);
                    }
                }
            }
        }
    }
}

} // namespace abm
} // namespace mio
