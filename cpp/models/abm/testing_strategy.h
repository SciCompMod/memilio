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
#ifndef EPI_ABM_TESTING_SCHEME_H
#define EPI_ABM_TESTING_SCHEME_H

#include "abm/parameters.h"
#include "abm/person.h"
#include "abm/location.h"
#include "abm/time.h"
#include "memilio/utils/random_number_generator.h"
#include <bitset>
#include <unordered_set>

namespace mio
{
namespace abm
{

/**
 * @brief TestingCriteria for TestingScheme.
 */
class TestingCriteria
{
public:
    /**
     * @brief Create a TestingCriteria where everyone is tested.
     */
    TestingCriteria() = default;

    /**
     * @brief Create a TestingCriteria.
     * @param[in] ages Vector of AgeGroup%s that are either allowed or required to be tested.
     * @param[in] infection_states Vector of #InfectionState%s that are either allowed or required to be tested.
     * An empty vector of ages or none bitset of #InfectionStates% means that no condition on the corresponding property
     * is set!
     */
    TestingCriteria(const std::vector<AgeGroup>& ages, const std::vector<InfectionState>& infection_states);

    /**
     * @brief Compares two TestingCriteria for functional equality.
     */
    bool operator==(const TestingCriteria& other) const;

    /**
     * @brief Add an AgeGroup to the set of AgeGroup%s that are either allowed or required to be tested.
     * @param[in] age_group AgeGroup to be added.
     */
    void add_age_group(const AgeGroup age_group);

    /**
     * @brief Remove an AgeGroup from the set of AgeGroup%s that are either allowed or required to be tested.
     * @param[in] age_group AgeGroup to be removed.
     */
    void remove_age_group(const AgeGroup age_group);

    /**
     * @brief Add an #InfectionState to the set of #InfectionState%s that are either allowed or required to be tested.
     * @param[in] infection_state #InfectionState to be added.
     */
    void add_infection_state(const InfectionState infection_state);

    /**
     * @brief Remove an #InfectionState from the set of #InfectionState%s that are either allowed or required to be
     * tested.
     * @param[in] infection_state #InfectionState to be removed.
     */
    void remove_infection_state(const InfectionState infection_state);

    /**
     * @brief Check if a Person and a Location meet all the required properties to get tested.
     * @param[in] p Person to be checked.
     * @param[in] t TimePoint when to evaluate the TestingCriteria.
     */
    bool evaluate(const Person& p, TimePoint t) const;

private:
    std::unordered_set<size_t> m_ages; ///< Set of #AgeGroup%s that are either allowed or required to be tested.
    std::bitset<(size_t)InfectionState::Count>
        m_infection_states; /**< BitSet of #InfectionState%s that are either allowed or required to
    be tested.*/
};

/**
 * @brief TestingScheme to regular test Person%s.
 */
class TestingScheme
{
public:
    /**
     * @brief Create a TestingScheme.
     * @param[in] testing_criteria Vector of TestingCriteria that are checked for testing.
     * @param[in] minimal_time_since_last_test TimeSpan of how often this scheme applies, i. e., when a new test is
     * performed after a Person's last test.
     * @param start_date Starting date of the scheme.
     * @param end_date Ending date of the scheme.
     * @param test_type The type of test to be performed.
     * @param probability Probability of the test to be performed if a testing rule applies.
     */
    TestingScheme(const TestingCriteria& testing_criteria, TimeSpan minimal_time_since_last_test, TimePoint start_date,
                  TimePoint end_date, const GenericTest& test_type, ScalarType probability);

    /**
     * @brief Compares two TestingScheme%s for functional equality.
     */
    bool operator==(const TestingScheme& other) const;

    /**
     * @brief Gets the activity status of the scheme.
     * @return Whether the TestingScheme is currently active.
     */
    bool is_active() const;

    /**
     * @brief Gets the activity status of the scheme at a given TimePoint.
     * @return Whether the TestingScheme is active at a given TimePoint.
     */
    bool is_active_at_time(TimePoint t) const;

    /**
     * @brief Checks if the scheme is active at a given time and updates activity status.
     * @param[in] t TimePoint to be updated at.
     */
    void update_activity_status(TimePoint t);

    /**
     * @brief Runs the TestingScheme and potentially tests a Person.
     * @param[inout] rng Person::RandomNumberGenerator for the Person being tested.
     * @param[in] person Person to check.
     * @param[in] t TimePoint when to run the scheme.
     * @return If the person is allowed to enter the Location by the scheme.
     */
    bool run_scheme(Person::RandomNumberGenerator& rng, Person& person, TimePoint t) const;

    /**
     * @brief Checks if the TestScheme is applicable for the given Person, trip TimePoint and current TimePoint
     * @param[in] person Person to check.
     * @param[in] trip_time TimePoint when the trip is schedules. 
     * @param[in] curr_time TimePoint of the current time.
     * @return If the TestScheme is applicable for the given Person, trip TimePoint and current TimePoint
     */
    bool is_applicable(const Person& person, TimePoint trip_time, TimePoint curr_time) const;

    /**
     * @brief Gets the type of the TestingScheme.
     * @return The type of the TestingScheme.
     */
    GenericTest get_type() const {
        return m_test_type;
    }

private:
    TestingCriteria m_testing_criteria; ///< TestingCriteria of the scheme.
    TimeSpan m_minimal_time_since_last_test; ///< Shortest period of time between two tests.
    TimePoint m_start_date; ///< Starting date of the scheme.
    TimePoint m_end_date; ///< Ending date of the scheme.
    GenericTest m_test_type; ///< Type of the test.
    ScalarType m_probability; ///< Probability of performing the test.
    bool m_is_active = false; ///< Whether the scheme is currently active.
};

/**
 * @brief Set of TestingSchemes that are checked for testing.
 */
class TestingStrategy
{
public:
    /**
     * @brief Create a TestingStrategy.
     * @param[in] testing_schemes Vector of TestingSchemes that are checked for testing.
     */
    TestingStrategy() = default;
    explicit TestingStrategy(const std::unordered_map<LocationId, std::vector<TestingScheme>>& location_to_schemes_map);

    /**
     * @brief Add a TestingScheme to the set of schemes that are checked for testing at a certain Location.
     * @param[in] loc_id LocationId key for TestingScheme to be added.
     * @param[in] scheme TestingScheme to be added.
     */
    void add_testing_scheme(const LocationId& loc_id, const TestingScheme& scheme);

    /**
     * @brief Add a TestingScheme to the set of schemes that are checked for testing at a certain LocationType.
     * A TestingScheme applies to all Location of the same type is store in 
     * LocationId{INVALID_LOCATION_INDEX, location_type} of m_location_to_schemes_map.
     * @param[in] loc_type LocationId key for TestingScheme to be added.
     * @param[in] scheme TestingScheme to be added.
     */
    void add_testing_scheme(const LocationType& loc_type, const TestingScheme& scheme)
    {
        add_testing_scheme(LocationId{INVALID_LOCATION_INDEX, loc_type}, scheme);
    }

    /**
     * @brief Remove a TestingScheme from the set of schemes that are checked for testing at a certain Location.
     * @param[in] loc_id LocationId key for TestingScheme to be remove.
     * @param[in] scheme TestingScheme to be removed.
     */
    void remove_testing_scheme(const LocationId& loc_id, const TestingScheme& scheme);

    /**
     * @brief Remove a TestingScheme from the set of schemes that are checked for testing at a certain Location.
     * A TestingScheme applies to all Location of the same type is store in 
     * LocationId{INVALID_LOCATION_INDEX, location_type} of m_location_to_schemes_map.
     * @param[in] loc_type LocationType key for TestingScheme to be remove.
     * @param[in] scheme TestingScheme to be removed.
     */
    void remove_testing_scheme(const LocationType& loc_type, const TestingScheme& scheme)
    {
        remove_testing_scheme(LocationId{INVALID_LOCATION_INDEX, loc_type}, scheme);
    }

    /**
     * @brief Checks if the given TimePoint is within the interval of start and end date of each TestingScheme and then
     * changes the activity status for each TestingScheme accordingly.
     * @param t TimePoint to check the activity status of each TestingScheme.
     */
    void update_activity_status(const TimePoint t);

    /**
     * @brief Runs the TestingStrategy and potentially tests a Person.
     * @param[inout] rng Person::RandomNumberGenerator for the Person being tested.
     * @param[in] person Person to check.
     * @param[in] location Location to check.
     * @param[in] t TimePoint when to run the strategy.
     * @return If the Person is allowed to enter the Location.
     */
    bool run_strategy(Person::RandomNumberGenerator& rng, Person& person, const Location& location, TimePoint t);

private:
    std::unordered_map<LocationId, std::vector<TestingScheme>>
        m_location_to_schemes_map; ///< Set of schemes that are checked for testing.
};

} // namespace abm
} // namespace mio

#endif
