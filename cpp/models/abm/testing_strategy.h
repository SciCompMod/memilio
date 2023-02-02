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
#ifndef EPI_ABM_TESTING_SCHEME_H
#define EPI_ABM_TESTING_SCHEME_H

#include "abm/parameters.h"
#include "abm/person.h"
#include "abm/location.h"

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
     * Create a testing criteria.
     * @param ages Vector of AgeGroup%s that are either allowed or required to be tested.
     * @param[in] location_types Vector of #LocationType%s that are either allowed or required to be tested.
     * @param[in] infection_states Vector of #InfectionState%s that are either allowed or required to be tested.
     * An empty vector of ages/#LocationType%s/#InfectionStates% means that no condition on the corresponding property
     * is set!
     */
    TestingCriteria() = default;
    TestingCriteria(const std::vector<AgeGroup>& ages, const std::vector<LocationType>& location_types,
                    const std::vector<InfectionState>& infection_states);

    /**
     * @brief Compares two TestingCriteria for functional equality.
     */
    bool operator==(TestingCriteria other) const;

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
     * @brief Add a #LocationType to the set of #LocationType%s that are either allowed or required to be tested.
     * @param[in] location_type #LocationType to be added.
     */
    void add_location_type(const LocationType location_type);
    /**
     * @brief Remove a #LocationType from the set of #LocationType%s that are either allowed or required to be tested.
     * @param[in] location_type #LocationType to be removed.
     */
    void remove_location_type(const LocationType location_type);
    /**
     * @brief Add an #InfectionState to the set of #InfectionState%s that are either allowed or required to be tested.
     * @param[in] infection_state #InfectionState to be added.
     */
    void add_infection_state(const InfectionState infection_state);
    /**
     * remove an infection state from the set of infection states that are either allowed or required to be tested
     * @param infection_state infection state to be removed
     */
    void remove_infection_state(const InfectionState infection_state);

    /**
     * check if a person and a location meet all the required properties to get tested
     * @param p person to be checked
     * @param l location to be checked
     */
    bool evaluate(const Person& p, const Location& l) const;

private:
    /**
     * check if a person has the required age to get tested
     * @param p person to be checked
     */
    bool has_requested_age(const Person& p) const;

    /**
     * check if a location is in the set of locations that are allowed for testing
     * @param l location to be checked
     */
    bool is_requested_location_type(const Location& l) const;

    /**
     * check if a person has the required infection state to get tested
     * @param p person to be checked
     */
    bool has_requested_infection_state(const Person& p) const;

    std::vector<AgeGroup> m_ages; ///< Set of ages that are either allowed or required to be tested.
    std::vector<LocationType> m_location_types; /**< Set of #LocationState%s that are either allowed or required to be 
    tested.*/
    std::vector<InfectionState> m_infection_states; /**< Set of #InfectionState%s that are either allowed or required to
    be tested.*/
};

/**
 * @brief TestingScheme to regular test Person%s.
 */
class TestingScheme
{
public:
    /**
     * Create a testing scheme.
     * @param testing_criteria vector of testing criteria that are checked for testing
     * @param minimal_time_since_last_test time length of how often this scheme applies, i. e., a new test is performed after a person's last test
     * @param start_date starting date of the scheme
     * @param end_date ending date of the scheme
     * @param probability probability of the test to be performed if a testing rule applies
     * @param test_type the type of test to be performed
     */
    TestingScheme(const std::vector<TestingCriteria>& testing_criteria, TimeSpan minimal_time_since_last_test,
                  TimePoint start_date, TimePoint end_date, const GenericTest& test_type, double probability);

    /**
     * @brief Compares two TestingScheme%s for functional equality.
     */
    bool operator==(const TestingScheme& other) const;

    /**
     * @brief Add a TestingCriteria to the set of age groups that are checked for testing.
     * @param[in] criteria TestingCriteria to be added.
     */
    void add_testing_criteria(const TestingCriteria criteria);

    /**
     * @brief Remove a TestingCriteria from the set of age groups that are checked for testing.
     * @param[in] criteria TestingCriteria to be removed.
     */
    void remove_testing_criteria(const TestingCriteria criteria);

    /**
     * @brief Get the activity status of the scheme.
     * @return Activity status of the scheme.
     */
    bool is_active() const;

    /**
     * checks if the scheme is active at a given time and updates activity status
     * @param t time to be updated at
     */
    void update_activity_status(const TimePoint t);

    /**
     * runs the testing scheme and tests a person if necessary
     * @return if the person is allowed to enter the location
     */
    bool run_scheme(Person& person, const Location& location) const;

private:
    std::vector<TestingCriteria> m_testing_criteria;
    TimeSpan m_minimal_time_since_last_test;
    TimePoint m_start_date;
    TimePoint m_end_date;
    GenericTest m_test_type;
    double m_probability;
    bool m_is_active = false;
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
    explicit TestingStrategy(const std::vector<TestingScheme>& testing_schemes);

    /**
     * @brief Add a TestingScheme to the set of schemes that are checked for testing.
     * @param[in] scheme TestingScheme to be added.
     */
    void add_testing_scheme(const TestingScheme& scheme);

    /**
     * @brief Remove a TestingScheme from the set of schemes that are checked for testing.
     * @param[in] scheme TestingScheme to be removed.
     */
    void remove_testing_scheme(const TestingScheme& scheme);

    /**
     * checks if the given time point t is within the interval of start and end date of each testing scheme and then changes the activity status for each testing scheme accordingly
     * @param t time point to check the activity status of each testing scheme
     */
    void update_activity_status(const TimePoint t);

    /**
     * run the testing strategy and tests a person if necessary
     * @return if the person is allowed to enter the location
     */
    bool run_strategy(Person& person, const Location& location) const;

private:
    std::vector<TestingScheme> m_testing_schemes; ///< Set of schemes that are checked for testing.
};

} // namespace abm
} // namespace mio

#endif
