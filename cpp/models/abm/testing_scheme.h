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
 * Testing Rule for Testing Scheme
 */
class TestingRule
{
public:
    /**
     * Create a testing rule.
     * @param ages vector of age groups that are allowed for testing
     * @param location_types vector of location types that are allowed for testing
     * @param infection_states vector of infection states that are allowed for testing
     * An empty vector of ages/location types/infection states means that no condition on the corresponding property is set!
     */
    TestingRule(const std::vector<AgeGroup> ages, const std::vector<LocationType> location_types,
                const std::vector<InfectionState> infection_states);

    /**
     * Compares two testing rules for equality.
     * Compare references. Still possible to clone rules.
     */
    bool operator==(const TestingRule& other) const
    {
        return ((this->m_ages == other.m_ages) && (this->m_infection_states == other.m_infection_states) &&
                (this->m_location_types == other.m_location_types));
    }

    /**
     * add an age group to the set of age groups that are allowed for testing
     * @param age_group age group to be added
     */
    void add_age_group(const AgeGroup age_group);

    /**
     * remove an age group from the set of age groups that are allowed for testing
     * @param age_group age group to be removed
     */
    void remove_age_group(const AgeGroup age_group);

    /**
     * add a location type to the set of location types that are allowed for testing
     * @param location_type location type to be added
     */
    void add_location_type(const LocationType location_type);

    /**
     * remove a location tpye from the set of location tpyes that are allowed for testing
     * @param location_type location type to be removed
     */
    void remove_location_type(const LocationType location_type);

    /**
     * add an infection state to the set of infection states that are allowed for testing
     * @param infection_state infection state to be added
     */
    void add_infection_state(const InfectionState infection_state);

    /**
     * remove an infection state from the set of infection states that are allowed for testing
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
    std::vector<AgeGroup> m_ages;
    std::vector<LocationType> m_location_types;
    std::vector<InfectionState> m_infection_states;
};

/**
 * Testing Scheme to regular test people
 */
class TestingScheme
{
public:
    /**
     * Create a testing scheme.
     * @param testing_rules vector of testing rules that are checked for testing
     * @param testing_frequency time length of how often this scheme applies, i. e., a new test is performed after a person's last test
     * @param start_date starting date of the scheme
     * @param end_date ending date of the scheme
     * @param probability probability of the test to be performed if a testing rule applies
     * @param test_type the type of test to be performed
     */
    TestingScheme(const std::vector<TestingRule> testing_rules, const TimeSpan testing_frequency, TimePoint start_date,
                  TimePoint end_date, const double probability, const GenericTest& test_type);

    /**
     * Compares two testing schemes for equality.
     * Compare references. Still possible to clone schemes.
     */
    bool operator==(const TestingScheme& other) const
    {
        return this->m_testing_rules == other.m_testing_rules &&
               this->m_testing_frequency == other.m_testing_frequency && this->m_start_date == other.m_start_date &&
               this->m_end_date == other.m_end_date && this->m_probability == other.m_probability &&
               this->m_test_type.get_default().sensitivity == other.m_test_type.get_default().sensitivity &&
               this->m_test_type.get_default().specificity == other.m_test_type.get_default().specificity;
        //To be adjusted and also TestType should be static.
    }

    /**
     * add a testing rule to the set of age groups that are checked for testing
     * @param rule testing rule to be added
     */
    void add_testing_rule(const TestingRule rule);

    /**
     * remove a testing rule from the set of age groups that are checked for testing
     * @param rule testing rule to be removed
     */
    void remove_testing_rule(const TestingRule rule);

    /**
     * @return activity status of the scheme
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
    std::vector<TestingRule> m_testing_rules;
    TimeSpan m_testing_frequency;
    double m_probability;
    TimePoint m_start_date;
    TimePoint m_end_date;
    GenericTest m_test_type;
    bool m_is_active = false;
};

} // namespace abm
} // namespace mio

#endif
