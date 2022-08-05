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
#include "abm/testing_rule.h"

namespace mio
{
namespace abm
{

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
        return this == &other;
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
