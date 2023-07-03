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

#include "abm/parameters.h" // IWYU pragma: keep
#include "abm/person.h"
#include "abm/location.h"
#include "abm/time.h"
#include "memilio/utils/random_number_generator.h" // IWYU pragma: keep

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
     * @brief Create a TestingCriteria.
     * @param[in] ages Vector of AgeGroup%s that are either allowed or required to be tested.
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
     * @brief Remove an #InfectionState from the set of #InfectionState%s that are either allowed or required to be
     * tested.
     * @param[in] infection_state #InfectionState to be removed.
     */
    void remove_infection_state(const InfectionState infection_state);

    /**
     * @brief Check if a Person and a Location meet all the required properties to get tested.
     * @param[in] p Person to be checked.
     * @param[in] l Location to be checked.
     * @param[in] t TimePoint when to evaluate the TestingCriteria.
     */
    template<typename FP=double>
    bool evaluate(const Person<FP>& p, const Location<FP>& l, TimePoint t) const
    {
        return has_requested_age(p) && is_requested_location_type(l) && has_requested_infection_state(p, t);
    }

private:
    /**
     * @brief Check if a Person has the required age to get tested.
     * @param[in] p Person to be checked.
     */
    template<typename FP=double>
    bool has_requested_age(const Person<FP>& p) const
    {
        if (m_ages.empty()) {
            return true; // no condition on the age
        }
        return std::find(m_ages.begin(), m_ages.end(), p.get_age()) != m_ages.end();
    }

    /**
     * @brief Check if a Location is in the set of Location%s that are allowed for testing.
     * @param[in] l Location to be checked.
     */
    template<typename FP=double>
    bool is_requested_location_type(const Location<FP>& l) const
    {
        if (m_location_types.empty()) {
            return true; // no condition on the location
        }
        return std::find(m_location_types.begin(), m_location_types.end(), l.get_type()) != m_location_types.end();
    }

    /**
     * @brief Check if a Person has the required InfectionState to get tested.
     * @param[in] p Person to be checked.
     * @param[in] t TimePoint when to check.
     */
    template<typename FP=double>
    bool has_requested_infection_state(const Person<FP>& p, TimePoint t) const
    {
        if (m_infection_states.empty()) {
            return true; // no condition on infection state
        }
        return std::find(m_infection_states.begin(), m_infection_states.end(), p.get_infection_state(t)) !=
               m_infection_states.end();
    }

    std::vector<AgeGroup> m_ages; ///< Set of #AgeGroup%s that are either allowed or required to be tested.
    std::vector<LocationType> m_location_types; /**< Set of #LocationState%s that are either allowed or required to be 
    tested.*/
    std::vector<InfectionState> m_infection_states; /**< Set of #InfectionState%s that are either allowed or required to
    be tested.*/
};

/**
 * @brief TestingScheme to regular test Person%s.
 */
template<typename FP=double>
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
    TestingScheme(const std::vector<TestingCriteria>& testing_criteria, TimeSpan minimal_time_since_last_test,
                  TimePoint start_date, TimePoint end_date, const GenericTest<FP>& test_type, ScalarType probability)
        : m_testing_criteria(testing_criteria)
        , m_minimal_time_since_last_test(minimal_time_since_last_test)
        , m_start_date(start_date)
        , m_end_date(end_date)
        , m_test_type(test_type)
        , m_probability(probability)
    {
    }


    /**
     * @brief Compares two TestingScheme%s for functional equality.
     */
    bool operator==(const TestingScheme& other) const
    {
        return this->m_testing_criteria == other.m_testing_criteria &&
               this->m_minimal_time_since_last_test == other.m_minimal_time_since_last_test &&
               this->m_start_date == other.m_start_date && this->m_end_date == other.m_end_date &&
               this->m_test_type.get_default().sensitivity == other.m_test_type.get_default().sensitivity &&
               this->m_test_type.get_default().specificity == other.m_test_type.get_default().specificity &&
               this->m_probability == other.m_probability;
        //To be adjusted and also TestType should be static.
    }

    /**
     * @brief Add a TestingCriteria to the set of TestingCriteria that are checked for testing.
     * @param[in] criteria TestingCriteria to be added.
     */
    void add_testing_criteria(const TestingCriteria criteria)
    {
        if (std::find(m_testing_criteria.begin(), m_testing_criteria.end(), criteria) == m_testing_criteria.end()) {
            m_testing_criteria.push_back(criteria);
        }
    }

    /**
     * @brief Remove a TestingCriteria from the set of TestingCriteria that are checked for testing.
     * @param[in] criteria TestingCriteria to be removed.
     */
    void remove_testing_criteria(const TestingCriteria criteria)
    {
        auto last = std::remove(m_testing_criteria.begin(), m_testing_criteria.end(), criteria);
        m_testing_criteria.erase(last, m_testing_criteria.end());
    }

    /**
     * @brief Get the activity status of the scheme.
     * @return Whether the TestingScheme is currently active.
     */
    bool is_active() const
    {
        return m_is_active;
    }

    /**
     * @brief Checks if the scheme is active at a given time and updates activity status.
     * @param[in] t TimePoint to be updated at.
     */
    void update_activity_status(TimePoint t)
    {
        m_is_active = (m_start_date <= t && t <= m_end_date);
    }

    /**
     * @brief Runs the TestingScheme and potentially tests a Person.
     * @param[in] person Person to check.
     * @param[in] location Location to check.
     * @param[in] t TimePoint when to run the scheme.
     * @return If the person is allowed to enter the Location by the scheme.
     */
    bool run_scheme(Person<FP>& person, const Location<FP>& location, TimePoint t) const
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

private:
    std::vector<TestingCriteria> m_testing_criteria; ///< Vector with all TestingCriteria of the scheme.
    TimeSpan m_minimal_time_since_last_test; ///< Shortest period of time between two tests.
    TimePoint m_start_date; ///< Starting date of the scheme.
    TimePoint m_end_date; ///< Ending date of the scheme.
    GenericTest<FP> m_test_type; ///< Type of the test.
    ScalarType m_probability; ///< Probability of performing the test.
    bool m_is_active = false; ///< Whether the scheme is currently active.
};

/**
 * @brief Set of TestingSchemes that are checked for testing.
 */
template<typename FP=double>
class TestingStrategy
{
public:
    /**
     * @brief Create a TestingStrategy.
     * @param[in] testing_schemes Vector of TestingSchemes that are checked for testing.
     */
    TestingStrategy() = default;
    explicit TestingStrategy(const std::vector<TestingScheme<FP>>& testing_schemes)
        : m_testing_schemes(testing_schemes)
    {
    }

    /**
     * @brief Add a TestingScheme to the set of schemes that are checked for testing.
     * @param[in] scheme TestingScheme to be added.
     */
    void add_testing_scheme(const TestingScheme<FP>& scheme)
    {
        if (std::find(m_testing_schemes.begin(), m_testing_schemes.end(), scheme) == m_testing_schemes.end()) {
            m_testing_schemes.push_back(scheme);
        }
    }


    /**
     * @brief Remove a TestingScheme from the set of schemes that are checked for testing.
     * @param[in] scheme TestingScheme to be removed.
     */
    void remove_testing_scheme(const TestingScheme<FP>& scheme)
    {
        auto last = std::remove(m_testing_schemes.begin(), m_testing_schemes.end(), scheme);
        m_testing_schemes.erase(last, m_testing_schemes.end());
    }

    /**
     * @brief Checks if the given TimePoint is within the interval of start and end date of each TestingScheme and then
     * changes the activity status for each TestingScheme accordingly.
     * @param t TimePoint to check the activity status of each TestingScheme.
     */
    void update_activity_status(const TimePoint t)
    {
        for (auto& ts : m_testing_schemes) {
            ts.update_activity_status(t);
        }
    }


    /**
     * @brief Runs the TestingStrategy and potentially tests a Person.
     * @param[in] person Person to check.
     * @param[in] location Location to check.
     * @param[in] t TimePoint when to run the strategy.
     * @return If the Person is allowed to enter the Location.
     */
    bool run_strategy(Person<FP>& person, const Location<FP>& location, TimePoint t) const
    {
        // Person who is in quarantine but not yet home should go home. Otherwise they can't because they test positive.
        if (location.get_type() == mio::abm::LocationType::Home && person.is_in_quarantine()) {
            return true;
        }
        return std::all_of(m_testing_schemes.begin(), m_testing_schemes.end(), [&person, location, t](TestingScheme<FP> ts) {
            if (ts.is_active()) {
                return ts.run_scheme(person, location, t);
            }
            return true;
        });
    }

private:
    std::vector<TestingScheme<FP>> m_testing_schemes; ///< Set of schemes that are checked for testing.
};

} // namespace abm
} // namespace mio

#endif
