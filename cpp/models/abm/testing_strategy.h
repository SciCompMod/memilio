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

#include "abm/config.h"
#include "abm/person.h"
#include "abm/time.h"
#include "memilio/utils/random_number_generator.h"
#include <bitset>
#include <unordered_set>
#include "memilio/utils/pointer_dereferencing_iterator.h"

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
     * @brief Add an #InfectionState to the set of #InfectionState%s that are either allowed or required to be tested.
     * @param[in] infection_state #InfectionState to be added.
     */
    void add_infection_state(const InfectionState infection_state);

    /**
     * @brief Check if a Person and a Location meet all the required properties to get tested.
     * @param[in] p Person to be checked.
     * @param[in] t TimePoint when to evaluate the TestingCriteria.
     */
    bool evaluate(const Person& p, TimePoint t) const;

private:
    std::bitset<MAX_NUM_AGE_GROUPS> m_ages; ///< Set of #AgeGroup%s that are either allowed or required to be tested.
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
     * @param test_parameters The parameters of test to be performed.
     * @param probability Probability of the test to be performed if a testing rule applies.
     */
    TestingScheme(const TestingCriteria& testing_criteria, TimeSpan minimal_time_since_last_test, TimePoint start_date,
                  TimePoint end_date, TestParameters test_parameters, ScalarType probability);

    /**
     * @brief Compares two TestingScheme%s for functional equality.
     */
    bool operator==(const TestingScheme& other) const;

    /**
     * @brief Get the activity status of the scheme.
     * @return Whether the TestingScheme is currently active.
     */
    bool is_active(const TimePoint t) const;

    /**
     * @brief Runs the TestingScheme and potentially tests a Person.
     * @param[inout] rng Person::RandomNumberGenerator for the Person being tested.
     * @param[in] person Person to check.
     * @param[in] t TimePoint when to run the scheme.
     * @return If the person is allowed to enter the Location by the scheme.
     */
    bool run_scheme(Person::RandomNumberGenerator& rng, Person& person, TimePoint t) const;

    mio::abm::TimePoint get_start_date() const
    {
        return m_start_date;
    }
    mio::abm::TimePoint get_end_date() const
    {
        return m_end_date;
    }

private:
    TestingCriteria m_testing_criteria; ///< TestingCriteria of the scheme.
    TimeSpan m_minimal_time_since_last_test; ///< Shortest period of time between two tests.
    TimePoint m_start_date; ///< Starting date of the scheme.
    TimePoint m_end_date; ///< Ending date of the scheme.
    TestParameters m_test_parameters; ///< Parameters of the test.
    ScalarType m_probability; ///< Probability of performing the test.
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
    explicit TestingStrategy(
        const std::unordered_map<LocationId, std::vector<mio::abm::TestingScheme>>& location_to_schemes_map);

    /**
     * @brief Add a TestingScheme to the set of schemes that are checked for testing at a certain Location.
     * @param[in] loc_id LocationId key for TestingScheme to be added.
     * @param[in] scheme TestingScheme to be added.
     */
    void add_testing_scheme(const LocationId& loc_id, const mio::abm::TestingScheme& scheme);

    /**
     * @brief Add a TestingScheme to the set of schemes that are checked for testing at a certain LocationType.
     * A TestingScheme applies to all Location of the same type is store in 
     * LocationId{INVALID_LOCATION_INDEX, location_type} of m_location_to_schemes_map.
     * @param[in] loc_type LocationId key for TestingScheme to be added.
     * @param[in] scheme TestingScheme to be added.
     */
    void add_testing_scheme(const LocationType& loc_type, const mio::abm::TestingScheme& scheme)
    {
        add_testing_scheme(LocationId{INVALID_LOCATION_INDEX, loc_type}, scheme);
    }

    /**
     * @brief Get the TestingSchemes that are checked for testing at a certain Location.
     * @return Vector of TestingScheme%s that are checked for testing at the Location.
     */
    const std::vector<std::pair<LocationId, std::vector<mio::abm::TestingScheme>>>& get_testing_schemes() const
    {
        return m_location_to_schemes_map;
    }
    using ConstLocationIterator = PointerDereferencingIterator<std::vector<std::unique_ptr<Location>>::const_iterator>;
    void update_location_testing_schemes(TimePoint t,
                                         mio::Range<std::pair<ConstLocationIterator, ConstLocationIterator>> locations);

    bool entry_allowed_testing_schemes(Person::RandomNumberGenerator& rng, Person& person, LocationId id,
                                       const mio::abm::TimePoint t);

private:
    std::vector<std::pair<LocationId, std::vector<mio::abm::TestingScheme>>>
        m_location_to_schemes_map; ///< Set of schemes that are checked for testing.
    std::vector<std::vector<mio::abm::TestingScheme>>
        m_testing_schemes_per_location; ///< List of TestingScheme%s that are checked for testing.
    std::vector<mio::abm::TimePoint> m_update_ts_scheduler; ///< List of TimePoint%s when to update the TestingScheme.
};

} // namespace abm
} // namespace mio

#endif
