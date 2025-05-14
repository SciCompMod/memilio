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
#ifndef MIO_ABM_TESTING_SCHEME_H
#define MIO_ABM_TESTING_SCHEME_H

#include "abm/config.h"
#include "abm/location_id.h"
#include "abm/location_type.h"
#include "abm/parameters.h"
#include "abm/person.h"
#include "abm/location.h"
#include "abm/time.h"
#include "memilio/io/default_serialize.h"

#include <bitset>
#include <vector>

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
     * @brief Check if a Person and a Location meet all the required properties to get tested.
     * @param[in] p Person to be checked.
     * @param[in] t TimePoint when to evaluate the TestingCriteria.
     */
    bool evaluate(const Person& p, TimePoint t) const;

    auto default_serialize()
    {
        return Members("TestingCriteria").add("ages", m_ages).add("infection_states", m_infection_states);
    }

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
     * @param validity_period The valid TimeSpan of the test. 
     * @param start_date Starting date of the scheme.
     * @param end_date Ending date of the scheme.
     * @param test_parameters The parameters of test to be performed.
     * @param probability Probability of the test to be performed if a testing rule applies.
     */
    TestingScheme(const TestingCriteria& testing_criteria, TimeSpan validity_period, TimePoint start_date,
                  TimePoint end_date, TestParameters test_parameters, ScalarType probability);

    /**
     * @brief Compares two TestingScheme%s for functional equality.
     */
    bool operator==(const TestingScheme& other) const;

    /**
     * @brief Gets the activity status of the scheme.
     * @return Whether the TestingScheme is currently active.
     */
    bool is_active(TimePoint t) const;

    /**
     * @brief Gets the start date of the scheme.
     * @return The start date of the scheme.
     */
    mio::abm::TimePoint get_start_date() const
    {
        return m_start_date;
    }

    /**
     * @brief Gets the end date of the scheme.
     * @return The end date of the scheme.
     */
    mio::abm::TimePoint get_end_date() const
    {
        return m_end_date;
    }


    /**
     * @brief Runs the TestingScheme and potentially tests a Person.
     * @param[inout] rng PersonalRandomNumberGenerator of the Person being tested.
     * @param[in] person Person to check.
     * @param[in] t TimePoint when to run the scheme.
     */
    void run_scheme(PersonalRandomNumberGenerator& rng, Person& person, TimePoint t) const;

    /// This method is used by the default serialization feature.
    auto default_serialize()
    {
        return Members("TestingScheme")
            .add("criteria", m_testing_criteria)
            .add("validity_period", m_validity_period)
            .add("start_date", m_start_date)
            .add("end_date", m_end_date)
            .add("test_params", m_test_parameters)
            .add("probability", m_probability);
    }

private:
    friend DefaultFactory<TestingScheme>;
    TestingScheme() = default;

    TestingCriteria m_testing_criteria; ///< TestingCriteria of the scheme.
    TimeSpan m_validity_period; ///< The valid TimeSpan of the test.
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
     * @brief List of testing schemes for a given LocationType and LocationId.
     * A LocalStrategy with id of value LocationId::invalid_id() is used for all Locations with LocationType type.
     */
    struct LocalStrategy {
        LocationType type;
        LocationId id;
        std::vector<TestingScheme> schemes;

        /// This method is used by the default serialization feature.
        auto default_serialize()
        {
            return Members("LocalStrategy").add("type", type).add("id", id).add("schemes", schemes);
        }
    };

    /**
     * @brief Create a TestingStrategy.
     * @param[in] testing_schemes Vector of TestingSchemes that are checked for testing.
     */
    TestingStrategy() = default;
    explicit TestingStrategy(const std::vector<LocalStrategy>& location_to_schemes_map);

    /**
     * @brief Add a TestingScheme to the set of schemes that are checked for testing at a certain Location.
     * A TestingScheme with loc_id of value LocationId::invalid_id() is used for all Locations of the given type.
     * @param[in] loc_type LocationType key for TestingScheme to be remove.
     * @param[in] loc_id LocationId key for TestingScheme to be added.
     * @param[in] scheme TestingScheme to be added.
     */
    void add_testing_scheme(const LocationType& loc_type, const LocationId& loc_id, const TestingScheme& scheme);

    /**
     * @brief Add a TestingScheme to the set of schemes that are checked for testing at a certain LocationType.
     * A TestingScheme with loc_id of value LocationId::invalid_id() is used for all Locations of the given type.
     * @param[in] loc_type LocationId key for TestingScheme to be added.
     * @param[in] scheme TestingScheme to be added.
     */
    void add_testing_scheme(const LocationType& loc_type, const TestingScheme& scheme)
    {
        add_testing_scheme(loc_type, LocationId::invalid_id(), scheme);
    }

    /**
     * @brief Runs the TestingStrategy and potentially tests a Person.
     * @param[inout] rng PersonalRandomNumberGenerator of the Person being tested.
     * @param[in] person Person to check.
     * @param[in] location Location to check.
     * @param[in] t TimePoint when to run the strategy.
     */
    void run_strategy(PersonalRandomNumberGenerator& rng, Person& person, const Location& location, TimePoint t);

    /// This method is used by the default serialization feature.
    auto default_serialize()
    {
        return Members("TestingStrategy").add("schemes", m_location_to_schemes_map);
    }

private:
    std::vector<LocalStrategy> m_location_to_schemes_map; ///< Set of schemes that are checked for testing.
    std::vector<mio::abm::TimePoint> m_testing_scheme_scheduler;
};

} // namespace abm
} // namespace mio

#endif
