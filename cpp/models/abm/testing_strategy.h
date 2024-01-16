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
#include "abm/parameters.h"
#include "abm/person.h"
#include "abm/location.h"
#include "abm/time.h"
#include "memilio/utils/random_number_generator.h"
#include <bitset>
#include <unordered_set> // IWYU pragma: keep

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
    template<typename FP=double>
    bool evaluate(const Person<FP>& p, TimePoint t) const
    {
        // An empty vector of ages or none bitset of #InfectionStates% means that no condition on the corresponding property is set.
        return (m_ages.none() || m_ages[static_cast<size_t>(p.get_age())]) &&
               (m_infection_states.none() || m_infection_states[static_cast<size_t>(p.get_infection_state(t))]);
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
    TestingScheme(const TestingCriteria& testing_criteria, TimeSpan minimal_time_since_last_test,
                  TimePoint start_date, TimePoint end_date, const GenericTest<FP>& test_type, double probability)
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
     * @param[inout] rng Person::RandomNumberGenerator for the Person being tested.
     * @param[in] person Person to check.
     * @param[in] t TimePoint when to run the scheme.
     * @return If the person is allowed to enter the Location by the scheme.
     */
    bool run_scheme(typename Person<FP>::RandomNumberGenerator& rng, Person<FP>& person,
                    TimePoint t) const
    {
        if (person.get_time_since_negative_test() > m_minimal_time_since_last_test) {
            if (m_testing_criteria.evaluate(person, t)) {
                double random = UniformDistribution<double>::get_instance()(rng);
                if (random < m_probability) {
                    return !person.get_tested(rng, t, m_test_type.get_default());
                }
            }
        }
        return true;
    }

private:
    TestingCriteria m_testing_criteria; ///< TestingCriteria of the scheme.
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
    explicit TestingStrategy(const std::unordered_map<LocationId, std::vector<TestingScheme<FP>>>& location_to_schemes_map)
        : m_location_to_schemes_map(location_to_schemes_map.begin(), location_to_schemes_map.end())
    {
    }

    /**
     * @brief Add a TestingScheme to the set of schemes that are checked for testing at a certain Location.
     * @param[in] loc_id LocationId key for TestingScheme to be added.
     * @param[in] scheme TestingScheme to be added.
     */
    void add_testing_scheme(const LocationId& loc_id, const TestingScheme<FP>& scheme)
    {
        auto iter_schemes =
            std::find_if(m_location_to_schemes_map.begin(), m_location_to_schemes_map.end(), [loc_id](auto& p) {
                return p.first == loc_id;
            });
        if (iter_schemes == m_location_to_schemes_map.end()) {
            //no schemes for this location yet, add a new list with one scheme
            m_location_to_schemes_map.emplace_back(loc_id, std::vector<TestingScheme<FP>>(1, scheme));
        }
        else {
            //add scheme to existing vector if the scheme doesn't exist yet
            auto& schemes = iter_schemes->second;
            if (std::find(schemes.begin(), schemes.end(), scheme) == schemes.end()) {
                schemes.push_back(scheme);
            }
        }
    }

    /**
     * @brief Add a TestingScheme to the set of schemes that are checked for testing at a certain LocationType.
     * A TestingScheme applies to all Location of the same type is store in 
     * LocationId{INVALID_LOCATION_INDEX, location_type} of m_location_to_schemes_map.
     * @param[in] loc_type LocationId key for TestingScheme to be added.
     * @param[in] scheme TestingScheme to be added.
     */
    void add_testing_scheme(const LocationType& loc_type, const TestingScheme<FP>& scheme)
    {
        add_testing_scheme(LocationId{INVALID_LOCATION_INDEX, loc_type}, scheme);
    }

    /**
     * @brief Remove a TestingScheme from the set of schemes that are checked for testing at a certain Location.
     * @param[in] loc_id LocationId key for TestingScheme to be remove.
     * @param[in] scheme TestingScheme to be removed.
     */
    void remove_testing_scheme(const LocationId& loc_id, const TestingScheme<FP>& scheme)
    {
        auto iter_schemes =
            std::find_if(m_location_to_schemes_map.begin(), m_location_to_schemes_map.end(), [loc_id](auto& p) {
                return p.first == loc_id;
            });
        if (iter_schemes != m_location_to_schemes_map.end()) {
            //remove the scheme from the list
            auto& schemes_vector = iter_schemes->second;
            auto last            = std::remove(schemes_vector.begin(), schemes_vector.end(), scheme);
            schemes_vector.erase(last, schemes_vector.end());
            //delete the list of schemes for this location if no schemes left
            if (schemes_vector.empty()) {
                m_location_to_schemes_map.erase(iter_schemes);
            }
        }
    }

    /**
     * @brief Remove a TestingScheme from the set of schemes that are checked for testing at a certain Location.
     * A TestingScheme applies to all Location of the same type is store in 
     * LocationId{INVALID_LOCATION_INDEX, location_type} of m_location_to_schemes_map.
     * @param[in] loc_type LocationType key for TestingScheme to be remove.
     * @param[in] scheme TestingScheme to be removed.
     */
    void remove_testing_scheme(const LocationType& loc_type, const TestingScheme<FP>& scheme)
    {
        remove_testing_scheme(LocationId{INVALID_LOCATION_INDEX, loc_type}, scheme);
    }

    /**
     * @brief Checks if the given TimePoint is within the interval of start and end date of each TestingScheme and then
     * changes the activity status for each TestingScheme accordingly.
     * @param t TimePoint to check the activity status of each TestingScheme.
     */
    void update_activity_status(const TimePoint t)
    {
        for (auto& [_, testing_schemes] : m_location_to_schemes_map) {
            for (auto& scheme : testing_schemes) {
                scheme.update_activity_status(t);
            }
        }
    }

    /**
     * @brief Runs the TestingStrategy and potentially tests a Person.
     * @param[inout] rng Person::RandomNumberGenerator for the Person being tested.
     * @param[in] person Person to check.
     * @param[in] location Location to check.
     * @param[in] t TimePoint when to run the strategy.
     * @return If the Person is allowed to enter the Location.
     */
    bool run_strategy(typename Person<FP>::RandomNumberGenerator& rng,
                      Person<FP>& person, const Location<FP>& location, TimePoint t)
    {
        // Person who is in quarantine but not yet home should go home. Otherwise they can't because they test positive.
        if (location.get_type() == mio::abm::LocationType::Home && person.is_in_quarantine()) {
            return true;
        }

        //lookup schemes for this specific location as well as the location type
        //lookup in std::vector instead of std::map should be much faster unless for large numbers of schemes
        for (auto loc_key : {LocationId{location.get_index(), location.get_type()},
                             LocationId{INVALID_LOCATION_INDEX, location.get_type()}}) {
            auto iter_schemes =
                std::find_if(m_location_to_schemes_map.begin(), m_location_to_schemes_map.end(), [loc_key](auto& p) {
                    return p.first == loc_key;
                });
            if (iter_schemes != m_location_to_schemes_map.end()) {
                //apply all testing schemes that are found
                auto& schemes = iter_schemes->second;
                if (!std::all_of(schemes.begin(), schemes.end(), [&rng, &person, t](TestingScheme<FP>& ts) {
                        return !ts.is_active() || ts.run_scheme(rng, person, t);
                    })) {
                    return false;
                }
            }
        }
        return true;
    }

private:
    std::vector<std::pair<LocationId, std::vector<TestingScheme<FP> > > >
        m_location_to_schemes_map; ///< Set of schemes that are checked for testing.
};

} // namespace abm
} // namespace mio

#endif
