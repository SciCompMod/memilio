/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Daniel Abele, Elisabeth Kluth, David Kerkmann, Sascha Korf, Martin J. Kuehn, Khoa Nguyen
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
#include "abm_helpers.h"
#include "random_number_test.h"

using TestTestingCriteria = RandomNumberTest;
/**
 * @brief Test for adding/removing age groups and infection states in TestingCriteria.
 */
TEST_F(TestTestingCriteria, addRemoveAndEvaluateTestCriteria)
{
    // Create test locations and a person in a specific infection state.
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location work(mio::abm::LocationType::Work, 0, num_age_groups);
    auto person =
        make_test_person(this->get_rng(), home, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms);

    mio::abm::TimePoint t{0};
    // Initialize testing criteria with no age group or infection state set.
    auto testing_criteria_empty = mio::abm::TestingCriteria();
    // Empty criteria should evaluate to true.
    EXPECT_EQ(testing_criteria_empty.evaluate(person, t), true);

    // Add infection states to the criteria.
    std::vector<mio::abm::InfectionState> test_infection_states = {mio::abm::InfectionState::InfectedSymptoms,
                                                                   mio::abm::InfectionState::InfectedNoSymptoms};
    std::vector<mio::AgeGroup> test_age_groups                  = {age_group_35_to_59};
    auto testing_criteria_false_ag = mio::abm::TestingCriteria(test_age_groups, test_infection_states);
    // Age group mismatch, should evaluate to false.
    EXPECT_EQ(testing_criteria_false_ag.evaluate(person, t), false);

    // Add age groups to the criteria.
    test_age_groups.push_back(age_group_15_to_34);
    auto testing_criteria_right_ag = mio::abm::TestingCriteria(test_age_groups, test_infection_states);
    // Now it should evaluate to true.
    EXPECT_EQ(testing_criteria_right_ag.evaluate(person, t), true);

    // Test inequality of testing criteria.
    EXPECT_EQ(testing_criteria_right_ag == testing_criteria_false_ag, false);
}

using TestTestingScheme = RandomNumberTest;

/**
 * @brief Test for checking TestingScheme's activity status and behavior during runtime.
 */
TEST_F(TestTestingScheme, runScheme)
{
    std::vector<mio::abm::InfectionState> test_infection_states1 = {mio::abm::InfectionState::InfectedSymptoms,
                                                                    mio::abm::InfectionState::InfectedNoSymptoms};
    std::vector<mio::abm::LocationType> test_location_types1     = {mio::abm::LocationType::Home,
                                                                    mio::abm::LocationType::Work};

    auto testing_criteria1                                   = mio::abm::TestingCriteria({}, test_infection_states1);
    std::vector<mio::abm::TestingCriteria> testing_criterias = {testing_criteria1};

    auto validity_period       = mio::abm::days(1);
    const auto start_date      = mio::abm::TimePoint(0) + mio::abm::seconds(1);
    const auto end_date        = mio::abm::TimePoint(0) + mio::abm::days(3);
    const auto probability     = 0.8;
    const auto test_params_pcr = mio::abm::TestParameters{0.9, 0.99, mio::abm::hours(48), mio::abm::TestType::PCR};

    std::vector<mio::abm::InfectionState> test_infection_states = {mio::abm::InfectionState::InfectedSymptoms,
                                                                   mio::abm::InfectionState::InfectedNoSymptoms};

    // Create a testing scheme based on the criteria and parameters.
    auto testing_scheme1 =
        mio::abm::TestingScheme(testing_criteria1, validity_period, start_date, end_date, test_params_pcr, probability);

    // Check the initial activity status.
    EXPECT_EQ(testing_scheme1.is_active(mio::abm::TimePoint(0)), false);
    EXPECT_EQ(testing_scheme1.is_active(mio::abm::TimePoint(10)), true);

    // Deactivate the scheme after the end date.
    EXPECT_EQ(testing_scheme1.is_active(end_date + mio::abm::seconds(200)), false);

    // Setup a second scheme with different infection states.
    std::vector<mio::abm::InfectionState> test_infection_states2 = {mio::abm::InfectionState::Recovered};
    auto testing_criteria2 = mio::abm::TestingCriteria({}, test_infection_states2);
    auto testing_scheme2 =
        mio::abm::TestingScheme(testing_criteria2, validity_period, start_date, end_date, test_params_pcr, probability);

    mio::abm::Location loc_home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location loc_work(mio::abm::LocationType::Work, 0, num_age_groups);
    // Since tests are performed before start_date, the InfectionState of all the Person have to take into account the test's required_time
    auto person1 =
        make_test_person(this->get_rng(), loc_home, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms,
                         start_date - test_params_pcr.required_time);
    auto rng_person1 = mio::abm::PersonalRandomNumberGenerator(person1);
    auto person2 = make_test_person(this->get_rng(), loc_home, age_group_15_to_34, mio::abm::InfectionState::Recovered,
                                    start_date - test_params_pcr.required_time);
    auto rng_person2 = mio::abm::PersonalRandomNumberGenerator(person2);

    // Mock uniform distribution to control random behavior in testing.
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly(5))
        .WillOnce(testing::Return(0.7)) // Person 1 got test
        .WillOnce(testing::Return(0.7)) // Test is positive
        .WillOnce(testing::Return(0.5)) // Person 1 complies to isolation
        .WillOnce(testing::Return(0.7)) // Person 2 got test
        .WillOnce(testing::Return(0.5)); // Person 2 tested negative and can enter

    EXPECT_EQ(testing_scheme1.run_and_test(rng_person1, person1, start_date),
              true); // Person tests and tests positive
    EXPECT_EQ(testing_scheme2.run_and_test(rng_person2, person2, start_date),
              false); // Person tests and tests negative
    EXPECT_EQ(testing_scheme1.run_and_test(rng_person1, person1, start_date),
              true); // Person doesn't test but used the last result (false to enter)
}

/**
 * @brief Test for TestingStrategy initialization and execution.
 */
TEST_F(TestTestingScheme, initAndRunTestingStrategy)
{
    auto rng                   = mio::RandomNumberGenerator();
    auto validity_period       = mio::abm::days(1);
    const auto start_date      = mio::abm::TimePoint(0);
    const auto end_date        = mio::abm::TimePoint(60 * 60 * 24 * 3);
    const auto probability     = 0.8;
    const auto test_params_pcr = mio::abm::TestParameters{0.9, 0.99, mio::abm::hours(48), mio::abm::TestType::PCR};

    std::vector<mio::abm::InfectionState> test_infection_states = {mio::abm::InfectionState::InfectedSymptoms,
                                                                   mio::abm::InfectionState::InfectedNoSymptoms};
    auto testing_criteria1                                      = mio::abm::TestingCriteria({}, test_infection_states);
    auto testing_scheme1 =
        mio::abm::TestingScheme(testing_criteria1, validity_period, start_date, end_date, test_params_pcr, probability);
    std::vector<mio::abm::InfectionState> test_infection_states2 = {mio::abm::InfectionState::Recovered};
    auto testing_criteria2 = mio::abm::TestingCriteria({}, test_infection_states2);
    auto testing_scheme2 =
        mio::abm::TestingScheme(testing_criteria2, validity_period, start_date, end_date, test_params_pcr, probability);
    mio::abm::Location loc_work(mio::abm::LocationType::Work, 0);
    // Since tests are performed before start_date, the InfectionState of all the Person have to take into account the test's required_time
    auto person1 =
        make_test_person(this->get_rng(), loc_work, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms,
                         start_date - test_params_pcr.required_time);
    auto rng_person1 = mio::abm::PersonalRandomNumberGenerator(person1);
    auto person2 = make_test_person(this->get_rng(), loc_work, age_group_15_to_34, mio::abm::InfectionState::Recovered,
                                    start_date - test_params_pcr.required_time);
    auto rng_person2 = mio::abm::PersonalRandomNumberGenerator(person2);

    // Mock uniform distribution to control random behavior in testing.
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly((5)))
        .WillOnce(testing::Return(0.7)) // Person 1 is tested for scheme 1
        .WillOnce(testing::Return(0.7)) // Test of Person 1 is positive
        .WillOnce(testing::Return(0.7)) // Person 1 complies to isolation
        .WillOnce(testing::Return(0.7)) // Person 2 is tested for scheme 2
        .WillOnce(testing::Return(0.5)); // Test of Person 2 is negative

    mio::abm::TestingStrategy test_strategy =
        mio::abm::TestingStrategy(std::vector<mio::abm::TestingStrategy::LocalStrategy>{},
                                  std::vector<mio::abm::TestingStrategy::LocalStrategy>{});
    test_strategy.add_scheme(mio::abm::LocationType::Work, testing_scheme1);
    test_strategy.add_scheme(mio::abm::LocationType::Work, testing_scheme2);
    EXPECT_EQ(test_strategy.run_and_check(rng_person1, person1, loc_work, start_date),
              false); // Person tests and tests positive
    EXPECT_EQ(test_strategy.run_and_check(rng_person2, person2, loc_work, start_date),
              true); // Person tests and tests negative
    EXPECT_EQ(test_strategy.run_and_check(rng_person1, person1, loc_work, start_date),
              false); // Person doesn't test but used the last result (false to enter)
}

/**
 * @brief Test for edge cases in TestingCriteria.
 */
TEST_F(TestTestingCriteria, testingCriteriaEdgeCases)
{
    // Create test locations and persons with different age groups and infection states
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);

    // Test with various infection states
    auto person_exposed =
        make_test_person(this->get_rng(), home, age_group_15_to_34, mio::abm::InfectionState::Exposed);
    auto person_symptoms =
        make_test_person(this->get_rng(), home, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms);
    auto person_no_symptoms =
        make_test_person(this->get_rng(), home, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms);
    auto person_recovered =
        make_test_person(this->get_rng(), home, age_group_15_to_34, mio::abm::InfectionState::Recovered);

    mio::abm::TimePoint t{0};

    // Test with only infection states criteria
    std::vector<mio::abm::InfectionState> test_infection_states = {mio::abm::InfectionState::InfectedSymptoms,
                                                                   mio::abm::InfectionState::InfectedNoSymptoms};
    auto testing_criteria_infection                             = mio::abm::TestingCriteria({}, test_infection_states);

    // Should match only infected persons
    EXPECT_EQ(testing_criteria_infection.evaluate(person_exposed, t), false);
    EXPECT_EQ(testing_criteria_infection.evaluate(person_symptoms, t), true);
    EXPECT_EQ(testing_criteria_infection.evaluate(person_no_symptoms, t), true);
    EXPECT_EQ(testing_criteria_infection.evaluate(person_recovered, t), false);

    // Test with only age group criteria
    std::vector<mio::AgeGroup> test_age_groups = {age_group_15_to_34, age_group_35_to_59};
    auto testing_criteria_age                  = mio::abm::TestingCriteria(test_age_groups, {});

    // Create persons with different age groups
    auto person_young =
        make_test_person(this->get_rng(), home, age_group_5_to_14, mio::abm::InfectionState::Susceptible);
    auto person_adult =
        make_test_person(this->get_rng(), home, age_group_15_to_34, mio::abm::InfectionState::Susceptible);
    auto person_older =
        make_test_person(this->get_rng(), home, age_group_35_to_59, mio::abm::InfectionState::Susceptible);
    auto person_senior =
        make_test_person(this->get_rng(), home, age_group_60_to_79, mio::abm::InfectionState::Susceptible);

    // Should match only specified age groups
    EXPECT_EQ(testing_criteria_age.evaluate(person_young, t), false);
    EXPECT_EQ(testing_criteria_age.evaluate(person_adult, t), true);
    EXPECT_EQ(testing_criteria_age.evaluate(person_older, t), true);
    EXPECT_EQ(testing_criteria_age.evaluate(person_senior, t), false);

    // Test with both age and infection state criteria
    auto testing_criteria_both = mio::abm::TestingCriteria(test_age_groups, test_infection_states);

    // Should match only when both criteria are met
    auto person_adult_infected =
        make_test_person(this->get_rng(), home, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms);
    auto person_young_infected =
        make_test_person(this->get_rng(), home, age_group_5_to_14, mio::abm::InfectionState::InfectedSymptoms);
    auto person_adult_recovered =
        make_test_person(this->get_rng(), home, age_group_15_to_34, mio::abm::InfectionState::Recovered);

    EXPECT_EQ(testing_criteria_both.evaluate(person_adult_infected, t), true);
    EXPECT_EQ(testing_criteria_both.evaluate(person_young_infected, t), false);
    EXPECT_EQ(testing_criteria_both.evaluate(person_adult_recovered, t), false);
}

/**
 * @brief Test for TestingScheme time-related functionality.
 */
TEST_F(TestTestingScheme, testingSchemeTimeValidity)
{
    auto testing_criteria  = mio::abm::TestingCriteria();
    auto validity_period   = mio::abm::days(1);
    const auto start_date  = mio::abm::TimePoint(100);
    const auto end_date    = mio::abm::TimePoint(500);
    const auto probability = 0.8;
    const auto test_params = mio::abm::TestParameters{0.9, 0.99, mio::abm::hours(48), mio::abm::TestType::PCR};

    auto testing_scheme =
        mio::abm::TestingScheme(testing_criteria, validity_period, start_date, end_date, test_params, probability);

    // Test boundary conditions of is_active method
    EXPECT_EQ(testing_scheme.is_active(mio::abm::TimePoint(99)), false); // Just before start
    EXPECT_EQ(testing_scheme.is_active(mio::abm::TimePoint(100)), true); // At start
    EXPECT_EQ(testing_scheme.is_active(mio::abm::TimePoint(300)), true); // Middle
    EXPECT_EQ(testing_scheme.is_active(mio::abm::TimePoint(499)), true); // Just before end
    EXPECT_EQ(testing_scheme.is_active(mio::abm::TimePoint(500)), false); // At end
    EXPECT_EQ(testing_scheme.is_active(mio::abm::TimePoint(501)), false); // After end
}

/**
 * @brief Test for TestingScheme result caching behavior.
 */
TEST_F(TestTestingScheme, testingSchemeResultCaching)
{
    auto validity_period       = mio::abm::days(1);
    const auto start_date      = mio::abm::TimePoint(0);
    const auto end_date        = mio::abm::TimePoint(0) + mio::abm::days(5);
    const auto probability     = 0.8;
    const auto test_params_pcr = mio::abm::TestParameters{0.9, 0.99, mio::abm::hours(2), mio::abm::TestType::PCR};

    auto testing_criteria = mio::abm::TestingCriteria({}, {});
    auto testing_scheme =
        mio::abm::TestingScheme(testing_criteria, validity_period, start_date, end_date, test_params_pcr, probability);

    // Create test person and location
    mio::abm::Location loc_home(mio::abm::LocationType::Home, 0, num_age_groups);
    auto person = make_test_person(this->get_rng(), loc_home, age_group_15_to_34,
                                   mio::abm::InfectionState::InfectedNoSymptoms, start_date);
    auto rng    = mio::abm::PersonalRandomNumberGenerator(person);

    // Mock uniform distribution to control test results
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .WillOnce(testing::Return(0.7)) // Test is performed
        .WillOnce(testing::Return(0.8)) // First test is positive
        .WillOnce(testing::Return(0.5)); // Isolation compliance

    // First test at t=10
    auto t1      = mio::abm::TimePoint(10);
    bool result1 = testing_scheme.run_and_test(rng, person, t1);
    EXPECT_EQ(result1, true); // We expect the test to be positive

    // Test result should be cached for validity period
    auto t2      = t1 + mio::abm::hours(12); // Within validity period
    bool result2 = testing_scheme.run_and_test(rng, person, t2);
    EXPECT_EQ(result2, true); // Should use cached result without calling RNG

    // After validity period expires, a new test should be performed
    // But we didn't mock additional RNG calls, so this would fail if it tries to run a new test
}

/**
 * @brief Test for different test types and parameters.
 */
TEST_F(TestTestingScheme, differentTestTypes)
{
    auto validity_period   = mio::abm::days(1);
    const auto start_date  = mio::abm::TimePoint(0);
    const auto end_date    = mio::abm::TimePoint(0) + mio::abm::days(5);
    const auto probability = 1.0; // Always test

    // Create PCR test parameters with high accuracy but long wait time
    const auto test_params_pcr = mio::abm::TestParameters{0.95, 0.99, mio::abm::hours(24), mio::abm::TestType::PCR};

    // Create rapid test parameters with lower accuracy (we need to set this to 24 hours to avoid the person getting healthy randomly)
    const auto test_params_rapid = mio::abm::TestParameters{0.8, 0.9, mio::abm::hours(24), mio::abm::TestType::Antigen};

    auto testing_criteria = mio::abm::TestingCriteria();

    // Create test schemes for both test types
    auto testing_scheme_pcr =
        mio::abm::TestingScheme(testing_criteria, validity_period, start_date, end_date, test_params_pcr, probability);

    auto testing_scheme_rapid = mio::abm::TestingScheme(testing_criteria, validity_period, start_date, end_date,
                                                        test_params_rapid, probability);

    // Create test persons with different infection states
    mio::abm::Location loc_home(mio::abm::LocationType::Home, 0, num_age_groups);
    auto person_infected =
        make_test_person(this->get_rng(), loc_home, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms,
                         start_date - test_params_pcr.required_time);
    auto rng_infected = mio::abm::PersonalRandomNumberGenerator(person_infected);

    auto person_healthy =
        make_test_person(this->get_rng(), loc_home, age_group_15_to_34, mio::abm::InfectionState::Susceptible,
                         start_date - test_params_pcr.required_time);
    auto rng_healthy = mio::abm::PersonalRandomNumberGenerator(person_healthy);

    // Mock uniform distribution to control test results
    // Mock uniform distribution to control test results for PCR test with infected person
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly(3))
        .WillOnce(testing::Return(0.7)) // PCR test is performed
        .WillOnce(testing::Return(0.05)) // PCR test correctly identifies infection (< sensitivity 0.95)
        .WillOnce(testing::Return(0.5)); // Person complies to isolation

    // Test PCR test with infected person
    bool pcr_infected_result = testing_scheme_pcr.run_and_test(rng_infected, person_infected, start_date);
    EXPECT_EQ(pcr_infected_result, true); // PCR should detect infection

    // Reset mock for PCR test with healthy person
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly(2))
        .WillOnce(testing::Return(0.7)) // PCR test is performed
        .WillOnce(testing::Return(0.98)); // PCR test correctly identifies no infection (< specificity 0.99)

    // Test PCR test with healthy person
    bool pcr_healthy_result = testing_scheme_pcr.run_and_test(rng_healthy, person_healthy, start_date);
    EXPECT_EQ(pcr_healthy_result, false); // PCR should correctly identify no infection

    // Reset mock for rapid test with infected person

    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly(3))
        .WillOnce(testing::Return(0.7)) // Rapid test is performed
        .WillOnce(testing::Return(0.0)) // Rapid test correctly identifies infection (< sensitivity 0.8)
        .WillOnce(testing::Return(0.1)); // Infected person complies for isolation

    // Test rapid test with infected person
    bool rapid_infected_result = testing_scheme_rapid.run_and_test(rng_infected, person_infected, start_date);
    EXPECT_EQ(rapid_infected_result, true); // Rapid test should detect infection
}

/**
 * @brief Test for combining multiple testing schemes in TestingStrategy.
 */
TEST_F(TestTestingScheme, multipleSchemesCombination)
{
    auto validity_period   = mio::abm::days(1);
    const auto start_date  = mio::abm::TimePoint(0);
    const auto end_date    = mio::abm::TimePoint(500);
    const auto probability = 0.8;
    const auto test_params = mio::abm::TestParameters{0.9, 0.99, mio::abm::hours(2), mio::abm::TestType::PCR};

    // Create testing schemes for different age groups
    std::vector<mio::AgeGroup> children_age_groups = {age_group_0_to_4, age_group_5_to_14};
    std::vector<mio::AgeGroup> adult_age_groups    = {age_group_15_to_34, age_group_35_to_59};
    std::vector<mio::AgeGroup> senior_age_groups   = {age_group_60_to_79, age_group_80_plus};

    auto testing_criteria_children = mio::abm::TestingCriteria(children_age_groups, {});
    auto testing_criteria_adults   = mio::abm::TestingCriteria(adult_age_groups, {});
    auto testing_criteria_seniors  = mio::abm::TestingCriteria(senior_age_groups, {});

    auto testing_scheme_children = mio::abm::TestingScheme(testing_criteria_children, validity_period, start_date,
                                                           end_date, test_params, probability);
    auto testing_scheme_adults = mio::abm::TestingScheme(testing_criteria_adults, validity_period, start_date, end_date,
                                                         test_params, probability);
    auto testing_scheme_seniors = mio::abm::TestingScheme(testing_criteria_seniors, validity_period, start_date,
                                                          end_date, test_params, probability);

    // Create TestingStrategy
    mio::abm::TestingStrategy test_strategy;

    // Add schemes to different location types
    test_strategy.add_scheme(mio::abm::LocationType::School, testing_scheme_children);
    test_strategy.add_scheme(mio::abm::LocationType::Work, testing_scheme_adults);
    test_strategy.add_scheme(mio::abm::LocationType::SocialEvent, testing_scheme_seniors);

    // Also add multiple location types at once
    std::vector<mio::abm::LocationType> public_locations = {mio::abm::LocationType::BasicsShop,
                                                            mio::abm::LocationType::SocialEvent};
    test_strategy.add_scheme(public_locations, testing_scheme_adults);

    // Create locations
    mio::abm::Location loc_home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location loc_school(mio::abm::LocationType::School, 1, num_age_groups);
    mio::abm::Location loc_work(mio::abm::LocationType::Work, 2, num_age_groups);
    mio::abm::Location loc_shop(mio::abm::LocationType::BasicsShop, 3, num_age_groups);

    // Create persons of different age groups
    auto child = make_test_person(this->get_rng(), loc_home, age_group_5_to_14, mio::abm::InfectionState::Susceptible,
                                  start_date);
    auto adult = make_test_person(this->get_rng(), loc_home, age_group_35_to_59, mio::abm::InfectionState::Susceptible,
                                  start_date);

    auto rng_child = mio::abm::PersonalRandomNumberGenerator(child);
    auto rng_adult = mio::abm::PersonalRandomNumberGenerator(adult);

    // Mock uniform distribution to control test results
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly(5))
        .WillOnce(testing::Return(0.7)) // Child gets tested at school
        .WillOnce(testing::Return(0.95)) // Child tests negative
        .WillOnce(testing::Return(0.7)) // Adult gets tested at work
        .WillOnce(testing::Return(0.999)) // Adult tests (false) positive
        .WillOnce(testing::Return(0.5)); // Adult complies to isolation

    // Test child at school - should pass the test
    bool school_result = test_strategy.run_and_check(rng_child, child, loc_school, start_date);
    EXPECT_EQ(school_result, true); // Child should be allowed to enter school

    // Test adult at work - should fail the test
    bool work_result = test_strategy.run_and_check(rng_adult, adult, loc_work, start_date);
    EXPECT_EQ(work_result, false); // Adult should not be allowed to enter work

    // Test adult at shop - should use cached result without testing again
    bool shop_result = test_strategy.run_and_check(rng_adult, adult, loc_shop, start_date);
    EXPECT_EQ(shop_result, false); // Adult should not be allowed to enter shop either
}

/**
 * @brief Test for TestingStrategy with location-specific schemes.
 */
TEST_F(TestTestingScheme, locationSpecificSchemes)
{
    auto validity_period   = mio::abm::days(1);
    const auto start_date  = mio::abm::TimePoint(0);
    const auto end_date    = mio::abm::TimePoint(500);
    const auto probability = 0.8;
    const auto test_params = mio::abm::TestParameters{0.9, 0.99, mio::abm::hours(2), mio::abm::TestType::PCR};

    // Create a testing scheme
    auto testing_criteria = mio::abm::TestingCriteria();
    auto testing_scheme =
        mio::abm::TestingScheme(testing_criteria, validity_period, start_date, end_date, test_params, probability);

    // Create TestingStrategy
    mio::abm::TestingStrategy test_strategy;

    // Add scheme to a specific location by ID
    mio::abm::LocationId specific_shop_id(42);
    test_strategy.add_scheme(specific_shop_id, testing_scheme);

    // Create locations with different IDs but same type
    mio::abm::Location shop1(mio::abm::LocationType::BasicsShop, 42, num_age_groups); // Has the specific ID
    mio::abm::Location shop2(mio::abm::LocationType::BasicsShop, 43, num_age_groups); // Different ID

    // Create a test person
    auto person =
        make_test_person(this->get_rng(), shop1, age_group_15_to_34, mio::abm::InfectionState::Susceptible, start_date);
    auto rng = mio::abm::PersonalRandomNumberGenerator(person);

    // Mock uniform distribution to control test results
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly(3))
        .WillOnce(testing::Return(0.7)) // Person gets tested at shop1
        .WillOnce(testing::Return(0.995)) // Test is positive
        .WillOnce(testing::Return(0.5)); // Person complies to isolation

    // Test at shop with specific ID - should run the scheme
    bool result1 = test_strategy.run_and_check(rng, person, shop1, start_date);
    EXPECT_EQ(result1, false); // Person should not be allowed to enter after positive test

    // Test at shop with different ID - no scheme should run
    // No need to mock RNG calls here as no test should be performed
    bool result2 = test_strategy.run_and_check(rng, person, shop2, start_date);
    EXPECT_EQ(result2, true); // Person should be allowed to enter as no test is required
}

TEST_F(TestTestingScheme, testCompliance)
{
    auto validity_period   = mio::abm::days(1);
    const auto start_date  = mio::abm::TimePoint(0);
    const auto end_date    = mio::abm::TimePoint(500);
    const auto probability = 0.8;
    const auto test_params = mio::abm::TestParameters{0.9, 0.99, mio::abm::hours(2), mio::abm::TestType::PCR};

    // Create a testing scheme
    auto testing_criteria = mio::abm::TestingCriteria();
    auto testing_scheme =
        mio::abm::TestingScheme(testing_criteria, validity_period, start_date, end_date, test_params, probability);

    // Create TestingStrategy
    mio::abm::TestingStrategy test_strategy;

    // Add scheme to a specific location by ID
    mio::abm::LocationId specific_shop_id(42);
    test_strategy.add_scheme(specific_shop_id, testing_scheme);

    // Create locations with different IDs but same type
    mio::abm::Location shop1(mio::abm::LocationType::BasicsShop, 42, num_age_groups); // Has the specific ID

    // Create a test person
    auto person =
        make_test_person(this->get_rng(), shop1, age_group_15_to_34, mio::abm::InfectionState::Susceptible, start_date);
    auto rng = mio::abm::PersonalRandomNumberGenerator(person);
    person.set_compliance(mio::abm::InterventionType::Testing, 0.1); // Set compliance for testing

    // Mock uniform distribution to control test results
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly(1))
        .WillOnce(testing::Return(0.2)); // Person is not compliant for testing

    // Test at shop with specific ID - should run the scheme
    bool result1 = test_strategy.run_and_check(rng, person, shop1, start_date);
    EXPECT_EQ(result1, false); // Person should not be allowed to enter after not complying to the test
}
