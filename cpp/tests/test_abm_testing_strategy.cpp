/* 
* Copyright (C) 2020-2024 MEmilio
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
#include "memilio/utils/random_number_generator.h"

TEST(TestTestingCriteria, addRemoveAndEvaluateTestCriteria)
{
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location work(mio::abm::LocationType::Work, 0, num_age_groups);
    auto person = make_test_person(home, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms);

    mio::abm::TimePoint t{0};
    auto testing_criteria = mio::abm::TestingCriteria();
    EXPECT_EQ(testing_criteria.evaluate(person, t), true);
    testing_criteria.add_infection_state(mio::abm::InfectionState::InfectedSymptoms);
    testing_criteria.add_infection_state(mio::abm::InfectionState::InfectedNoSymptoms);

    testing_criteria.add_age_group(age_group_35_to_59);
    EXPECT_EQ(testing_criteria.evaluate(person, t),
              false); // now it isn't empty and get's evaluated against age group
    testing_criteria.remove_age_group(age_group_35_to_59);
    EXPECT_EQ(testing_criteria.evaluate(person, t), true);

    testing_criteria.remove_infection_state(mio::abm::InfectionState::InfectedSymptoms);
    EXPECT_EQ(testing_criteria.evaluate(person, t), false);
    testing_criteria.add_infection_state(mio::abm::InfectionState::InfectedSymptoms);

    auto testing_criteria_manual = mio::abm::TestingCriteria(
        std::vector<mio::AgeGroup>({age_group_15_to_34}),
        std::vector<mio::abm::InfectionState>({mio::abm::InfectionState::InfectedNoSymptoms}));
    EXPECT_EQ(testing_criteria == testing_criteria_manual, false);
    testing_criteria_manual.add_infection_state(mio::abm::InfectionState::InfectedSymptoms);
    testing_criteria_manual.remove_age_group(age_group_15_to_34);
    EXPECT_EQ(testing_criteria == testing_criteria_manual, true);
}

TEST(TestTestingScheme, runScheme)
{
    auto rng = mio::RandomNumberGenerator();

    std::vector<mio::abm::InfectionState> test_infection_states1 = {mio::abm::InfectionState::InfectedSymptoms,
                                                                    mio::abm::InfectionState::InfectedNoSymptoms};
    std::vector<mio::abm::LocationType> test_location_types1     = {mio::abm::LocationType::Home,
                                                                    mio::abm::LocationType::Work};

    auto testing_criteria1                                   = mio::abm::TestingCriteria({}, test_infection_states1);
    std::vector<mio::abm::TestingCriteria> testing_criterias = {testing_criteria1};

    auto validity_period       = mio::abm::days(1);
    const auto start_date      = mio::abm::TimePoint(0);
    const auto end_date        = mio::abm::TimePoint(60 * 60 * 24 * 3);
    const auto probability     = 0.8;
    const auto test_params_pcr = mio::abm::TestParameters{0.9, 0.99, mio::abm::hours(48), mio::abm::TestType::PCR};

    std::vector<mio::abm::InfectionState> test_infection_states = {mio::abm::InfectionState::InfectedSymptoms,
                                                                   mio::abm::InfectionState::InfectedNoSymptoms};

    auto testing_scheme1 =
        mio::abm::TestingScheme(testing_criteria1, validity_period, start_date, end_date, test_params_pcr, probability);

    EXPECT_EQ(testing_scheme1.is_active(), false);
    testing_scheme1.update_activity_status(mio::abm::TimePoint(10));
    EXPECT_EQ(testing_scheme1.is_active(), true);
    testing_scheme1.update_activity_status(mio::abm::TimePoint(60 * 60 * 24 * 3 + 200));
    EXPECT_EQ(testing_scheme1.is_active(), false);
    testing_scheme1.update_activity_status(mio::abm::TimePoint(0));

    std::vector<mio::abm::InfectionState> test_infection_states2 = {mio::abm::InfectionState::Recovered};
    auto testing_criteria2 = mio::abm::TestingCriteria({}, test_infection_states2);
    auto testing_scheme2 =
        mio::abm::TestingScheme(testing_criteria2, validity_period, start_date, end_date, test_params_pcr, probability);

    mio::abm::Location loc_home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location loc_work(mio::abm::LocationType::Work, 0, num_age_groups);
    // Since tests are performed before start_date, the InfectionState of all the Person have to take into account the test's required_time
    auto person1     = make_test_person(loc_home, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms,
                                        start_date - test_params_pcr.required_time);
    auto rng_person1 = mio::abm::PersonalRandomNumberGenerator(rng, person1);
    auto person2     = make_test_person(loc_home, age_group_15_to_34, mio::abm::InfectionState::Recovered,
                                        start_date - test_params_pcr.required_time);
    auto rng_person2 = mio::abm::PersonalRandomNumberGenerator(rng, person2);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly(5))
        .WillOnce(testing::Return(0.7)) // Person 1 got test
        .WillOnce(testing::Return(0.7)) // Test is positive
        .WillOnce(testing::Return(0.5)) // Person 1 complies to isolation
        .WillOnce(testing::Return(0.7)) // Person 2 got test
        .WillOnce(testing::Return(0.5)); // Person 2 tested negative and can enter

    EXPECT_EQ(testing_scheme1.run_scheme(rng_person1, person1, start_date),
              false); // Person tests and tests positive
    EXPECT_EQ(testing_scheme2.run_scheme(rng_person2, person2, start_date),
              true); // Person tests and tests negative
    EXPECT_EQ(testing_scheme1.run_scheme(rng_person1, person1, start_date),
              false); // Person doesn't test but used the last result (false to enter)
}

TEST(TestTestingScheme, initAndRunTestingStrategy)
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
    testing_scheme1.update_activity_status(mio::abm::TimePoint(0));
    std::vector<mio::abm::InfectionState> test_infection_states2 = {mio::abm::InfectionState::Recovered};
    auto testing_criteria2 = mio::abm::TestingCriteria({}, test_infection_states2);
    auto testing_scheme2 =
        mio::abm::TestingScheme(testing_criteria2, validity_period, start_date, end_date, test_params_pcr, probability);
    testing_scheme2.update_activity_status(mio::abm::TimePoint(0));
    mio::abm::Location loc_work(mio::abm::LocationType::Work, 0);
    // Since tests are performed before start_date, the InfectionState of all the Person have to take into account the test's required_time
    auto person1     = make_test_person(loc_work, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms,
                                        start_date - test_params_pcr.required_time);
    auto rng_person1 = mio::abm::PersonalRandomNumberGenerator(rng, person1);
    auto person2     = make_test_person(loc_work, age_group_15_to_34, mio::abm::InfectionState::Recovered,
                                        start_date - test_params_pcr.required_time);
    auto rng_person2 = mio::abm::PersonalRandomNumberGenerator(rng, person2);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly((8)))
        .WillOnce(testing::Return(0.7)) // Person 1 complies to testing
        .WillOnce(testing::Return(0.7)) // Person 1 is tested for scheme 1
        .WillOnce(testing::Return(0.7)) // Test of Person 1 is positive
        .WillOnce(testing::Return(0.7)) // Person 1 complies to isolation
        .WillOnce(testing::Return(0.7)) // Person 2 complies to testing
        .WillOnce(testing::Return(0.7)) // Person 2 is tested for scheme 2
        .WillOnce(testing::Return(0.5)) // Test of Person 2 is negative
        .WillOnce(testing::Return(0.7)); // Person 1 complies to testing

    mio::abm::TestingStrategy test_strategy =
        mio::abm::TestingStrategy(std::vector<mio::abm::TestingStrategy::LocalStrategy>{});
    test_strategy.add_testing_scheme(mio::abm::LocationType::Work, testing_scheme1);
    test_strategy.add_testing_scheme(mio::abm::LocationType::Work, testing_scheme2);
    EXPECT_EQ(test_strategy.run_strategy(rng_person1, person1, loc_work, start_date),
              false); // Person tests and tests positive
    EXPECT_EQ(test_strategy.run_strategy(rng_person2, person2, loc_work, start_date),
              true); // Person tests and tests negative
    EXPECT_EQ(test_strategy.run_strategy(rng_person1, person1, loc_work, start_date),
              false); // Person doesn't test but used the last result (false to enter)
}
