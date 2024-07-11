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
    ASSERT_EQ(testing_criteria.evaluate(person, t), true);
    testing_criteria.add_infection_state(mio::abm::InfectionState::InfectedSymptoms);
    testing_criteria.add_infection_state(mio::abm::InfectionState::InfectedNoSymptoms);

    testing_criteria.add_age_group(age_group_35_to_59);
    ASSERT_EQ(testing_criteria.evaluate(person, t),
              false); // now it isn't empty and get's evaluated against age group
    testing_criteria.remove_age_group(age_group_35_to_59);
    ASSERT_EQ(testing_criteria.evaluate(person, t), true);

    testing_criteria.remove_infection_state(mio::abm::InfectionState::InfectedSymptoms);
    ASSERT_EQ(testing_criteria.evaluate(person, t), false);
    testing_criteria.add_infection_state(mio::abm::InfectionState::InfectedSymptoms);

    auto testing_criteria_manual = mio::abm::TestingCriteria(
        std::vector<mio::AgeGroup>({age_group_15_to_34}),
        std::vector<mio::abm::InfectionState>({mio::abm::InfectionState::InfectedNoSymptoms}));
    ASSERT_EQ(testing_criteria == testing_criteria_manual, false);
    testing_criteria_manual.add_infection_state(mio::abm::InfectionState::InfectedSymptoms);
    testing_criteria_manual.remove_age_group(age_group_15_to_34);
    ASSERT_EQ(testing_criteria == testing_criteria_manual, true);
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

    const auto testing_min_time = mio::abm::days(1);
    const auto start_date       = mio::abm::TimePoint(0);
    const auto end_date         = mio::abm::TimePoint(60 * 60 * 24 * 3);
    const auto probability      = 0.8;
    const auto test_type        = mio::abm::PCRTest();

    std::vector<mio::abm::InfectionState> test_infection_states = {mio::abm::InfectionState::InfectedSymptoms,
                                                                   mio::abm::InfectionState::InfectedNoSymptoms};

    auto testing_scheme1 =
        mio::abm::TestingScheme(testing_criteria1, testing_min_time, start_date, end_date, test_type, probability);

    ASSERT_EQ(testing_scheme1.is_active(), false);
    testing_scheme1.update_activity_status(mio::abm::TimePoint(10));
    ASSERT_EQ(testing_scheme1.is_active(), true);
    testing_scheme1.update_activity_status(mio::abm::TimePoint(60 * 60 * 24 * 3 + 200));
    ASSERT_EQ(testing_scheme1.is_active(), false);
    testing_scheme1.update_activity_status(mio::abm::TimePoint(0));

    std::vector<mio::abm::InfectionState> test_infection_states2 = {mio::abm::InfectionState::Recovered};
    auto testing_criteria2 = mio::abm::TestingCriteria({}, test_infection_states2);
    auto testing_scheme2 =
        mio::abm::TestingScheme(testing_criteria2, testing_min_time, start_date, end_date, test_type, probability);

    mio::abm::Location loc_home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location loc_work(mio::abm::LocationType::Work, 0, num_age_groups);
    auto person1     = make_test_person(loc_home, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms);
    auto rng_person1 = mio::abm::PersonalRandomNumberGenerator(rng, person1);
    auto person2     = make_test_person(loc_home, age_group_15_to_34, mio::abm::InfectionState::Recovered);
    auto rng_person2 = mio::abm::PersonalRandomNumberGenerator(rng, person2);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly(4))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.5))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.5));
    ASSERT_EQ(testing_scheme1.run_scheme(rng_person1, person1, start_date), false); // Person tests and tests positive
    ASSERT_EQ(testing_scheme2.run_scheme(rng_person2, person2, start_date), true); // Person tests and tests negative
    ASSERT_EQ(testing_scheme1.run_scheme(rng_person1, person1, start_date),
              true); // Person doesn't test
}

TEST(TestTestingScheme, initAndRunTestingStrategy)
{
    auto rng                    = mio::RandomNumberGenerator();
    const auto testing_min_time = mio::abm::days(1);
    const auto start_date       = mio::abm::TimePoint(0);
    const auto end_date         = mio::abm::TimePoint(60 * 60 * 24 * 3);
    const auto probability      = 0.8;
    const auto test_type        = mio::abm::PCRTest();

    std::vector<mio::abm::InfectionState> test_infection_states = {mio::abm::InfectionState::InfectedSymptoms,
                                                                   mio::abm::InfectionState::InfectedNoSymptoms};
    auto testing_criteria1                                      = mio::abm::TestingCriteria({}, test_infection_states);
    auto testing_scheme1 =
        mio::abm::TestingScheme(testing_criteria1, testing_min_time, start_date, end_date, test_type, probability);
    testing_scheme1.update_activity_status(mio::abm::TimePoint(0));
    std::vector<mio::abm::InfectionState> test_infection_states2 = {mio::abm::InfectionState::Recovered};
    auto testing_criteria2 = mio::abm::TestingCriteria({}, test_infection_states2);
    auto testing_scheme2 =
        mio::abm::TestingScheme(testing_criteria2, testing_min_time, start_date, end_date, test_type, probability);

    mio::abm::Location loc_work(mio::abm::LocationType::Work, 0);
    auto person1     = make_test_person(loc_work, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms);
    auto rng_person1 = mio::abm::PersonalRandomNumberGenerator(rng, person1);
    auto person2     = make_test_person(loc_work, age_group_15_to_34, mio::abm::InfectionState::Recovered);
    auto rng_person2 = mio::abm::PersonalRandomNumberGenerator(rng, person2);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly(2)) //only sampled twice, testing criteria don't apply to third person
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.5));

    mio::abm::TestingStrategy test_strategy =
        mio::abm::TestingStrategy(std::vector<mio::abm::TestingStrategy::LocalStrategy>{});
    test_strategy.add_testing_scheme(mio::abm::LocationType::Work, testing_scheme1);
    test_strategy.add_testing_scheme(mio::abm::LocationType::Work, testing_scheme2);
    ASSERT_EQ(test_strategy.run_strategy(rng_person1, person1, loc_work, start_date),
              false); // Person tests and tests positive
    ASSERT_EQ(test_strategy.run_strategy(rng_person2, person2, loc_work, start_date),
              true); // Person tests and tests negative
    ASSERT_EQ(test_strategy.run_strategy(rng_person1, person1, loc_work, start_date), true); // Person doesn't test
}
