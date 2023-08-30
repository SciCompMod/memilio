/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#include "abm_helpers.h"

TEST(TestTestingCriteria, addRemoveAndEvaluateTestCriteria)
{
    auto home   = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto work   = mio::abm::Location(mio::abm::LocationType::Work, 0);
    auto person = make_test_person(home, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSymptoms);

    mio::abm::TimePoint t{0};
    auto testing_criteria = mio::abm::TestingCriteria();
    ASSERT_EQ(testing_criteria.evaluate(person, work, t), true);
    testing_criteria.add_infection_state(mio::abm::InfectionState::InfectedSymptoms);
    testing_criteria.add_infection_state(mio::abm::InfectionState::InfectedNoSymptoms);
    testing_criteria.add_location_type(mio::abm::LocationType::Home);
    testing_criteria.add_location_type(mio::abm::LocationType::Work);

    ASSERT_EQ(testing_criteria.evaluate(person, work, t), true);
    ASSERT_EQ(testing_criteria.evaluate(person, home, t), true);

    testing_criteria.add_age_group(mio::abm::AgeGroup::Age35to59);
    ASSERT_EQ(testing_criteria.evaluate(person, home, t),
              false); // now it isn't empty and get's evaluated against age group
    testing_criteria.remove_age_group(mio::abm::AgeGroup::Age35to59);
    ASSERT_EQ(testing_criteria.evaluate(person, home, t), true);

    testing_criteria.remove_infection_state(mio::abm::InfectionState::InfectedSymptoms);
    ASSERT_EQ(testing_criteria.evaluate(person, home, t), false);

    testing_criteria.add_infection_state(mio::abm::InfectionState::InfectedSymptoms);
    testing_criteria.remove_location_type(mio::abm::LocationType::Home);
    ASSERT_EQ(testing_criteria.evaluate(person, home, t), false);

    auto testing_criteria_manual =
        mio::abm::TestingCriteria({}, std::vector<mio::abm::LocationType>({mio::abm::LocationType::Work}),
                                  std::vector<mio::abm::InfectionState>({mio::abm::InfectionState::InfectedNoSymptoms,
                                                                         mio::abm::InfectionState::InfectedSymptoms}));
    ASSERT_EQ(testing_criteria == testing_criteria_manual, true);
    testing_criteria_manual.remove_infection_state(mio::abm::InfectionState::InfectedSymptoms);
    ASSERT_EQ(testing_criteria == testing_criteria_manual, false);
}

TEST(TestTestingScheme, runScheme)
{
    std::vector<mio::abm::InfectionState> test_infection_states1 = {mio::abm::InfectionState::InfectedSymptoms,
                                                                    mio::abm::InfectionState::InfectedNoSymptoms};
    std::vector<mio::abm::LocationType> test_location_types1     = {mio::abm::LocationType::Home,
                                                                    mio::abm::LocationType::Work};

    auto testing_criteria1 = mio::abm::TestingCriteria({}, test_location_types1, test_infection_states1);
    std::vector<mio::abm::TestingCriteria> testing_criterias = {testing_criteria1};

    const auto testing_min_time = mio::abm::days(1);
    const auto start_date       = mio::abm::TimePoint(0);
    const auto end_date         = mio::abm::TimePoint(60 * 60 * 24 * 3);
    const auto probability      = 0.8;
    const auto test_type        = mio::abm::PCRTest();

    auto testing_scheme =
        mio::abm::TestingScheme(testing_criterias, testing_min_time, start_date, end_date, test_type, probability);

    ASSERT_EQ(testing_scheme.is_active(), false);
    testing_scheme.update_activity_status(mio::abm::TimePoint(10));
    ASSERT_EQ(testing_scheme.is_active(), true);
    testing_scheme.update_activity_status(mio::abm::TimePoint(60 * 60 * 24 * 3 + 200));
    ASSERT_EQ(testing_scheme.is_active(), false);
    testing_scheme.update_activity_status(mio::abm::TimePoint(0));

    std::vector<mio::abm::InfectionState> test_infection_states2 = {mio::abm::InfectionState::Recovered};
    std::vector<mio::abm::LocationType> test_location_types2     = {mio::abm::LocationType::Home};
    auto testing_criteria2 = mio::abm::TestingCriteria({}, test_location_types2, test_infection_states2);
    testing_scheme.add_testing_criteria(testing_criteria2);

    auto loc_home = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto loc_work = mio::abm::Location(mio::abm::LocationType::Work, 0);
    auto person1  = make_test_person(loc_home, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedNoSymptoms);
    auto person2 =
        make_test_person(loc_home, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::Recovered);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly(5))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.5))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.5))
        .WillOnce(testing::Return(0.9));
    ASSERT_EQ(testing_scheme.run_scheme(person1, loc_home, start_date), false); // Person tests and tests positive
    ASSERT_EQ(testing_scheme.run_scheme(person2, loc_work, start_date), true); // Person tests and  tests negative
    ASSERT_EQ(testing_scheme.run_scheme(person1, loc_home, start_date),
              true); // Person is in quarantine and wants to go home -> can do so
    ASSERT_EQ(testing_scheme.run_scheme(person1, loc_work, start_date), true); // Person doesn't test

    testing_scheme.add_testing_criteria(testing_criteria1);
    testing_scheme.remove_testing_criteria(testing_criteria1);
    ASSERT_EQ(testing_scheme.run_scheme(person1, loc_home, start_date), true);
}
