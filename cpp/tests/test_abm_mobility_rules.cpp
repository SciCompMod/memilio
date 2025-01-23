/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "abm/location_type.h"
#include "abm/mobility_rules.h"
#include "abm/person.h"
#include "abm_helpers.h"
#include "random_number_test.h"

using TestMobilityRules = RandomNumberTest;

/**
 * @brief Test random mobility transitions based on lockdown status and time.
 */
TEST_F(TestMobilityRules, random_mobility)
{
    int t = 0, dt = 1;
    auto default_type = mio::abm::LocationType::Cemetery;
    auto person       = mio::abm::Person(this->get_rng(), default_type, 0, 0, age_group_15_to_34);
    auto p_rng        = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), person);
    auto params       = mio::abm::Parameters(num_age_groups);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>> mock_exp_dist;
    EXPECT_CALL(mock_exp_dist.get_mock(), invoke)
        .Times(testing::Exactly(2))
        // Values for v in mio::abm::random_transition
        .WillOnce(testing::Return((t + dt) / 2.))
        .WillOnce(testing::Return(t + 2. * dt));

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_disc_dist;
    EXPECT_CALL(mock_disc_dist.get_mock(), invoke)
        .Times(testing::Exactly(1))
        // Arbitrary value for random_idx in mio::abm::random_transition
        .WillOnce(testing::Return(2));

    const auto random_mobility = [&]() {
        return mio::abm::random_mobility(p_rng, person, mio::abm::TimePoint{t}, mio::abm::days(dt), params);
    };

    params.set<mio::abm::LockdownDate>(mio::abm::TimePoint{t + 2 * dt});

    // Ensure random transitions occur
    const auto dest0 = random_mobility();
    EXPECT_NE(dest0, default_type) << "should return a new location type (via random_transition)";

    const auto dest1 = random_mobility();
    EXPECT_EQ(dest1, default_type) << "should return current location type (via random_transition)";

    // Once lockdown is active, ensure no mobility
    params.set<mio::abm::LockdownDate>(mio::abm::TimePoint{t});

    const auto dest2 = random_mobility();
    EXPECT_EQ(dest2, default_type) << "should return current location type";
}

/**
 * @brief Test school transitions.
 */
TEST_F(TestMobilityRules, student_goes_to_school)
{
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillRepeatedly(testing::Return(1.0));

    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    auto p_child =
        mio::abm::Person(this->get_rng(), home.get_type(), home.get_id(), home.get_model_id(), age_group_5_to_14);
    auto p_adult =
        mio::abm::Person(this->get_rng(), home.get_type(), home.get_id(), home.get_model_id(), age_group_15_to_34);

    auto t_morning              = mio::abm::TimePoint(0) + mio::abm::hours(7);
    auto t_weekend              = mio::abm::TimePoint(0) + mio::abm::days(5) + mio::abm::hours(7);
    auto dt                     = mio::abm::hours(1);
    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);
    auto child_rng              = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_child);
    auto adult_rng              = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_child);
    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    params.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    params.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 or 35-59)
    params.get<mio::abm::AgeGroupGotoWork>()                     = false;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_15_to_34] = true;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_35_to_59] = true;

    // Test that child goes to school
    EXPECT_EQ(mio::abm::go_to_school(child_rng, p_child, t_morning, dt, params), mio::abm::LocationType::School);
    // Test that adult does not go to school
    EXPECT_EQ(mio::abm::go_to_school(adult_rng, p_adult, t_morning, dt, params), mio::abm::LocationType::Home);
    // Test that child goes back home after school
    EXPECT_EQ(mio::abm::go_to_school(child_rng, p_child, t_weekend, dt, params), mio::abm::LocationType::Home);
}

/**
 * @brief Test children go to school at different times.
 */
TEST_F(TestMobilityRules, students_go_to_school_in_different_times)
{
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        // Mocking the random values will define at what time the student should go to school, i.e:
        // random is in [0,1/3] -> goes to school at 6
        // random is in [1/3,2/3] -> goes to school at 7
        // random is in [2/3,1.] -> goes to school at 8
        .WillOnce(testing::Return(0.0))
        .WillOnce(testing::Return(0.0))
        .WillOnce(testing::Return(0.0))
        .WillOnce(testing::Return(0.0))
        .WillOnce(testing::Return(0.8))
        .WillOnce(testing::Return(0.8))
        .WillOnce(testing::Return(0.8))
        .WillOnce(testing::Return(0.8))
        .WillRepeatedly(testing::Return(1.0));

    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    auto p_child_goes_to_school_at_6 =
        mio::abm::Person(this->get_rng(), home.get_type(), home.get_id(), home.get_model_id(), age_group_5_to_14);
    auto rng_child_goes_to_school_at_6 =
        mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_child_goes_to_school_at_6);
    auto p_child_goes_to_school_at_8 =
        mio::abm::Person(this->get_rng(), home.get_type(), home.get_id(), home.get_model_id(), age_group_5_to_14);
    auto rng_child_goes_to_school_at_8 =
        mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_child_goes_to_school_at_8);

    auto t_morning_6 = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto t_morning_8 = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto dt          = mio::abm::hours(1);

    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);
    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    params.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    params.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 or 35-59)
    params.get<mio::abm::AgeGroupGotoWork>()                     = false;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_15_to_34] = true;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_35_to_59] = true;

    // Mock randomness ensures children leave home at various school start times.
    EXPECT_EQ(
        mio::abm::go_to_school(rng_child_goes_to_school_at_6, p_child_goes_to_school_at_6, t_morning_6, dt, params),
        mio::abm::LocationType::School);
    EXPECT_EQ(
        mio::abm::go_to_school(rng_child_goes_to_school_at_6, p_child_goes_to_school_at_6, t_morning_8, dt, params),
        mio::abm::LocationType::Home);
    EXPECT_EQ(
        mio::abm::go_to_school(rng_child_goes_to_school_at_8, p_child_goes_to_school_at_8, t_morning_6, dt, params),
        mio::abm::LocationType::Home);
    EXPECT_EQ(
        mio::abm::go_to_school(rng_child_goes_to_school_at_8, p_child_goes_to_school_at_8, t_morning_8, dt, params),
        mio::abm::LocationType::School);
}

/**
 * @brief Test students going to school at different times with smaller time steps.
 * This test checks that students leave for school at various times in the morning,
 * with smaller time steps, based on different random values.
 */
TEST_F(TestMobilityRules, students_go_to_school_in_different_times_with_smaller_time_steps)
{
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        //Mocking the random values will define at what time the student should go to school, i.e:
        // random is in [0,1/6] -> goes to school at 6
        // random is in [1/6,2/6] -> goes to school at 6:30
        // random is in [2/6,3/6] -> goes to school at 7:00
        // random is in [3/6,4/6] -> goes to school at 7:30
        // random is in [4/6,5/6] -> goes to school at 8:00
        // random is in [5/6,6/6] -> goes to school at 8:30
        .WillOnce(testing::Return(0.0))
        .WillOnce(testing::Return(0.0))
        .WillOnce(testing::Return(0.0))
        .WillOnce(testing::Return(0.0))
        .WillOnce(testing::Return(0.9))
        .WillOnce(testing::Return(0.9))
        .WillOnce(testing::Return(0.9))
        .WillOnce(testing::Return(0.9))
        .WillRepeatedly(testing::Return(1.0));

    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    // First student goes to school at 6:00 AM
    auto p_child_goes_to_school_at_6 =
        mio::abm::Person(this->get_rng(), home.get_type(), home.get_id(), home.get_model_id(), age_group_5_to_14);
    auto rng_child_goes_to_school_at_6 =
        mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_child_goes_to_school_at_6);

    // Second student goes to school at 8:30 AM
    auto p_child_goes_to_school_at_8_30 =
        mio::abm::Person(this->get_rng(), home.get_type(), home.get_id(), home.get_model_id(), age_group_5_to_14);
    auto rng_child_goes_to_school_at_8_30 =
        mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_child_goes_to_school_at_8_30);

    // Time points for the morning scenarios
    auto t_morning_6            = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto t_morning_8_30         = mio::abm::TimePoint(0) + mio::abm::hours(8) + mio::abm::seconds(1800);
    auto dt                     = mio::abm::seconds(1800);
    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);
    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    params.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    params.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 or 35-59)
    params.get<mio::abm::AgeGroupGotoWork>()                     = false;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_15_to_34] = true;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_35_to_59] = true;

    // Check that the first student goes to school at 6:00 AM
    EXPECT_EQ(
        mio::abm::go_to_school(rng_child_goes_to_school_at_6, p_child_goes_to_school_at_6, t_morning_6, dt, params),
        mio::abm::LocationType::School);
    EXPECT_EQ(
        mio::abm::go_to_school(rng_child_goes_to_school_at_6, p_child_goes_to_school_at_6, t_morning_8_30, dt, params),
        mio::abm::LocationType::Home);
    // Check that the second student goes to school at 8:30 AM
    EXPECT_EQ(mio::abm::go_to_school(rng_child_goes_to_school_at_8_30, p_child_goes_to_school_at_8_30, t_morning_6, dt,
                                     params),
              mio::abm::LocationType::Home);
    EXPECT_EQ(mio::abm::go_to_school(rng_child_goes_to_school_at_8_30, p_child_goes_to_school_at_8_30, t_morning_8_30,
                                     dt, params),
              mio::abm::LocationType::School);
}
/**
 * @brief Test return home from school.
 */
TEST_F(TestMobilityRules, school_return)
{
    mio::abm::Location school(mio::abm::LocationType::School, 0, num_age_groups);
    auto p_child =
        mio::abm::Person(this->get_rng(), school.get_type(), school.get_id(), school.get_model_id(), age_group_5_to_14);
    auto rng_child = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_child);

    // Simulate a time point after school hours
    auto t  = mio::abm::TimePoint(0) + mio::abm::hours(15);
    auto dt = mio::abm::hours(1);

    // Ensure that the child returns home after school is over
    EXPECT_EQ(mio::abm::go_to_school(rng_child, p_child, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Home);
}

/**
 * @brief Test work transitions.
 */
TEST_F(TestMobilityRules, worker_goes_to_work)
{
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    // Mock the uniform distribution to control the randomness for workers' decisions.
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillRepeatedly(testing::Return(1.0));

    auto p_retiree =
        mio::abm::Person(this->get_rng(), home.get_type(), home.get_id(), home.get_model_id(), age_group_60_to_79);
    auto rng_retiree = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_retiree);
    auto p_adult =
        mio::abm::Person(this->get_rng(), home.get_type(), home.get_id(), home.get_model_id(), age_group_15_to_34);
    auto rng_adult = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_adult);

    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto t_night   = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(4);
    auto dt        = mio::abm::hours(1);

    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);
    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    params.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    params.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 or 35-59)
    params.get<mio::abm::AgeGroupGotoWork>()                     = false;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_15_to_34] = true;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_35_to_59] = true;

    // Check that the retiree (age group 60-79) should stay home and not go to work.
    EXPECT_EQ(mio::abm::go_to_work(rng_retiree, p_retiree, t_morning, dt, params), mio::abm::LocationType::Home);
    // Check that the adult (age group 15-34) should go to work at 8:00 AM.
    EXPECT_EQ(mio::abm::go_to_work(rng_adult, p_adult, t_morning, dt, params), mio::abm::LocationType::Home);
    // Check that during the night (4:00 AM), the adult should stay home and not go to work.
    EXPECT_EQ(mio::abm::go_to_work(rng_adult, p_adult, t_night, dt, params), mio::abm::LocationType::Home);
}

/**
 * @brief This test checks if workers go to work based on their schedule when the time step is not evenly divisible.
 */
TEST_F(TestMobilityRules, worker_goes_to_work_with_non_dividable_timespan)
{
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);

    // Mock the uniform distribution to control the randomness for workers' decisions.
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillRepeatedly(testing::Return(1.0));

    // Set up two people: one retiree and one working adult.
    auto p_retiree =
        mio::abm::Person(this->get_rng(), home.get_type(), home.get_id(), home.get_model_id(), age_group_60_to_79);
    auto rng_retiree = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_retiree);
    auto p_adult =
        mio::abm::Person(this->get_rng(), home.get_type(), home.get_id(), home.get_model_id(), age_group_15_to_34);
    auto rng_adult = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_adult);

    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto t_night   = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(4);
    auto dt        = mio::abm::minutes(53);

    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);
    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    params.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    params.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 or 35-59)
    params.get<mio::abm::AgeGroupGotoWork>()                     = false;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_15_to_34] = true;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_35_to_59] = true;

    // Check that the retiree (age group 60-79) should stay home and not go to work even with non-dividable time step.
    EXPECT_EQ(mio::abm::go_to_work(rng_retiree, p_retiree, t_morning, dt, params), mio::abm::LocationType::Home);
    // Check that the adult (age group 15-34) should still go to work at 8:00 AM even with the non-dividable time step.
    EXPECT_EQ(mio::abm::go_to_work(rng_adult, p_adult, t_morning, dt, params), mio::abm::LocationType::Home);
    // Check that during the night (4:00 AM), the adult should stay home and not go to work even with the non-dividable time step.
    EXPECT_EQ(mio::abm::go_to_work(rng_adult, p_adult, t_night, dt, params), mio::abm::LocationType::Home);
}

/**
 * @brief Test workers going to work at different times and follow their schedules accordingly.
 */
TEST_F(TestMobilityRules, workers_go_to_work_in_different_times)
{
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);

    // Mock the uniform distribution to control the randomness of the workers' go-to-work times.
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        // Mocked random values define when workers go to work
        // Values:
        // [0, 1/3] -> Worker goes to work at 6 AM
        // [1/3, 2/3] -> Worker goes to work at 7 AM
        // [2/3, 1] -> Worker goes to work at 8 AM
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.9))
        .WillOnce(testing::Return(0.9))
        .WillOnce(testing::Return(0.9))
        .WillOnce(testing::Return(0.9))
        .WillRepeatedly(testing::Return(1.0));

    // Create two workers: one goes to work at 6 AM and the other at 8 AM.
    auto p_adult_goes_to_work_at_6 =
        mio::abm::Person(this->get_rng(), home.get_type(), home.get_id(), home.get_model_id(), age_group_15_to_34);
    auto rng_adult_goes_to_work_at_6 =
        mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_adult_goes_to_work_at_6);
    auto p_adult_goes_to_work_at_8 =
        mio::abm::Person(this->get_rng(), home.get_type(), home.get_id(), home.get_model_id(), age_group_15_to_34);
    auto rng_adult_goes_to_work_at_8 =
        mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_adult_goes_to_work_at_8);

    auto t_morning_6            = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto t_morning_8            = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto t_night                = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(4);
    auto dt                     = mio::abm::hours(1);
    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);
    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    params.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    params.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 or 35-59)
    params.get<mio::abm::AgeGroupGotoWork>()                     = false;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_15_to_34] = true;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_35_to_59] = true;

    // Check that the worker going to work at 6 AM goes to work
    EXPECT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_6, p_adult_goes_to_work_at_6, t_morning_6, dt, params),
              mio::abm::LocationType::Work);
    // Check that the worker going to work at 6 AM stays home at 8 AM (since they are already at work)
    EXPECT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_6, p_adult_goes_to_work_at_6, t_morning_8, dt, params),
              mio::abm::LocationType::Home);
    // Check that the worker going to work at 6 AM returns home at night
    EXPECT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_6, p_adult_goes_to_work_at_6, t_night, dt, params),
              mio::abm::LocationType::Home);
    // Check that the worker going to work at 8 AM stays home at 6 AM
    EXPECT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_8, p_adult_goes_to_work_at_8, t_morning_6, dt, params),
              mio::abm::LocationType::Home);
    // Check that the worker going to work at 8 AM goes to work at 8 AM
    EXPECT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_8, p_adult_goes_to_work_at_8, t_morning_8, dt, params),
              mio::abm::LocationType::Work);
    // Check that the worker going to work at 8 AM returns home at night
    EXPECT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_8, p_adult_goes_to_work_at_8, t_night, dt, params),
              mio::abm::LocationType::Home);
}

/**
 * @brief Test that workers return home after finishing their workday.
 */
TEST_F(TestMobilityRules, work_return)
{
    mio::abm::Location work(mio::abm::LocationType::Work, 0, num_age_groups);
    // Set up a random number generator and a worker who is currently at work
    auto p_adult =
        mio::abm::Person(this->get_rng(), work.get_type(), work.get_id(), work.get_model_id(), age_group_35_to_59);
    auto rng_adult = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_adult);
    // Set the time to 5 PM (17:00) when the worker should return home
    auto t  = mio::abm::TimePoint(0) + mio::abm::hours(17);
    auto dt = mio::abm::hours(1);
    // Test that the worker, who is currently at work, goes home after 5 PM
    EXPECT_EQ(mio::abm::go_to_work(rng_adult, p_adult, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Home);
}

/**
 * @brief Test quarantine behavior.
 */
TEST_F(TestMobilityRules, quarantine)
{
    auto t           = mio::abm::TimePoint(12346);
    auto dt          = mio::abm::hours(1);
    auto test_params = mio::abm::TestParameters{1.0, 1.0, mio::abm::minutes(30), mio::abm::TestType::Generic};
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location work(mio::abm::LocationType::Work, 0, num_age_groups);
    mio::abm::Location hospital(mio::abm::LocationType::Hospital, 0, num_age_groups);

    auto p_inf1 =
        make_test_person(this->get_rng(), work, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms, t);
    auto rng_inf1 = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_inf1);
    p_inf1.get_tested(rng_inf1, t, test_params);
    // Check detected infected person quarantines at home
    EXPECT_EQ(mio::abm::go_to_quarantine(rng_inf1, p_inf1, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Home);

    auto p_inf2 =
        make_test_person(this->get_rng(), work, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms, t);
    auto rng_inf2 = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_inf2);
    // Check that undetected infected person does not quaratine
    EXPECT_EQ(mio::abm::go_to_quarantine(rng_inf2, p_inf2, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Work);

    auto p_inf3 =
        make_test_person(this->get_rng(), hospital, age_group_15_to_34, mio::abm::InfectionState::InfectedSevere, t);
    auto rng_inf3 = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_inf3);
    p_inf1.get_tested(rng_inf3, t, test_params);
    // Check that detected infected person does not leave hospital to quarantine
    EXPECT_EQ(mio::abm::go_to_quarantine(rng_inf3, p_inf3, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Hospital);
}

/**
 * @brief Test hospital transition.
 */
TEST_F(TestMobilityRules, hospital)
{
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    auto t  = mio::abm::TimePoint(12346);
    auto dt = mio::abm::hours(1);
    auto p_inf =
        make_test_person(this->get_rng(), home, age_group_15_to_34, mio::abm::InfectionState::InfectedSevere, t);
    auto rng_inf = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_inf);

    // Ensure person goes to the hospital when severely infected
    EXPECT_EQ(mio::abm::go_to_hospital(rng_inf, p_inf, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Hospital);

    auto p_car =
        make_test_person(this->get_rng(), home, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms);
    auto rng_car = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_car);
    // Ensure person has infection symptoms still stay at home
    EXPECT_EQ(mio::abm::go_to_hospital(rng_car, p_car, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Home);
}

/**
 * @brief Test the mobility of people when going shopping.
 */
TEST_F(TestMobilityRules, go_shopping)
{
    mio::abm::Location hospital(mio::abm::LocationType::Hospital, 0, num_age_groups);
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);

    auto t_weekday = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(9);
    auto t_sunday  = mio::abm::TimePoint(0) + mio::abm::days(6) + mio::abm::hours(9);
    auto t_night   = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(1);
    auto dt        = mio::abm::hours(1);

    // Create an infected child in the hospital
    auto p_hosp   = make_test_person(this->get_rng(), hospital, age_group_0_to_4,
                                     mio::abm::InfectionState::InfectedSymptoms, t_weekday);
    auto rng_hosp = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_hosp);
    // Create a healthy elderly person at home
    auto p_home =
        mio::abm::Person(this->get_rng(), home.get_type(), home.get_id(), home.get_model_id(), age_group_60_to_79);
    auto rng_home = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_home);

    // Check that an infected person stays in the hospital and doesn't go shopping
    EXPECT_EQ(mio::abm::go_to_shop(rng_hosp, p_hosp, t_weekday, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Hospital);
    // Check that a person at home doesn't go shopping on a Sunday
    EXPECT_EQ(mio::abm::go_to_shop(rng_home, p_home, t_sunday, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Home);
    // Check that a person at home doesn't go shopping on a Sunday
    EXPECT_EQ(mio::abm::go_to_shop(rng_home, p_home, t_night, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Home);

    // Mocking the random distribution to simulate a person going shopping on a weekday
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    // Test that a person goes to a basic shop on a weekday at 9 AM
    EXPECT_EQ(mio::abm::go_to_shop(rng_home, p_home, t_weekday, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::BasicsShop);
}

/**
 * @brief Test the behavior of a person returning home after spending time at a shop.
 */
TEST_F(TestMobilityRules, shop_return)
{
    auto params = mio::abm::Parameters(num_age_groups);
    auto t      = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(9);
    auto dt     = mio::abm::hours(1);

    // Create a person at a basic shop who is asymptomatically infected
    mio::abm::Location shop(mio::abm::LocationType::BasicsShop, 0, num_age_groups);
    auto p =
        make_test_person(this->get_rng(), shop, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms, t);
    auto rng_p = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p);
    // Simulate the person spending 1 hour at the shop
    p.add_time_at_location(dt);

    // After spending sufficient time at the shop, the person should return home
    EXPECT_EQ(mio::abm::go_to_shop(rng_p, p, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Home);
}

/**
 * @brief Test the behavior of people going to a social event based on their current location, time, and random factors.
 */
TEST_F(TestMobilityRules, go_event)
{
    // Initialize two people, one at work and one at home
    mio::abm::Location work(mio::abm::LocationType::Work, 0, num_age_groups);
    auto p_work =
        mio::abm::Person(this->get_rng(), work.get_type(), work.get_id(), work.get_model_id(), age_group_35_to_59);
    auto rng_work = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_work);
    mio::abm::Location home(mio::abm::LocationType::Home, 1, num_age_groups);
    auto p_home =
        mio::abm::Person(this->get_rng(), home.get_type(), home.get_id(), home.get_model_id(), age_group_60_to_79);
    auto rng_home = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_home);

    auto t_weekday  = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(20);
    auto t_saturday = mio::abm::TimePoint(0) + mio::abm::days(5) + mio::abm::hours(10);
    auto t_night    = mio::abm::TimePoint(0) + mio::abm::days(5) + mio::abm::hours(1);
    auto dt         = mio::abm::hours(1);

    // Check that person at work should not go to an event during a weekday evening, so they stay at work
    EXPECT_EQ(mio::abm::go_to_event(rng_work, p_work, t_weekday, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Work);
    // Check that person at home during nighttime should not go to an event, so they stay at home
    EXPECT_EQ(mio::abm::go_to_event(rng_home, p_home, t_night, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Home);
    // Mock the random distribution to simulate event participation based on random draw
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    // Check that person at home goes to a social event on a weekday evening based on a low random value
    EXPECT_EQ(mio::abm::go_to_event(rng_home, p_home, t_weekday, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::SocialEvent);

    // Check that person at home goes to a social event on a Saturday morning based on a low random value
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    EXPECT_EQ(mio::abm::go_to_event(rng_home, p_home, t_saturday, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::SocialEvent);
}

/**
 * @brief Test the behavior of a person returning home after spending time at a social event.
 */
TEST_F(TestMobilityRules, event_return)
{
    auto params = mio::abm::Parameters(num_age_groups);
    auto t      = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(21);
    auto dt     = mio::abm::hours(3);

    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location social_event(mio::abm::LocationType::SocialEvent, 1, num_age_groups);
    // Initialize the person at the social event location
    auto p     = mio::abm::Person(this->get_rng(), social_event.get_type(), social_event.get_id(),
                                  social_event.get_model_id(), age_group_15_to_34);
    auto rng_p = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p);
    // Simulate the person spending 3 hours at the social event
    p.add_time_at_location(dt);
    // After spending the time at the social event, the person should return home
    EXPECT_EQ(mio::abm::go_to_event(rng_p, p, t, dt, params), mio::abm::LocationType::Home);
}

/**
 * @brief Test ICU transition.
 */
TEST_F(TestMobilityRules, icu)
{
    mio::abm::Location hospital(mio::abm::LocationType::Hospital, 0, num_age_groups);
    auto t  = mio::abm::TimePoint(12346);
    auto dt = mio::abm::hours(1);
    auto p_hosp =
        make_test_person(this->get_rng(), hospital, age_group_15_to_34, mio::abm::InfectionState::InfectedCritical, t);
    auto rng_hosp = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_hosp);

    // Ensure critically infected person goes to the ICU
    EXPECT_EQ(mio::abm::go_to_icu(rng_hosp, p_hosp, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::ICU);

    mio::abm::Location work(mio::abm::LocationType::Work, 1, num_age_groups);
    auto p_work =
        make_test_person(this->get_rng(), work, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms, t);
    auto rng_work = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_work);
    // Ensure infected with symptions person can still go to work
    EXPECT_EQ(mio::abm::go_to_icu(rng_work, p_work, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Work);
}

/**
 * @brief Test recovery transition.
 */
TEST_F(TestMobilityRules, recover)
{
    mio::abm::Location hospital(mio::abm::LocationType::Hospital, 0);
    auto t  = mio::abm::TimePoint(12346);
    auto dt = mio::abm::hours(1);
    auto p_rec =
        make_test_person(this->get_rng(), hospital, age_group_60_to_79, mio::abm::InfectionState::Recovered, t);
    auto rng_rec = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_rec);
    auto p_inf =
        make_test_person(this->get_rng(), hospital, age_group_60_to_79, mio::abm::InfectionState::InfectedSevere, t);
    auto rng_inf = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_inf);
    // Ensure recovered person returns home and infected severe person stay in hospital
    EXPECT_EQ(mio::abm::return_home_when_recovered(rng_rec, p_rec, t, dt, {num_age_groups}),
              mio::abm::LocationType::Home);
    EXPECT_EQ(mio::abm::return_home_when_recovered(rng_inf, p_inf, t, dt, {num_age_groups}),
              mio::abm::LocationType::Hospital);
}

/**
 * @brief Test the ability to go to cemetery for deceased individuals.
 */
TEST_F(TestMobilityRules, dead)
{
    mio::abm::Location icu(mio::abm::LocationType::ICU, 0);
    auto t      = mio::abm::TimePoint(12346);
    auto dt     = mio::abm::hours(1);
    auto p_dead = make_test_person(this->get_rng(), icu, age_group_60_to_79, mio::abm::InfectionState::Dead, t);
    auto p_rng  = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p_dead);

    EXPECT_EQ(mio::abm::get_buried(p_rng, p_dead, t, dt, {num_age_groups}), mio::abm::LocationType::Cemetery);
}
