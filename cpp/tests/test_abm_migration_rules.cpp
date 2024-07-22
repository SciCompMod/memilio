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
#include "abm/location_type.h"
#include "abm/migration_rules.h"
#include "abm/person.h"
#include "abm_helpers.h"
#include "memilio/utils/random_number_generator.h"

TEST(TestMigrationRules, random_migration)
{
    int t = 0, dt = 1;
    auto rng          = mio::RandomNumberGenerator();
    auto default_type = mio::abm::LocationType::Cemetery;
    auto person       = mio::abm::Person(rng, default_type, 0, age_group_15_to_34);
    auto p_rng        = mio::abm::PersonalRandomNumberGenerator(rng, person);
    auto params       = mio::abm::Parameters(num_age_groups);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>> mock_exp_dist;
    EXPECT_CALL(mock_exp_dist.get_mock(), invoke)
        .Times(testing::Exactly(2))
        // values for v in mio::abm::random_transition
        .WillOnce(testing::Return((t + dt) / 2.))
        .WillOnce(testing::Return(t + 2. * dt));

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_disc_dist;
    EXPECT_CALL(mock_disc_dist.get_mock(), invoke)
        .Times(testing::Exactly(1))
        // arbitrary value for random_idx in mio::abm::random_transition
        .WillOnce(testing::Return(2));

    const auto random_migration = [&]() {
        return mio::abm::random_migration(p_rng, person, mio::abm::TimePoint{t}, mio::abm::days(dt), params);
    };

    params.set<mio::abm::LockdownDate>(mio::abm::TimePoint{t + 2 * dt});

    const auto dest0 = random_migration();
    EXPECT_NE(dest0, default_type) << "should return a new location type (via random_transition)";

    const auto dest1 = random_migration();
    EXPECT_EQ(dest1, default_type) << "should return current location type (via random_transition)";

    params.set<mio::abm::LockdownDate>(mio::abm::TimePoint{t});

    const auto dest2 = random_migration();
    EXPECT_EQ(dest2, default_type) << "should return current location type";
}

TEST(TestMigrationRules, student_goes_to_school)
{
    auto rng = mio::RandomNumberGenerator();
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillRepeatedly(testing::Return(1.0));

    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    auto p_child = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_5_to_14);
    auto p_adult = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_15_to_34);

    auto t_morning              = mio::abm::TimePoint(0) + mio::abm::hours(7);
    auto t_weekend              = mio::abm::TimePoint(0) + mio::abm::days(5) + mio::abm::hours(7);
    auto dt                     = mio::abm::hours(1);
    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);
    auto child_rng              = mio::abm::PersonalRandomNumberGenerator(rng, p_child);
    auto adult_rng              = mio::abm::PersonalRandomNumberGenerator(rng, p_child);
    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    params.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    params.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 or 35-59)
    params.get<mio::abm::AgeGroupGotoWork>()                     = false;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_15_to_34] = true;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_35_to_59] = true;

    EXPECT_EQ(mio::abm::go_to_school(child_rng, p_child, t_morning, dt, params), mio::abm::LocationType::School);
    EXPECT_EQ(mio::abm::go_to_school(adult_rng, p_adult, t_morning, dt, params), mio::abm::LocationType::Home);
    EXPECT_EQ(mio::abm::go_to_school(child_rng, p_child, t_weekend, dt, params), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, students_go_to_school_in_different_times)
{
    auto rng = mio::RandomNumberGenerator();
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        //Mocking the random values will define at what time the student should go to school, i.e:
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
    auto p_child_goes_to_school_at_6   = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_5_to_14);
    auto rng_child_goes_to_school_at_6 = mio::abm::PersonalRandomNumberGenerator(rng, p_child_goes_to_school_at_6);
    auto p_child_goes_to_school_at_8   = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_5_to_14);
    auto rng_child_goes_to_school_at_8 = mio::abm::PersonalRandomNumberGenerator(rng, p_child_goes_to_school_at_8);

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

TEST(TestMigrationRules, students_go_to_school_in_different_times_with_smaller_time_steps)
{
    auto rng = mio::RandomNumberGenerator();
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
    auto p_child_goes_to_school_at_6    = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_5_to_14);
    auto rng_child_goes_to_school_at_6  = mio::abm::PersonalRandomNumberGenerator(rng, p_child_goes_to_school_at_6);
    auto p_child_goes_to_school_at_8_30 = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_5_to_14);
    auto rng_child_goes_to_school_at_8_30 =
        mio::abm::PersonalRandomNumberGenerator(rng, p_child_goes_to_school_at_8_30);

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

    EXPECT_EQ(
        mio::abm::go_to_school(rng_child_goes_to_school_at_6, p_child_goes_to_school_at_6, t_morning_6, dt, params),
        mio::abm::LocationType::School);
    EXPECT_EQ(
        mio::abm::go_to_school(rng_child_goes_to_school_at_6, p_child_goes_to_school_at_6, t_morning_8_30, dt, params),
        mio::abm::LocationType::Home);
    EXPECT_EQ(mio::abm::go_to_school(rng_child_goes_to_school_at_8_30, p_child_goes_to_school_at_8_30, t_morning_6, dt,
                                     params),
              mio::abm::LocationType::Home);
    EXPECT_EQ(mio::abm::go_to_school(rng_child_goes_to_school_at_8_30, p_child_goes_to_school_at_8_30, t_morning_8_30,
                                     dt, params),
              mio::abm::LocationType::School);
}

TEST(TestMigrationRules, school_return)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location school(mio::abm::LocationType::School, 0, num_age_groups);
    auto p_child   = mio::abm::Person(rng, school.get_type(), school.get_id(), age_group_5_to_14);
    auto rng_child = mio::abm::PersonalRandomNumberGenerator(rng, p_child);

    auto t  = mio::abm::TimePoint(0) + mio::abm::hours(15);
    auto dt = mio::abm::hours(1);

    EXPECT_EQ(mio::abm::go_to_school(rng_child, p_child, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, worker_goes_to_work)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
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

    auto p_retiree   = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_60_to_79);
    auto rng_retiree = mio::abm::PersonalRandomNumberGenerator(rng, p_retiree);
    auto p_adult     = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_15_to_34);
    auto rng_adult   = mio::abm::PersonalRandomNumberGenerator(rng, p_adult);

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

    EXPECT_EQ(mio::abm::go_to_work(rng_retiree, p_retiree, t_morning, dt, params), mio::abm::LocationType::Home);
    EXPECT_EQ(mio::abm::go_to_work(rng_adult, p_adult, t_morning, dt, params), mio::abm::LocationType::Home);
    EXPECT_EQ(mio::abm::go_to_work(rng_adult, p_adult, t_night, dt, params), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, worker_goes_to_work_with_non_dividable_timespan)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
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

    auto p_retiree   = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_60_to_79);
    auto rng_retiree = mio::abm::PersonalRandomNumberGenerator(rng, p_retiree);
    auto p_adult     = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_15_to_34);
    auto rng_adult   = mio::abm::PersonalRandomNumberGenerator(rng, p_adult);

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

    EXPECT_EQ(mio::abm::go_to_work(rng_retiree, p_retiree, t_morning, dt, params), mio::abm::LocationType::Home);
    EXPECT_EQ(mio::abm::go_to_work(rng_adult, p_adult, t_morning, dt, params), mio::abm::LocationType::Home);
    EXPECT_EQ(mio::abm::go_to_work(rng_adult, p_adult, t_night, dt, params), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, workers_go_to_work_in_different_times)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.9))
        .WillOnce(testing::Return(0.9))
        .WillOnce(testing::Return(0.9))
        .WillOnce(testing::Return(0.9))

        .WillRepeatedly(testing::Return(1.0));

    auto p_adult_goes_to_work_at_6   = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_15_to_34);
    auto rng_adult_goes_to_work_at_6 = mio::abm::PersonalRandomNumberGenerator(rng, p_adult_goes_to_work_at_6);
    auto p_adult_goes_to_work_at_8   = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_15_to_34);
    auto rng_adult_goes_to_work_at_8 = mio::abm::PersonalRandomNumberGenerator(rng, p_adult_goes_to_work_at_8);

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

    EXPECT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_6, p_adult_goes_to_work_at_6, t_morning_6, dt, params),
              mio::abm::LocationType::Work);
    EXPECT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_6, p_adult_goes_to_work_at_6, t_morning_8, dt, params),
              mio::abm::LocationType::Home);
    EXPECT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_6, p_adult_goes_to_work_at_6, t_night, dt, params),
              mio::abm::LocationType::Home);
    EXPECT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_8, p_adult_goes_to_work_at_8, t_morning_6, dt, params),
              mio::abm::LocationType::Home);
    EXPECT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_8, p_adult_goes_to_work_at_8, t_morning_8, dt, params),
              mio::abm::LocationType::Work);
    EXPECT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_8, p_adult_goes_to_work_at_8, t_night, dt, params),
              mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, work_return)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location work(mio::abm::LocationType::Work, 0, num_age_groups);
    auto p_adult   = mio::abm::Person(rng, work.get_type(), work.get_id(), age_group_35_to_59);
    auto rng_adult = mio::abm::PersonalRandomNumberGenerator(rng, p_adult);
    auto t         = mio::abm::TimePoint(0) + mio::abm::hours(17);
    auto dt        = mio::abm::hours(1);
    EXPECT_EQ(mio::abm::go_to_work(rng_adult, p_adult, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, quarantine)
{
    auto rng         = mio::RandomNumberGenerator();
    auto t           = mio::abm::TimePoint(12346);
    auto dt          = mio::abm::hours(1);
    auto test_params = mio::abm::TestParameters{1.0, 1.0};

    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location work(mio::abm::LocationType::Work, 0, num_age_groups);
    mio::abm::Location hospital(mio::abm::LocationType::Hospital, 0, num_age_groups);

    auto p_inf1   = make_test_person(work, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms, t);
    auto rng_inf1 = mio::abm::PersonalRandomNumberGenerator(rng, p_inf1);
    p_inf1.get_tested(rng_inf1, t, test_params);
    EXPECT_EQ(mio::abm::go_to_quarantine(rng_inf1, p_inf1, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Home); //detected infected person quarantines at home

    auto p_inf2   = make_test_person(work, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms, t);
    auto rng_inf2 = mio::abm::PersonalRandomNumberGenerator(rng, p_inf2);
    EXPECT_EQ(mio::abm::go_to_quarantine(rng_inf2, p_inf2, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Work); //undetected infected person does not quaratine

    auto p_inf3   = make_test_person(hospital, age_group_15_to_34, mio::abm::InfectionState::InfectedSevere, t);
    auto rng_inf3 = mio::abm::PersonalRandomNumberGenerator(rng, p_inf3);
    p_inf1.get_tested(rng_inf3, t, test_params);
    EXPECT_EQ(mio::abm::go_to_quarantine(rng_inf3, p_inf3, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Hospital); //detected infected person does not leave hospital to quarantine
}

TEST(TestMigrationRules, hospital)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    auto t       = mio::abm::TimePoint(12346);
    auto dt      = mio::abm::hours(1);
    auto p_inf   = make_test_person(home, age_group_15_to_34, mio::abm::InfectionState::InfectedSevere, t);
    auto rng_inf = mio::abm::PersonalRandomNumberGenerator(rng, p_inf);

    EXPECT_EQ(mio::abm::go_to_hospital(rng_inf, p_inf, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Hospital);

    auto p_car   = make_test_person(home, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms);
    auto rng_car = mio::abm::PersonalRandomNumberGenerator(rng, p_car);
    EXPECT_EQ(mio::abm::go_to_hospital(rng_car, p_car, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, go_shopping)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location hospital(mio::abm::LocationType::Hospital, 0, num_age_groups);
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);

    auto t_weekday = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(9);
    auto t_sunday  = mio::abm::TimePoint(0) + mio::abm::days(6) + mio::abm::hours(9);
    auto t_night   = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(1);
    auto dt        = mio::abm::hours(1);

    auto p_hosp   = make_test_person(hospital, age_group_0_to_4, mio::abm::InfectionState::InfectedSymptoms, t_weekday);
    auto rng_hosp = mio::abm::PersonalRandomNumberGenerator(rng, p_hosp);
    auto p_home   = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_60_to_79);
    auto rng_home = mio::abm::PersonalRandomNumberGenerator(rng, p_home);

    EXPECT_EQ(mio::abm::go_to_shop(rng_hosp, p_hosp, t_weekday, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Hospital);
    EXPECT_EQ(mio::abm::go_to_shop(rng_home, p_home, t_sunday, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Home);
    EXPECT_EQ(mio::abm::go_to_shop(rng_home, p_home, t_night, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Home);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    EXPECT_EQ(mio::abm::go_to_shop(rng_home, p_home, t_weekday, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::BasicsShop);
}

TEST(TestMigrationRules, shop_return)
{
    auto rng    = mio::RandomNumberGenerator();
    auto params = mio::abm::Parameters(num_age_groups);
    auto t      = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(9);
    auto dt     = mio::abm::hours(1);

    mio::abm::Location shop(mio::abm::LocationType::BasicsShop, 0, num_age_groups);
    auto p     = make_test_person(shop, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms, t);
    auto rng_p = mio::abm::PersonalRandomNumberGenerator(rng, p);

    p.add_time_at_location(dt);

    EXPECT_EQ(mio::abm::go_to_shop(rng_p, p, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, go_event)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location work(mio::abm::LocationType::Work, 0, num_age_groups);
    auto p_work   = mio::abm::Person(rng, work.get_type(), work.get_id(), age_group_35_to_59);
    auto rng_work = mio::abm::PersonalRandomNumberGenerator(rng, p_work);
    mio::abm::Location home(mio::abm::LocationType::Home, 1, num_age_groups);
    auto p_home   = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_60_to_79);
    auto rng_home = mio::abm::PersonalRandomNumberGenerator(rng, p_home);

    auto t_weekday  = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(20);
    auto t_saturday = mio::abm::TimePoint(0) + mio::abm::days(5) + mio::abm::hours(10);
    auto t_night    = mio::abm::TimePoint(0) + mio::abm::days(5) + mio::abm::hours(1);
    auto dt         = mio::abm::hours(1);

    EXPECT_EQ(mio::abm::go_to_event(rng_work, p_work, t_weekday, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Work);
    EXPECT_EQ(mio::abm::go_to_event(rng_home, p_home, t_night, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Home);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    EXPECT_EQ(mio::abm::go_to_event(rng_home, p_home, t_weekday, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::SocialEvent);

    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    EXPECT_EQ(mio::abm::go_to_event(rng_home, p_home, t_saturday, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::SocialEvent);
}

TEST(TestMigrationRules, event_return)
{
    auto rng    = mio::RandomNumberGenerator();
    auto params = mio::abm::Parameters(num_age_groups);
    auto t      = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(21);
    auto dt     = mio::abm::hours(3);

    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location social_event(mio::abm::LocationType::SocialEvent, 1, num_age_groups);
    auto p     = mio::abm::Person(rng, social_event.get_type(), social_event.get_id(), age_group_15_to_34);
    auto rng_p = mio::abm::PersonalRandomNumberGenerator(rng, p);

    p.add_time_at_location(dt);

    EXPECT_EQ(mio::abm::go_to_event(rng_p, p, t, dt, params), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, icu)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location hospital(mio::abm::LocationType::Hospital, 0, num_age_groups);
    auto t        = mio::abm::TimePoint(12346);
    auto dt       = mio::abm::hours(1);
    auto p_hosp   = make_test_person(hospital, age_group_15_to_34, mio::abm::InfectionState::InfectedCritical, t);
    auto rng_hosp = mio::abm::PersonalRandomNumberGenerator(rng, p_hosp);

    EXPECT_EQ(mio::abm::go_to_icu(rng_hosp, p_hosp, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::ICU);

    mio::abm::Location work(mio::abm::LocationType::Work, 1, num_age_groups);
    auto p_work   = make_test_person(work, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms, t);
    auto rng_work = mio::abm::PersonalRandomNumberGenerator(rng, p_work);
    EXPECT_EQ(mio::abm::go_to_icu(rng_work, p_work, t, dt, mio::abm::Parameters(num_age_groups)),
              mio::abm::LocationType::Work);
}

TEST(TestMigrationRules, recover)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location hospital(mio::abm::LocationType::Hospital, 0);
    auto t       = mio::abm::TimePoint(12346);
    auto dt      = mio::abm::hours(1);
    auto p_rec   = make_test_person(hospital, age_group_60_to_79, mio::abm::InfectionState::Recovered, t);
    auto rng_rec = mio::abm::PersonalRandomNumberGenerator(rng, p_rec);
    auto p_inf   = make_test_person(hospital, age_group_60_to_79, mio::abm::InfectionState::InfectedSevere, t);
    auto rng_inf = mio::abm::PersonalRandomNumberGenerator(rng, p_inf);
    EXPECT_EQ(mio::abm::return_home_when_recovered(rng_rec, p_rec, t, dt, {num_age_groups}),
              mio::abm::LocationType::Home);
    EXPECT_EQ(mio::abm::return_home_when_recovered(rng_inf, p_inf, t, dt, {num_age_groups}),
              mio::abm::LocationType::Hospital);
}

TEST(TestMigrationRules, dead)
{
    auto rng = mio::RandomNumberGenerator();

    mio::abm::Location icu(mio::abm::LocationType::ICU, 0);
    auto t      = mio::abm::TimePoint(12346);
    auto dt     = mio::abm::hours(1);
    auto p_dead = make_test_person(icu, age_group_60_to_79, mio::abm::InfectionState::Dead, t);
    auto p_rng  = mio::abm::PersonalRandomNumberGenerator(rng, p_dead);

    EXPECT_EQ(mio::abm::get_buried(p_rng, p_dead, t, dt, {num_age_groups}), mio::abm::LocationType::Cemetery);
}
