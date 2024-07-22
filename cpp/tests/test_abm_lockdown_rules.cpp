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
#include "abm/lockdown_rules.h"
#include "abm/migration_rules.h"
#include "abm/person.h"
#include "abm_helpers.h"
#include "memilio/utils/random_number_generator.h"

TEST(TestLockdownRules, school_closure)
{
    auto rng       = mio::RandomNumberGenerator();
    auto t         = mio::abm::TimePoint(0);
    auto dt        = mio::abm::hours(1);
    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(6);

    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location school(mio::abm::LocationType::School, 1, num_age_groups);

    //setup rng mock so one person is home schooled and the other goes to school
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        .WillOnce(testing::Return(0.4))
        .WillOnce(testing::Return(0.4))
        .WillOnce(testing::Return(0.4))
        .WillOnce(testing::Return(0.4))
        .WillOnce(testing::Return(0.2))
        .WillOnce(testing::Return(0.2))
        .WillOnce(testing::Return(0.2))
        .WillOnce(testing::Return(0.2))
        .WillRepeatedly(testing::Return(1.0));

    auto p1 = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_5_to_14);
    p1.set_assigned_location(home.get_type(), home.get_id());
    p1.set_assigned_location(school.get_type(), school.get_id());
    auto p2 = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_5_to_14);
    p2.set_assigned_location(home.get_type(), home.get_id());
    p2.set_assigned_location(school.get_type(), school.get_id());
    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);
    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    params.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    params.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 or 35-59)
    params.get<mio::abm::AgeGroupGotoWork>()                     = false;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_15_to_34] = true;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_35_to_59] = true;
    mio::abm::set_school_closure(t, 0.7, params);

    auto p1_rng = mio::abm::PersonalRandomNumberGenerator(rng, p1);
    ASSERT_EQ(mio::abm::go_to_school(p1_rng, p1, t_morning, dt, params), mio::abm::LocationType::Home);
    auto p2_rng = mio::abm::PersonalRandomNumberGenerator(rng, p2);
    ASSERT_EQ(mio::abm::go_to_school(p2_rng, p2, t_morning, dt, params), mio::abm::LocationType::School);
}

TEST(TestLockdownRules, school_opening)
{
    auto rng       = mio::RandomNumberGenerator();
    auto t_closing = mio::abm::TimePoint(0);
    auto t_opening = mio::abm::TimePoint(0) + mio::abm::days(1);
    auto dt        = mio::abm::hours(1);
    auto t_morning = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(7);

    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location school(mio::abm::LocationType::School, 1, num_age_groups);
    //setup rng mock so the person is homeschooled in case of lockdown
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillRepeatedly(testing::Return(1.0));
    auto p = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_5_to_14);
    p.set_assigned_location(home.get_type(), home.get_id());
    p.set_assigned_location(school.get_type(), school.get_id());
    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);
    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    params.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    params.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 or 35-59)
    params.get<mio::abm::AgeGroupGotoWork>()                     = false;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_15_to_34] = true;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_35_to_59] = true;
    mio::abm::set_school_closure(t_closing, 1., params);
    mio::abm::set_school_closure(t_opening, 0., params);

    auto p_rng = mio::abm::PersonalRandomNumberGenerator(rng, p);
    ASSERT_EQ(mio::abm::go_to_school(p_rng, p, t_morning, dt, params), mio::abm::LocationType::School);
}

TEST(TestLockdownRules, home_office)
{
    auto rng       = mio::RandomNumberGenerator();
    auto t         = mio::abm::TimePoint(0);
    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto dt        = mio::abm::hours(1);

    mio::abm::Location home(mio::abm::LocationType::Home, 0);
    mio::abm::Location work(mio::abm::LocationType::Work, 1);
    mio::abm::Parameters params(num_age_groups);

    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    params.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    params.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 or 35-59)
    params.get<mio::abm::AgeGroupGotoWork>()                     = false;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_15_to_34] = true;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_35_to_59] = true;

    mio::abm::set_home_office(t, 0.4, params);

    //setup rng mock so one person goes to work and the other works at home
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(4))
        .WillOnce(testing::Return(0.5))
        .WillOnce(testing::Return(0.5))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.7))
        .WillRepeatedly(testing::Return(1.0));

    auto person1 = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_15_to_34);
    auto person2 = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_15_to_34);
    person1.set_assigned_location(home.get_type(), home.get_id());
    person1.set_assigned_location(work.get_type(), work.get_id());
    person2.set_assigned_location(home.get_type(), home.get_id());
    person2.set_assigned_location(work.get_type(), work.get_id());

    auto p1_rng = mio::abm::PersonalRandomNumberGenerator(rng, person1);
    ASSERT_EQ(mio::abm::go_to_work(p1_rng, person1, t_morning, dt, params), mio::abm::LocationType::Work);
    auto p2_rng = mio::abm::PersonalRandomNumberGenerator(rng, person2);
    ASSERT_EQ(mio::abm::go_to_work(p2_rng, person2, t_morning, dt, params), mio::abm::LocationType::Home);
}

TEST(TestLockdownRules, no_home_office)
{
    auto rng       = mio::RandomNumberGenerator();
    auto t_closing = mio::abm::TimePoint(0);
    auto t_opening = mio::abm::TimePoint(0) + mio::abm::days(1);
    auto dt        = mio::abm::hours(1);
    auto t_morning = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(8);

    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location work(mio::abm::LocationType::Work, 0, num_age_groups);

    //setup rng mock so the person works in home office
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.7))
        .WillRepeatedly(testing::Return(1.0));

    auto p = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_15_to_34);
    p.set_assigned_location(home.get_type(), home.get_id());
    p.set_assigned_location(work.get_type(), work.get_id());
    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);
    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    params.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    params.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 or 35-59)
    params.get<mio::abm::AgeGroupGotoWork>()                     = false;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_15_to_34] = true;
    params.get<mio::abm::AgeGroupGotoWork>()[age_group_35_to_59] = true;

    mio::abm::set_home_office(t_closing, 0.5, params);
    mio::abm::set_home_office(t_opening, 0., params);

    auto p_rng = mio::abm::PersonalRandomNumberGenerator(rng, p);
    ASSERT_EQ(mio::abm::go_to_work(p_rng, p, t_morning, dt, params), mio::abm::LocationType::Work);
}

TEST(TestLockdownRules, social_event_closure)
{
    auto rng       = mio::RandomNumberGenerator();
    auto t         = mio::abm::TimePoint(0);
    auto dt        = mio::abm::hours(1);
    auto t_evening = mio::abm::TimePoint(0) + mio::abm::hours(19);

    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location event(mio::abm::LocationType::SocialEvent, 1, num_age_groups);
    auto p = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_5_to_14);
    p.set_assigned_location(home.get_type(), home.get_id());
    p.set_assigned_location(event.get_type(), event.get_id());
    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);

    mio::abm::close_social_events(t, 1, params);

    auto p_rng = mio::abm::PersonalRandomNumberGenerator(rng, p);
    ASSERT_EQ(mio::abm::go_to_event(p_rng, p, t_evening, dt, params), mio::abm::LocationType::Home);
}

TEST(TestLockdownRules, social_events_opening)
{
    auto rng       = mio::RandomNumberGenerator();
    auto t_closing = mio::abm::TimePoint(0);
    auto t_opening = mio::abm::TimePoint(0) + mio::abm::days(1);
    auto dt        = mio::abm::hours(1);
    auto t_evening = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(19);

    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location event(mio::abm::LocationType::SocialEvent, 1, num_age_groups);
    auto p = mio::abm::Person(rng, home.get_type(), home.get_id(), age_group_5_to_14);
    p.set_assigned_location(event.get_type(), event.get_id());
    p.set_assigned_location(home.get_type(), home.get_id());
    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);

    mio::abm::close_social_events(t_closing, 1, params);
    mio::abm::close_social_events(t_opening, 0, params);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    auto p_rng = mio::abm::PersonalRandomNumberGenerator(rng, p);
    ASSERT_EQ(mio::abm::go_to_event(p_rng, p, t_evening, dt, params), mio::abm::LocationType::SocialEvent);
}
