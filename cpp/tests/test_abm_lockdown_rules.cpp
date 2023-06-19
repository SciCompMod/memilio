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

TEST(TestLockdownRules, school_closure)
{
    auto t         = mio::abm::TimePoint(0);
    auto dt        = mio::abm::hours(1);
    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto school    = mio::abm::Location(mio::abm::LocationType::School, 0);

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

    auto p1 = mio::abm::Person(home, mio::abm::AgeGroup::Age5to14);
    p1.set_assigned_location(home);
    p1.set_assigned_location(school);
    auto p2 = mio::abm::Person(home, mio::abm::AgeGroup::Age5to14);
    p2.set_assigned_location(home);
    p2.set_assigned_location(school);
    mio::abm::MigrationParameters params;

    mio::abm::set_school_closure(t, 0.7, params);

    ASSERT_EQ(mio::abm::go_to_school(p1, t_morning, dt, params), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_school(p2, t_morning, dt, params), mio::abm::LocationType::School);
}

TEST(TestLockdownRules, school_opening)
{
    auto t_closing = mio::abm::TimePoint(0);
    auto t_opening = mio::abm::TimePoint(0) + mio::abm::days(1);
    auto dt        = mio::abm::hours(1);
    auto t_morning = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(7);
    auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto school    = mio::abm::Location(mio::abm::LocationType::School, 0);
    //setup rng mock so the person is homeschooled in case of lockdown
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillRepeatedly(testing::Return(1.0));
    auto p = mio::abm::Person(home, mio::abm::AgeGroup::Age5to14);
    p.set_assigned_location(home);
    p.set_assigned_location(school);
    mio::abm::MigrationParameters params;

    mio::abm::set_school_closure(t_closing, 1., params);
    mio::abm::set_school_closure(t_opening, 0., params);

    ASSERT_EQ(mio::abm::go_to_school(p, t_morning, dt, params), mio::abm::LocationType::School);
}

TEST(TestLockdownRules, home_office)
{
    auto t         = mio::abm::TimePoint(0);
    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto dt        = mio::abm::hours(1);
    auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto work      = mio::abm::Location(mio::abm::LocationType::Work, 0);
    mio::abm::MigrationParameters params;

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

    auto person1 = mio::abm::Person(home, mio::abm::AgeGroup::Age15to34);
    auto person2 = mio::abm::Person(home, mio::abm::AgeGroup::Age15to34);
    person1.set_assigned_location(home);
    person1.set_assigned_location(work);
    person2.set_assigned_location(home);
    person2.set_assigned_location(work);

    ASSERT_EQ(mio::abm::go_to_work(person1, t_morning, dt, params), mio::abm::LocationType::Work);
    ASSERT_EQ(mio::abm::go_to_work(person2, t_morning, dt, params), mio::abm::LocationType::Home);
}

TEST(TestLockdownRules, no_home_office)
{
    auto t_closing = mio::abm::TimePoint(0);
    auto t_opening = mio::abm::TimePoint(0) + mio::abm::days(1);
    auto dt        = mio::abm::hours(1);
    auto t_morning = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(8);
    auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto work      = mio::abm::Location(mio::abm::LocationType::Work, 0);

    //setup rng mock so the person works in home office
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.7))
        .WillRepeatedly(testing::Return(1.0));

    auto p = mio::abm::Person(home, mio::abm::AgeGroup::Age15to34);
    p.set_assigned_location(home);
    p.set_assigned_location(work);
    mio::abm::MigrationParameters params;

    mio::abm::set_home_office(t_closing, 0.5, params);
    mio::abm::set_home_office(t_opening, 0., params);

    ASSERT_EQ(mio::abm::go_to_work(p, t_morning, dt, params), mio::abm::LocationType::Work);
}

TEST(TestLockdownRules, social_event_closure)
{
    auto t         = mio::abm::TimePoint(0);
    auto dt        = mio::abm::hours(1);
    auto t_evening = mio::abm::TimePoint(0) + mio::abm::hours(19);
    auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto event     = mio::abm::Location(mio::abm::LocationType::SocialEvent, 0);
    auto p         = mio::abm::Person(home, mio::abm::AgeGroup::Age5to14);
    p.set_assigned_location(home);
    p.set_assigned_location(event);
    mio::abm::MigrationParameters params;

    mio::abm::close_social_events(t, 1, params);

    ASSERT_EQ(mio::abm::go_to_event(p, t_evening, dt, params), mio::abm::LocationType::Home);
}

TEST(TestLockdownRules, social_events_opening)
{
    auto t_closing = mio::abm::TimePoint(0);
    auto t_opening = mio::abm::TimePoint(0) + mio::abm::days(1);
    auto dt        = mio::abm::hours(1);
    auto t_evening = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(19);
    auto home      = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto event     = mio::abm::Location(mio::abm::LocationType::SocialEvent, 0);
    auto p         = mio::abm::Person(home, mio::abm::AgeGroup::Age5to14);
    p.set_assigned_location(event);
    p.set_assigned_location(home);
    mio::abm::MigrationParameters params;

    mio::abm::close_social_events(t_closing, 1, params);
    mio::abm::close_social_events(t_opening, 0, params);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(mio::abm::go_to_event(p, t_evening, dt, params), mio::abm::LocationType::SocialEvent);
}
