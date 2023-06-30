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

TEST(TestMigrationRules, student_goes_to_school)
{
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(8))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillRepeatedly(testing::Return(1.0));

    auto home    = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto p_child = mio::abm::Person(home, mio::abm::AgeGroup::Age5to14);
    auto p_adult = mio::abm::Person(home, mio::abm::AgeGroup::Age15to34);

    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(7);
    auto t_weekend = mio::abm::TimePoint(0) + mio::abm::days(5) + mio::abm::hours(7);
    auto dt        = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_school(p_child, t_morning, dt, {}), mio::abm::LocationType::School);
    ASSERT_EQ(mio::abm::go_to_school(p_adult, t_morning, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_school(p_child, t_weekend, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, students_go_to_school_in_different_times)
{
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

    auto home                        = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto p_child_goes_to_school_at_6 = mio::abm::Person(home, mio::abm::AgeGroup::Age5to14);
    auto p_child_goes_to_school_at_8 = mio::abm::Person(home, mio::abm::AgeGroup::Age5to14);

    auto t_morning_6 = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto t_morning_8 = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto dt          = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_school(p_child_goes_to_school_at_6, t_morning_6, dt, {}), mio::abm::LocationType::School);
    ASSERT_EQ(mio::abm::go_to_school(p_child_goes_to_school_at_6, t_morning_8, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_school(p_child_goes_to_school_at_8, t_morning_6, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_school(p_child_goes_to_school_at_8, t_morning_8, dt, {}), mio::abm::LocationType::School);
}

TEST(TestMigrationRules, students_go_to_school_in_different_times_with_smaller_time_steps)
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

    auto home                           = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto p_child_goes_to_school_at_6    = mio::abm::Person(home, mio::abm::AgeGroup::Age5to14);
    auto p_child_goes_to_school_at_8_30 = mio::abm::Person(home, mio::abm::AgeGroup::Age5to14);

    auto t_morning_6    = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto t_morning_8_30 = mio::abm::TimePoint(0) + mio::abm::hours(8) + mio::abm::seconds(1800);
    auto dt             = mio::abm::seconds(1800);

    ASSERT_EQ(mio::abm::go_to_school(p_child_goes_to_school_at_6, t_morning_6, dt, {}), mio::abm::LocationType::School);
    ASSERT_EQ(mio::abm::go_to_school(p_child_goes_to_school_at_6, t_morning_8_30, dt, {}),
              mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_school(p_child_goes_to_school_at_8_30, t_morning_6, dt, {}),
              mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_school(p_child_goes_to_school_at_8_30, t_morning_8_30, dt, {}),
              mio::abm::LocationType::School);
}

TEST(TestMigrationRules, school_return)
{
    auto school  = mio::abm::Location(mio::abm::LocationType::School, 0);
    auto p_child = mio::abm::Person(school, mio::abm::AgeGroup::Age5to14);

    auto t  = mio::abm::TimePoint(0) + mio::abm::hours(15);
    auto dt = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_school(p_child, t, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, worker_goes_to_work)
{
    auto home = mio::abm::Location(mio::abm::LocationType::Home, 0);
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

    auto p_retiree = mio::abm::Person(home, mio::abm::AgeGroup::Age60to79);
    auto p_adult   = mio::abm::Person(home, mio::abm::AgeGroup::Age15to34);

    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto t_night   = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(4);
    auto dt        = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_work(p_retiree, t_morning, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(p_adult, t_morning, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(p_adult, t_night, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, worker_goes_to_work_with_non_dividable_timespan)
{
    auto home = mio::abm::Location(mio::abm::LocationType::Home, 0);
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

    auto p_retiree = mio::abm::Person(home, mio::abm::AgeGroup::Age60to79);
    auto p_adult   = mio::abm::Person(home, mio::abm::AgeGroup::Age15to34);

    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto t_night   = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(4);
    auto dt        = mio::abm::minutes(53);

    ASSERT_EQ(mio::abm::go_to_work(p_retiree, t_morning, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(p_adult, t_morning, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(p_adult, t_night, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, workers_go_to_work_in_different_times)
{
    auto home = mio::abm::Location(mio::abm::LocationType::Home, 0);
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

    auto p_adult_goes_to_work_at_6 = mio::abm::Person(home, mio::abm::AgeGroup::Age15to34);
    auto p_adult_goes_to_work_at_8 = mio::abm::Person(home, mio::abm::AgeGroup::Age15to34);

    auto t_morning_6 = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto t_morning_8 = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto t_night     = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(4);
    auto dt          = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_work(p_adult_goes_to_work_at_6, t_morning_6, dt, {}), mio::abm::LocationType::Work);
    ASSERT_EQ(mio::abm::go_to_work(p_adult_goes_to_work_at_6, t_morning_8, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(p_adult_goes_to_work_at_6, t_night, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(p_adult_goes_to_work_at_8, t_morning_6, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(p_adult_goes_to_work_at_8, t_morning_8, dt, {}), mio::abm::LocationType::Work);
    ASSERT_EQ(mio::abm::go_to_work(p_adult_goes_to_work_at_8, t_night, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, work_return)
{
    auto work    = mio::abm::Location(mio::abm::LocationType::Work, 0);
    auto p_adult = mio::abm::Person(work, mio::abm::AgeGroup::Age35to59);
    auto t       = mio::abm::TimePoint(0) + mio::abm::hours(17);
    auto dt      = mio::abm::hours(1);
    ASSERT_EQ(mio::abm::go_to_work(p_adult, t, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, quarantine)
{
    auto t  = mio::abm::TimePoint(12346);
    auto dt = mio::abm::hours(1);

    auto home     = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto work     = mio::abm::Location(mio::abm::LocationType::Work, 0);
    auto hospital = mio::abm::Location(mio::abm::LocationType::Hospital, 0);

    auto p_inf1 = make_test_person(work, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSymptoms, t);
    p_inf1.detect_infection(t);
    ASSERT_EQ(mio::abm::go_to_quarantine(p_inf1, t, dt, {}),
              mio::abm::LocationType::Home); //detected infected person quarantines at home

    auto p_inf2 = make_test_person(work, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSymptoms, t);
    ASSERT_EQ(mio::abm::go_to_quarantine(p_inf2, t, dt, {}),
              mio::abm::LocationType::Work); //undetected infected person does not quaratine

    auto p_inf3 =
        make_test_person(hospital, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSevere, t);
    p_inf3.detect_infection(t);
    ASSERT_EQ(mio::abm::go_to_quarantine(p_inf3, t, dt, {}),
              mio::abm::LocationType::Hospital); //detected infected person does not leave hospital to quarantine
}

TEST(TestMigrationRules, hospital)
{
    auto home  = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto t     = mio::abm::TimePoint(12346);
    auto dt    = mio::abm::hours(1);
    auto p_inf = make_test_person(home, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSevere, t);

    ASSERT_EQ(mio::abm::go_to_hospital(p_inf, t, dt, {}), mio::abm::LocationType::Hospital);

    auto p_car = make_test_person(home, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSymptoms);
    ASSERT_EQ(mio::abm::go_to_hospital(p_car, t, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, go_shopping)
{
    auto hospital = mio::abm::Location(mio::abm::LocationType::Hospital, 0);
    auto home     = mio::abm::Location(mio::abm::LocationType::Home, 0);

    auto t_weekday = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(9);
    auto t_sunday  = mio::abm::TimePoint(0) + mio::abm::days(6) + mio::abm::hours(9);
    auto t_night   = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(1);
    auto dt        = mio::abm::hours(1);

    auto p_hosp =
        make_test_person(hospital, mio::abm::AgeGroup::Age0to4, mio::abm::InfectionState::InfectedSymptoms, t_weekday);
    auto p_home = mio::abm::Person(home, mio::abm::AgeGroup::Age60to79);

    ASSERT_EQ(mio::abm::go_to_shop(p_hosp, t_weekday, dt, {}), mio::abm::LocationType::Hospital);
    ASSERT_EQ(mio::abm::go_to_shop(p_home, t_sunday, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_shop(p_home, t_night, dt, {}), mio::abm::LocationType::Home);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(mio::abm::go_to_shop(p_home, t_weekday, dt, {}), mio::abm::LocationType::BasicsShop);
}

TEST(TestMigrationRules, shop_return)
{
    auto params = mio::abm::GlobalInfectionParameters{};
    auto t      = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(9);
    auto dt     = mio::abm::hours(1);

    auto home = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto shop = mio::abm::Location(mio::abm::LocationType::BasicsShop, 0);
    auto p    = make_test_person(home, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedNoSymptoms, t);
    home.add_person(p);
    p.migrate_to(shop);
    p.interact(t, dt, params); //person only returns home after some time passed

    ASSERT_EQ(mio::abm::go_to_shop(p, t, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, go_event)
{
    auto work   = mio::abm::Location(mio::abm::LocationType::Work, 0);
    auto p_work = mio::abm::Person(work, mio::abm::AgeGroup::Age35to59);
    auto home   = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto p_home = mio::abm::Person(home, mio::abm::AgeGroup::Age60to79);

    auto t_weekday  = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(20);
    auto t_saturday = mio::abm::TimePoint(0) + mio::abm::days(5) + mio::abm::hours(10);
    auto t_night    = mio::abm::TimePoint(0) + mio::abm::days(5) + mio::abm::hours(1);
    auto dt         = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_event(p_work, t_weekday, dt, {}), mio::abm::LocationType::Work);
    ASSERT_EQ(mio::abm::go_to_event(p_home, t_night, dt, {}), mio::abm::LocationType::Home);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(mio::abm::go_to_event(p_home, t_weekday, dt, {}), mio::abm::LocationType::SocialEvent);

    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(mio::abm::go_to_event(p_home, t_saturday, dt, {}), mio::abm::LocationType::SocialEvent);
}

TEST(TestMigrationRules, event_return)
{
    auto params = mio::abm::GlobalInfectionParameters{};
    auto t      = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(21);
    auto dt     = mio::abm::hours(3);

    auto home = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto shop = mio::abm::Location(mio::abm::LocationType::SocialEvent, 0);
    auto p    = mio::abm::Person(home, mio::abm::AgeGroup::Age15to34);
    home.add_person(p);
    p.migrate_to(shop);
    p.interact(t, dt, params);

    ASSERT_EQ(mio::abm::go_to_event(p, t, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, icu)
{
    auto hospital = mio::abm::Location(mio::abm::LocationType::Hospital, 0);
    auto t        = mio::abm::TimePoint(12346);
    auto dt       = mio::abm::hours(1);
    auto p_hosp =
        make_test_person(hospital, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedCritical, t);

    ASSERT_EQ(mio::abm::go_to_icu(p_hosp, t, dt, {}), mio::abm::LocationType::ICU);

    auto work   = mio::abm::Location(mio::abm::LocationType::Work, 0);
    auto p_work = make_test_person(work, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSymptoms, t);
    ASSERT_EQ(mio::abm::go_to_icu(p_work, t, dt, {}), mio::abm::LocationType::Work);
}

TEST(TestMigrationRules, recover)
{
    auto hospital = mio::abm::Location(mio::abm::LocationType::Hospital, 0);
    auto t        = mio::abm::TimePoint(12346);
    auto dt       = mio::abm::hours(1);
    auto p_rec =
        make_test_person(hospital, mio::abm::AgeGroup::Age60to79, mio::abm::InfectionState::Recovered, t);
    auto p_inf =
        make_test_person(hospital, mio::abm::AgeGroup::Age60to79, mio::abm::InfectionState::InfectedSevere, t);

    ASSERT_EQ(mio::abm::return_home_when_recovered(p_rec, t, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::return_home_when_recovered(p_inf, t, dt, {}), mio::abm::LocationType::Hospital);
}

TEST(TestMigrationRules, dead)
{
    auto icu    = mio::abm::Location(mio::abm::LocationType::ICU, 0);
    auto t      = mio::abm::TimePoint(12346);
    auto dt       = mio::abm::hours(1);
    auto p_dead = make_test_person(icu, mio::abm::AgeGroup::Age60to79, mio::abm::InfectionState::Dead, t);

    ASSERT_EQ(mio::abm::get_buried(p_dead, t, dt, {}), mio::abm::LocationType::Cemetery);
}
