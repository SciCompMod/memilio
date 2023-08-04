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
#include "abm/person.h"
#include "abm_helpers.h"
#include "memilio/utils/random_number_generator.h"

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

    mio::abm::Location home(mio::abm::LocationType::Home, 0);
    auto p_child = mio::abm::Person(rng, home, mio::abm::AgeGroup::Age5to14);
    auto p_adult = mio::abm::Person(rng, home, mio::abm::AgeGroup::Age15to34);

    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(7);
    auto t_weekend = mio::abm::TimePoint(0) + mio::abm::days(5) + mio::abm::hours(7);
    auto dt        = mio::abm::hours(1);

    auto child_rng = mio::abm::Person::RandomNumberGenerator(rng, p_child);
    auto adult_rng = mio::abm::Person::RandomNumberGenerator(rng, p_child);
    ASSERT_EQ(mio::abm::go_to_school(child_rng, p_child, t_morning, dt, {}), mio::abm::LocationType::School);
    ASSERT_EQ(mio::abm::go_to_school(adult_rng, p_adult, t_morning, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_school(child_rng, p_child, t_weekend, dt, {}), mio::abm::LocationType::Home);
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

    mio::abm::Location home(mio::abm::LocationType::Home, 0);
    auto p_child_goes_to_school_at_6   = mio::abm::Person(rng, home, mio::abm::AgeGroup::Age5to14);
    auto rng_child_goes_to_school_at_6 = mio::abm::Person::RandomNumberGenerator(rng, p_child_goes_to_school_at_6);
    auto p_child_goes_to_school_at_8   = mio::abm::Person(rng, home, mio::abm::AgeGroup::Age5to14);
    auto rng_child_goes_to_school_at_8 = mio::abm::Person::RandomNumberGenerator(rng, p_child_goes_to_school_at_8);

    auto t_morning_6 = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto t_morning_8 = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto dt          = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_school(rng_child_goes_to_school_at_6, p_child_goes_to_school_at_6, t_morning_6, dt, {}),
              mio::abm::LocationType::School);
    ASSERT_EQ(mio::abm::go_to_school(rng_child_goes_to_school_at_6, p_child_goes_to_school_at_6, t_morning_8, dt, {}),
              mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_school(rng_child_goes_to_school_at_8, p_child_goes_to_school_at_8, t_morning_6, dt, {}),
              mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_school(rng_child_goes_to_school_at_8, p_child_goes_to_school_at_8, t_morning_8, dt, {}),
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

    mio::abm::Location home(mio::abm::LocationType::Home, 0);
    auto p_child_goes_to_school_at_6    = mio::abm::Person(rng, home, mio::abm::AgeGroup::Age5to14);
    auto rng_child_goes_to_school_at_6  = mio::abm::Person::RandomNumberGenerator(rng, p_child_goes_to_school_at_6);
    auto p_child_goes_to_school_at_8_30 = mio::abm::Person(rng, home, mio::abm::AgeGroup::Age5to14);
    auto rng_child_goes_to_school_at_8_30 =
        mio::abm::Person::RandomNumberGenerator(rng, p_child_goes_to_school_at_8_30);

    auto t_morning_6    = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto t_morning_8_30 = mio::abm::TimePoint(0) + mio::abm::hours(8) + mio::abm::seconds(1800);
    auto dt             = mio::abm::seconds(1800);

    ASSERT_EQ(mio::abm::go_to_school(rng_child_goes_to_school_at_6, p_child_goes_to_school_at_6, t_morning_6, dt, {}),
              mio::abm::LocationType::School);
    ASSERT_EQ(
        mio::abm::go_to_school(rng_child_goes_to_school_at_6, p_child_goes_to_school_at_6, t_morning_8_30, dt, {}),
        mio::abm::LocationType::Home);
    ASSERT_EQ(
        mio::abm::go_to_school(rng_child_goes_to_school_at_8_30, p_child_goes_to_school_at_8_30, t_morning_6, dt, {}),
        mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_school(rng_child_goes_to_school_at_8_30, p_child_goes_to_school_at_8_30, t_morning_8_30,
                                     dt, {}),
              mio::abm::LocationType::School);
}

TEST(TestMigrationRules, school_return)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location school(mio::abm::LocationType::School, 0);
    auto p_child   = mio::abm::Person(rng, school, mio::abm::AgeGroup::Age5to14);
    auto rng_child = mio::abm::Person::RandomNumberGenerator(rng, p_child);

    auto t  = mio::abm::TimePoint(0) + mio::abm::hours(15);
    auto dt = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_school(rng_child, p_child, t, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, worker_goes_to_work)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location home(mio::abm::LocationType::Home, 0);
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

    auto p_retiree   = mio::abm::Person(rng, home, mio::abm::AgeGroup::Age60to79);
    auto rng_retiree = mio::abm::Person::RandomNumberGenerator(rng, p_retiree);
    auto p_adult     = mio::abm::Person(rng, home, mio::abm::AgeGroup::Age15to34);
    auto rng_adult   = mio::abm::Person::RandomNumberGenerator(rng, p_adult);

    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto t_night   = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(4);
    auto dt        = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_work(rng_retiree, p_retiree, t_morning, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(rng_adult, p_adult, t_morning, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(rng_adult, p_adult, t_night, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, worker_goes_to_work_with_non_dividable_timespan)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location home(mio::abm::LocationType::Home, 0);
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

    auto p_retiree   = mio::abm::Person(rng, home, mio::abm::AgeGroup::Age60to79);
    auto rng_retiree = mio::abm::Person::RandomNumberGenerator(rng, p_retiree);
    auto p_adult     = mio::abm::Person(rng, home, mio::abm::AgeGroup::Age15to34);
    auto rng_adult   = mio::abm::Person::RandomNumberGenerator(rng, p_adult);

    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto t_night   = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(4);
    auto dt        = mio::abm::minutes(53);

    ASSERT_EQ(mio::abm::go_to_work(rng_retiree, p_retiree, t_morning, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(rng_adult, p_adult, t_morning, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(rng_adult, p_adult, t_night, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, workers_go_to_work_in_different_times)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location home(mio::abm::LocationType::Home, 0);
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

    auto p_adult_goes_to_work_at_6   = mio::abm::Person(rng, home, mio::abm::AgeGroup::Age15to34);
    auto rng_adult_goes_to_work_at_6 = mio::abm::Person::RandomNumberGenerator(rng, p_adult_goes_to_work_at_6);
    auto p_adult_goes_to_work_at_8   = mio::abm::Person(rng, home, mio::abm::AgeGroup::Age15to34);
    auto rng_adult_goes_to_work_at_8 = mio::abm::Person::RandomNumberGenerator(rng, p_adult_goes_to_work_at_8);

    auto t_morning_6 = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto t_morning_8 = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto t_night     = mio::abm::TimePoint(0) + mio::abm::days(1) + mio::abm::hours(4);
    auto dt          = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_6, p_adult_goes_to_work_at_6, t_morning_6, dt, {}),
              mio::abm::LocationType::Work);
    ASSERT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_6, p_adult_goes_to_work_at_6, t_morning_8, dt, {}),
              mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_6, p_adult_goes_to_work_at_6, t_night, dt, {}),
              mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_8, p_adult_goes_to_work_at_8, t_morning_6, dt, {}),
              mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_8, p_adult_goes_to_work_at_8, t_morning_8, dt, {}),
              mio::abm::LocationType::Work);
    ASSERT_EQ(mio::abm::go_to_work(rng_adult_goes_to_work_at_8, p_adult_goes_to_work_at_8, t_night, dt, {}),
              mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, work_return)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location work(mio::abm::LocationType::Work, 0);
    auto p_adult   = mio::abm::Person(rng, work, mio::abm::AgeGroup::Age35to59);
    auto rng_adult = mio::abm::Person::RandomNumberGenerator(rng, p_adult);
    auto t         = mio::abm::TimePoint(0) + mio::abm::hours(17);
    auto dt        = mio::abm::hours(1);
    ASSERT_EQ(mio::abm::go_to_work(rng_adult, p_adult, t, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, quarantine)
{
    auto rng = mio::RandomNumberGenerator();
    auto t   = mio::abm::TimePoint(12346);
    auto dt  = mio::abm::hours(1);

    mio::abm::Location home(mio::abm::LocationType::Home, 0);
    mio::abm::Location work(mio::abm::LocationType::Work, 0);
    mio::abm::Location hospital(mio::abm::LocationType::Hospital, 0);

    auto p_inf1   = make_test_person(work, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSymptoms, t);
    auto rng_inf1 = mio::abm::Person::RandomNumberGenerator(rng, p_inf1);
    p_inf1.detect_infection(t);
    ASSERT_EQ(mio::abm::go_to_quarantine(rng_inf1, p_inf1, t, dt, {}),
              mio::abm::LocationType::Home); //detected infected person quarantines at home

    auto p_inf2   = make_test_person(work, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSymptoms, t);
    auto rng_inf2 = mio::abm::Person::RandomNumberGenerator(rng, p_inf2);
    ASSERT_EQ(mio::abm::go_to_quarantine(rng_inf2, p_inf2, t, dt, {}),
              mio::abm::LocationType::Work); //undetected infected person does not quaratine

    auto p_inf3 =
        make_test_person(hospital, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSevere, t);
    auto rng_inf3 = mio::abm::Person::RandomNumberGenerator(rng, p_inf3);
    p_inf3.detect_infection(t);
    ASSERT_EQ(mio::abm::go_to_quarantine(rng_inf3, p_inf3, t, dt, {}),
              mio::abm::LocationType::Hospital); //detected infected person does not leave hospital to quarantine
}

TEST(TestMigrationRules, hospital)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location home(mio::abm::LocationType::Home, 0);
    auto t       = mio::abm::TimePoint(12346);
    auto dt      = mio::abm::hours(1);
    auto p_inf   = make_test_person(home, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSevere, t);
    auto rng_inf = mio::abm::Person::RandomNumberGenerator(rng, p_inf);

    ASSERT_EQ(mio::abm::go_to_hospital(rng_inf, p_inf, t, dt, {}), mio::abm::LocationType::Hospital);

    auto p_car   = make_test_person(home, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSymptoms);
    auto rng_car = mio::abm::Person::RandomNumberGenerator(rng, p_car);
    ASSERT_EQ(mio::abm::go_to_hospital(rng_car, p_car, t, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, go_shopping)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location hospital(mio::abm::LocationType::Hospital, 0);
    mio::abm::Location home(mio::abm::LocationType::Home, 0);

    auto t_weekday = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(9);
    auto t_sunday  = mio::abm::TimePoint(0) + mio::abm::days(6) + mio::abm::hours(9);
    auto t_night   = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(1);
    auto dt        = mio::abm::hours(1);

    auto p_hosp =
        make_test_person(hospital, mio::abm::AgeGroup::Age0to4, mio::abm::InfectionState::InfectedSymptoms, t_weekday);
    auto rng_hosp = mio::abm::Person::RandomNumberGenerator(rng, p_hosp);
    auto p_home   = mio::abm::Person(rng, home, mio::abm::AgeGroup::Age60to79);
    auto rng_home = mio::abm::Person::RandomNumberGenerator(rng, p_home);

    ASSERT_EQ(mio::abm::go_to_shop(rng_hosp, p_hosp, t_weekday, dt, {}), mio::abm::LocationType::Hospital);
    ASSERT_EQ(mio::abm::go_to_shop(rng_home, p_home, t_sunday, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::go_to_shop(rng_home, p_home, t_night, dt, {}), mio::abm::LocationType::Home);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(mio::abm::go_to_shop(rng_home, p_home, t_weekday, dt, {}), mio::abm::LocationType::BasicsShop);
}

TEST(TestMigrationRules, shop_return)
{
    auto rng    = mio::RandomNumberGenerator();
    auto params = mio::abm::GlobalInfectionParameters{};
    auto t      = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(9);
    auto dt     = mio::abm::hours(1);

    mio::abm::Location home(mio::abm::LocationType::Home, 0);
    mio::abm::Location shop(mio::abm::LocationType::BasicsShop, 0);
    auto p     = make_test_person(home, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedNoSymptoms, t);
    auto rng_p = mio::abm::Person::RandomNumberGenerator(rng, p);
    home.add_person(p);
    p.migrate_to(shop);
    p.interact(rng_p, t, dt, params); //person only returns home after some time passed

    ASSERT_EQ(mio::abm::go_to_shop(rng_p, p, t, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, go_event)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location work(mio::abm::LocationType::Work, 0);
    auto p_work   = mio::abm::Person(rng, work, mio::abm::AgeGroup::Age35to59);
    auto rng_work = mio::abm::Person::RandomNumberGenerator(rng, p_work);
    mio::abm::Location home(mio::abm::LocationType::Home, 0);
    auto p_home   = mio::abm::Person(rng, home, mio::abm::AgeGroup::Age60to79);
    auto rng_home = mio::abm::Person::RandomNumberGenerator(rng, p_home);

    auto t_weekday  = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(20);
    auto t_saturday = mio::abm::TimePoint(0) + mio::abm::days(5) + mio::abm::hours(10);
    auto t_night    = mio::abm::TimePoint(0) + mio::abm::days(5) + mio::abm::hours(1);
    auto dt         = mio::abm::hours(1);

    ASSERT_EQ(mio::abm::go_to_event(rng_work, p_work, t_weekday, dt, {}), mio::abm::LocationType::Work);
    ASSERT_EQ(mio::abm::go_to_event(rng_home, p_home, t_night, dt, {}), mio::abm::LocationType::Home);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(mio::abm::go_to_event(rng_home, p_home, t_weekday, dt, {}), mio::abm::LocationType::SocialEvent);

    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(testing::Return(0.01));
    ASSERT_EQ(mio::abm::go_to_event(rng_home, p_home, t_saturday, dt, {}), mio::abm::LocationType::SocialEvent);
}

TEST(TestMigrationRules, event_return)
{
    auto rng    = mio::RandomNumberGenerator();
    auto params = mio::abm::GlobalInfectionParameters{};
    auto t      = mio::abm::TimePoint(0) + mio::abm::days(4) + mio::abm::hours(21);
    auto dt     = mio::abm::hours(3);

    mio::abm::Location home(mio::abm::LocationType::Home, 0);
    mio::abm::Location social_event(mio::abm::LocationType::SocialEvent, 0);
    auto p     = mio::abm::Person(rng, home, mio::abm::AgeGroup::Age15to34);
    auto rng_p = mio::abm::Person::RandomNumberGenerator(rng, p);
    home.add_person(p);
    p.migrate_to(social_event);
    p.interact(rng_p, t, dt, params);

    ASSERT_EQ(mio::abm::go_to_event(rng_p, p, t, dt, {}), mio::abm::LocationType::Home);
}

TEST(TestMigrationRules, icu)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location hospital(mio::abm::LocationType::Hospital, 0);
    auto t  = mio::abm::TimePoint(12346);
    auto dt = mio::abm::hours(1);
    auto p_hosp =
        make_test_person(hospital, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedCritical, t);
    auto rng_hosp = mio::abm::Person::RandomNumberGenerator(rng, p_hosp);

    ASSERT_EQ(mio::abm::go_to_icu(rng_hosp, p_hosp, t, dt, {}), mio::abm::LocationType::ICU);

    mio::abm::Location work(mio::abm::LocationType::Work, 0);
    auto p_work   = make_test_person(work, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSymptoms, t);
    auto rng_work = mio::abm::Person::RandomNumberGenerator(rng, p_work);
    ASSERT_EQ(mio::abm::go_to_icu(rng_work, p_work, t, dt, {}), mio::abm::LocationType::Work);
}

TEST(TestMigrationRules, recover)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location hospital(mio::abm::LocationType::Hospital, 0);
    auto t  = mio::abm::TimePoint(12346);
    auto dt = mio::abm::hours(1);
    auto p_rec =
        make_test_person(hospital, mio::abm::AgeGroup::Age60to79, mio::abm::InfectionState::Recovered, t);
    auto rng_rec = mio::abm::Person::RandomNumberGenerator(rng, p_rec);
    auto p_inf =
        make_test_person(hospital, mio::abm::AgeGroup::Age60to79, mio::abm::InfectionState::InfectedSevere, t);
    auto rng_inf = mio::abm::Person::RandomNumberGenerator(rng, p_inf);

    ASSERT_EQ(mio::abm::return_home_when_recovered(rng_rec, p_rec, t, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(mio::abm::return_home_when_recovered(rng_inf, p_inf, t, dt, {}), mio::abm::LocationType::Hospital);
}

TEST(TestMigrationRules, dead)
{
    auto rng = mio::RandomNumberGenerator();

    mio::abm::Location icu(mio::abm::LocationType::ICU, 0);
    auto t      = mio::abm::TimePoint(12346);
    auto dt     = mio::abm::hours(1);
    auto p_dead = make_test_person(icu, mio::abm::AgeGroup::Age60to79, mio::abm::InfectionState::Dead, t);
    auto p_rng  = mio::abm::Person::RandomNumberGenerator(rng, p_dead);

    ASSERT_EQ(mio::abm::get_buried(p_rng, p_dead, t, dt, {}), mio::abm::LocationType::Cemetery);
}
