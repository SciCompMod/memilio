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
#include "abm_helpers.h"
#include "abm/common_abm_loggers.h"
#include "memilio/io/history.h"

TEST(TestSimulation, advance_random)
{
    auto world     = mio::abm::World(NUM_AGE_GROUPS);
    auto location1 = world.add_location(mio::abm::LocationType::School);
    auto location2 = world.add_location(mio::abm::LocationType::School);
    auto& p1       = world.add_person(location1, AGE_GROUP_5_TO_14);
    auto& p2       = world.add_person(location1, AGE_GROUP_5_TO_14);
    auto& p3       = world.add_person(location2, AGE_GROUP_5_TO_14);
    auto& p4       = world.add_person(location2, AGE_GROUP_5_TO_14);
    p1.set_assigned_location(location1);
    p2.set_assigned_location(location1);
    p3.set_assigned_location(location2);
    p4.set_assigned_location(location2);

    auto sim = mio::abm::Simulation(mio::abm::TimePoint(0), std::move(world));

    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
        Eigen::Index(mio::abm::InfectionState::Count)};

    sim.advance(mio::abm::TimePoint(0) + mio::abm::hours(50), historyTimeSeries);
    auto log = std::get<0>(historyTimeSeries.get_log());
    ASSERT_EQ(log.get_num_time_points(), 51);
    ASSERT_THAT(log.get_times(), ElementsAreLinspace(0.0, 50.0 / 24.0, 51));
    for (auto&& v : log) {
        ASSERT_EQ(v.sum(), 4);
    }
}

TEST(TestSimulation, advance_subpopulation)
{
    auto world       = mio::abm::World(NUM_AGE_GROUPS);
    auto location_id = world.add_location(mio::abm::LocationType::School);
    auto& school     = world.get_individualized_location(location_id);
    auto& person1 = add_test_person(world, location_id, AGE_GROUP_5_TO_14, mio::abm::InfectionState::InfectedSymptoms);
    auto& person2 = add_test_person(world, location_id, AGE_GROUP_15_TO_34, mio::abm::InfectionState::InfectedSymptoms);
    auto& person3 = add_test_person(world, location_id, AGE_GROUP_35_TO_59, mio::abm::InfectionState::Exposed);
    person1.set_assigned_location(location_id);
    person2.set_assigned_location(location_id);
    person3.set_assigned_location(location_id);

    auto sim = mio::abm::Simulation(mio::abm::TimePoint(0), std::move(world));
    sim.advance(mio::abm::TimePoint(0) + mio::abm::hours(50));

    for (size_t i = 0; i < 51; i++) {
        auto v = school.get_subpopulations().get_value(i);
        // Check whether the number of persons in infected state at the location is consistent
        ASSERT_LE(v[size_t(mio::abm::InfectionState::InfectedSymptoms)], 3);
        // Check the time evolution is correct
        ASSERT_EQ(school.get_subpopulations().get_time(i), ScalarType(i) / 24);
    }
}

TEST(TestSimulation, initializeSubpopulation)
{
    auto world  = mio::abm::World(NUM_AGE_GROUPS);
    auto loc_id = world.add_location(mio::abm::LocationType::PublicTransport, 3);
    auto& loc   = world.get_individualized_location(loc_id);
    ASSERT_EQ(loc.get_subpopulations().get_num_time_points(), 0);

    auto t   = mio::abm::TimePoint(0);
    auto sim = mio::abm::Simulation(t + mio::abm::days(7), std::move(world));

    ASSERT_EQ(sim.get_world().get_individualized_location(loc_id).get_subpopulations().get_time(0), 7);
}

TEST(TestSimulation, getWorldAndTimeConst)
{

    auto t     = mio::abm::TimePoint(0);
    auto world = mio::abm::World(NUM_AGE_GROUPS);
    auto sim   = mio::abm::Simulation(t + mio::abm::days(7), std::move(world));

    auto t_test = mio::abm::days(7);
    ASSERT_EQ(sim.get_time(), mio::abm::TimePoint(t_test.seconds()));

    const mio::abm::World world_test{std::move(sim.get_world())};
    EXPECT_EQ(world_test.get_locations().size(), 1);
}

TEST(TestSimulation, advanceWithHistory)
{
    auto world       = mio::abm::World(NUM_AGE_GROUPS);
    auto home_id     = world.add_location(mio::abm::LocationType::Home);
    auto work_id     = world.add_location(mio::abm::LocationType::Work);
    auto icu_id      = world.add_location(mio::abm::LocationType::ICU);
    auto hospital_id = world.add_location(mio::abm::LocationType::Hospital);
    auto school_id   = world.add_location(mio::abm::LocationType::School);
    auto social_id   = world.add_location(mio::abm::LocationType::SocialEvent);
    auto basics_id   = world.add_location(mio::abm::LocationType::BasicsShop);
    auto public_id   = world.add_location(mio::abm::LocationType::PublicTransport);

    auto& person1 = add_test_person(world, home_id, AGE_GROUP_5_TO_14, mio::abm::InfectionState::Exposed);
    auto& person2 = add_test_person(world, home_id, AGE_GROUP_15_TO_34, mio::abm::InfectionState::Exposed);
    auto& person3 = add_test_person(world, home_id, AGE_GROUP_35_TO_59, mio::abm::InfectionState::Dead);
    person1.set_assigned_location(home_id);
    person2.set_assigned_location(home_id);
    person3.set_assigned_location(home_id);
    person1.set_assigned_location(school_id);
    person2.set_assigned_location(work_id);
    person2.set_assigned_location(icu_id);
    person2.set_assigned_location(hospital_id);
    person1.set_assigned_location(social_id);
    person2.set_assigned_location(social_id);
    person3.set_assigned_location(social_id);
    person1.set_assigned_location(basics_id);
    person2.set_assigned_location(basics_id);
    person3.set_assigned_location(basics_id);
    person2.set_assigned_location(public_id);

    auto sim = mio::abm::Simulation(mio::abm::TimePoint(0), std::move(world));
    mio::History<mio::DataWriterToMemory, mio::abm::LogLocationInformation, mio::abm::LogPersonInformation,
                 mio::abm::LogDataForMovement>
        historyPersonInf;
    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
        Eigen::Index(mio::abm::InfectionState::Count)};
    mio::History<mio::abm::DataWriterToMemoryDelta, mio::abm::LogDataForMovement> historyPersonInfDelta;
    mio::abm::TripList& trip_list = world.get_trip_list();

    mio::abm::Trip trip1(person2.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(2), work_id);
    mio::abm::Trip trip2(person2.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(3), icu_id);
    mio::abm::Trip trip3(person2.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(4), hospital_id);
    mio::abm::Trip trip4(person2.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(5), social_id);
    mio::abm::Trip trip5(person2.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(6), basics_id);
    mio::abm::Trip trip6(person2.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(7), public_id);
    mio::abm::Trip trip7(person2.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(8), home_id);

    trip_list.add_trip(trip1);
    trip_list.add_trip(trip2);
    trip_list.add_trip(trip3);
    trip_list.add_trip(trip4);
    trip_list.add_trip(trip5);
    trip_list.add_trip(trip6);
    trip_list.add_trip(trip7);

    sim.advance(mio::abm::TimePoint(0) + mio::abm::hours(24), historyPersonInf, historyTimeSeries,
                historyPersonInfDelta);

    auto logLocationInfo      = std::get<0>(historyPersonInf.get_log());
    auto logPersonInfo        = std::get<1>(historyPersonInf.get_log());
    auto logMovementInfo      = std::get<2>(historyPersonInf.get_log());
    auto logTimeSeries        = std::get<0>(historyTimeSeries.get_log());
    auto logMovementInfoDelta = std::get<0>(historyPersonInfDelta.get_log());

    ASSERT_EQ(logLocationInfo[0].size(), 9);
    ASSERT_EQ(logPersonInfo[0].size(), 3);
    ASSERT_EQ(logMovementInfo.size(), 25);
    ASSERT_EQ(logTimeSeries.get_num_time_points(), 25);
    ASSERT_EQ(logMovementInfoDelta.size(), 26);
    ASSERT_EQ(logMovementInfoDelta[0].size(), 3);
    ASSERT_EQ(logMovementInfoDelta[1].size(), 3);
    ASSERT_EQ(logMovementInfoDelta[2].size(), 1);
    ASSERT_EQ(logMovementInfoDelta[17].size(), 1); //everyone returns from school at 15 o 'clock
}
