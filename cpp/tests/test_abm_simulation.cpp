/* 
* Copyright (C) 2020-2026 MEmilio
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
#include "abm/simulation.h"
#include "abm/result_simulation.h"
#include "abm/time.h"
#include "abm_helpers.h"
#include "abm/common_abm_loggers.h"
#include "matchers.h"
#include "memilio/io/history.h"
#include <cstddef>
#include <cstdint>

TEST(TestSimulation, advance_random)
{
    auto model     = mio::abm::Model(num_age_groups);
    auto location1 = model.add_location(mio::abm::LocationType::School);
    auto location2 = model.add_location(mio::abm::LocationType::School);
    auto p1        = model.add_person(location1, age_group_5_to_14);
    auto p2        = model.add_person(location1, age_group_5_to_14);
    auto p3        = model.add_person(location2, age_group_5_to_14);
    auto p4        = model.add_person(location2, age_group_5_to_14);

    model.assign_location(p1, location1);
    model.assign_location(p2, location1);
    model.assign_location(p3, location2);
    model.assign_location(p4, location2);

    auto sim = mio::abm::Simulation(mio::abm::TimePoint(0), std::move(model));

    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
        Eigen::Index(mio::abm::InfectionState::Count)};

    sim.advance(mio::abm::TimePoint(0) + mio::abm::hours(50), historyTimeSeries);
    auto log = std::get<0>(historyTimeSeries.get_log());
    EXPECT_EQ(log.get_num_time_points(), 51);
    EXPECT_THAT(log.get_times(), ElementsAreLinspace(0.0, 50.0 / 24.0, 51));
    for (auto&& v : log) {
        EXPECT_EQ(v.sum(), 4);
    }
}

TEST(TestSimulation, getModelAndTimeConst)
{
    auto t     = mio::abm::TimePoint(0);
    auto model = mio::abm::Model(num_age_groups);
    auto sim   = mio::abm::Simulation(t + mio::abm::days(7), std::move(model));

    auto t_test = mio::abm::days(7);
    EXPECT_EQ(sim.get_time(), mio::abm::TimePoint(t_test.seconds()));

    const mio::abm::Model model_test{std::move(sim.get_model())};
    EXPECT_EQ(model_test.get_locations().size(), 1);
}

TEST(TestSimulation, advanceWithCommonHistory)
{
    auto model       = mio::abm::Model(num_age_groups);
    auto home_id     = model.add_location(mio::abm::LocationType::Home);
    auto work_id     = model.add_location(mio::abm::LocationType::Work);
    auto icu_id      = model.add_location(mio::abm::LocationType::ICU);
    auto hospital_id = model.add_location(mio::abm::LocationType::Hospital);
    auto school_id   = model.add_location(mio::abm::LocationType::School);
    auto social_id   = model.add_location(mio::abm::LocationType::SocialEvent);
    auto basics_id   = model.add_location(mio::abm::LocationType::BasicsShop);
    auto public_id   = model.add_location(mio::abm::LocationType::PublicTransport);

    auto person1 = add_test_person(model, home_id, age_group_5_to_14, mio::abm::InfectionState::Exposed);
    auto person2 = add_test_person(model, home_id, age_group_15_to_34, mio::abm::InfectionState::Exposed);
    auto person3 = add_test_person(model, home_id, age_group_35_to_59, mio::abm::InfectionState::Dead);

    model.assign_location(person1, home_id);
    model.assign_location(person2, home_id);
    model.assign_location(person3, home_id);
    model.assign_location(person1, school_id);
    model.assign_location(person2, work_id);
    model.assign_location(person2, icu_id);
    model.assign_location(person2, hospital_id);
    model.assign_location(person1, social_id);
    model.assign_location(person2, social_id);
    model.assign_location(person3, social_id);
    model.assign_location(person1, basics_id);
    model.assign_location(person2, basics_id);
    model.assign_location(person3, basics_id);
    model.assign_location(person2, public_id);

    mio::abm::TripList& trip_list = model.get_trip_list();

    // We add trips for person two to test the history and if it is working correctly
    mio::abm::Trip trip1(static_cast<uint64_t>(person2.get()), mio::abm::TimePoint(0) + mio::abm::hours(2), home_id);
    mio::abm::Trip trip2(static_cast<uint64_t>(person2.get()), mio::abm::TimePoint(0) + mio::abm::hours(3), home_id);
    mio::abm::Trip trip3(static_cast<uint64_t>(person2.get()), mio::abm::TimePoint(0) + mio::abm::hours(4), home_id);
    mio::abm::Trip trip4(static_cast<uint64_t>(person2.get()), mio::abm::TimePoint(0) + mio::abm::hours(5), home_id);
    mio::abm::Trip trip5(static_cast<uint64_t>(person2.get()), mio::abm::TimePoint(0) + mio::abm::hours(6), home_id);
    mio::abm::Trip trip6(static_cast<uint64_t>(person2.get()), mio::abm::TimePoint(0) + mio::abm::hours(7), home_id);
    mio::abm::Trip trip7(static_cast<uint64_t>(person2.get()), mio::abm::TimePoint(0) + mio::abm::hours(8), home_id);

    // Add to one vector
    auto trips = std::vector<mio::abm::Trip>{trip1, trip2, trip3, trip4, trip5, trip6, trip7};

    trip_list.add_trips(trips);

    mio::History<mio::DataWriterToMemory, mio::abm::LogLocationInformation, mio::abm::LogPersonInformation,
                 mio::abm::LogDataForMobility>
        historyPersonInf;
    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
        Eigen::Index(mio::abm::InfectionState::Count)};
    mio::History<mio::abm::DataWriterToMemoryDelta, mio::abm::LogDataForMobility> historyPersonInfDelta;
    auto sim = mio::abm::Simulation(mio::abm::TimePoint(0), std::move(model));
    sim.advance(mio::abm::TimePoint(0) + mio::abm::hours(24), historyPersonInf, historyTimeSeries,
                historyPersonInfDelta);

    auto logLocationInfo      = std::get<0>(historyPersonInf.get_log());
    auto logPersonInfo        = std::get<1>(historyPersonInf.get_log());
    auto logMobilityInfo      = std::get<2>(historyPersonInf.get_log());
    auto logTimeSeries        = std::get<0>(historyTimeSeries.get_log());
    auto logMobilityInfoDelta = std::get<0>(historyPersonInfDelta.get_log());

    EXPECT_EQ(logLocationInfo[0].size(), 9); // Check if all locations are in the log, 9 locations
    EXPECT_EQ(logPersonInfo[0].size(), 3); // Check if all persons are in the log, 3 persons
    EXPECT_EQ(
        logMobilityInfo.size(),
        25); // Check if for all time points Mobility data is in the log, 25 time points (24 hours + 1 for the initial state)
    EXPECT_EQ(logTimeSeries.get_num_time_points(),
              25); // Check if all time points are in the log, 25 time points (24 hours + 1 for the initial state)
    EXPECT_EQ(
        logMobilityInfoDelta.size(),
        26); // Check if for all time points Mobility data is in the log, 26 time points (24 hours + 1 for the initial state + 1 helper entry for calculating the delta)
    EXPECT_EQ(logMobilityInfoDelta[0].size(),
              3); // Check if all persons are in the delta-logger Mobility helper entry 0, 3 persons
    EXPECT_EQ(logMobilityInfoDelta[1].size(), 3); // Check if all persons are in the delta-log first entry, 3 persons
}

TEST(TestSimulation, ResultSimulation)
{
    // run a ResultSimulation on a minimal setup
    auto model    = mio::abm::Model(num_age_groups);
    auto location = model.add_location(mio::abm::LocationType::Home);
    auto person   = model.add_person(location, age_group_15_to_34);

    model.assign_location(person, location);

    const auto t0   = mio::abm::TimePoint(0) + mio::abm::hours(100);
    const auto tmax = t0 + mio::abm::hours(50);
    auto sim        = mio::abm::ResultSimulation(std::move(model), t0);

    // run simulation. expect one timepoint per day, but nothing to change in the results
    sim.advance(tmax);
    const size_t N = (size_t)(tmax - t0).hours() + 1;
    ASSERT_EQ(sim.get_result().get_num_time_points(), N);
    EXPECT_THAT(sim.get_result().get_times(), ElementsAreLinspace(t0.days(), tmax.days(), N));
    for (const auto& tp : sim.get_result()) {
        EXPECT_EQ(tp.sum(), 1.0);
    }
    EXPECT_EQ(sim.get_result().get_value(0)[(Eigen::Index)mio::abm::InfectionState::Susceptible], 1.0);
    EXPECT_EQ(sim.get_result().get_value(N - 1)[(Eigen::Index)mio::abm::InfectionState::Susceptible], 1.0);
}
