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
#include "memilio/io/history.h"

struct LogTimePoint : mio::LogAlways {
    using Type = double;
    static Type log(const mio::abm::Simulation& sim)
    {
        return sim.get_time().hours();
    }
};

TEST(TestSimulation, advance_random)
{
    auto world     = mio::abm::World();
    auto location1 = world.add_location(mio::abm::LocationType::School);
    auto location2 = world.add_location(mio::abm::LocationType::School);
    auto& p1       = world.add_person(location1, mio::abm::AgeGroup::Age5to14);
    auto& p2       = world.add_person(location1, mio::abm::AgeGroup::Age5to14);
    auto& p3       = world.add_person(location2, mio::abm::AgeGroup::Age5to14);
    auto& p4       = world.add_person(location2, mio::abm::AgeGroup::Age5to14);
    p1.set_assigned_location(location1);
    p2.set_assigned_location(location1);
    p3.set_assigned_location(location2);
    p4.set_assigned_location(location2);

    auto sim = mio::abm::Simulation(mio::abm::TimePoint(0), std::move(world));

    sim.advance(mio::abm::TimePoint(0) + mio::abm::hours(50));
    ASSERT_EQ(sim.get_result().get_num_time_points(), 51);
    ASSERT_THAT(sim.get_result().get_times(), ElementsAreLinspace(0.0, 50.0 / 24.0, 51));
    for (auto&& v : sim.get_result()) {
        ASSERT_EQ(v.sum(), 4);
    }
}

TEST(TestDiscreteDistribution, generate)
{
    using namespace mio;
    auto distribution = mio::DiscreteDistribution<size_t>();

    std::vector<double> weights;
    for (size_t i = 0; i < 50; i++) {
        weights = {};
        ASSERT_EQ(distribution(weights), 0);

        weights = {0.5};
        ASSERT_EQ(distribution(weights), 0);

        weights = {0.5, 1.3, 0.1, 0.4, 0.3};
        auto d  = distribution(weights);
        ASSERT_GE(d, 0);
        ASSERT_LE(d, 4);
    }
}

TEST(TestSimulation, advance_subpopulation)
{
    auto world       = mio::abm::World();
    auto location_id = world.add_location(mio::abm::LocationType::School);
    auto& school     = world.get_individualized_location(location_id);
    auto& person1 =
        add_test_person(world, location_id, mio::abm::AgeGroup::Age5to14, mio::abm::InfectionState::InfectedSymptoms);
    auto& person2 =
        add_test_person(world, location_id, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSymptoms);
    auto& person3 =
        add_test_person(world, location_id, mio::abm::AgeGroup::Age35to59, mio::abm::InfectionState::Exposed);
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
    auto world  = mio::abm::World();
    auto loc_id = world.add_location(mio::abm::LocationType::PublicTransport, 3);
    auto loc    = world.get_individualized_location(loc_id);
    ASSERT_EQ(loc.get_subpopulations().get_num_time_points(), 0);

    auto t   = mio::abm::TimePoint(0);
    auto sim = mio::abm::Simulation(t + mio::abm::days(7), std::move(world));

    ASSERT_EQ(sim.get_world().get_individualized_location(loc_id).get_subpopulations().get_time(0), 7);
}

TEST(TestSimulation, getWorldAndTimeConst)
{

    auto t     = mio::abm::TimePoint(0);
    auto world = mio::abm::World();
    auto sim   = mio::abm::Simulation(t + mio::abm::days(7), std::move(world));

    auto t_test = mio::abm::days(7);
    ASSERT_EQ(sim.get_time(), mio::abm::TimePoint(t_test.seconds()));

    const mio::abm::World world_test{std::move(sim.get_world())};
    EXPECT_EQ(world_test.get_locations().size(), 1);
}

TEST(TestSimulation, advanceWithHistory)
{

    auto world = mio::abm::World();
    auto sim   = mio::abm::Simulation(mio::abm::TimePoint(0), std::move(world));
    mio::HistoryWithMemoryWriter<LogTimePoint> history;

    sim.advance(mio::abm::TimePoint(0) + mio::abm::hours(2), history);
    ASSERT_EQ(std::get<0>(history.get_log()).size(), 3);
    ASSERT_NEAR(std::get<0>(history.get_log())[0], 0.0, 1e-14);
    ASSERT_NEAR(std::get<0>(history.get_log())[1], 1.0, 1e-14);
    auto test_get_templ_log = history.get_log<LogTimePoint>();
    ASSERT_NEAR(test_get_templ_log[2], 2.0, 1e-14);
}
