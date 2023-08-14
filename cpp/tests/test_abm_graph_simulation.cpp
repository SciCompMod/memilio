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
#include "models/graph_abm/graph_world.h"
#include "models/abm/world.h"
#include "models/abm/person.h"
#include "models/graph_abm/graph_simulation.h"
#include "models/graph_abm/graph_simulation.cpp"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "abm_helpers.h"

TEST(TestGraphSimulation, advance_random)
{
    mio::abm::GlobalInfectionParameters infection_params;
    auto world     = mio::graph_abm::GraphWorld(infection_params, 0);
    auto location1 = world.add_location(mio::abm::LocationType::School);
    auto location2 = world.add_location(mio::abm::LocationType::School);
    auto home      = world.add_location(mio::abm::LocationType::Home);
    auto& p1       = world.add_person(location1, mio::abm::AgeGroup::Age5to14);
    auto& p2       = world.add_person(location1, mio::abm::AgeGroup::Age5to14);
    auto& p3       = world.add_person(location2, mio::abm::AgeGroup::Age5to14);
    auto& p4       = world.add_person(location2, mio::abm::AgeGroup::Age5to14);
    p1.set_assigned_location(location1);
    p2.set_assigned_location(location1);
    p3.set_assigned_location(location2);
    p4.set_assigned_location(location2);
    p1.set_assigned_location(home);
    p2.set_assigned_location(home);
    p3.set_assigned_location(home);
    p4.set_assigned_location(home);

    auto sim = mio::graph_abm::GraphSimulation(mio::abm::TimePoint(0), std::move(world));

    sim.advance(mio::abm::TimePoint(0) + mio::abm::hours(18));
    ASSERT_EQ(sim.get_result().get_num_time_points(), 19);
    ASSERT_THAT(sim.get_result().get_times(), ElementsAreLinspace(0.0, 18.0 / 24.0, 19));
    for (auto&& v : sim.get_result()) {
        ASSERT_EQ(v.sum(), 4);
    }
}

TEST(TestGraphSimulation, advance_subpopulation)
{
    mio::abm::GlobalInfectionParameters infection_params;
    auto world       = mio::graph_abm::GraphWorld(infection_params, 0);
    auto location_id = world.add_location(mio::abm::LocationType::School);
    auto& school     = world.get_individualized_location(location_id);
    auto person1 =
        add_test_person(world, location_id, mio::abm::AgeGroup::Age5to14, mio::abm::InfectionState::InfectedSymptoms);
    auto person2 =
        add_test_person(world, location_id, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSymptoms);
    auto person3 =
        add_test_person(world, location_id, mio::abm::AgeGroup::Age35to59, mio::abm::InfectionState::Exposed);

    auto sim = mio::graph_abm::GraphSimulation(mio::abm::TimePoint(0), std::move(world));
    sim.advance(mio::abm::TimePoint(0) + mio::abm::hours(50));

    for (size_t i = 0; i < 51; i++) {
        auto v = school.get_subpopulations().get_value(i);
        // Check whether the number of persons in infected state at the location is consistent
        ASSERT_LE(v[size_t(mio::abm::InfectionState::InfectedSymptoms)], 3);
        // Check the time evolution is correct
        ASSERT_EQ(school.get_subpopulations().get_time(i), ScalarType(i) / 24);
    }
}

TEST(TestGraphSimulation, initializeSubpopulation)
{
    mio::abm::GlobalInfectionParameters infection_params;
    auto world  = mio::graph_abm::GraphWorld(infection_params, 0);
    auto loc_id = world.add_location(mio::abm::LocationType::PublicTransport);
    auto loc    = world.get_individualized_location(loc_id);
    ASSERT_EQ(loc.get_subpopulations().get_num_time_points(), 0);

    auto t   = mio::abm::TimePoint(0);
    auto sim = mio::abm::Simulation(t + mio::abm::days(7), std::move(world));

    ASSERT_EQ(sim.get_world().get_individualized_location(loc_id).get_subpopulations().get_time(0), 7);
}

TEST(TestGraphSimulation, stopsAtTmaxGraphABM)
{
    mio::abm::GlobalInfectionParameters infection_params;
    auto world1 = mio::graph_abm::GraphWorld(infection_params, 0);
    auto world2 = mio::graph_abm::GraphWorld(infection_params, 1);

    auto t0   = mio::abm::TimePoint(0);
    auto dt   = mio::abm::hours(12);
    auto tmax = mio::abm::TimePoint(0) + mio::abm::days(5);

    mio::Graph<mio::SimulationNode<mio::graph_abm::GraphSimulation>, mio::MigrationEdgeABM> g;
    g.add_node(0, t0, std::move(world1));
    g.add_node(1, t0, std::move(world2));
    g.add_edge(0, 1);

    auto sim = mio::make_migration_sim(t0, dt, std::move(g));

    sim.advance(tmax, true);

    EXPECT_EQ(sim.get_t(), tmax);
}
