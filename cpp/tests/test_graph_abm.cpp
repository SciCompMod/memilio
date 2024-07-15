/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Julia Bicker
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

#include "abm/world.h"
#include "abm/location_type.h"
#include "abm/time.h"
#include "graph_abm/graph_abm_mobility.h"
#include "graph_abm/mobility_rules.h"
#include "memilio/epidemiology/age_group.h"
#include <gtest/gtest.h>

TEST(TestGraphAbm, test_activessness)
{
    auto world = mio::abm::World(size_t(1));
    world.parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple({mio::AgeGroup(0)}, true);
    auto work_id = world.add_location(mio::abm::LocationType::Work);
    auto home_id = world.add_location(mio::abm::LocationType::Home);
    auto& p1     = world.add_person(home_id, mio::AgeGroup(0));
    auto& p2     = world.add_person(home_id, mio::AgeGroup(0));
    p1.set_assigned_location(work_id);
    p2.set_assigned_location(work_id);
    p1.set_assigned_location(home_id);
    p2.set_assigned_location(home_id);

    auto& home = world.get_individualized_location(home_id);
    auto& work = world.get_individualized_location(work_id);

    EXPECT_EQ(p1.get_location(), home);
    EXPECT_EQ(p2.get_location(), home);
    EXPECT_EQ(world.get_activeness_statuses().size(), 2);

    world.change_activeness(p1.get_person_id());
    EXPECT_EQ(world.get_activeness_statuses()[p1.get_person_id()], false);
    EXPECT_EQ(world.get_activeness_statuses()[p2.get_person_id()], true);

    auto t  = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto dt = mio::abm::hours(3);

    world.evolve(t, dt);

    //inactive persons do not move
    EXPECT_EQ(p1.get_location(), home);
    EXPECT_EQ(p2.get_location(), work);
}

struct MockHistory {

    template <class T>
    void log(const T& t)
    {
        mio::unused(t);
    }
};

TEST(TestGraphAbm, test_evolve_node)
{
    auto t  = mio::abm::TimePoint(0);
    auto dt = mio::abm::hours(2);
    mio::ABMSimulationNode<MockHistory> node(MockHistory{}, t, size_t(1));
    node.evolve(t, dt);

    EXPECT_EQ(node.get_simulation().get_time(), mio::abm::TimePoint(dt.seconds()));
}

TEST(TestGraphAbm, test_apply_mobility)
{
    auto world_1   = mio::abm::World(size_t(1), 1);
    world_1.parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple({mio::AgeGroup(0)}, true);
    auto work_id_1 = world_1.add_location(mio::abm::LocationType::Work);
    world_1.get_individualized_location(work_id_1).set_world_id(world_1.get_id());
    auto work_id_2 = world_1.add_location(mio::abm::LocationType::Work);
    auto& work_2 = world_1.get_individualized_location(work_id_2);
    work_2.set_world_id(2);
    auto home_id = world_1.add_location(mio::abm::LocationType::Home);
    world_1.get_individualized_location(home_id).set_world_id(world_1.get_id());
    auto& p1 = world_1.add_person(home_id, mio::AgeGroup(0));
    auto& p2 = world_1.add_person(home_id, mio::AgeGroup(0));
    p1.set_assigned_location(work_id_1);
    p2.set_assigned_location(work_id_2);
    p1.set_assigned_location(home_id);
    p2.set_assigned_location(home_id);
    auto world_2 = mio::abm::World(world_1, 2);
    // Deactivate persons in world 2
    world_2.change_activeness(p1.get_person_id());
    world_2.change_activeness(p2.get_person_id());

    auto t0 = mio::abm::TimePoint(0) + mio::abm::hours(6);
    mio::ABMSimulationNode<MockHistory> node1(MockHistory{}, t0, std::move(world_1));
    mio::ABMSimulationNode<MockHistory> node2(MockHistory{}, t0, std::move(world_2));

    mio::ABMMobilityEdge<MockHistory> edge({p2.get_person_id()}, {&mio::apply_commuting});
    edge.apply_mobility(node1, node2, t0+mio::abm::hours(2));

    EXPECT_EQ(work_2.get_number_persons(), 1);
    EXPECT_EQ(node1.get_simulation().get_world().get_activeness_statuses()[p1.get_person_id()], true);
    EXPECT_EQ(node1.get_simulation().get_world().get_activeness_statuses()[p2.get_person_id()], false);
    EXPECT_EQ(node2.get_simulation().get_world().get_activeness_statuses()[p1.get_person_id()], false);
    EXPECT_EQ(node2.get_simulation().get_world().get_activeness_statuses()[p2.get_person_id()], true);
}