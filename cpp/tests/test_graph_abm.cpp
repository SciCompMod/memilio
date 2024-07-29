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
#include <cstddef>
#include <gtest/gtest.h>

TEST(TestGraphAbm, test_activessness)
{
    auto world                                                           = mio::abm::World(size_t(1));
    world.parameters.get<mio::abm::AgeGroupGotoWork>()[mio::AgeGroup(0)] = true;
    auto work_id = world.add_location(mio::abm::LocationType::Work);
    auto home_id = world.add_location(mio::abm::LocationType::Home);
    auto& home   = world.get_location(home_id);
    auto& work   = world.get_location(work_id);
    auto p1_id   = world.add_person(home_id, mio::AgeGroup(0));
    auto p2_id   = world.add_person(home_id, mio::AgeGroup(0));
    auto& p1     = world.get_person(p1_id);
    auto& p2     = world.get_person(p2_id);
    p1.set_assigned_location(work.get_type(), work.get_id(), work.get_world_id());
    p2.set_assigned_location(work.get_type(), work.get_id(), work.get_world_id());
    p1.set_assigned_location(home.get_type(), home.get_id(), home.get_world_id());
    p2.set_assigned_location(home.get_type(), home.get_id(), home.get_world_id());

    EXPECT_EQ(p1.get_location(), home_id);
    EXPECT_EQ(p2.get_location(), home_id);
    EXPECT_EQ(p1.get_location_world_id(), world.get_id());
    EXPECT_EQ(p2.get_location_world_id(), world.get_id());
    EXPECT_EQ(world.get_activeness_statuses().size(), 2);

    world.change_activeness(p1.get_id());
    EXPECT_EQ(world.get_activeness_statuses()[p1.get_id().get()], false);
    EXPECT_EQ(world.get_activeness_statuses()[p2.get_id().get()], true);

    auto t  = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto dt = mio::abm::hours(3);

    world.evolve(t, dt);

    //inactive persons do not move
    EXPECT_EQ(p1.get_location(), home_id);
    EXPECT_EQ(p2.get_location(), work_id);
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
    auto world_1                                                           = mio::abm::World(size_t(1), 1);
    auto world_2                                                           = mio::abm::World(size_t(1), 2);
    world_1.parameters.get<mio::abm::AgeGroupGotoWork>()[mio::AgeGroup(0)] = true;
    world_2.parameters.get<mio::abm::AgeGroupGotoWork>()[mio::AgeGroup(0)] = true;
    auto work_id_1 = world_1.add_location(mio::abm::LocationType::Work);
    auto home_id   = world_1.add_location(mio::abm::LocationType::Home);
    auto work_id_2 = world_2.add_location(mio::abm::LocationType::Work);
    auto& work_1   = world_1.get_location(work_id_1);
    auto& work_2   = world_2.get_location(work_id_2);
    auto& home     = world_1.get_location(home_id);

    EXPECT_EQ(work_1.get_world_id(), 1);
    EXPECT_EQ(work_2.get_world_id(), 2);

    auto p1_id = world_1.add_person(home_id, mio::AgeGroup(0));
    auto p2_id = world_1.add_person(home_id, mio::AgeGroup(0));
    auto& p1   = world_1.get_person(p1_id);
    auto& p2   = world_1.get_person(p2_id);
    p1.set_assigned_location(work_1.get_type(), work_1.get_id(), work_1.get_world_id());
    p2.set_assigned_location(work_2.get_type(), work_2.get_id(), work_2.get_world_id());
    p1.set_assigned_location(home.get_type(), home.get_id(), home.get_world_id());
    p2.set_assigned_location(home.get_type(), home.get_id(), home.get_world_id());

    //copy persons to world 2
    world_2.copy_persons_from_other_world(world_1);

    auto t = mio::abm::TimePoint(0) + mio::abm::hours(6);
    mio::ABMSimulationNode<MockHistory> node1(MockHistory{}, t, std::move(world_1));
    mio::ABMSimulationNode<MockHistory> node2(MockHistory{}, t, std::move(world_2));

    mio::ABMMobilityEdge<MockHistory> edge({p2.get_id().get()}, {&mio::apply_commuting});
    edge.apply_mobility(node1, node2, t);

    EXPECT_EQ(node2.get_simulation().get_world().get_number_persons(work_id_2), 1);
    EXPECT_EQ(node1.get_simulation().get_world().get_activeness_statuses()[p1.get_id().get()], true);
    EXPECT_EQ(node1.get_simulation().get_world().get_activeness_statuses()[p2.get_id().get()], false);
    EXPECT_EQ(node2.get_simulation().get_world().get_activeness_statuses()[p1.get_id().get()], false);
    EXPECT_EQ(node2.get_simulation().get_world().get_activeness_statuses()[p2.get_id().get()], true);

    //return home
    t += mio::abm::hours(12);
    edge.apply_mobility(node1, node2, t);
    EXPECT_EQ(node2.get_simulation().get_world().get_number_persons(work_id_2), 0);
    EXPECT_EQ(node1.get_simulation().get_world().get_activeness_statuses()[p1.get_id().get()], true);
    EXPECT_EQ(node1.get_simulation().get_world().get_activeness_statuses()[p2.get_id().get()], true);
    EXPECT_EQ(node2.get_simulation().get_world().get_activeness_statuses()[p1.get_id().get()], false);
    EXPECT_EQ(node2.get_simulation().get_world().get_activeness_statuses()[p2.get_id().get()], false);
}