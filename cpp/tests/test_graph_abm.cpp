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

#include "abm/model.h"
#include "abm/location_type.h"
#include "abm/time.h"
#include "graph_abm/graph_abm_mobility.h"
#include "graph_abm/mobility_rules.h"
#include "memilio/epidemiology/age_group.h"
#include <cstddef>
#include <gtest/gtest.h>

TEST(TestGraphAbm, test_activessness)
{
    auto model                                                           = mio::abm::Model(size_t(1));
    model.parameters.get<mio::abm::AgeGroupGotoWork>()[mio::AgeGroup(0)] = true;
    auto work_id = model.add_location(mio::abm::LocationType::Work);
    auto home_id = model.add_location(mio::abm::LocationType::Home);
    auto& home   = model.get_location(home_id);
    auto& work   = model.get_location(work_id);
    auto p1_id   = model.add_person(home_id, mio::AgeGroup(0));
    auto p2_id   = model.add_person(home_id, mio::AgeGroup(0));
    auto& p1     = model.get_person(p1_id);
    auto& p2     = model.get_person(p2_id);
    p1.set_assigned_location(work.get_type(), work.get_id(), work.get_model_id());
    p2.set_assigned_location(work.get_type(), work.get_id(), work.get_model_id());
    p1.set_assigned_location(home.get_type(), home.get_id(), home.get_model_id());
    p2.set_assigned_location(home.get_type(), home.get_id(), home.get_model_id());

    EXPECT_EQ(p1.get_location(), home_id);
    EXPECT_EQ(p2.get_location(), home_id);
    EXPECT_EQ(p1.get_location_model_id(), model.get_id());
    EXPECT_EQ(p2.get_location_model_id(), model.get_id());
    EXPECT_EQ(model.get_activeness_statuses().size(), 2);

    model.change_activeness(p1.get_id());
    EXPECT_EQ(model.get_activeness_statuses()[p1.get_id().get()], false);
    EXPECT_EQ(model.get_activeness_statuses()[p2.get_id().get()], true);

    auto t  = mio::abm::TimePoint(0) + mio::abm::hours(6);
    auto dt = mio::abm::hours(3);

    model.evolve(t, dt);

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
    auto model_1                                                           = mio::abm::Model(size_t(1), 1);
    auto model_2                                                           = mio::abm::Model(size_t(1), 2);
    model_1.parameters.get<mio::abm::AgeGroupGotoWork>()[mio::AgeGroup(0)] = true;
    model_2.parameters.get<mio::abm::AgeGroupGotoWork>()[mio::AgeGroup(0)] = true;
    auto work_id_1 = model_1.add_location(mio::abm::LocationType::Work);
    auto home_id   = model_1.add_location(mio::abm::LocationType::Home);
    auto work_id_2 = model_2.add_location(mio::abm::LocationType::Work);
    auto& work_1   = model_1.get_location(work_id_1);
    auto& work_2   = model_2.get_location(work_id_2);
    auto& home     = model_1.get_location(home_id);

    EXPECT_EQ(work_1.get_model_id(), 1);
    EXPECT_EQ(work_2.get_model_id(), 2);

    auto p1_id = model_1.add_person(home_id, mio::AgeGroup(0));
    auto p2_id = model_1.add_person(home_id, mio::AgeGroup(0));
    auto& p1   = model_1.get_person(p1_id);
    auto& p2   = model_1.get_person(p2_id);
    p1.set_assigned_location(work_1.get_type(), work_1.get_id(), work_1.get_model_id());
    p2.set_assigned_location(work_2.get_type(), work_2.get_id(), work_2.get_model_id());
    p1.set_assigned_location(home.get_type(), home.get_id(), home.get_model_id());
    p2.set_assigned_location(home.get_type(), home.get_id(), home.get_model_id());

    //copy persons to model 2
    model_2.copy_persons_from_other_model(model_1);

    auto t = mio::abm::TimePoint(0) + mio::abm::hours(6);
    mio::ABMSimulationNode<MockHistory> node1(MockHistory{}, t, std::move(model_1));
    mio::ABMSimulationNode<MockHistory> node2(MockHistory{}, t, std::move(model_2));

    mio::ABMMobilityEdge<MockHistory> edge({p2.get_id().get()}, {&mio::apply_commuting});
    edge.apply_mobility(node1, node2, t);

    EXPECT_EQ(node2.get_simulation().get_model().get_number_persons(work_id_2), 1);
    EXPECT_EQ(node1.get_simulation().get_model().get_activeness_statuses()[p1.get_id().get()], true);
    EXPECT_EQ(node1.get_simulation().get_model().get_activeness_statuses()[p2.get_id().get()], false);
    EXPECT_EQ(node2.get_simulation().get_model().get_activeness_statuses()[p1.get_id().get()], false);
    EXPECT_EQ(node2.get_simulation().get_model().get_activeness_statuses()[p2.get_id().get()], true);

    //return home
    t += mio::abm::hours(12);
    edge.apply_mobility(node1, node2, t);
    EXPECT_EQ(node2.get_simulation().get_model().get_number_persons(work_id_2), 0);
    EXPECT_EQ(node1.get_simulation().get_model().get_activeness_statuses()[p1.get_id().get()], true);
    EXPECT_EQ(node1.get_simulation().get_model().get_activeness_statuses()[p2.get_id().get()], true);
    EXPECT_EQ(node2.get_simulation().get_model().get_activeness_statuses()[p1.get_id().get()], false);
    EXPECT_EQ(node2.get_simulation().get_model().get_activeness_statuses()[p2.get_id().get()], false);
}