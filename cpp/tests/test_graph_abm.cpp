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

#include "abm/location.h"
#include "abm/location_id.h"
#include "abm/mobility_data.h"
#include "abm/model.h"
#include "graph_abm/model_wrapper.h"
#include "abm/location_type.h"
#include "abm/time.h"
#include "graph_abm/graph_abm_mobility.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/miompi.h"
#include "memilio/mobility/graph.h"
#include <algorithm>
#include <cstddef>
#include <gtest/gtest.h>

struct MockHistory {

    template <class T>
    void log(const T& t)
    {
        mio::unused(t);
    }
};

TEST(TestGraphAbm, test_evolve_node)
{
    auto t                                                               = mio::abm::TimePoint(0);
    auto dt                                                              = mio::abm::hours(10);
    auto model                                                           = mio::ModelWrapper(size_t(1), 1);
    model.parameters.get<mio::abm::AgeGroupGotoWork>()[mio::AgeGroup(0)] = true;
    auto home_id = model.add_location(mio::abm::LocationType::Home);
    auto& home   = model.get_location(home_id);
    auto work    = mio::abm::Location(mio::abm::LocationType::Work, mio::abm::LocationId(0), size_t(1), 1);
    auto pid     = model.add_person(home_id, mio::AgeGroup(0));
    auto index   = model.get_person_index(pid);
    auto& p      = model.get_person(pid);
    p.set_assigned_location(home.get_type(), home.get_id(), home.get_model_id());
    p.set_assigned_location(work.get_type(), work.get_id(), 2);
    mio::ABMSimulationNode<MockHistory> node(MockHistory{}, t, std::move(model));
    EXPECT_EQ(node.get_simulation().get_model().get_activeness_statuses()[index], true);
    node.evolve(t, dt);

    EXPECT_EQ(node.get_simulation().get_time(), mio::abm::TimePoint(dt.seconds()));
    EXPECT_EQ(node.get_simulation().get_model().get_activeness_statuses()[index], false);
    EXPECT_EQ(node.get_simulation().get_model().get_person_buffer().size(), 1);
    EXPECT_EQ(node.get_simulation().get_model().get_person_buffer()[0], index);
}

TEST(TestGraphAbm, test_apply_mobility)
{
    auto model1                                                           = mio::ModelWrapper(size_t(2), 1);
    auto model2                                                           = mio::ModelWrapper(size_t(2), 2);
    model1.parameters.get<mio::abm::AgeGroupGotoWork>()[mio::AgeGroup(0)] = true;
    model2.parameters.get<mio::abm::AgeGroupGotoWork>()[mio::AgeGroup(0)] = true;
    auto work_id_1  = model1.add_location(mio::abm::LocationType::Work);
    auto home_id    = model1.add_location(mio::abm::LocationType::Home);
    auto work_id_2  = model2.add_location(mio::abm::LocationType::Work);
    auto event_id_1 = model1.add_location(mio::abm::LocationType::SocialEvent);
    auto event_id_2 = model2.add_location(mio::abm::LocationType::SocialEvent);
    auto& work_1    = model1.get_location(work_id_1);
    auto& work_2    = model2.get_location(work_id_2);
    auto& home      = model1.get_location(home_id);
    auto& event_1   = model1.get_location(event_id_1);
    auto& event_2   = model2.get_location(event_id_2);

    EXPECT_EQ(work_1.get_model_id(), 1);
    EXPECT_EQ(work_2.get_model_id(), 2);

    auto p1_id    = model1.add_person(home_id, mio::AgeGroup(0));
    auto p2_id    = model1.add_person(home_id, mio::AgeGroup(0));
    auto p2_index = model1.get_person_index(p2_id);
    auto p3_id    = model1.add_person(home_id, mio::AgeGroup(1));
    auto p4_id    = model1.add_person(home_id, mio::AgeGroup(1));
    //auto p4_index = model1.get_person_index(p4_id);
    auto& p1 = model1.get_person(p1_id);
    auto& p2 = model1.get_person(p2_id);
    auto& p3 = model1.get_person(p3_id);
    auto& p4 = model1.get_person(p4_id);
    p1.set_assigned_location(work_1.get_type(), work_1.get_id(), work_1.get_model_id());
    p2.set_assigned_location(work_2.get_type(), work_2.get_id(), work_2.get_model_id());
    p1.set_assigned_location(home.get_type(), home.get_id(), home.get_model_id());
    p2.set_assigned_location(home.get_type(), home.get_id(), home.get_model_id());
    p3.set_assigned_location(home.get_type(), home.get_id(), home.get_model_id());
    p4.set_assigned_location(home.get_type(), home.get_id(), home.get_model_id());
    p3.set_assigned_location(event_1.get_type(), event_1.get_id(), event_1.get_model_id());
    p4.set_assigned_location(event_2.get_type(), event_2.get_id(), event_2.get_model_id());

    mio::abm::TripList& trips = model1.get_trip_list();
    mio::abm::Trip trip1(p3.get_id(), mio::abm::TimePoint(0) + mio::abm::hours(8), event_id_1, model1.get_id(), home_id,
                         model1.get_id(), mio::abm::TransportMode::Unknown, mio::abm::LocationType::SocialEvent);
    mio::abm::Trip trip2(p4.get_id(), mio::abm::TimePoint(0) + mio::abm::hours(8), event_id_2, model2.get_id(), home_id,
                         model1.get_id(), mio::abm::TransportMode::Unknown, mio::abm::LocationType::SocialEvent);

    trips.add_trip(trip1);
    trips.add_trip(trip2);

    auto t  = mio::abm::TimePoint(0);
    auto dt = mio::abm::hours(12);
    mio::ABMSimulationNode<MockHistory> node1(MockHistory{}, t, std::move(model1));
    mio::ABMSimulationNode<MockHistory> node2(MockHistory{}, t, std::move(model2));

    node1.evolve(t, dt);
    node2.evolve(t, dt);

    EXPECT_EQ(node2.get_simulation().get_model().get_persons().size(), 0);
    EXPECT_EQ(node1.get_simulation().get_model().get_persons().size(), 4);
    EXPECT_EQ(node1.get_simulation().get_model().get_activeness_statuses()[p2_index], false);
    //EXPECT_EQ(node1.get_simulation().get_model().get_activeness_statuses()[p4_index], false);

    // mio::ABMMobilityEdge<MockHistory> edge;
    // edge.apply_mobility(node1, node2, t);

    // EXPECT_EQ(node2.get_simulation().get_model().get_persons().size(), 2);
    // EXPECT_EQ(node2.get_simulation().get_model().get_persons().size(), 2);
    // EXPECT_EQ(node1.get_simulation().get_model().get_person_buffer().size(), 0);
}

TEST(TestGraphABM, test_graph_simulation)
{
    auto model1 = mio::ModelWrapper(size_t(1), 1);
    auto model2 = mio::ModelWrapper(size_t(1), 2);

    mio::abm::TimePoint t0   = mio::abm::TimePoint(0);
    mio::abm::TimePoint tmax = t0 + mio::abm::days(5);

    mio::Graph<mio::ABMSimulationNode<MockHistory>, mio::ABMMobilityEdge<MockHistory>> graph;
    graph.add_node(model1.get_id(), MockHistory{}, t0, std::move(model1));
    graph.add_node(model2.get_id(), MockHistory{}, t0, std::move(model2));
    graph.add_edge(0, 1);
    graph.add_edge(1, 0);

    auto sim = mio::make_abm_graph_sim<MockHistory>(t0, mio::abm::hours(12), std::move(graph));
    sim.advance(tmax);

    EXPECT_EQ(sim.get_t(), tmax);
}
