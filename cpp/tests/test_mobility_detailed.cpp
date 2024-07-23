/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Henrik Zunker
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
#include "memilio/mobility/metapopulation_mobility_detailed.h"
#include "matchers.h"

#include <gtest/gtest.h>

// Mock classes for Graph, Node, and Edge
struct NodeProperty {
    double stay_duration;
};

struct EdgeProperty {
    double travel_time;
    std::vector<size_t> path;
};

struct Node {
    NodeProperty property;
};

struct Edge {
    size_t start_node_idx;
    size_t end_node_idx;
    EdgeProperty property;
};

class MockGraph
{
public:
    const std::vector<Node>& nodes() const
    {
        return nodes_;
    }

    const std::vector<Edge>& edges() const
    {
        return edges_;
    }

    void add_node(const Node& node)
    {
        nodes_.push_back(node);
    }

    void add_edge(const Edge& edge)
    {
        edges_.push_back(edge);
    }

private:
    std::vector<Node> nodes_;
    std::vector<Edge> edges_;
};

class ScheduleManagerTest : public ::testing::Test
{
protected:
    ScheduleManagerTest()
        : manager(100)
    {
    }

    mio::ScheduleManager manager;
};

// Test case for computing the schedule with a simple graph
TEST_F(ScheduleManagerTest, ComputeSimpleSchedule)
{
    MockGraph graph;

    Node node1 = {NodeProperty{0.4}};
    Node node2 = {NodeProperty{0.4}};
    graph.add_node(node1);
    graph.add_node(node2);

    Edge edge1 = {0, 1, EdgeProperty{0.1, {0, 1}}};
    Edge edge2 = {1, 0, EdgeProperty{0.1, {1, 0}}};
    graph.add_edge(edge1);
    graph.add_edge(edge2);

    auto schedule = manager.compute_schedule(graph);

    // Check that the schedule edges and mobility schedule edges are correctly filled
    EXPECT_EQ(schedule.schedule_edges.size(), 2);
    EXPECT_EQ(schedule.mobility_schedule_edges.size(), 2);

    // Check the contents of the schedule
    for (const auto& edge_schedule : schedule.schedule_edges) {
        EXPECT_EQ(edge_schedule.size(), 100); // timesteps
    }

    for (const auto& mobility_schedule : schedule.mobility_schedule_edges) {
        EXPECT_EQ(mobility_schedule.size(), 100); // timesteps
    }
}

// Test case for computing the schedule with no edges
TEST_F(ScheduleManagerTest, ComputeEmptySchedule)
{
    MockGraph graph;

    Node node1 = {NodeProperty{1.0}};
    Node node2 = {NodeProperty{1.0}};
    graph.add_node(node1);
    graph.add_node(node2);

    auto schedule = manager.compute_schedule(graph);

    // Check that the schedule edges and mobility schedule edges are empty
    EXPECT_TRUE(schedule.schedule_edges.empty());
    EXPECT_TRUE(schedule.mobility_schedule_edges.empty());
}

// Test case for verifying first_mobility vector
TEST_F(ScheduleManagerTest, VerifyFirstMobility)
{
    MockGraph graph;

    Node node1 = {NodeProperty{0.5}};
    Node node2 = {NodeProperty{0.5}};
    graph.add_node(node1);
    graph.add_node(node2);

    Edge edge1 = {0, 1, EdgeProperty{0.2, {0, 1}}};
    graph.add_edge(edge1);

    auto schedule = manager.compute_schedule(graph);

    // First mobility should be at timestep 10
    EXPECT_EQ(schedule.first_mobility.size(), 1);
    EXPECT_EQ(schedule.first_mobility[0], 10);
}

// Test case for verifying local and mobility integration schedules
TEST_F(ScheduleManagerTest, VerifyIntegrationSchedules)
{
    MockGraph graph;

    Node node1 = {NodeProperty{0.4}};
    Node node2 = {NodeProperty{0.4}};
    graph.add_node(node1);
    graph.add_node(node2);

    Edge edge1 = {0, 1, EdgeProperty{0.1, {0, 1}}};
    graph.add_edge(edge1);

    auto schedule = manager.compute_schedule(graph);

    // Check that local and mobility integration schedules are filled
    EXPECT_EQ(schedule.local_int_schedule.size(), 2);
    EXPECT_EQ(schedule.mobility_int_schedule.size(), 2);

    // the local_int_schedule for the first node should contain only the start, end time and the start time of the mobility.
    EXPECT_EQ(schedule.local_int_schedule[0].size(), 3);
    EXPECT_EQ(schedule.local_int_schedule[0][0], 0);
    EXPECT_EQ(schedule.local_int_schedule[0][1], 40);
    EXPECT_EQ(schedule.local_int_schedule[0][2], 99);

    // the second node should have the start and end time and additionally the arrival and departure time of the
    // individual from the first node.
    EXPECT_EQ(schedule.local_int_schedule[1].size(), 4);
    EXPECT_EQ(schedule.local_int_schedule[1][0], 0);
    EXPECT_EQ(schedule.local_int_schedule[1][1], 50);
    EXPECT_EQ(schedule.local_int_schedule[1][2], 90);
    EXPECT_EQ(schedule.local_int_schedule[1][3], 99);

    // since the travel time is 0.1, the mobility_int_schedule for the first node should contain the start time of the mobility,
    // the time where the individual arrives at the second mobility model and the time where the individual arrives again before getting back.
    // the last time step isnt considered since we do this by default at the end of the day.
    EXPECT_EQ(schedule.mobility_int_schedule[0].size(), 3);
    EXPECT_EQ(schedule.mobility_int_schedule[0][0], 40);
    EXPECT_EQ(schedule.mobility_int_schedule[0][1], 45);
    EXPECT_EQ(schedule.mobility_int_schedule[0][2], 95);

    // for the second node, all times should be shifted by 5 and add the time point where the individual start their home trip
    EXPECT_EQ(schedule.mobility_int_schedule[1].size(), 4);
    EXPECT_EQ(schedule.mobility_int_schedule[1][0], 45);
    EXPECT_EQ(schedule.mobility_int_schedule[1][1], 50);
    EXPECT_EQ(schedule.mobility_int_schedule[1][2], 90);
    EXPECT_EQ(schedule.mobility_int_schedule[1][3], 95);
}

// Test case where travel time per region lower than allowed minimum of 0.01
TEST_F(ScheduleManagerTest, MinimalTravelTimePerRegion)
{
    MockGraph graph;

    Node node1 = {NodeProperty{0.4}};
    Node node2 = {NodeProperty{0.4}};
    graph.add_node(node1);
    graph.add_node(node2);

    Edge edge1 = {0, 1, EdgeProperty{0.001, {0, 1}}}; // Minimal travel time per region
    graph.add_edge(edge1);

    auto schedule = manager.compute_schedule(graph);

    // Check that the schedule edges and mobility schedule edges are correctly filled
    EXPECT_EQ(schedule.schedule_edges.size(), 1);
    EXPECT_EQ(schedule.mobility_schedule_edges.size(), 1);

    // Check that the travel time was correctly set to the minimum value
    EXPECT_EQ(schedule.mobility_schedule_edges[0][56], true);
    EXPECT_EQ(schedule.mobility_schedule_edges[0][57], true);
    EXPECT_EQ(schedule.mobility_schedule_edges[0][98], true);
    EXPECT_EQ(schedule.mobility_schedule_edges[0][99], true);

    // check that the number of true values is 4
    size_t count = 0;
    for (const auto& mobility_schedule : schedule.mobility_schedule_edges) {
        for (const auto& mobility : mobility_schedule) {
            if (mobility) {
                count++;
            }
        }
    }
    EXPECT_EQ(count, 4);
}