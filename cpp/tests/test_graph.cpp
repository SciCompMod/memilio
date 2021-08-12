/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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
#include <epidemiology/utils/graph.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

TEST(TestGraph, creation)
{
    epi::Graph<int, int> g;
    g.add_node(0, 6);
    g.add_node(1, 4);
    g.add_node(2, 8);
    g.add_edge(0, 1, 1);
    g.add_edge(2, 1, 3);
    g.add_edge(1, 2, 2);

    std::vector<epi::Node<int>> n = {{0, 6}, {1, 4}, {2, 8}};
    EXPECT_THAT(g.nodes(), testing::ElementsAreArray(n));

    std::vector<epi::Edge<int>> v = {{0, 1, 1}, {1, 2, 2}, {2, 1, 3}};
    EXPECT_THAT(g.edges(), testing::ElementsAreArray(v));
}

TEST(TestGraph, duplicate_edge)
{
    epi::Graph<int, int> g;
    g.add_node(6);
    g.add_node(4);
    g.add_node(8);
    g.add_edge(0, 1, 1);
    g.add_edge(2, 1, 3);
    g.add_edge(1, 2, 2);
    g.add_edge(1, 2, 3);

    EXPECT_EQ(g.edges()[1], (epi::Edge<int>{1, 2, 3}));
}

TEST(TestGraph, ot_edges)
{
    epi::Graph<int, int> g;
    g.add_node(0);
    g.add_node(1);
    g.add_node(2);
    g.add_node(3);
    g.add_edge(0, 1, 0);
    g.add_edge(1, 2, 1);
    g.add_edge(0, 2, 2);
    g.add_edge(3, 0, 3);

    std::vector<epi::Edge<int>> v0 = {{0, 1, 0}, {0, 2, 2}};
    EXPECT_THAT(g.out_edges(0), testing::ElementsAreArray(v0));

    std::vector<epi::Edge<int>> v1 = {{1, 2, 1}};
    EXPECT_THAT(g.out_edges(1), testing::ElementsAreArray(v1));
}
