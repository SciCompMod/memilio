#include <epidemiology/graph.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

TEST(TestGraph, creation)
{
    epi::Graph<int, int> g;
    g.add_node(6);
    g.add_node(4);
    g.add_node(8);
    g.add_edge(0, 1, 1);
    g.add_edge(2, 1, 3);
    g.add_edge(1, 2, 2);

    EXPECT_THAT(g.nodes(), testing::ElementsAre(6, 4, 8));

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
    
    EXPECT_EQ(g.edges()[1], (epi::Edge<int>{ 1, 2, 3 }));
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
    
    std::vector<epi::Edge<int>> v0 = { { 0, 1, 0 }, { 0, 2, 2 } };
    EXPECT_THAT(g.out_edges(0), testing::ElementsAreArray(v0));

    std::vector<epi::Edge<int>> v1 = { { 1, 2, 1 } };
    EXPECT_THAT(g.out_edges(1), testing::ElementsAreArray(v1));
}

class MockNodeFunc
{
public:
    MOCK_METHOD(void, invoke, (double t, double dt, int i), ());
    void operator()(double t, double dt, int& i) { invoke(t, dt, i); };
};

class MockEdgeFunc
{
public:
    MOCK_METHOD(void, invoke, (double t, double dt, int e, int n1, int n2), ());
    void operator()(double t, double dt, int& e, int& n1, int& n2) { invoke(t, dt, e, n1, n2); };
};

TEST(TestGraph, simulate)
{
    using testing::_;

    epi::Graph<int, int> g;
    g.add_node(0);
    g.add_node(1);
    g.add_node(2);
    g.add_node(3);
    g.add_edge(0, 1, 0);
    g.add_edge(1, 2, 1);
    g.add_edge(0, 2, 2);
    g.add_edge(3, 0, 3);

    MockEdgeFunc edge_func;
    MockNodeFunc node_func;

    const auto t0 = 1;
    const auto tmax = 3;
    const auto dt = 1;

    for (int t = t0; t < tmax; ++t)
    {
        EXPECT_CALL(edge_func, invoke(t + 1, 1, 0, 0, 1)).Times(1);
        EXPECT_CALL(edge_func, invoke(t + 1, 1, 1, 1, 2)).Times(1);
        EXPECT_CALL(edge_func, invoke(t + 1, 1, 2, 0, 2)).Times(1);
        EXPECT_CALL(edge_func, invoke(t + 1, 1, 3, 3, 0)).Times(1);

        EXPECT_CALL(node_func, invoke(t, 1, 0)).Times(1);
        EXPECT_CALL(node_func, invoke(t, 1, 1)).Times(1);
        EXPECT_CALL(node_func, invoke(t, 1, 2)).Times(1);
        EXPECT_CALL(node_func, invoke(t, 1, 3)).Times(1);
    }

    simulate_graph(t0, tmax, dt, g, node_func, edge_func);
}

TEST(TestGraph, persistentChangesDuringSimulation)
{
    epi::Graph<int, int> g;
    g.add_node(6);
    g.add_node(4);
    g.add_node(8);
    g.add_edge(0, 1, 1);
    g.add_edge(0, 2, 2);
    g.add_edge(1, 2, 3);

    auto node_func = [] (auto&& t, auto&& dt, auto&& n) { ++n; };
    auto edge_func = [] (auto&& t, auto&& dt, auto&& e, auto&& n1, auto&& n2) { ++e; ++n2; };

    int num_steps = 2;
    simulate_graph(0, num_steps, 1, g, node_func, edge_func);
    
    EXPECT_THAT(g.nodes(), testing::ElementsAre(6 + num_steps, 4 + num_steps + num_steps, 8 + num_steps + 2 * num_steps));
    std::vector<epi::Edge<int>> v = {{0, 1, 1 + num_steps}, {0, 2, 2 + num_steps}, {1, 2, 3 + num_steps}};
    EXPECT_THAT(g.edges(), testing::ElementsAreArray(v));
}