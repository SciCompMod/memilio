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

namespace
{

struct MoveOnly {
    MoveOnly();
    MoveOnly(const MoveOnly&) = delete;
    MoveOnly& operator=(const MoveOnly&) = delete;
    MoveOnly(MoveOnly&&)                 = default;
    MoveOnly& operator=(MoveOnly&&) = default;
};
using MoveOnlyGraph = epi::Graph<MoveOnly, MoveOnly>;

template <class G>
using add_node_expr_t = decltype(std::declval<G>().add_node(int()));
template <class G>
using add_edge_expr_t = decltype(std::declval<G>().add_edge(size_t(), size_t()));

} // namespace

static_assert(std::is_constructible<MoveOnlyGraph>::value, "Graph should support move-only node and edge properties.");
static_assert(std::is_move_constructible<MoveOnlyGraph>::value && std::is_move_assignable<MoveOnlyGraph>::value,
              "Graph should support move-only node and edge properties.");
static_assert(epi::is_expression_valid<add_node_expr_t, MoveOnlyGraph>::value,
              "Graph should support move-only node and edge properties.");
static_assert(epi::is_expression_valid<add_edge_expr_t, MoveOnlyGraph>::value,
              "Graph should support move-only node and edge properties.");
