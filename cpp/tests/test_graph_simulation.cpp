/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/mobility/metapopulation_mobility_stochastic.h"
#include "memilio/compartments/simulation.h"
#include "ode_seir/model.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

class MockNodeFunc
{
public:
    MOCK_METHOD(void, invoke, (double t, double dt, int& i), ());
    void operator()(double t, double dt, int& i)
    {
        invoke(t, dt, i);
    };
};

class MockEdgeFunc
{
public:
    MOCK_METHOD(void, invoke, (double t, double dt, int& e, int& n1, int& n2), ());
    void operator()(double t, double dt, int& e, int& n1, int& n2)
    {
        invoke(t, dt, e, n1, n2);
    };
};

TEST(TestGraphSimulation, simulate)
{
    using testing::_;
    using testing::Eq;

    mio::Graph<int, int> g;
    g.add_node(6, 0);
    g.add_node(8, 1);
    g.add_node(4, 2);
    g.add_node(2, 3);
    g.add_edge(0, 1, 0);
    g.add_edge(1, 2, 1);
    g.add_edge(0, 2, 2);
    g.add_edge(3, 0, 3);

    MockEdgeFunc edge_func;
    MockNodeFunc node_func;

    const auto t0   = 1;
    const auto tmax = 3.0;
    const auto dt   = 1.0;

    testing::ExpectationSet node_func_calls;

    node_func_calls += EXPECT_CALL(node_func, invoke(1, 1, Eq(0))).Times(1);
    node_func_calls += EXPECT_CALL(node_func, invoke(1, 1, Eq(1))).Times(1);
    node_func_calls += EXPECT_CALL(node_func, invoke(1, 1, Eq(2))).Times(1);
    node_func_calls += EXPECT_CALL(node_func, invoke(1, 1, Eq(3))).Times(1);

    EXPECT_CALL(edge_func, invoke(2, 1, Eq(0), Eq(0), Eq(1))).Times(1).After(node_func_calls);
    EXPECT_CALL(edge_func, invoke(2, 1, Eq(2), Eq(0), Eq(2))).Times(1).After(node_func_calls);
    EXPECT_CALL(edge_func, invoke(2, 1, Eq(1), Eq(1), Eq(2))).Times(1).After(node_func_calls);
    EXPECT_CALL(edge_func, invoke(2, 1, Eq(3), Eq(3), Eq(0))).Times(1).After(node_func_calls);

    node_func_calls += EXPECT_CALL(node_func, invoke(2, 1, Eq(0))).Times(1);
    node_func_calls += EXPECT_CALL(node_func, invoke(2, 1, Eq(1))).Times(1);
    node_func_calls += EXPECT_CALL(node_func, invoke(2, 1, Eq(2))).Times(1);
    node_func_calls += EXPECT_CALL(node_func, invoke(2, 1, Eq(3))).Times(1);

    EXPECT_CALL(edge_func, invoke(3, 1, Eq(0), Eq(0), Eq(1))).Times(1).After(node_func_calls);
    EXPECT_CALL(edge_func, invoke(3, 1, Eq(2), Eq(0), Eq(2))).Times(1).After(node_func_calls);
    EXPECT_CALL(edge_func, invoke(3, 1, Eq(1), Eq(1), Eq(2))).Times(1).After(node_func_calls);
    EXPECT_CALL(edge_func, invoke(3, 1, Eq(3), Eq(3), Eq(0))).Times(1).After(node_func_calls);

    auto sim = mio::make_graph_sim(
        t0, dt, g,
        [&node_func](auto&& t, auto&& dt_, auto&& n) {
            node_func(t, dt_, n);
        },
        [&edge_func](auto&& t, auto&& dt_, auto&& e, auto&& n1, auto&& n2) {
            edge_func(t, dt_, e, n1, n2);
        });

    sim.advance(tmax);

    EXPECT_NEAR(sim.get_t(), tmax, 1e-15);
}

TEST(TestGraphSimulation, stopsAtTmax)
{
    using testing::_;
    using testing::Eq;

    mio::Graph<int, int> g;
    g.add_node(6, 0);
    g.add_node(8, 1);
    g.add_edge(0, 1, 0);

    const auto t0   = 1.0;
    const auto tmax = 3.123;
    const auto dt   = 0.076;

    auto sim = mio::make_graph_sim(
        t0, dt, g, [](auto&&, auto&&, auto&&) {}, [](auto&&, auto&&, auto&&, auto&&, auto&&) {});

    sim.advance(tmax);

    EXPECT_NEAR(sim.get_t(), tmax, 1e-15);
}

TEST(TestGraphSimulation, stopsAtTmaxStochastic)
{
    using testing::_;
    using testing::Eq;

    const auto t0   = 1.0;
    const auto tmax = 5.;
    const auto dt   = 0.076;

    mio::oseir::Model model;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Susceptible)}] = 0.9;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}]     = 0.1;
    model.populations.set_total(1000);

    mio::Graph<mio::SimulationNode<mio::Simulation<mio::oseir::Model>>, mio::MigrationEdgeStochastic> g;
    g.add_node(0, model, t0);
    g.add_node(1, model, t0);
    g.add_edge(0, 1, Eigen::VectorXd::Constant(4, 0.001));

    auto sim = mio::make_migration_sim(t0, dt, std::move(g));

    sim.advance(tmax);

    EXPECT_NEAR(sim.get_t(), tmax, 1e-15);
}

TEST(TestGraphSimulation, persistentChangesDuringSimulation)
{
    mio::Graph<int, int> g;
    g.add_node(0, 6);
    g.add_node(1, 4);
    g.add_node(2, 8);
    g.add_edge(0, 1, 1);
    g.add_edge(0, 2, 2);
    g.add_edge(1, 2, 3);

    auto node_func = [](auto&& /*t*/, auto&& /*dt*/, auto&& n) {
        ++n;
    };
    auto edge_func = [](auto&& /*t*/, auto&& /*dt*/, auto&& e, auto&& /*n1*/, auto&& n2) {
        ++e;
        ++n2;
    };

    auto t0       = 0;
    auto dt       = 1;
    auto sim      = mio::make_graph_sim(t0, dt, g, node_func, edge_func);
    int num_steps = 2;
    sim.advance(t0 + num_steps * dt);

    std::vector<mio::Node<int>> v_n = {{0, 6 + num_steps}, {1, 4 + 2 * num_steps}, {2, 8 + 3 * num_steps}};
    EXPECT_THAT(sim.get_graph().nodes(), testing::ElementsAreArray(v_n));
    std::vector<mio::Edge<int>> v_e = {{0, 1, 1 + num_steps}, {0, 2, 2 + num_steps}, {1, 2, 3 + num_steps}};
    EXPECT_THAT(sim.get_graph().edges(), testing::ElementsAreArray(v_e));
}

TEST(TestGraphSimulation, consistencyStochasticMobility)
{
    using testing::_;
    using testing::Eq;

    //set seeds
    mio::thread_local_rng().seed({114381446, 2427727386, 806223567, 832414962, 4121923627, 1581162203});

    const auto t0   = 0.0;
    const auto tmax = 20.;
    const auto dt   = 0.076;

    mio::oseir::Model model;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Susceptible)}] = 0.7;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}]     = 0.3;
    model.populations.set_total(1000);

    mio::Graph<mio::SimulationNode<mio::Simulation<mio::oseir::Model>>, mio::MigrationEdgeStochastic> g;
    g.add_node(0, model, t0);
    g.add_node(1, model, t0);
    g.add_edge(0, 1, Eigen::VectorXd::Constant(4, 0.001));

    auto sim = mio::make_migration_sim(t0, dt, std::move(g));

    sim.advance(tmax);

    auto result_n0 = sim.get_graph().nodes()[0].property.get_result().get_last_value();
    auto result_n1 = sim.get_graph().nodes()[1].property.get_result().get_last_value();

    auto expected_values_n0 = std::vector<double>{691.0, 6.3518514260209971, 31.303976182729517, 257.34417239124963};
    auto actual_values_n0   = std::vector<double>{result_n0[0], result_n0[1], result_n0[2], result_n0[3]};
    auto expected_values_n1 = std::vector<double>{709.0, 6.4651647420642382, 33.101208735481720, 265.43362652245412};
    auto actual_values_n1   = std::vector<double>{result_n1[0], result_n1[1], result_n1[2], result_n1[3]};

    for (size_t i = 0; i < expected_values_n0.size(); ++i) {
        EXPECT_THAT(expected_values_n0[i], testing::DoubleNear(actual_values_n0[i], 1e-9));
        EXPECT_THAT(expected_values_n1[i], testing::DoubleNear(actual_values_n1[i], 1e-9));
    }
}

namespace
{

struct MoveOnly {
    MoveOnly();
    MoveOnly(const MoveOnly&)            = delete;
    MoveOnly& operator=(const MoveOnly&) = delete;
    MoveOnly(MoveOnly&&)                 = default;
    MoveOnly& operator=(MoveOnly&&)      = default;
};
using MoveOnlyGraph    = mio::Graph<MoveOnly, MoveOnly>;
using MoveOnlyGraphSim = mio::GraphSimulation<MoveOnlyGraph>;

} // namespace

static_assert(std::is_constructible<MoveOnlyGraphSim, double, double, MoveOnlyGraph&&, MoveOnlyGraphSim::node_function,
                                    MoveOnlyGraphSim::edge_function>::value,
              "GraphSimulation should support move-only graphs.");
static_assert(std::is_move_constructible<MoveOnlyGraphSim>::value && std::is_move_assignable<MoveOnlyGraphSim>::value,
              "GraphSimulation should support move-only graphs.");
