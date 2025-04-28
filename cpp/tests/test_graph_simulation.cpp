/*
* Copyright (C) 2020-2025 MEmilio
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
#include "abm_helpers.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/mobility/metapopulation_mobility_stochastic.h"
#include "memilio/compartments/simulation.h"
#include "ode_seir/model.h"
#include "gtest/gtest.h"
#include "load_test_data.h"
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

    const auto t0   = 1.0;
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

    auto sim =
        mio::make_graph_sim(t0, dt, g, [](auto&&, auto&&, auto&&) {}, [](auto&&, auto&&, auto&&, auto&&, auto&&) {});

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

    mio::oseir::Model<double> model(1);
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 0.9;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = 0.1;
    model.populations.set_total(1000);

    mio::Graph<mio::SimulationNode<mio::Simulation<double, mio::oseir::Model<double>>>, mio::MobilityEdgeStochastic> g;
    g.add_node(0, model, t0);
    g.add_node(1, model, t0);
    g.add_edge(0, 1, Eigen::VectorXd::Constant(4, 0.001));

    auto sim = mio::make_mobility_sim(t0, dt, std::move(g));

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

    const auto t0   = 0.0;
    const auto tmax = 10.;
    const auto dt   = 0.076;

    mio::oseir::Model<double> model(1);
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 0.7;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = 0.3;
    model.populations.set_total(1000);

    mio::Graph<mio::SimulationNode<mio::Simulation<double, mio::oseir::Model<double>>>, mio::MobilityEdgeStochastic> g;
    g.add_node(0, model, t0);
    g.add_node(1, model, t0);
    g.add_edge(0, 1, Eigen::VectorXd::Constant(4, 0.001));

    auto sim = mio::make_mobility_sim(t0, dt, std::move(g));

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<ScalarType>>>>
        mock_exponential_dist;
    // use pregenerated exp(1) random values
    // all values are used to set normalized_waiting_time in GraphSimulationStochastic<...>::advance,
    // the first value is used at the function start, all others later during the while loop
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke)
        .Times(testing::Exactly(10))
        .WillOnce(testing::Return(0.446415))
        .WillOnce(testing::Return(1.04048))
        .WillOnce(testing::Return(0.136687))
        .WillOnce(testing::Return(2.50697))
        .WillOnce(testing::Return(1.61943))
        .WillOnce(testing::Return(0.267578))
        .WillOnce(testing::Return(1.03696))
        .WillOnce(testing::Return(0.58395))
        .WillOnce(testing::Return(0.113943))
        .WillOnce(testing::Return(1.204045));

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    // these values determine which transition event should occur in GraphSimulationStochastic<...>::advance
    // during this short sim, the chance of event==0 is ~70% every time
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke)
        .Times(testing::Exactly(9))
        .WillOnce(testing::Return(0))
        .WillOnce(testing::Return(1))
        .WillRepeatedly(testing::Return(0));

    sim.advance(tmax);

    auto result_n0 = sim.get_graph().nodes()[0].property.get_result().get_last_value();
    auto result_n1 = sim.get_graph().nodes()[1].property.get_result().get_last_value();

    auto expected_values_n0 = std::vector<double>{692.0, 43.630772796677256, 95.750528156188381, 159.61869904713436};
    auto actual_values_n0   = std::vector<double>{result_n0[0], result_n0[1], result_n0[2], result_n0[3]};
    auto expected_values_n1 = std::vector<double>{708.0, 44.063147085799322, 96.485223892060375, 160.45162902214025};
    auto actual_values_n1   = std::vector<double>{result_n1[0], result_n1[1], result_n1[2], result_n1[3]};

    for (size_t i = 0; i < expected_values_n0.size(); ++i) {
        EXPECT_THAT(expected_values_n0[i], testing::DoubleNear(actual_values_n0[i], 1e-7));
        EXPECT_THAT(expected_values_n1[i], testing::DoubleNear(actual_values_n1[i], 1e-7));
    }
}

template <typename Graph>
mio::GraphSimulation<Graph> create_simulation(Graph&& g, mio::oseir::Model<double>& model, double t0, double tmax,
                                              double dt)
{
    g.add_node(0, model, t0);
    g.add_node(1, model, t0);
    g.add_node(2, model, t0);
    for (size_t county_idx_i = 0; county_idx_i < g.nodes().size(); ++county_idx_i) {
        for (size_t county_idx_j = 0; county_idx_j < g.nodes().size(); ++county_idx_j) {
            if (county_idx_i == county_idx_j)
                continue;
            g.add_edge(county_idx_i, county_idx_j, Eigen::VectorXd::Constant(4, 0.001));
        }
    }

    auto sim = mio::make_mobility_sim(t0, dt, std::move(g));

    sim.advance(tmax);

    return sim;
}

TEST(TestGraphSimulation, consistencyFlowMobility)
{
    double t0   = 0;
    double tmax = 1;
    double dt   = 0.001;

    mio::oseir::Model<double> model(1);
    double total_population                                                      = 10000;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]   = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]  = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}] = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] =
        total_population - model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}];
    model.parameters.set<mio::oseir::TimeExposed<double>>(5.2);
    model.parameters.set<mio::oseir::TimeInfected<double>>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(0.04);
    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::oseir::ContactPatterns<double>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(10);

    model.check_constraints();

    auto sim_no_flows =
        create_simulation(mio::Graph<mio::SimulationNode<mio::Simulation<double, mio::oseir::Model<double>>>,
                                     mio::MobilityEdge<double>>(),
                          model, t0, tmax, dt);

    auto sim_flows =
        create_simulation(mio::Graph<mio::SimulationNode<mio::FlowSimulation<double, mio::oseir::Model<double>>>,
                                     mio::MobilityEdge<double>>(),
                          model, t0, tmax, dt);

    //test if all results of both simulations are equal for all nodes
    for (size_t node_id = 0; node_id < sim_no_flows.get_graph().nodes().size(); ++node_id) {
        auto& results_no_flows = sim_no_flows.get_graph().nodes()[node_id].property.get_result();
        auto& results_flows    = sim_flows.get_graph().nodes()[node_id].property.get_result();
        EXPECT_EQ((size_t)results_no_flows.get_num_time_points(), (size_t)results_flows.get_num_time_points());
        for (size_t t_indx = 0; t_indx < (size_t)results_no_flows.get_num_time_points(); t_indx++) {
            EXPECT_NEAR(results_no_flows.get_time((Eigen::Index)t_indx), results_flows.get_time((Eigen::Index)t_indx),
                        1e-10);
            auto tmp_sol_no_flows = results_no_flows.get_value((Eigen::Index)t_indx);
            auto tmp_sol_flows    = results_flows.get_value((Eigen::Index)t_indx);
            EXPECT_NEAR(tmp_sol_no_flows[0], tmp_sol_flows[0], 1e-10);
            EXPECT_NEAR(tmp_sol_no_flows[1], tmp_sol_flows[1], 1e-10);
            EXPECT_NEAR(tmp_sol_no_flows[2], tmp_sol_flows[2], 1e-10);
        }
    }

    // test all values from one node to the provided reference data for both simulations
    const auto& res_sim = sim_flows.get_graph().nodes()[0].property.get_result();
    const auto compare  = load_test_data_csv<ScalarType>("graphsimulation-compare.csv");
    EXPECT_EQ((size_t)compare.size(), (size_t)res_sim.get_num_time_points());
    for (size_t t_indx = 0; t_indx < (size_t)res_sim.get_num_time_points(); t_indx++) {
        EXPECT_NEAR(compare[t_indx][0], res_sim.get_time((Eigen::Index)t_indx), 1e-10);
        auto temp_sol = res_sim.get_value((Eigen::Index)t_indx);
        EXPECT_NEAR(compare[t_indx][1], temp_sol[0], 1e-10);
        EXPECT_NEAR(compare[t_indx][2], temp_sol[1], 1e-10);
        EXPECT_NEAR(compare[t_indx][3], temp_sol[2], 1e-10);
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
