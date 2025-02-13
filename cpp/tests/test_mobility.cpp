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
#define _USE_MATH_DEFINES

#include "matchers.h"
#include "memilio/compartments/simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "ode_seir/model.h"
#include "ode_secir/model.h"

#include "gtest/gtest.h"

TEST(TestMobility, compareNoMobilityWithSingleIntegration)
{
    auto t0   = 0.0;
    auto tmax = 5.0;
    auto dt   = 0.5;

    mio::oseir::Model<double> model1(1);
    model1.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 0.9;
    model1.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = 0.1;
    model1.populations.set_total(1000);
    model1.parameters.get<mio::oseir::ContactPatterns<double>>().get_cont_freq_mat()[0].get_baseline().setConstant(10);
    model1.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(0.4);
    model1.parameters.set<mio::oseir::TimeExposed<double>>(4);
    model1.parameters.set<mio::oseir::TimeInfected<double>>(10);

    auto model2                                                                     = model1;
    model2.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 1.;
    model2.populations.set_total(500);

    auto graph_sim =
        mio::make_mobility_sim(t0, dt,
                               mio::Graph<mio::SimulationNode<mio::Simulation<double, mio::oseir::Model<double>>>,
                                          mio::MobilityEdge<double>>());
    auto& g = graph_sim.get_graph();
    g.add_node(0, model1, t0);
    g.add_node(1, model2, t0);

    g.nodes()[0].property.get_simulation().set_integrator(std::make_shared<mio::EulerIntegratorCore<double>>());
    g.nodes()[1].property.get_simulation().set_integrator(std::make_shared<mio::EulerIntegratorCore<double>>());

    g.add_edge(0, 1, Eigen::VectorXd::Constant(4, 0)); //no mobility along this edge
    g.add_edge(1, 0, Eigen::VectorXd::Constant(4, 0));

    auto single_sim1 = mio::Simulation<double, mio::oseir::Model<double>>(model1, t0);
    auto single_sim2 = mio::Simulation<double, mio::oseir::Model<double>>(model2, t0);
    single_sim1.set_integrator(std::make_shared<mio::EulerIntegratorCore<double>>());
    single_sim2.set_integrator(std::make_shared<mio::EulerIntegratorCore<double>>());

    graph_sim.advance(tmax);
    single_sim1.advance(tmax);
    single_sim2.advance(tmax);

    EXPECT_DOUBLE_EQ(g.nodes()[0].property.get_result().get_last_time(), single_sim1.get_result().get_last_time());
    EXPECT_DOUBLE_EQ(g.nodes()[1].property.get_result().get_last_time(), single_sim2.get_result().get_last_time());

    //graph may have different time steps, so we can't expect high accuracy here
    EXPECT_NEAR(
        (g.nodes()[0].property.get_result().get_last_value() - single_sim1.get_result().get_last_value()).norm(), 0.0,
        1e-6);
    EXPECT_NEAR(
        (g.nodes()[1].property.get_result().get_last_value() - single_sim2.get_result().get_last_value()).norm(), 0.0,
        1e-6);
}

TEST(TestMobility, nodeAdvance)
{
    using Model = mio::osecir::Model<double>;
    Model model(1);
    auto& params = model.parameters;

    auto& cm = static_cast<mio::ContactMatrixGroup&>(model.parameters.get<mio::osecir::ContactPatterns<double>>());
    cm[0].get_minimum()(0, 0) = 5.0;

    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}] = 100;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}, 1000);
    params.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0]            = 1.;
    params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0] = 1.;
    params.apply_constraints();

    double t0 = 2.835;
    double dt = 0.5;

    mio::SimulationNode<mio::Simulation<double, Model>> node(model, t0);
    node.advance(t0, dt);
    ASSERT_DOUBLE_EQ(node.get_result().get_last_time(), t0 + dt);
    ASSERT_EQ(print_wrap(node.get_result().get_last_value()), print_wrap(node.get_last_state()));
}

TEST(TestMobility, edgeApplyMobility)
{
    using Model = mio::osecir::Model<double>;

    //setup nodes
    Model model(1);
    auto& params = model.parameters;
    auto& cm     = static_cast<mio::ContactMatrixGroup&>(model.parameters.get<mio::osecir::ContactPatterns<double>>());
    cm[0].get_baseline()(0, 0) = 5.0;

    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]            = 10;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}, 1000);
    params.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0] = 1.;
    params.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0]   = 1.;
    params.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0]   = 1.;
    params.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0]        = 0.5;
    params.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0]                      = 1.;
    params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]           = 1.;
    params.apply_constraints();
    double t = 3.125;
    mio::SimulationNode<mio::osecir::Simulation<>> node1(model, t);
    mio::SimulationNode<mio::osecir::Simulation<>> node2(model, t);

    //setup edge
    mio::MobilityEdge<double> edge(Eigen::VectorXd::Constant(10, 0.1));

    //forward mobility
    edge.apply_mobility(t, 0.5, node1, node2);
    EXPECT_EQ(print_wrap(node1.get_result().get_last_value()),
              print_wrap((Eigen::VectorXd(10) << 990 - 99, 0, 0, 0, 10 - 1, 0, 0, 0, 0, 0).finished()));
    EXPECT_EQ(print_wrap(node2.get_result().get_last_value()),
              print_wrap((Eigen::VectorXd(10) << 990 + 99, 0, 0, 0, 10 + 1, 0, 0, 0, 0, 0).finished()));

    //returns
    node1.advance(t, 0.5);
    node2.advance(t, 0.5);
    t += 0.5;
    edge.apply_mobility(t, 0.5, node1, node2);
    auto v = node1.get_result().get_last_value();
    EXPECT_DOUBLE_EQ(v.sum(), 1000);
    EXPECT_LT(v[0], 990);
    EXPECT_NEAR(v[0], 990, 50);
    EXPECT_GT(v[1], 0);
    EXPECT_NEAR(v[1], 0, 50.);
    EXPECT_GT(v[2], 0);
    EXPECT_NEAR(v[2], 0, 5.);
    EXPECT_NEAR(v[4], 10, 5.);
    EXPECT_GT(v[6], 0);
    EXPECT_NEAR(v[6], 0, 5.);
    EXPECT_DOUBLE_EQ(node2.get_result().get_last_value().sum(), 1000);

    //change node again
    node1.advance(t, 0.5);
    node2.advance(t, 0.5);
    t += 0.5;
    edge.apply_mobility(t, 0.5, node1, node2);
    EXPECT_DOUBLE_EQ(node1.get_result().get_last_value().sum(), 900);
    EXPECT_DOUBLE_EQ(node2.get_result().get_last_value().sum(), 1100);
}

TEST(TestMobility, add_mobility_result_time_point)
{
    using Model = mio::osecir::Model<double>;

    //setup nodes
    const size_t num_groups = 1;
    Model model(num_groups);
    auto& params = model.parameters;
    auto& cm     = static_cast<mio::ContactMatrixGroup&>(model.parameters.get<mio::osecir::ContactPatterns<double>>());
    cm[0].get_baseline()(0, 0) = 5.0;

    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]          = 10;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]            = 20;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}, 1000);
    params.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0] = 1.;
    params.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0]   = 1.;
    params.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0]   = 1.;
    params.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0]        = 0.5;
    params.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0]                      = 1.;
    params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]           = 1.;
    params.apply_constraints();

    // get indices of INS and ISy compartments.
    std::vector<std::vector<size_t>> indices_save_edges(2);

    // Reserve Space. The multiplication by 2 is necessary because we have the
    // base and the confirmed compartments for each age group.
    for (auto& vec : indices_save_edges) {
        vec.reserve(2 * num_groups);
    }

    // get indices and write them to the vector
    for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(num_groups); ++i) {
        indices_save_edges[0].emplace_back(
            model.populations.get_flat_index({i, mio::osecir::InfectionState::InfectedNoSymptoms}));
        indices_save_edges[0].emplace_back(
            model.populations.get_flat_index({i, mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}));
        indices_save_edges[1].emplace_back(
            model.populations.get_flat_index({i, mio::osecir::InfectionState::InfectedSymptoms}));
        indices_save_edges[1].emplace_back(
            model.populations.get_flat_index({i, mio::osecir::InfectionState::InfectedSymptomsConfirmed}));
    }

    //setup different edges
    double t = 0.;
    mio::SimulationNode<mio::osecir::Simulation<>> node1(model, t);
    mio::SimulationNode<mio::osecir::Simulation<>> node2(model, t);
    mio::MobilityEdge edge1(Eigen::VectorXd::Constant(10, 0.1), indices_save_edges);
    edge1.apply_mobility(t, 0.0, node1, node2);
    auto mobility = edge1.get_mobility_results().get_last_value();
    EXPECT_NEAR(mobility[0], 1.0, 1e-12);
    EXPECT_NEAR(mobility[1], 2.0, 1e-12);
    EXPECT_NEAR(mobility[2], 100.0, 1e-12);

    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 100;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 30;
    mio::SimulationNode<mio::osecir::Simulation<>> node3(model, t);
    mio::SimulationNode<mio::osecir::Simulation<>> node4(model, t);
    mio::MobilityEdge edge2(Eigen::VectorXd::Constant(10, 0.1), indices_save_edges);
    edge2.apply_mobility(t, 0.5, node3, node4);
    mobility = edge2.get_mobility_results().get_last_value();
    EXPECT_NEAR(mobility[0], 11.0, 1e-12);
    EXPECT_NEAR(mobility[1], 5.0, 1e-12);
    EXPECT_NEAR(mobility[2], 113.0, 1e-12);
}
