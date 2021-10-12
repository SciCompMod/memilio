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
#define _USE_MATH_DEFINES

#include "memilio/mobility/mobility.h"
#include "seir/seir.h"
#include "secir/secir.h"
#include "memilio/math/eigen_util.h"
#include "memilio/math/eigen.h"
#include "memilio/compartments/simulation.h"
#include "matchers.h"
#include "gtest/gtest.h"

#include <cmath>

TEST(TestMobility, compareNoMigrationWithSingleIntegration)
{
    auto t0   = 0;
    auto tmax = 5;
    auto dt   = 0.5;

    mio::SeirModel model1;
    model1.populations[{mio::Index<mio::SeirInfType>(mio::SeirInfType::S)}] = 0.9;
    model1.populations[{mio::Index<mio::SeirInfType>(mio::SeirInfType::E)}] = 0.1;
    model1.populations.set_total(1000);
    model1.parameters.get<mio::ContactFrequency>().get_baseline()(0, 0) = 10;
    model1.parameters.set<mio::TransmissionRisk>(0.4);
    model1.parameters.set<mio::StageTimeIncubationInv>(1./4);
    model1.parameters.set<mio::StageTimeInfectiousInv>(1./10);

    auto model2 = model1;
    model2.populations[{mio::Index<mio::SeirInfType>(mio::SeirInfType::S)}] = 1.;
    model2.populations.set_total(500);

    auto graph_sim = mio::make_migration_sim(
        t0, dt,
        mio::Graph<mio::SimulationNode<mio::Simulation<mio::SeirModel>>, mio::MigrationEdge>());
    auto& g = graph_sim.get_graph();
    g.add_node(0, model1, t0);
    g.add_node(1, model2, t0);

    g.nodes()[0].property.get_simulation().set_integrator(std::make_shared<mio::EulerIntegratorCore>());
    g.nodes()[1].property.get_simulation().set_integrator(std::make_shared<mio::EulerIntegratorCore>());

    g.add_edge(0, 1, Eigen::VectorXd::Constant(4, 0)); //no migration along this edge
    g.add_edge(1, 0, Eigen::VectorXd::Constant(4, 0));

    auto single_sim1 = mio::Simulation<mio::SeirModel>(model1, t0);
    auto single_sim2 = mio::Simulation<mio::SeirModel>(model2, t0);
    single_sim1.set_integrator(std::make_shared<mio::EulerIntegratorCore>());
    single_sim2.set_integrator(std::make_shared<mio::EulerIntegratorCore>());

    graph_sim.advance(tmax);
    single_sim1.advance(tmax);
    single_sim2.advance(tmax);

    EXPECT_DOUBLE_EQ(g.nodes()[0].property.get_result().get_last_time(), single_sim1.get_result().get_last_time());
    EXPECT_DOUBLE_EQ(g.nodes()[1].property.get_result().get_last_time(), single_sim2.get_result().get_last_time());

    //graph may have different time steps, so we can't expect high accuracy here
    EXPECT_NEAR((g.nodes()[0].property.get_result().get_last_value() - single_sim1.get_result().get_last_value()).norm(),
                0.0, 1e-6);
    EXPECT_NEAR((g.nodes()[1].property.get_result().get_last_value() - single_sim2.get_result().get_last_value()).norm(),
                0.0, 1e-6);
}

TEST(TestMobility, nodeEvolve)
{
    using Model = mio::SecirModel;
    Model model(1);
    auto& params = model.parameters;

    auto& cm = static_cast<mio::ContactMatrixGroup&>(model.parameters.get<mio::ContactPatterns>());
    cm[0].get_minimum()(0, 0) = 5.0;

    model.populations[{mio::AgeGroup(0), mio::InfectionState::Exposed}] = 100;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::InfectionState::Susceptible}, 1000);
    params.get<mio::SerialInterval>()[(mio::AgeGroup)0] = 1.5;
    params.get<mio::IncubationTime>()[(mio::AgeGroup)0] = 2.;
    params.apply_constraints();

    double t0 = 2.835;
    double dt = 0.5;

    mio::SimulationNode<mio::Simulation<Model>> node(model, t0);
    node.evolve(t0, dt);
    ASSERT_DOUBLE_EQ(node.get_result().get_last_time(), t0 + dt);
    ASSERT_EQ(print_wrap(node.get_result().get_last_value()), print_wrap(node.get_last_state()));
}

TEST(TestMobility, edgeApplyMigration)
{
    using Model = mio::SecirModel;

    //setup nodes
    Model model(1);
    auto& params = model.parameters;
    auto& cm = static_cast<mio::ContactMatrixGroup&>(model.parameters.get<mio::ContactPatterns>());
    cm[0].get_baseline()(0, 0) = 5.0;

    model.populations[{mio::AgeGroup(0), mio::InfectionState::Infected}] = 10;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::InfectionState::Susceptible}, 1000);
    params.get<mio::InfectionProbabilityFromContact>()[(mio::AgeGroup)0] = 1.;
    params.get<mio::RiskOfInfectionFromSympomatic>()[(mio::AgeGroup)0] = 1.;
    params.get<mio::RelativeCarrierInfectability>()[(mio::AgeGroup)0] = 1.;
    params.get<mio::HospitalizedCasesPerInfectious>()[(mio::AgeGroup)0] = 0.5;
    params.get<mio::SerialInterval>()[(mio::AgeGroup)0] = 1.5;
    params.get<mio::IncubationTime>()[(mio::AgeGroup)0] = 2.;
    params.apply_constraints();
    double t = 3.125;
    mio::SimulationNode<mio::Simulation<Model>> node1(model, t);
    mio::SimulationNode<mio::Simulation<Model>> node2(model, t);

    //setup edge
    mio::MigrationEdge edge(Eigen::VectorXd::Constant(8, 0.1));

    //forward migration
    edge.apply_migration(t, 0.5, node1, node2);
    EXPECT_EQ(print_wrap(node1.get_result().get_last_value()), print_wrap((Eigen::VectorXd(8) << 990-99, 0, 0, 10-1, 0, 0, 0, 0).finished()));
    EXPECT_EQ(print_wrap(node2.get_result().get_last_value()), print_wrap((Eigen::VectorXd(8) << 990+99, 0, 0, 10+1, 0, 0, 0, 0).finished()));

    //returns
    node1.evolve(t, 0.5);
    node2.evolve(t, 0.5);
    t += 0.5;
    edge.apply_migration(t, 0.5, node1, node2);
    auto v = node1.get_result().get_last_value();
    EXPECT_DOUBLE_EQ(v.sum(), 1000);
    EXPECT_LT(v[0], 990);
    EXPECT_NEAR(v[0], 990, 50);
    EXPECT_GT(v[1], 0);
    EXPECT_NEAR(v[1], 0, 50.);
    EXPECT_GT(v[2], 0);
    EXPECT_NEAR(v[2], 0, 5.);
    EXPECT_NEAR(v[3], 10, 5.);
    EXPECT_GT(v[4], 0);
    EXPECT_NEAR(v[4], 0, 5.);
    EXPECT_DOUBLE_EQ(node2.get_result().get_last_value().sum(), 1000);
    
    
    //migrate again
    node1.evolve(t, 0.5);
    node2.evolve(t, 0.5);
    t += 0.5;
    edge.apply_migration(t, 0.5, node1, node2);
    EXPECT_DOUBLE_EQ(node1.get_result().get_last_value().sum(), 900);
    EXPECT_DOUBLE_EQ(node2.get_result().get_last_value().sum(), 1100);
}
