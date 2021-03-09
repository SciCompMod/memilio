#define _USE_MATH_DEFINES

#include "epidemiology/migration/migration.h"
#include "epidemiology/secir/seir.h"
#include "epidemiology/secir/secir.h"
#include "epidemiology/utils/eigen_util.h"
#include "epidemiology/utils/eigen.h"
#include "epidemiology/model/simulation.h"
#include "matchers.h"
#include "gtest/gtest.h"

#include <cmath>

TEST(TestMigration, compareNoMigrationWithSingleIntegration)
{
    auto t0   = 0;
    auto tmax = 5;
    auto dt   = 0.5;

    epi::SeirModel model1;
    model1.populations[{epi::SeirInfType::S}] = 0.90000000000000002;
    model1.populations[{epi::SeirInfType::E}] = 0.10000000000000001;
    model1.populations.set_total(1000);
    model1.parameters.get<epi::ContactFrequency>().get_baseline()(0, 0) = 10;
    model1.parameters.set<epi::TransmissionRisk>(0.4);
    model1.parameters.set<epi::StageTimeIncubationInv>(1./4);
    model1.parameters.set<epi::StageTimeInfectiousInv>(1./10);

    auto model2 = model1;
    model2.populations[{epi::SeirInfType::S}] = 1.;
    model2.populations.set_total(500);

    auto graph_sim = epi::make_migration_sim(
        t0, dt,
        epi::Graph<epi::ModelNode<epi::Simulation<epi::SeirModel>>, epi::MigrationEdge>());
    auto& g = graph_sim.get_graph();
    g.add_node(0, model1, t0);
    g.add_node(1, model2, t0);

    g.nodes()[0].property.get_simulation().set_integrator(std::make_shared<epi::EulerIntegratorCore>());
    g.nodes()[1].property.get_simulation().set_integrator(std::make_shared<epi::EulerIntegratorCore>());

    g.add_edge(0, 1, Eigen::VectorXd::Constant(4, 0)); //no migration along this edge
    g.add_edge(1, 0, Eigen::VectorXd::Constant(4, 0));

    auto single_sim1 = epi::Simulation<epi::SeirModel>(model1, t0);
    auto single_sim2 = epi::Simulation<epi::SeirModel>(model2, t0);
    single_sim1.set_integrator(std::make_shared<epi::EulerIntegratorCore>());
    single_sim2.set_integrator(std::make_shared<epi::EulerIntegratorCore>());

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

TEST(TestMigration, nodeEvolve)
{
    using Model = epi::SecirModel<epi::AgeGroup1>;
    Model model;
    auto& params = model.parameters;

    auto& cm = static_cast<epi::ContactMatrixGroup&>(model.parameters.get_contact_patterns());
    cm[0].get_minimum()(0, 0) = 5.0;

    model.populations[{(epi::AgeGroup1)0, epi::InfectionType::E}] = 100;
    model.populations.set_difference_from_total(1000, (epi::AgeGroup1)0, epi::InfectionType::S);
    params.times[0].set_serialinterval(1.5);
    params.times[0].set_incubation(2.0);
    params.apply_constraints();

    double t0 = 2.835;
    double dt = 0.5;

    epi::ModelNode<epi::Simulation<Model>> node(model, t0);
    node.evolve(t0, dt);
    ASSERT_DOUBLE_EQ(node.get_result().get_last_time(), t0 + dt);
    ASSERT_EQ(print_wrap(node.get_result().get_last_value()), print_wrap(node.get_last_state()));
}

TEST(TestMigration, edgeApplyMigration)
{
    using Model = epi::SecirModel<epi::AgeGroup1>;

    //setup nodes
    Model model;
    auto& params = model.parameters;
    auto& cm = static_cast<epi::ContactMatrixGroup&>(model.parameters.get_contact_patterns());
    cm[0].get_baseline()(0, 0) = 5.0;

    model.populations[{(epi::AgeGroup1)0, epi::InfectionType::I}] = 10;
    model.populations.set_difference_from_total(1000, (epi::AgeGroup1)0, epi::InfectionType::S);
    params.probabilities[0].set_infection_from_contact(1.0);
    params.probabilities[0].set_risk_from_symptomatic(1.0);
    params.probabilities[0].set_carrier_infectability(1.0);
    params.probabilities[0].set_hospitalized_per_infectious(0.5);
    params.times[0].set_serialinterval(1.5);
    params.times[0].set_incubation(2.0);
    params.apply_constraints();
    double t = 3.125;
    epi::ModelNode<epi::Simulation<Model>> node1(model, t);
    epi::ModelNode<epi::Simulation<Model>> node2(model, t);

    //setup edge
    epi::MigrationEdge edge(Eigen::VectorXd::Constant(8, 0.1));

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
