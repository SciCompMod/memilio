#define _USE_MATH_DEFINES

#include "epidemiology/migration/migration.h"
#include "epidemiology/secir/seir.h"
#include "epidemiology/secir/secir.h"
#include "epidemiology/utils/eigen_util.h"
#include "epidemiology/utils/eigen.h"
#include "matchers.h"
#include "gtest/gtest.h"

#include <cmath>

TEST(TestMigration, compareNoMigrationWithSingleIntegration)
{
    auto t0   = 0;
    auto tmax = 5;
    auto dt = 0.5;
    
    auto params1 = epi::SeirParams();
    params1.populations.set({epi::SeirCompartments::SeirS}, 0.9);
    params1.populations.set({epi::SeirCompartments::SeirE}, 0.1);
    params1.populations.set_total(1000);
    params1.contact_frequency.get_baseline()(0, 0) = 2.5;
    params1.times.set_incubation(4);
    params1.times.set_infectious(10);

    auto params2 = params1;
    params2.populations.set({epi::SeirCompartments::SeirS}, 1.);
    params2.populations.set_total(500);

    auto graph_sim =
        epi::make_migration_sim(t0, dt, epi::Graph<epi::ModelNode<epi::SeirSimulation>, epi::MigrationEdge>());
    auto& g = graph_sim.get_graph();
    g.add_node(params1, t0);
    g.add_node(params2, t0);
    g.add_edge(0, 1, Eigen::VectorXd::Constant(4, 0)); //no migration along this edge
    g.add_edge(1, 0, Eigen::VectorXd::Constant(4, 0));

    auto single_sim1 = epi::SeirSimulation(params1, t0);
    auto single_sim2 = epi::SeirSimulation(params2, t0);

    graph_sim.advance(tmax);
    single_sim1.advance(tmax);
    single_sim2.advance(tmax);

    EXPECT_DOUBLE_EQ(g.nodes()[0].get_result().get_last_time(), single_sim1.get_result().get_last_time());
    EXPECT_DOUBLE_EQ(g.nodes()[1].get_result().get_last_time(), single_sim2.get_result().get_last_time());

    //graph may have different time steps, so we can't expect high accuracy here
    EXPECT_NEAR((g.nodes()[0].get_result().get_last_value() - single_sim1.get_result().get_last_value()).norm(),
                0.0, 1e-6);
    EXPECT_NEAR((g.nodes()[1].get_result().get_last_value() - single_sim2.get_result().get_last_value()).norm(),
                0.0, 1e-6);
}

TEST(TestMigration, nodeEvolve)
{
    auto cm = epi::ContactMatrixGroup(1, 1);
    cm[0].get_minimum()(0, 0) = 5.0;
    epi::SecirParams params(cm);
    params.populations.set({0, epi::E}, 100);
    params.populations.set_difference_from_total({0, epi::S}, 1000);
    params.times[0].set_serialinterval(1.5);
    params.times[0].set_incubation(2.0);
    params.apply_constraints();

    double t0 = 2.835;
    double dt = 0.5;

    epi::ModelNode<epi::SecirSimulation> node(params, t0);
    node.evolve(t0, dt);
    ASSERT_DOUBLE_EQ(node.get_result().get_last_time(), t0 + dt);
    ASSERT_EQ(print_wrap(node.get_result().get_last_value()), print_wrap(node.get_last_state()));
}

TEST(TestMigration, edgeApplyMigration)
{
    //setup nodes
    auto cm = epi::ContactMatrixGroup(1, 1);
    cm[0].get_baseline()(0, 0) = 5.0;
    epi::SecirParams params(cm);
    params.populations.set({0, epi::I}, 10);
    params.populations.set_difference_from_total({0, epi::S}, 1000);
    params.probabilities[0].set_infection_from_contact(1.0);
    params.probabilities[0].set_risk_from_symptomatic(1.0);
    params.probabilities[0].set_carrier_infectability(1.0);
    params.probabilities[0].set_hospitalized_per_infectious(0.5);
    params.times[0].set_serialinterval(1.5);
    params.times[0].set_incubation(2.0);
    params.apply_constraints();
    double t = 3.125;
    epi::ModelNode<epi::SecirSimulation> node1(params, t);
    epi::ModelNode<epi::SecirSimulation> node2(params, t);

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
