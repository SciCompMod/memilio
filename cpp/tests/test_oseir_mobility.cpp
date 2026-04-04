/*
* Copyright (C) 2020-2026 MEmilio
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
#include "memilio/compartments/simulation.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/math/stepper_wrapper.h"
#include "ode_seir/model.h"
#include "ode_seir/mobility.h"
#include <gtest/gtest.h>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <cmath>

namespace
{

mio::oseir::Model<double> create_seir_model(size_t num_agegroups, double total_pop = 10000.0)
{
    mio::oseir::Model<double> model(static_cast<int>(num_agegroups));

    for (size_t g = 0; g < num_agegroups; ++g) {
        auto ag                                                        = mio::AgeGroup(static_cast<int>(g));
        model.populations[{ag, mio::oseir::InfectionState::Exposed}]   = 50.0;
        model.populations[{ag, mio::oseir::InfectionState::Infected}]  = 30.0;
        model.populations[{ag, mio::oseir::InfectionState::Recovered}] = 20.0;
        model.populations[{ag, mio::oseir::InfectionState::Susceptible}] =
            total_pop / num_agegroups - 50.0 - 30.0 - 20.0;
    }

    model.parameters.set<mio::oseir::TimeExposed<double>>(5.2);
    model.parameters.set<mio::oseir::TimeInfected<double>>(6.0);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(0.1);

    auto& cm = model.parameters.get<mio::oseir::ContactPatterns<double>>();
    cm.get_cont_freq_mat()[0].get_baseline().setConstant(2.7);

    return model;
}

} // namespace

TEST(TestOseirMobility, noMobilityMatchesSingleSim)
{
    const double t0   = 0.0;
    const double tmax = 10.0;
    const double dt   = 0.5;

    auto model = create_seir_model(1);

    // Graph simulation with zero mobility coefficients
    mio::Graph<mio::SimulationNode<double, mio::Simulation<double, mio::oseir::Model<double>>>,
               mio::MobilityEdge<double>>
        g;
    g.add_node(0, model, t0);
    g.add_node(1, model, t0);
    g.add_edge(0, 1, Eigen::VectorXd::Constant(4, 0.0));
    g.add_edge(1, 0, Eigen::VectorXd::Constant(4, 0.0));

    auto graph_sim = mio::oseir::make_seir_mobility_sim<double>(t0, dt, std::move(g));
    graph_sim.advance(tmax);

    // Single simulation with SeirRK4IntegratorCore (same as used in graph)
    auto single_sim = mio::Simulation<double, mio::oseir::Model<double>>(model, t0);
    single_sim.set_integrator_core(std::make_unique<mio::oseir::SeirRK4IntegratorCore<double>>(single_sim.get_model()));
    single_sim.get_dt() = dt;
    single_sim.advance(tmax);

    auto diff = (graph_sim.get_graph().nodes()[0].property.get_result().get_last_value() -
                 single_sim.get_result().get_last_value())
                    .norm();
    EXPECT_NEAR(diff, 0.0, 1e-10);
}

TEST(TestOseirMobility, mobilityReturnsPreservation)
{
    const double t0   = 0.0;
    const double tmax = 5.0;
    const double dt   = 0.5;

    auto model1 = create_seir_model(1, 10000.0);
    auto model2 = create_seir_model(1, 5000.0);

    mio::Graph<mio::SimulationNode<double, mio::Simulation<double, mio::oseir::Model<double>>>,
               mio::MobilityEdge<double>>
        g;
    g.add_node(0, model1, t0);
    g.add_node(1, model2, t0);
    g.add_edge(0, 1, Eigen::VectorXd::Constant(4, 0.05));
    g.add_edge(1, 0, Eigen::VectorXd::Constant(4, 0.03));

    auto graph_sim = mio::oseir::make_seir_mobility_sim<double>(t0, dt, std::move(g));
    graph_sim.advance(tmax);

    double total_pop = 0.0;
    for (auto& n : graph_sim.get_graph().nodes()) {
        total_pop += n.property.get_result().get_last_value().sum();
    }
    EXPECT_NEAR(total_pop, 15000.0, 1e-6);
}

TEST(TestOseirMobility, seirIntegratorMatchesBoostRK4)
{
    const double t0   = 0.0;
    const double tmax = 10.0;
    const double dt   = 0.5;

    auto model = create_seir_model(1);

    // Run 1: SeirRK4IntegratorCore
    auto sim1 = mio::Simulation<double, mio::oseir::Model<double>>(model, t0);
    sim1.set_integrator_core(std::make_unique<mio::oseir::SeirRK4IntegratorCore<double>>(sim1.get_model()));
    sim1.get_dt() = dt;
    sim1.advance(tmax);

    // Run 2: Standard ExplicitStepperWrapper<runge_kutta4>
    auto sim2 = mio::Simulation<double, mio::oseir::Model<double>>(model, t0);
    sim2.set_integrator_core(
        std::make_unique<mio::ExplicitStepperWrapper<double, boost::numeric::odeint::runge_kutta4>>());
    sim2.get_dt() = dt;
    sim2.advance(tmax);

    auto diff = (sim1.get_result().get_last_value() - sim2.get_result().get_last_value()).norm();
    EXPECT_NEAR(diff, 0.0, 1e-10) << "SeirRK4IntegratorCore must match boost runge_kutta4.\n"
                                  << "Custom: " << sim1.get_result().get_last_value().transpose() << "\n"
                                  << "Boost:  " << sim2.get_result().get_last_value().transpose();
}

TEST(TestOseirMobility, multiAgeGroupMobilityReturns)
{
    const double t0   = 0.0;
    const double tmax = 5.0;
    const double dt   = 0.5;
    const size_t NG   = 3;

    auto model = create_seir_model(NG, 30000.0);

    mio::Graph<mio::SimulationNode<double, mio::Simulation<double, mio::oseir::Model<double>>>,
               mio::MobilityEdge<double>>
        g;
    g.add_node(0, model, t0);
    g.add_node(1, model, t0);

    const size_t NC = 4 * NG;
    g.add_edge(0, 1, Eigen::VectorXd::Constant(NC, 0.05));
    g.add_edge(1, 0, Eigen::VectorXd::Constant(NC, 0.05));

    auto graph_sim = mio::oseir::make_seir_mobility_sim<double>(t0, dt, std::move(g));
    graph_sim.advance(tmax);

    double total_pop = 0.0;
    for (auto& n : graph_sim.get_graph().nodes()) {
        total_pop += n.property.get_result().get_last_value().sum();
    }
    EXPECT_NEAR(total_pop, 60000.0, 1e-4);

    // Check non-negativity
    for (auto& n : graph_sim.get_graph().nodes()) {
        auto v = n.property.get_result().get_last_value();
        for (int i = 0; i < v.size(); ++i) {
            EXPECT_GE(v[i], -1e-12) << "Negative compartment at index " << i;
        }
    }
}

TEST(TestOseirMobility, evaluatorMatchesGetDerivatives)
{
    auto model     = create_seir_model(2, 20000.0);
    const double t = 3.5;

    Eigen::VectorXd z(8);
    z << 8000, 500, 300, 1200, 7000, 600, 400, 2000;

    Eigen::VectorXd dydt_ref(8);
    model.get_derivatives(z, z, t, dydt_ref);

    mio::oseir::SeirCoefficientEvaluator<double> eval(model);
    std::vector<double> lambda(2);
    eval.compute_foi(z, t, lambda);
    Eigen::VectorXd dydt_eval(8);
    eval.apply_rhs(lambda, z, dydt_eval);

    EXPECT_NEAR((dydt_ref - dydt_eval).norm(), 0.0, 1e-12)
        << "SeirCoefficientEvaluator must match model.get_derivatives for totals.\n"
        << "Reference: " << dydt_ref.transpose() << "\n"
        << "Evaluator: " << dydt_eval.transpose();
}

TEST(TestOseirMobility, evaluatorMatchesSubpopDerivatives)
{
    auto model     = create_seir_model(2, 20000.0);
    const double t = 1.0;

    Eigen::VectorXd z(8);
    z << 8000, 500, 300, 1200, 7000, 600, 400, 2000;

    Eigen::VectorXd c(8);
    c << 200, 15, 10, 30, 180, 20, 12, 50;

    Eigen::VectorXd dydt_ref(8);
    model.get_derivatives(z, c, t, dydt_ref);

    mio::oseir::SeirCoefficientEvaluator<double> eval(model);
    std::vector<double> lambda(2);
    eval.compute_foi(z, t, lambda);
    Eigen::VectorXd dydt_eval(8);
    eval.apply_rhs(lambda, c, dydt_eval);

    EXPECT_NEAR((dydt_ref - dydt_eval).norm(), 0.0, 1e-12)
        << "SeirCoefficientEvaluator must match model.get_derivatives for subpopulations.\n"
        << "Reference: " << dydt_ref.transpose() << "\n"
        << "Evaluator: " << dydt_eval.transpose();
}

TEST(TestOseirMobility, cachedReturnsMatchFallback)
{
    const double t0 = 0.0;
    const double dt = 0.5;

    auto model                                                                     = create_seir_model(1);
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]    = 200.0;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 10000.0 - 50.0 - 200.0 - 20.0;

    Eigen::VectorXd total(4);
    total << 9730, 50, 200, 20;

    Eigen::VectorXd commuters(4);
    commuters << 500, 3, 10, 2;

    // (a) Cached path: step() fills cache, then integrate_commuters
    mio::oseir::SeirRK4IntegratorCore<double> integrator(model);
    Eigen::VectorXd ytp1(4);
    double t_step  = t0;
    double dt_step = dt;
    integrator.step(mio::DerivFunction<double>{}, total, t_step, dt_step, ytp1);
    ASSERT_TRUE(integrator.has_valid_cache());

    Eigen::VectorXd c_cached = commuters;
    integrator.integrate_commuters(c_cached, dt);

    // (b) Fallback path: full co-integration
    mio::oseir::SeirCoefficientEvaluator<double> eval(model);
    Eigen::VectorXd c_fallback = commuters;
    mio::oseir::detail::calculate_mobility_returns_seir<double>(eval, c_fallback, total, t0, dt);

    EXPECT_NEAR((c_cached - c_fallback).norm(), 0.0, 1e-12)
        << "Cached path must exactly match fallback co-integration.\n"
        << "Cached:   " << c_cached.transpose() << "\n"
        << "Fallback: " << c_fallback.transpose();
}

TEST(TestOseirMobility, fullGraphSimNonNegative)
{
    const double t0   = 0.0;
    const double tmax = 15.0;
    const double dt   = 0.5;

    auto model                                                                     = create_seir_model(1);
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]    = 200.0;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 10000.0 - 50.0 - 200.0 - 20.0;

    mio::Graph<mio::SimulationNode<double, mio::Simulation<double, mio::oseir::Model<double>>>,
               mio::MobilityEdge<double>>
        g;
    g.add_node(0, model, t0);
    g.add_node(1, model, t0);
    g.add_edge(0, 1, Eigen::VectorXd::Constant(4, 0.1));
    g.add_edge(1, 0, Eigen::VectorXd::Constant(4, 0.1));

    auto graph_sim = mio::oseir::make_seir_mobility_sim<double>(t0, dt, std::move(g));
    graph_sim.advance(tmax);

    double total_pop = 0.0;
    for (auto& n : graph_sim.get_graph().nodes()) {
        total_pop += n.property.get_result().get_last_value().sum();
    }
    EXPECT_NEAR(total_pop, 20000.0, 1e-4);

    for (auto& n : graph_sim.get_graph().nodes()) {
        auto v = n.property.get_result().get_last_value();
        for (int i = 0; i < v.size(); ++i) {
            EXPECT_GE(v[i], -1e-6) << "Negative compartment at index " << i;
        }
    }
}
