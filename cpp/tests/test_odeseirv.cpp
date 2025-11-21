/*
* Copyright (C) 2020-2025 MEmilio
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

#include "memilio/config.h"
#include "memilio/math/integrator.h"
#include "memilio/math/euler.h"
#include "memilio/utils/time_series.h"
#include "memilio/compartments/simulation.h"
#include "memilio/compartments/flow_simulation.h"

#include "ode_seirv/model.h"
#include "ode_seirv/infection_state.h"
#include "ode_seirv/parameters.h"

#include <gtest/gtest.h>

#include <memory>
#include <cmath>

TEST(TestOdeSeirv, simulateDefault)
{
    double t0   = 0.0;
    double tmax = 1.0;
    double dt   = 0.1;

    mio::oseirv::Model<double> model(1);
    auto result = simulate(t0, tmax, dt, model); // generic template simulate

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-12);
    EXPECT_EQ(result.get_num_elements(), (Eigen::Index)mio::oseirv::InfectionState::Count);
}

TEST(TestOdeSeirv, populationConservation)
{
    mio::oseirv::Model<double> model(1);
    const double S  = 300;
    const double SV = 200;
    const double E  = 50;
    const double EV = 25;
    const double I  = 40;
    const double IV = 10;
    const double R  = 0;
    const double RV = 0;
    double total    = S + SV + E + EV + I + IV + R + RV;

    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::Susceptible}]           = S;
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::SusceptibleVaccinated}] = SV;
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::Exposed}]               = E;
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::ExposedVaccinated}]     = EV;
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::Infected}]              = I;
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::InfectedVaccinated}]    = IV;
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::Recovered}]             = R;
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::RecoveredVaccinated}]   = RV;

    mio::ContactMatrixGroup<ScalarType>& cm_h = model.parameters.get<mio::oseirv::ContactPatternsHealthy<double>>();
    cm_h[0]                                   = mio::ContactMatrix<double>(Eigen::MatrixXd::Constant(1, 1, 1));
    mio::ContactMatrixGroup<ScalarType>& cm_s = model.parameters.get<mio::oseirv::ContactPatternsSick<double>>();
    cm_s[0]                                   = mio::ContactMatrix<double>(Eigen::MatrixXd::Constant(1, 1, 1));

    model.parameters.set<mio::oseirv::BaselineTransmissibility<double>>(1.0);
    model.parameters.set<mio::oseirv::TimeExposed<double>>(1.0);
    model.parameters.set<mio::oseirv::TimeInfected<double>>(1.0);

    double t0 = 0.0, tmax = 5.0, dt = 0.5;
    auto sim  = simulate(t0, tmax, dt, model);
    auto last = sim.get_last_value();
    EXPECT_NEAR(last.sum(), total, 1e-8);
}

TEST(TestOdeSeirv, applyConstraints)
{
    mio::oseirv::Model<double> model(1);
    // First: defaults should need no correction
    EXPECT_FALSE(model.parameters.apply_constraints());

    model.parameters.set<mio::oseirv::TimeExposed<double>>(0.0); // invalid, will be set to minimum time
    model.parameters.set<mio::oseirv::TimeInfected<double>>(-1.0); // invalid, will be set to minimum time
    model.parameters.set<mio::oseirv::ClusteringExponent<double>>(-0.5); // invalid, must become >0 (1.0)
    model.parameters.set<mio::oseirv::BaselineTransmissibility<double>>(-2.); // invalid -> 0
    model.parameters.set<mio::oseirv::OutsideFoI<double>>(-0.5); // invalid -> 0
    EXPECT_TRUE(model.parameters.apply_constraints());
    EXPECT_DOUBLE_EQ((double)model.parameters.get<mio::oseirv::TimeExposed<double>>(), 1e-1);
    EXPECT_DOUBLE_EQ((double)model.parameters.get<mio::oseirv::TimeInfected<double>>(), 1e-1);
    EXPECT_EQ(model.parameters.get<mio::oseirv::ClusteringExponent<double>>(), 1.0);
    EXPECT_EQ(model.parameters.get<mio::oseirv::BaselineTransmissibility<double>>(), 0.0);
    EXPECT_EQ(model.parameters.get<mio::oseirv::OutsideFoI<double>>(), 0.0);
}

TEST(TestOdeSeirv, checkConstraints)
{
    mio::oseirv::Model<double> model(1);
    EXPECT_FALSE(model.parameters.check_constraints());

    model.parameters.set<mio::oseirv::TimeExposed<double>>(0.05);
    EXPECT_TRUE(model.parameters.check_constraints());
    model.parameters.set<mio::oseirv::TimeExposed<double>>(0.5);
    EXPECT_FALSE(model.parameters.check_constraints());

    model.parameters.set<mio::oseirv::TimeInfected<double>>(0.05);
    EXPECT_TRUE(model.parameters.check_constraints());
    model.parameters.set<mio::oseirv::TimeInfected<double>>(0.5);
    EXPECT_FALSE(model.parameters.check_constraints());

    model.parameters.set<mio::oseirv::ClusteringExponent<double>>(0.0);
    EXPECT_TRUE(model.parameters.check_constraints());
    model.parameters.set<mio::oseirv::ClusteringExponent<double>>(1.0);
    EXPECT_FALSE(model.parameters.check_constraints());

    model.parameters.set<mio::oseirv::BaselineTransmissibility<double>>(-1.0);
    EXPECT_TRUE(model.parameters.check_constraints());
    model.parameters.set<mio::oseirv::BaselineTransmissibility<double>>(0.5);
    EXPECT_FALSE(model.parameters.check_constraints());

    model.parameters.set<mio::oseirv::OutsideFoI<double>>(-0.5);
    EXPECT_TRUE(model.parameters.check_constraints());
}

TEST(TestOdeSeirv, flowsSingleAgeGroup)
{
    mio::oseirv::Model<double> model(1);

    // Populations
    const double S = 300, SV = 200, E = 50, EV = 25, I = 40, IV = 10; // R, RV = 0
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::Susceptible}]           = S;
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::SusceptibleVaccinated}] = SV;
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::Exposed}]               = E;
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::ExposedVaccinated}]     = EV;
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::Infected}]              = I;
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::InfectedVaccinated}]    = IV;

    // Parameters: identity contacts, TimeInfected=1, ClusteringExponent=1, BaselineTransmissibility=1, others zero ⇒ λ = (I+IV)/N
    mio::ContactMatrixGroup<ScalarType>& cm_h = model.parameters.get<mio::oseirv::ContactPatternsHealthy<double>>();
    cm_h[0]                                   = mio::ContactMatrix<double>(Eigen::MatrixXd::Constant(1, 1, 1));
    mio::ContactMatrixGroup<ScalarType>& cm_s = model.parameters.get<mio::oseirv::ContactPatternsSick<double>>();
    cm_s[0]                                   = mio::ContactMatrix<double>(Eigen::MatrixXd::Constant(1, 1, 1));
    model.parameters.set<mio::oseirv::BaselineTransmissibility<double>>(1.0);
    model.parameters.set<mio::oseirv::TimeExposed<double>>(1.0);
    model.parameters.set<mio::oseirv::TimeInfected<double>>(1.0);
    model.parameters.set<mio::oseirv::ClusteringExponent<double>>(1.0);
    model.parameters.set<mio::oseirv::OutsideFoI<double>>(0.0);
    model.parameters.set<mio::oseirv::SeasonalityAmplitude<double>>(0.0);

    auto y0  = model.get_initial_values();
    auto pop = y0; // same vector for signature
    Eigen::VectorXd flows(6); // number of flows for one age group
    flows.setZero();
    model.get_flows(pop, y0, 0.0, flows);

    const double inv_time_exposed  = 1.0 / model.parameters.get<mio::oseirv::TimeExposed<double>>();
    const double inv_time_infected = 1.0 / model.parameters.get<mio::oseirv::TimeInfected<double>>();
    const double N                 = S + SV + E + EV + I + IV; // (R,RV = 0)
    const double lambda            = (I + IV) / N; // expected force of infection
    const double f_SE              = S * lambda;
    const double f_SV_EV           = SV * lambda;
    const double f_E_I             = inv_time_exposed * E;
    const double f_EV_IV           = inv_time_exposed * EV;
    const double f_I_R             = inv_time_infected * I;
    const double f_IV_RV           = inv_time_infected * IV;

    auto idx_SE = model.template get_flat_flow_index<mio::oseirv::InfectionState::Susceptible,
                                                     mio::oseirv::InfectionState::Exposed>(mio::AgeGroup(0));
    auto idx_SV = model.template get_flat_flow_index<mio::oseirv::InfectionState::SusceptibleVaccinated,
                                                     mio::oseirv::InfectionState::ExposedVaccinated>(mio::AgeGroup(0));
    auto idx_EI =
        model.template get_flat_flow_index<mio::oseirv::InfectionState::Exposed, mio::oseirv::InfectionState::Infected>(
            mio::AgeGroup(0));
    auto idx_EIV =
        model.template get_flat_flow_index<mio::oseirv::InfectionState::ExposedVaccinated,
                                           mio::oseirv::InfectionState::InfectedVaccinated>(mio::AgeGroup(0));
    auto idx_IR = model.template get_flat_flow_index<mio::oseirv::InfectionState::Infected,
                                                     mio::oseirv::InfectionState::Recovered>(mio::AgeGroup(0));
    auto idx_IVR =
        model.template get_flat_flow_index<mio::oseirv::InfectionState::InfectedVaccinated,
                                           mio::oseirv::InfectionState::RecoveredVaccinated>(mio::AgeGroup(0));

    EXPECT_NEAR(flows[idx_SE], f_SE, 1e-12);
    EXPECT_NEAR(flows[idx_SV], f_SV_EV, 1e-12);
    EXPECT_NEAR(flows[idx_EI], f_E_I, 1e-12);
    EXPECT_NEAR(flows[idx_EIV], f_EV_IV, 1e-12);
    EXPECT_NEAR(flows[idx_IR], f_I_R, 1e-12);
    EXPECT_NEAR(flows[idx_IVR], f_IV_RV, 1e-12);
}

TEST(TestOdeSeirv, flowsTwoAgeGroupsIdentityContacts)
{
    mio::oseirv::Model<double> model(2);
    mio::ContactMatrixGroup<ScalarType>& cm_h = model.parameters.get<mio::oseirv::ContactPatternsHealthy<double>>();
    // Let each group have only contacts with itself and set other parameters, so that λ_i = I_i / N_i
    Eigen::MatrixXd Id2                       = Eigen::MatrixXd::Identity(2, 2);
    cm_h[0]                                   = mio::ContactMatrix<double>(Id2);
    mio::ContactMatrixGroup<ScalarType>& cm_s = model.parameters.get<mio::oseirv::ContactPatternsSick<double>>();
    cm_s[0]                                   = mio::ContactMatrix<double>(Eigen::MatrixXd::Zero(2, 2));
    model.parameters.set<mio::oseirv::BaselineTransmissibility<double>>(1.0);
    model.parameters.set<mio::oseirv::TimeExposed<double>>(1.0);
    model.parameters.set<mio::oseirv::TimeInfected<double>>(1.0);
    model.parameters.set<mio::oseirv::ClusteringExponent<double>>(1.0);
    model.parameters.set<mio::oseirv::SeasonalityAmplitude<double>>(0.0);
    model.parameters.set<mio::oseirv::OutsideFoI<double>>(0.0);

    // Only population in the non-vaccinated susceptible and infected compartments
    // Group 0
    double S0 = 100, I0 = 20; // others zero
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::Susceptible}] = S0;
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::Infected}]    = I0;
    // Group 1
    double S1 = 80, I1 = 10;
    model.populations[{mio::AgeGroup(1), mio::oseirv::InfectionState::Susceptible}] = S1;
    model.populations[{mio::AgeGroup(1), mio::oseirv::InfectionState::Infected}]    = I1;

    auto y0  = model.get_initial_values();
    auto pop = y0;
    Eigen::VectorXd flows(12); // 6 flows * 2 age groups
    flows.setZero();
    model.get_flows(pop, y0, 0.0, flows);

    double N0                      = S0 + I0;
    double N1                      = S1 + I1;
    double lambda0                 = I0 / N0; // identity contacts => only own group contributes
    double lambda1                 = I1 / N1;
    const double inv_time_infected = 1.0 / model.parameters.get<mio::oseirv::TimeInfected<double>>();

    auto idx_SE_0 = model.template get_flat_flow_index<mio::oseirv::InfectionState::Susceptible,
                                                       mio::oseirv::InfectionState::Exposed>(mio::AgeGroup(0));
    auto idx_SE_1 = model.template get_flat_flow_index<mio::oseirv::InfectionState::Susceptible,
                                                       mio::oseirv::InfectionState::Exposed>(mio::AgeGroup(1));
    auto idx_EI_0 =
        model.template get_flat_flow_index<mio::oseirv::InfectionState::Exposed, mio::oseirv::InfectionState::Infected>(
            mio::AgeGroup(0));
    auto idx_EI_1 =
        model.template get_flat_flow_index<mio::oseirv::InfectionState::Exposed, mio::oseirv::InfectionState::Infected>(
            mio::AgeGroup(1));
    auto idx_IR_0 = model.template get_flat_flow_index<mio::oseirv::InfectionState::Infected,
                                                       mio::oseirv::InfectionState::Recovered>(mio::AgeGroup(0));
    auto idx_IR_1 = model.template get_flat_flow_index<mio::oseirv::InfectionState::Infected,
                                                       mio::oseirv::InfectionState::Recovered>(mio::AgeGroup(1));

    EXPECT_NEAR(flows[idx_SE_0], S0 * lambda0, 1e-12);
    EXPECT_NEAR(flows[idx_SE_1], S1 * lambda1, 1e-12);
    // Exposed are zero => progressions must be zero
    EXPECT_NEAR(flows[idx_EI_0], 0.0, 1e-12);
    EXPECT_NEAR(flows[idx_EI_1], 0.0, 1e-12);
    EXPECT_NEAR(flows[idx_IR_0], inv_time_infected * I0, 1e-12);
    EXPECT_NEAR(flows[idx_IR_1], inv_time_infected * I1, 1e-12);
}

TEST(TestOdeSeirv, zeroPopulationNoNan)
{
    mio::oseirv::Model<double> model(1);
    model.populations.set_total(0.0);
    auto y0  = model.get_initial_values();
    auto pop = y0;
    Eigen::VectorXd flows(6);
    flows.setZero();
    model.get_flows(pop, y0, 0.0, flows);
    for (int i = 0; i < flows.size(); ++i) {
        EXPECT_FALSE(std::isnan(flows[i]));
        EXPECT_EQ(flows[i], 0.0);
    }
}

TEST(TestOdeSeirv, simulationEuler)
{
    mio::oseirv::Model<double> model(1);
    // Simple initial values
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::Exposed}]  = 10;
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::Infected}] = 5;
    model.populations.set_difference_from_group_total<mio::AgeGroup>(
        {mio::AgeGroup(0), mio::oseirv::InfectionState::Susceptible}, 100.0);

    // Identity contacts
    mio::ContactMatrixGroup<ScalarType>& cm_h = model.parameters.get<mio::oseirv::ContactPatternsHealthy<double>>();
    cm_h[0]                                   = mio::ContactMatrix<double>(Eigen::MatrixXd::Constant(1, 1, 1));
    mio::ContactMatrixGroup<ScalarType>& cm_s = model.parameters.get<mio::oseirv::ContactPatternsSick<double>>();
    cm_s[0]                                   = mio::ContactMatrix<double>(Eigen::MatrixXd::Constant(1, 1, 0.));
    model.parameters.set<mio::oseirv::BaselineTransmissibility<double>>(1.0);
    model.parameters.set<mio::oseirv::TimeExposed<double>>(1.0);
    model.parameters.set<mio::oseirv::TimeInfected<double>>(1.0);

    double t0 = 0.0, tmax = 1.0, dt = 0.5;
    std::unique_ptr<mio::OdeIntegratorCore<double>> integrator = std::make_unique<mio::EulerIntegratorCore<double>>();
    auto sim                                                   = simulate(t0, tmax, dt, model, std::move(integrator));
    EXPECT_EQ(sim.get_num_time_points(), 3); // t=0,0.5,1.0
    // Sanity: all values non-negative
    for (Eigen::Index i = 0; i < sim.get_last_value().size(); ++i) {
        EXPECT_GE(sim.get_last_value()[i], 0.0);
    }
}

TEST(TestOdeSeirv, flowSimulationEuler)
{
    mio::oseirv::Model<double> model(1);
    model.populations[{mio::AgeGroup(0), mio::oseirv::InfectionState::Infected}] = 5;
    model.populations.set_difference_from_group_total<mio::AgeGroup>(
        {mio::AgeGroup(0), mio::oseirv::InfectionState::Susceptible}, 100.0);

    mio::ContactMatrixGroup<ScalarType>& cm_h = model.parameters.get<mio::oseirv::ContactPatternsHealthy<double>>();
    cm_h[0]                                   = mio::ContactMatrix<double>(Eigen::MatrixXd::Constant(1, 1, 1));
    model.parameters.set<mio::oseirv::BaselineTransmissibility<double>>(1.0);
    model.parameters.set<mio::oseirv::TimeExposed<double>>(1.0);
    model.parameters.set<mio::oseirv::TimeInfected<double>>(1.0);

    double t0 = 0.0, tmax = 0.5, dt = 0.5;
    std::unique_ptr<mio::OdeIntegratorCore<double>> integrator = std::make_unique<mio::EulerIntegratorCore<double>>();
    auto sim = simulate_flows(t0, tmax, dt, model, std::move(integrator));
    EXPECT_EQ(sim[0].get_num_time_points(), 2);
    EXPECT_EQ(sim[1].get_num_time_points(), 2);
    // Flow vector size should be 6 for one age group
    EXPECT_EQ(sim[1].get_last_value().size(), 6);
}
