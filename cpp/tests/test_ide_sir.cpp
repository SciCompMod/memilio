/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Anna Wendler
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
#include "load_test_data.h"
#include "ide_sir/infection_state.h"
#include "ide_sir/model.h"
#include "ide_sir/model_simplified.h"
#include "ide_sir/parameters.h"
#include "ide_sir/simulation.h"
#include "ide_sir/simulation_simplified.h"
#include "ide_sir/gregory_weights.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/utils/time_series.h"
#include <gtest/gtest.h>
#include <vector>
#include <cmath>

// Check that the Gregory weights are set correctly for orders 1, 2 and 3.
TEST(IdeSir, checkGregoryWeights)
{

    // Define t_max. Below, we will test the Gregory weights for n=gregory_order, ..., n_max.
    size_t n_max = 6;

    // Define Gregory orders that we test.
    std::vector<size_t> gregory_orders = {1, 2, 3};

    for (size_t gregory_order : gregory_orders) {

        // Here, we define the matrices for Sigma and the vector corresponding to Omega explicitly for n up to n_max. We
        // will then compare this to our implementation in get_gregoryweights(), where we implemented the matrices and
        // vectors in a reduced form independent of n_max.
        std::vector<std::vector<ScalarType>> gregoryWeights_sigma_expected;
        std::vector<ScalarType> gregoryWeights_omega_expected;
        switch (gregory_order) {
        case 1:
            gregoryWeights_sigma_expected = {{1. / 2.}, {1. / 2.}, {1. / 2.}, {1. / 2.}, {1. / 2.}, {1. / 2.}};
            gregoryWeights_omega_expected = {1. / 2., 1., 1., 1., 1., 1.};
            break;
        case 2:
            gregoryWeights_sigma_expected = {{5. / 12., 14. / 12.},
                                             {5. / 12., 13. / 12.},
                                             {5. / 12., 13. / 12.},
                                             {5. / 12., 13. / 12.},
                                             {5. / 12., 13. / 12.}};
            gregoryWeights_omega_expected = {5. / 12., 13. / 12., 12. / 12., 12. / 12., 12. / 12.};
            break;
        case 3:
            gregoryWeights_sigma_expected = {{9. / 24., 27. / 24., 27. / 24.},
                                             {9. / 24., 28. / 24., 22. / 24.},
                                             {9. / 24., 28. / 24., 23. / 24.},
                                             {9. / 24., 28. / 24., 23. / 24.}};
            gregoryWeights_omega_expected = {9. / 24., 28. / 24., 23. / 24., 24. / 24.};
            break;
        }
        // Get Gregory weights from implementation in get_gregoryweights().
        std::vector<Eigen::MatrixX<ScalarType>> vec_gregoryweights = mio::isir::get_gregoryweights(gregory_order);
        Eigen::MatrixX<ScalarType> gregoryWeights_sigma            = vec_gregoryweights[0];
        Eigen::MatrixX<ScalarType> gregoryWeights_omega            = vec_gregoryweights[1];

        size_t row_index, column_index;
        for (size_t n = gregory_order; n <= n_max; n++) {
            for (size_t j = 0; j < gregory_order; j++) {

                // Check sigma.

                // Compute row_index and column_index according to n and j as in the implementation of sum_part1.
                // If n <= gregory_order + gregory_order - 2, then the row_index corresponds to n - gregory_order.
                if (n <= gregory_order + gregory_order - 2) {
                    row_index = n - gregory_order;
                }
                // Else, for n >= m_gregory_order - 1, the entries in gregoryWeights_sigma do not change anymore and the
                // corresponding row_index is given by gregory_order - 1.
                else {
                    row_index = gregory_order - 1;
                }
                // The column index only depends on the current index of the sum j.
                column_index = j;

                EXPECT_EQ(gregoryWeights_sigma(row_index, column_index),
                          gregoryWeights_sigma_expected[n - gregory_order][j]);
            }
        }

        size_t weight_index;
        for (size_t n = gregory_order; n <= n_max; n++) {
            for (size_t j = gregory_order; j <= n; j++) {
                // Check omega.

                // Compute weight_index according to n and j as in the implementation of sum_part2.
                // Depending on the Gregory order and the current time step, we determine the required weight index of gregoryWeights_omega.
                // This is necessary because we implemented gregoryWeights_omega in a reduced way.
                if (n - j <= gregory_order) {
                    weight_index = n - j;
                }
                else {
                    weight_index = gregory_order;
                }

                EXPECT_EQ(gregoryWeights_omega(weight_index), gregoryWeights_omega_expected[n - j]);
            }

            // Define weight_indices that we want to check the values for in gregoryWeights_omega.
            std::vector<size_t> weight_indices = {0, 1, 2, 3};
        }
    }
}

// Test that ModelMessinaExtended simplifies to ModelMessina when choosing the contact rate appropriately.
TEST(IdeSir, compareModelMessinaAndModelMessinaExtended)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;
    // Define simulation parameters.
    ScalarType t0   = 0.;
    ScalarType tmax = 1.;
    ScalarType dt   = 0.1;

    ScalarType S0               = 90.;
    ScalarType I0               = 10.;
    ScalarType R0               = 0.;
    ScalarType total_population = S0 + I0 + R0;

    // We want to check in particular that choosing beta based on cont_freq and the total_population leads to the same results.

    ScalarType cont_freq = 0.1;
    ScalarType beta      = cont_freq / total_population;

    size_t gregory_order = 2;

    // Define init_populations.
    mio::TimeSeries<ScalarType> init_populations((size_t)mio::isir::InfectionState::Count);
    mio::TimeSeries<ScalarType> init_populations_extended((size_t)mio::isir::InfectionState::Count);
    // Initialize first (gregory_order-1) time points with constant values.
    // Only set S because this is the only compartment we consider at the moment.
    Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));
    vec_init[(size_t)mio::isir::InfectionState::Susceptible] = S0;
    // Add time points S_0, S_1, S_{n0-1} to init_populations as these values are assumed to be known in the groundtruth.
    init_populations.add_time_point(t0, vec_init);
    init_populations_extended.add_time_point(t0, vec_init);
    while (init_populations.get_last_time() < (gregory_order - 1) * dt - 1e-10) {
        init_populations.add_time_point(init_populations.get_last_time() + dt, vec_init);
        init_populations_extended.add_time_point(init_populations_extended.get_last_time() + dt, vec_init);
    }

    // Set up ModelMessina and ModelMessinaextended.
    mio::isir::ModelMessina model(std::move(init_populations), total_population, gregory_order);
    mio::isir::ModelMessinaExtended model_extended(std::move(init_populations_extended), total_population,
                                                   gregory_order);

    mio::NormalDistributionDensity normaldensity(0.4, 0.6);
    mio::StateAgeFunctionWrapper dist(normaldensity);
    std::vector<mio::StateAgeFunctionWrapper> vec_dist((size_t)mio::isir::InfectionTransition::Count, dist);
    model.parameters.get<mio::isir::TransitionDistributions>()          = vec_dist;
    model_extended.parameters.get<mio::isir::TransitionDistributions>() = vec_dist;

    mio::ConstantFunction transmissiononcontact(1.5);
    mio::StateAgeFunctionWrapper transmissiononcontact_wrapper(transmissiononcontact);
    model.parameters.get<mio::isir::TransmissionProbabilityOnContact>()          = transmissiononcontact_wrapper;
    model_extended.parameters.get<mio::isir::TransmissionProbabilityOnContact>() = transmissiononcontact_wrapper;

    mio::ConstantFunction riskofinfection(1.);
    mio::StateAgeFunctionWrapper riskofinfection_wrapper(riskofinfection);
    model.parameters.get<mio::isir::RiskOfInfectionFromSymptomatic>()          = riskofinfection_wrapper;
    model_extended.parameters.get<mio::isir::RiskOfInfectionFromSymptomatic>() = riskofinfection_wrapper;

    // In ModelMessina we use beta fpr the contact rate, in ModelMessinaextended we define the contacts via the contact
    // matrix and the total population.
    model.parameters.get<mio::isir::beta>() = beta;

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    model_extended.parameters.get<mio::isir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    // Simulate.
    mio::isir::SimulationMessina sim(model, dt);
    sim.advance_messina(tmax);
    mio::TimeSeries<ScalarType> results = sim.get_result();

    mio::isir::SimulationMessinaExtended sim_extended(model_extended, dt);
    sim_extended.advance_messina(tmax);
    mio::TimeSeries<ScalarType> results_extended = sim_extended.get_result();

    // Compare results.
    for (size_t i = 0; i < (size_t)results.get_num_time_points(); i++) {
        EXPECT_NEAR(results[i][(Eigen::Index)mio::isir::InfectionState::Susceptible],
                    results_extended[i][(Eigen::Index)mio::isir::InfectionState::Susceptible], 1e-6);
    }
}

// Test that ModelMessinaExtendedDetailedInit simplifies to ModelMessinaExtended when choosing the initial conditions appropriately.
TEST(IdeSir, compareModelMessinaExtendedAndModelMessinaExtendedDetailedInit)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;
    // Define simulation parameters.
    size_t gregory_order           = 2;
    size_t finite_difference_order = 1;

    ScalarType dt = 0.1;

    ScalarType t0 = 0.;

    ScalarType tmax = 1.;

    ScalarType S0               = 90.;
    ScalarType I0               = 10.;
    ScalarType R0               = 0.;
    ScalarType total_population = S0 + I0 + R0;

    ScalarType cont_freq = 0.1;

    // Define init_populations.
    mio::TimeSeries<ScalarType> init_populations_extended((size_t)mio::isir::InfectionState::Count);
    mio::TimeSeries<ScalarType> init_populations_extended_detailed_init((size_t)mio::isir::InfectionState::Count);
    // Initialize first (gregory_order-1) time points with constant values.
    // Only set S because this is the only compartment we consider at the moment.
    Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));
    vec_init[(size_t)mio::isir::InfectionState::Susceptible] = S0;
    vec_init[(size_t)mio::isir::InfectionState::Infected]    = I0;
    vec_init[(size_t)mio::isir::InfectionState::Recovered]   = R0;

    // Add time points t0 to init_populations.
    init_populations_extended.add_time_point(t0, vec_init);
    init_populations_extended_detailed_init.add_time_point(t0, vec_init);

    // Set up ModelMessinaExtended and ModelMessinaExtendedDetailedInit.
    mio::isir::ModelMessinaExtended model_extended(std::move(init_populations_extended), total_population,
                                                   gregory_order, finite_difference_order);
    mio::isir::ModelMessinaExtendedDetailedInit model_extended_detailed_init(
        std::move(init_populations_extended_detailed_init), total_population, gregory_order, finite_difference_order);

    mio::ExponentialSurvivalFunction exp(2.);
    mio::StateAgeFunctionWrapper dist(exp);
    std::vector<mio::StateAgeFunctionWrapper> vec_dist((size_t)mio::isir::InfectionTransition::Count, dist);
    model_extended.parameters.get<mio::isir::TransitionDistributions>()               = vec_dist;
    model_extended_detailed_init.parameters.get<mio::isir::TransitionDistributions>() = vec_dist;

    mio::ConstantFunction transmissiononcontact(1.5);
    mio::StateAgeFunctionWrapper transmissiononcontact_wrapper(transmissiononcontact);
    model_extended.parameters.get<mio::isir::TransmissionProbabilityOnContact>() = transmissiononcontact_wrapper;
    model_extended_detailed_init.parameters.get<mio::isir::TransmissionProbabilityOnContact>() =
        transmissiononcontact_wrapper;

    mio::ConstantFunction riskofinfection(1.);
    mio::StateAgeFunctionWrapper riskofinfection_wrapper(riskofinfection);
    model_extended.parameters.get<mio::isir::RiskOfInfectionFromSymptomatic>()               = riskofinfection_wrapper;
    model_extended_detailed_init.parameters.get<mio::isir::RiskOfInfectionFromSymptomatic>() = riskofinfection_wrapper;

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    model_extended.parameters.get<mio::isir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);
    model_extended_detailed_init.parameters.get<mio::isir::ContactPatterns>() =
        mio::UncertainContactMatrix(contact_matrix);

    // Simulate.
    mio::isir::SimulationMessinaExtended sim_extended(model_extended, dt);
    sim_extended.advance_messina(tmax);
    mio::TimeSeries<ScalarType> results_extended = sim_extended.get_result();
    mio::TimeSeries<ScalarType> flows_extended   = sim_extended.get_flows();

    mio::isir::SimulationMessinaExtendedDetailedInit sim_extended_detailed_init(model_extended_detailed_init, dt);
    sim_extended_detailed_init.advance(tmax);
    mio::TimeSeries<ScalarType> results_extended_detailed_init = sim_extended_detailed_init.get_result();
    mio::TimeSeries<ScalarType> flows_extended_detailed_init   = sim_extended_detailed_init.get_flows();

    // Compare results.
    for (size_t i = 0; i < (size_t)results_extended.get_num_time_points(); i++) {
        for (size_t compartment = 0; compartment < (size_t)mio::isir::InfectionState::Count; compartment++) {
            EXPECT_NEAR(results_extended[i][compartment], results_extended_detailed_init[i][compartment], 1e-6);
        }
    }
}

TEST(IdeSir, testFiniteDifferenceApproximation)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType t0   = 0.;
    ScalarType tmax = 2.;

    size_t gregory_order        = 1;
    ScalarType total_population = 2.;

    std::vector<size_t> finite_difference_orders = {1, 2, 4};

    for (size_t finite_difference_order : finite_difference_orders) {
        // std::cout << "Finite diff order: " << finite_difference_order << std::endl;

        // Set values of S to sin(x) on interval from 0 to 2.

        std::vector<ScalarType> dt_exponents = {0, 1, 2, 3};

        std::vector<ScalarType> errors = {};

        for (size_t dt_exponent : dt_exponents) {

            ScalarType dt = pow(10, -(ScalarType)dt_exponent);

            mio::TimeSeries<ScalarType> init_populations((size_t)mio::isir::InfectionState::Count);

            for (size_t i = t0 / dt; i <= std::round(tmax / dt); i++) {
                Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));
                vec_init[(size_t)mio::isir::InfectionState::Susceptible] = sin(i * dt);
                vec_init[(size_t)mio::isir::InfectionState::Infected] =
                    total_population - vec_init[(size_t)mio::isir::InfectionState::Susceptible];
                // vec_init[(size_t)mio::isir::InfectionState::Recovered]   = 0.;

                // Add time points t0 to init_populations.
                init_populations.add_time_point(i * dt, vec_init);
            }

            mio::isir::ModelMessinaExtendedDetailedInit model(std::move(init_populations), total_population,
                                                              gregory_order, finite_difference_order);

            model.flows.add_time_point(
                0., mio::TimeSeries<ScalarType>::Vector::Constant((size_t)mio::isir::InfectionTransition::Count, 0.));
            model.flows.get_value(0)[(Eigen::Index)mio::isir::InfectionTransition::SusceptibleToInfected] =
                -(model.populations.get_value(1)[(Eigen::Index)mio::isir::InfectionState::Susceptible] -
                  model.populations.get_value(0)[(Eigen::Index)mio::isir::InfectionState::Susceptible]) /
                dt;
            for (size_t i = 1; i < (size_t)model.populations.get_num_time_points(); i++) {
                model.flows.add_time_point(i * dt, mio::TimeSeries<ScalarType>::Vector::Constant(
                                                       (size_t)mio::isir::InfectionTransition::Count, 0.));
                model.compute_S_deriv(dt, i);
            }
            // Append error.
            errors.push_back(
                abs(-model.flows.get_last_value()[(Eigen::Index)mio::isir::InfectionTransition::SusceptibleToInfected] -
                    cos(tmax)));
        }

        for (size_t i = 0; i < errors.size(); i++) {
            // std::cout << "error: " << errors[i] << std::endl;
        }

        // Compute order of convergence.
        for (size_t i = 0; i < errors.size() - 1; i++) {
            ScalarType order =
                log(errors[i + 1] / errors[i]) / log(pow(10, -dt_exponents[i + 1]) / pow(10, -dt_exponents[i]));
            // std::cout << "Order: " << order << std::endl;
            EXPECT_NEAR(order, finite_difference_order, 0.3);
        }
    }
}
