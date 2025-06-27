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
#include "ide_sir/parameters.h"
#include "ide_sir/simulation.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include <gtest/gtest.h>
#include <vector>

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
        std::cout << "gregory order: " << gregory_order << std::endl;
        switch (gregory_order) {
        case 1:
            gregoryWeights_sigma_expected = {{1. / 2.}, {1. / 2.}, {1. / 2.}, {1. / 2.}, {1. / 2.}, {1. / 2.}};
            gregoryWeights_omega_expected = {1. / 2., 1., 1., 1., 1., 1.};
            break;
        case 2:
            std::cout << "case 2 \n";
            gregoryWeights_sigma_expected = {{5. / 12., 14. / 12.},
                                             {5. / 12., 13. / 12.},
                                             {5. / 12., 13. / 12.},
                                             {5. / 12., 13. / 12.},
                                             {5. / 12., 13. / 12.}};
            gregoryWeights_omega_expected = {5. / 12., 13. / 12., 12. / 12., 12. / 12., 12. / 12.};
            break;
        case 3:
            std::cout << "case 3 \n";
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
            // std::cout << "n: " << n << std::endl;
            for (size_t j = 0; j < gregory_order; j++) {
                // std::cout << "j: " << j << std::endl;

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

                // std::cout << "lhs: " << gregoryWeights_sigma(row_index, column_index) << std::endl;
                // std::cout << "rhs: " << gregoryWeights_sigma_expected[n - gregory_order][j] << std::endl;
                EXPECT_EQ(gregoryWeights_sigma(row_index, column_index),
                          gregoryWeights_sigma_expected[n - gregory_order][j]);
            }
        }

        size_t weight_index;
        for (size_t n = gregory_order; n <= n_max; n++) {
            // std::cout << "n: " << n << std::endl;
            for (size_t j = gregory_order; j <= n; j++) {
                // std::cout << "j: " << j << std::endl;
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

                mio::unused(weight_index);

                // std::cout << "n-j: " << n - j << std::endl;

                EXPECT_EQ(gregoryWeights_omega(weight_index), gregoryWeights_omega_expected[n - j]);
            }

            // Define weight_indices that we want to check the values for in gregoryWeights_omega.
            std::vector<size_t> weight_indices = {0, 1, 2, 3};
        }
    }
}

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
    vec_init[(size_t)mio::isir::InfectionState::Susceptible] = 90.;
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
    model.parameters.get<mio::isir::TransmissionProbabilityOnContact>()        = riskofinfection_wrapper;
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
