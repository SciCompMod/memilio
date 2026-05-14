/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Lena Ploetzke
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

#include "lct_secir/model.h"
#include "lct_secir/infection_state.h"
#include "lct_secir/parameters.h"
#include "ode_secir/model.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/math/eigen.h"
#include "memilio/compartments/simulation.h"
#include "load_test_data.h"

#include <iostream>
#include <vector>

#include <gtest/gtest.h>
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"

// Test confirms that default construction of an LCT model works.
TEST(TestLCTSecir, simulateDefault)
{
    using InfState  = mio::lsecir::InfectionState;
    using LctState  = mio::LctInfectionState<ScalarType, InfState, 1, 1, 1, 1, 1, 1, 1, 1>;
    using Model     = mio::lsecir::Model<ScalarType, LctState>;
    ScalarType t0   = 0;
    ScalarType tmax = 1;
    ScalarType dt   = 0.1;

    Eigen::VectorX<ScalarType> init = Eigen::VectorX<ScalarType>::Constant((Eigen::Index)InfState::Count, 15);
    init[0]                         = 200;
    init[3]                         = 50;
    init[5]                         = 30;

    Model model;
    for (size_t i = 0; i < LctState::Count; i++) {
        model.populations[i] = init[i];
    }

    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
    ScalarType sum_pop = init.sum();
    for (Eigen::Index i = 0; i < result.get_num_time_points(); i++) {
        EXPECT_NEAR(sum_pop, result[i].sum(), 1e-5);
    }
}

/* Test compares the result for an LCT SECIR model with one single subcompartment for each infection state
    with the result of the equivalent ODE SECIR model. */
TEST(TestLCTSecir, compareWithOdeSecir)
{
    using InfState  = mio::lsecir::InfectionState;
    using LctState  = mio::LctInfectionState<ScalarType, InfState, 1, 1, 1, 1, 1, 1, 1, 1>;
    using Model     = mio::lsecir::Model<ScalarType, LctState>;
    ScalarType t0   = 0;
    ScalarType tmax = 5;
    ScalarType dt   = 0.1;

    // Initialization vector for both models.
    Eigen::VectorX<ScalarType> init = Eigen::VectorX<ScalarType>::Constant((Eigen::Index)InfState::Count, 15);
    init[0]                         = 200;
    init[3]                         = 50;
    init[5]                         = 30;

    // Define LCT model.
    Model model_lct;
    //Set initial values
    for (size_t i = 0; i < LctState::Count; i++) {
        model_lct.populations[i] = init[i];
    }

    // Set Parameters.
    model_lct.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[0]            = 3.2;
    model_lct.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[0] = 2;
    model_lct.parameters.template get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[0]   = 5.8;
    model_lct.parameters.template get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[0]     = 9.5;
    model_lct.parameters.template get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[0]   = 7.1;

    model_lct.parameters.template get<mio::lsecir::TransmissionProbabilityOnContact<ScalarType>>()[0] = 0.05;

    mio::ContactMatrixGroup<ScalarType>& contact_matrix_lct =
        model_lct.parameters.template get<mio::lsecir::ContactPatterns<ScalarType>>();
    contact_matrix_lct[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(1, 1, 10));
    contact_matrix_lct[0].add_damping(0.7, mio::SimulationTime<ScalarType>(2.));

    model_lct.parameters.template get<mio::lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[0] = 0.7;
    model_lct.parameters.template get<mio::lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[0] = 0.25;
    model_lct.parameters.template get<mio::lsecir::StartDay<ScalarType>>()                          = 50;
    model_lct.parameters.template get<mio::lsecir::Seasonality<ScalarType>>()                       = 0.1;
    model_lct.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[0] = 0.09;
    model_lct.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[0]      = 0.2;
    model_lct.parameters.template get<mio::lsecir::CriticalPerSevere<ScalarType>>()[0]              = 0.25;
    model_lct.parameters.template get<mio::lsecir::DeathsPerCritical<ScalarType>>()[0]              = 0.3;

    // Simulate.
    mio::TimeSeries<ScalarType> result_lct = mio::simulate<ScalarType, Model>(
        t0, tmax, dt, model_lct,
        std::make_unique<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>());

    // Initialize ODE model with one age group.
    mio::osecir::Model<ScalarType> model_ode(1);
    // Set initial distribution of the population.
    model_ode.populations[{(mio::AgeGroup)0, mio::osecir::InfectionState::Exposed}] =
        init[(Eigen::Index)InfState::Exposed];
    model_ode.populations[{(mio::AgeGroup)0, mio::osecir::InfectionState::InfectedNoSymptoms}] =
        init[(Eigen::Index)InfState::InfectedNoSymptoms];
    model_ode.populations[{(mio::AgeGroup)0, mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model_ode.populations[{(mio::AgeGroup)0, mio::osecir::InfectionState::InfectedSymptoms}] =
        init[(Eigen::Index)InfState::InfectedSymptoms];
    model_ode.populations[{(mio::AgeGroup)0, mio::osecir::InfectionState::InfectedSymptomsConfirmed}] = 0;
    model_ode.populations[{(mio::AgeGroup)0, mio::osecir::InfectionState::InfectedSevere}] =
        init[(Eigen::Index)InfState::InfectedSevere];
    model_ode.populations[{(mio::AgeGroup)0, mio::osecir::InfectionState::InfectedCritical}] =
        init[(Eigen::Index)InfState::InfectedCritical];
    model_ode.populations[{(mio::AgeGroup)0, mio::osecir::InfectionState::Recovered}] =
        init[(Eigen::Index)InfState::Recovered];
    model_ode.populations[{(mio::AgeGroup)0, mio::osecir::InfectionState::Dead}] = init[(Eigen::Index)InfState::Dead];
    model_ode.populations.set_difference_from_total({(mio::AgeGroup)0, mio::osecir::InfectionState::Susceptible},
                                                    init.sum());

    // Set parameters according to the parameters of the LCT model.
    // No restrictions by additional parameters.
    model_ode.parameters.template get<mio::osecir::TestAndTraceCapacity<ScalarType>>() =
        std::numeric_limits<ScalarType>::max();
    model_ode.parameters.template get<mio::osecir::ICUCapacity<ScalarType>>() = std::numeric_limits<ScalarType>::max();

    model_ode.parameters.set<mio::osecir::StartDay<ScalarType>>(50);
    model_ode.parameters.set<mio::osecir::Seasonality<ScalarType>>(0.1);
    model_ode.parameters.template get<mio::osecir::TimeExposed<ScalarType>>()[(mio::AgeGroup)0]            = 3.2;
    model_ode.parameters.template get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0] = 2.0;
    model_ode.parameters.template get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0]   = 5.8;
    model_ode.parameters.template get<mio::osecir::TimeInfectedSevere<ScalarType>>()[(mio::AgeGroup)0]     = 9.5;
    model_ode.parameters.template get<mio::osecir::TimeInfectedCritical<ScalarType>>()[(mio::AgeGroup)0]   = 7.1;

    mio::ContactMatrixGroup<ScalarType>& contact_matrix_ode =
        model_ode.parameters.template get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix_ode[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(1, 1, 10));
    contact_matrix_ode[0].add_damping(0.7, mio::SimulationTime<ScalarType>(2.));

    model_ode.parameters.template get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()[(mio::AgeGroup)0] =
        0.05;
    model_ode.parameters.template get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        0.7;
    model_ode.parameters.template get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        0.09;
    model_ode.parameters.template get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[(mio::AgeGroup)0] =
        0.25;
    model_ode.parameters.template get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0] = 0.2;
    model_ode.parameters.template get<mio::osecir::CriticalPerSevere<ScalarType>>()[(mio::AgeGroup)0]         = 0.25;
    model_ode.parameters.template get<mio::osecir::DeathsPerCritical<ScalarType>>()[(mio::AgeGroup)0]         = 0.3;

    // Simulate.
    mio::TimeSeries<ScalarType> result_ode = mio::osecir::simulate<ScalarType>(
        t0, tmax, dt, model_ode,
        std::make_unique<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>());

    // Simulation results should be equal.
    ASSERT_EQ(result_lct.get_num_time_points(), result_ode.get_num_time_points());
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(result_lct.get_time(i), result_ode.get_time(i), 1e-5);

        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::Susceptible],
                    result_ode[i][(Eigen::Index)mio::osecir::InfectionState::Susceptible], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::Exposed],
                    result_ode[i][(Eigen::Index)mio::osecir::InfectionState::Exposed], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::InfectedNoSymptoms],
                    result_ode[i][(Eigen::Index)mio::osecir::InfectionState::InfectedNoSymptoms], 1e-5);
        EXPECT_NEAR(0, result_ode[i][(Eigen::Index)mio::osecir::InfectionState::InfectedNoSymptomsConfirmed], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::InfectedSymptoms],
                    result_ode[i][(Eigen::Index)mio::osecir::InfectionState::InfectedSymptoms], 1e-5);
        EXPECT_NEAR(0, result_ode[i][(Eigen::Index)mio::osecir::InfectionState::InfectedSymptomsConfirmed], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::InfectedCritical],
                    result_ode[i][(Eigen::Index)mio::osecir::InfectionState::InfectedCritical], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::InfectedSevere],
                    result_ode[i][(Eigen::Index)mio::osecir::InfectionState::InfectedSevere], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::Recovered],
                    result_ode[i][(Eigen::Index)mio::osecir::InfectionState::Recovered], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::Dead],
                    result_ode[i][(Eigen::Index)mio::osecir::InfectionState::Dead], 1e-5);
    }
}

// Test if the function get_derivatives() is working using a hand calculated result.
TEST(TestLCTSecir, testEvalRightHandSide)
{
    // Define model.
    // Chose more than one subcompartment for all compartments except S, R, D so that the function is correct for all selections.
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<ScalarType, InfState, 1, 2, 3, 2, 2, 2, 1, 1>;
    using Model    = mio::lsecir::Model<ScalarType, LctState>;
    Model model;

    // Define initial population distribution in infection states, one entry per subcompartment.
    std::vector<std::vector<ScalarType>> initial_populations = {{750},    {30, 20}, {20, 10, 10}, {30, 20},
                                                                {40, 10}, {10, 20}, {20},         {10}};

    // Transfer the initial values in initial_populations to the model.
    std::vector<ScalarType> flat_initial_populations;
    for (auto&& vec : initial_populations) {
        flat_initial_populations.insert(flat_initial_populations.end(), vec.begin(), vec.end());
    }
    for (size_t i = 0; i < (size_t)LctState::Count; i++) {
        model.populations[i] = flat_initial_populations[i];
    }

    // Set parameters.
    model.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[0]            = 3.2;
    model.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[0] = 2;
    model.parameters.template get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[0]   = 5.8;
    model.parameters.template get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[0]     = 9.5;
    model.parameters.template get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[0]   = 7.1;

    model.parameters.template get<mio::lsecir::TransmissionProbabilityOnContact<ScalarType>>()[0] = 0.05;

    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        model.parameters.template get<mio::lsecir::ContactPatterns<ScalarType>>();
    contact_matrix[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(1, 1, 10));

    model.parameters.template get<mio::lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[0] = 0.7;
    model.parameters.template get<mio::lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[0] = 0.25;
    model.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[0] = 0.09;
    model.parameters.set<mio::lsecir::Seasonality<ScalarType>>(0.);
    model.parameters.set<mio::lsecir::StartDay<ScalarType>>(0);
    model.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[0] = 0.2;
    model.parameters.template get<mio::lsecir::CriticalPerSevere<ScalarType>>()[0]         = 0.25;
    model.parameters.template get<mio::lsecir::DeathsPerCritical<ScalarType>>()[0]         = 0.3;

    // Compare the result of get_derivatives() with a hand calculated result.
    size_t num_subcompartments = LctState::Count;
    Eigen::VectorX<ScalarType> dydt(num_subcompartments);
    model.get_derivatives(model.get_initial_values(), model.get_initial_values(), 0.0, dydt);

    Eigen::VectorX<ScalarType> compare(num_subcompartments);
    compare << -15.3409, -3.4091, 6.25, -17.5, 15, 0, 3.3052, 3.4483, -7.0417, 6.3158, -2.2906, -2.8169, 12.3899,
        1.6901;

    for (size_t i = 0; i < num_subcompartments; i++) {
        EXPECT_NEAR(compare[i], dydt[i], 1e-3) << " at index " << i << ".\n";
    }
}

// Test if the function get_derivatives() is working with different groups using a result calculated by hand.
TEST(TestLCTSecir, testEvalRightHandSideGroups)
{
    // Define model.
    // Choose more than one subcompartment for all compartments except S, R, D so that the function is correct for all selections.
    using InfState  = mio::lsecir::InfectionState;
    using LctState1 = mio::LctInfectionState<ScalarType, InfState, 1, 2, 3, 2, 2, 2, 1, 1>;
    using LctState2 = mio::LctInfectionState<ScalarType, InfState, 1, 1, 1, 1, 1, 1, 1, 1>;
    using Model     = mio::lsecir::Model<ScalarType, LctState1, LctState2>;
    Model model;
    size_t num_groups          = Model::num_groups;
    size_t num_subcompartments = model.populations.get_num_compartments();

    // Define initial population distribution in infection states, one entry per subcompartment.
    std::vector<std::vector<ScalarType>> initial_populations1 = {{750},    {30, 20}, {20, 10, 10}, {30, 20},
                                                                 {40, 10}, {10, 20}, {20},         {10}};
    std::vector<std::vector<ScalarType>> initial_populations2 = {{750}, {10}, {50}, {1}, {9}, {0}, {30}, {100}};

    // Transfer the initial values in initial_populations to the model.
    std::vector<ScalarType> flat_initial_populations;
    for (auto&& vec : initial_populations1) {
        flat_initial_populations.insert(flat_initial_populations.end(), vec.begin(), vec.end());
    }
    for (auto&& vec : initial_populations2) {
        flat_initial_populations.insert(flat_initial_populations.end(), vec.begin(), vec.end());
    }
    for (size_t i = 0; i < num_subcompartments; i++) {
        model.populations[i] = flat_initial_populations[i];
    }

    // Set parameters.
    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        model.parameters.template get<mio::lsecir::ContactPatterns<ScalarType>>();
    contact_matrix[0] =
        mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(num_groups, num_groups, 10));
    // Group 1.
    model.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[0]                      = 5.;
    model.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[0]           = 5.;
    model.parameters.template get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[0]             = 5.;
    model.parameters.template get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[0]               = 5.;
    model.parameters.template get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[0]             = 5.;
    model.parameters.template get<mio::lsecir::TransmissionProbabilityOnContact<ScalarType>>()[0] = 0.1;
    model.parameters.template get<mio::lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[0]   = 1.;
    model.parameters.template get<mio::lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[0]   = 1.;
    model.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[0]   = 0.1;
    model.parameters.set<mio::lsecir::Seasonality<ScalarType>>(0.);
    model.parameters.set<mio::lsecir::StartDay<ScalarType>>(0);
    model.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[0] = 0.1;
    model.parameters.template get<mio::lsecir::CriticalPerSevere<ScalarType>>()[0]         = 0.1;
    model.parameters.template get<mio::lsecir::DeathsPerCritical<ScalarType>>()[0]         = 0.1;
    // Group 2.
    model.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[1]                      = 2.;
    model.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[1]           = 2.;
    model.parameters.template get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[1]             = 2.;
    model.parameters.template get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[1]               = 2.;
    model.parameters.template get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[1]             = 2.;
    model.parameters.template get<mio::lsecir::TransmissionProbabilityOnContact<ScalarType>>()[1] = 0.1;
    model.parameters.template get<mio::lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[1]   = 0.5;
    model.parameters.template get<mio::lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[1]   = 0.5;
    model.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[1]   = 0.5;
    model.parameters.set<mio::lsecir::Seasonality<ScalarType>>(0.);
    model.parameters.set<mio::lsecir::StartDay<ScalarType>>(0);
    model.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[1] = 0.5;
    model.parameters.template get<mio::lsecir::CriticalPerSevere<ScalarType>>()[1]         = 0.5;
    model.parameters.template get<mio::lsecir::DeathsPerCritical<ScalarType>>()[1]         = 0.5;

    // Compare the result of get_derivatives() with a hand calculated result.
    Eigen::VectorX<ScalarType> dydt(num_subcompartments);
    model.get_derivatives(model.get_initial_values(), model.get_initial_values(), 0.0, dydt);

    Eigen::VectorX<ScalarType> compare(num_subcompartments);
    compare << -90.6818, 78.6818, 4, -4, 6, 0, -6.6, 4, -15.2, 12, -3.6, -4, 18.6, 0.8, -90.6818, 85.6818, -20, 12,
        -4.25, 2.25, 15, 0;

    for (size_t i = 0; i < num_subcompartments; i++) {
        EXPECT_NEAR(compare[i], dydt[i], 1e-3) << " at index " << i << ".\n";
    }
    // Also test function compute_compartments with this setup.
    mio::TimeSeries<ScalarType> result(num_subcompartments);
    // Define TimeSeries as input for the function.
    result.add_time_point(0, model.get_initial_values());
    result.add_time_point(1, model.get_initial_values() + dydt);
    mio::TimeSeries<ScalarType> population = model.calculate_compartments(result);
    // Sum of subcompartments in initial_population.
    Eigen::VectorX<ScalarType> compare_population0(2 * (size_t)InfState::Count);
    compare_population0 << 750, 50, 40, 50, 50, 30, 20, 10, 750, 10, 50, 1, 9, 0, 30, 100;
    ASSERT_EQ(compare_population0.size(), static_cast<size_t>(population.get_num_elements()));
    EXPECT_NEAR(result.get_time(0), population.get_time(0), 1e-7);
    for (size_t i = 0; i < 2 * (size_t)InfState::Count; i++) {
        EXPECT_NEAR(compare_population0[i], population.get_value(0)[i], 1e-3) << " at index " << i << ".\n";
    }
    // Sum of subcompartments in compare vector.
    Eigen::VectorX<ScalarType> compare_population1(2 * (size_t)InfState::Count);
    compare_population1 << -90.6818, 82.6818, 2, -2.6, -3.2, -7.6, 18.6, 0.8, -90.6818, 85.6818, -20, 12, -4.25, 2.25,
        15, 0;
    compare_population1 = compare_population0 + compare_population1;
    EXPECT_NEAR(result.get_time(1), population.get_time(1), 1e-7);
    for (size_t i = 0; i < 2 * (size_t)InfState::Count; i++) {
        EXPECT_NEAR(compare_population1[i], population.get_value(1)[i], 1e-3) << " at index " << i << ".\n";
    }
}

/* Test if the function get_derivatives() is working for more than one group. 
 * The parameters and LctStates are chosen, such that the sum of the groups should produce the 
 * same output as the function with one group. This is tested. */
TEST(TestLCTSecir, testEvalRightHandSideThreeGroupsEqual)
{
    // Define model.
    // Chose more than one subcompartment for all compartments except S, R, D so that the function is correct for all selections.
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<ScalarType, InfState, 1, 2, 3, 2, 2, 2, 1, 1>;
    using Model    = mio::lsecir::Model<ScalarType, LctState, LctState, LctState>;
    Model model;
    size_t num_groups          = Model::num_groups;
    size_t num_subcompartments = LctState::Count;

    // Define initial population distribution in infection states, one entry per subcompartment.
    std::vector<std::vector<ScalarType>> initial_populations = {{750},    {30, 20}, {20, 10, 10}, {30, 20},
                                                                {40, 10}, {10, 20}, {20},         {10}};

    // Transfer the initial values in initial_populations to the model.
    std::vector<ScalarType> flat_initial_populations;
    for (auto&& vec : initial_populations) {
        flat_initial_populations.insert(flat_initial_populations.end(), vec.begin(), vec.end());
    }
    for (size_t group = 0; group < num_groups; group++) {
        for (size_t i = 0; i < num_subcompartments; i++) {
            model.populations[num_subcompartments * group + i] = flat_initial_populations[i] / (ScalarType)num_groups;
        }
    }
    // Set parameters.
    for (size_t group = 0; group < num_groups; group++) {
        model.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[group]                      = 3.2;
        model.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[group]           = 2;
        model.parameters.template get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[group]             = 5.8;
        model.parameters.template get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[group]               = 9.5;
        model.parameters.template get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[group]             = 7.1;
        model.parameters.template get<mio::lsecir::TransmissionProbabilityOnContact<ScalarType>>()[group] = 0.05;

        model.parameters.template get<mio::lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[group] = 0.7;
        model.parameters.template get<mio::lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[group] = 0.25;
        model.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[group] = 0.09;
        model.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[group]      = 0.2;
        model.parameters.template get<mio::lsecir::CriticalPerSevere<ScalarType>>()[group]              = 0.25;
        model.parameters.template get<mio::lsecir::DeathsPerCritical<ScalarType>>()[group]              = 0.3;
    }
    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        model.parameters.template get<mio::lsecir::ContactPatterns<ScalarType>>();
    contact_matrix[0] = mio::ContactMatrix<ScalarType>(
        Eigen::MatrixX<ScalarType>::Constant(num_groups, num_groups, 10. / (ScalarType)num_groups));
    model.parameters.set<mio::lsecir::Seasonality<ScalarType>>(0.);
    model.parameters.set<mio::lsecir::StartDay<ScalarType>>(0);

    // Compare the result of get_derivatives() with a hand calculated result.

    Eigen::VectorX<ScalarType> dydt(num_groups * num_subcompartments);

    model.get_derivatives(model.get_initial_values(), model.get_initial_values(), 0.0, dydt);

    Eigen::VectorX<ScalarType> compare(num_subcompartments);
    compare << -15.3409, -3.4091, 6.25, -17.5, 15, 0, 3.3052, 3.4483, -7.0417, 6.3158, -2.2906, -2.8169, 12.3899,
        1.6901;
    ScalarType sunum_groups = 0;
    for (size_t i = 0; i < num_subcompartments; i++) {
        for (size_t group = 0; group < num_groups; group++) {
            sunum_groups += dydt[group * num_subcompartments + i];
        }
        EXPECT_NEAR(sunum_groups, compare[i], 1e-4) << "at subcompartment number " << i;
        sunum_groups = 0;
    }
}

// Model setup to compare result with a previous output.
class ModelTestLCTSecir : public testing::Test
{
public:
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<ScalarType, InfState, 1, 2, 3, 1, 1, 5, 1, 1>;
    using Model    = mio::lsecir::Model<ScalarType, LctState>;

protected:
    virtual void SetUp()
    {
        // Define initial distribution of the population in the subcompartments.
        std::vector<std::vector<ScalarType>> initial_populations = {{750}, {30, 20},          {20, 10, 10}, {50},
                                                                    {50},  {10, 10, 5, 3, 2}, {20},         {10}};
        model                                                    = new Model();
        // Transfer the initial values in initial_populations to the model.
        std::vector<ScalarType> flat_initial_populations;
        for (auto&& vec : initial_populations) {
            flat_initial_populations.insert(flat_initial_populations.end(), vec.begin(), vec.end());
        }
        for (size_t i = 0; i < LctState::Count; i++) {
            model->populations[i] = flat_initial_populations[i];
        }

        // Set parameters.
        model->parameters.get<mio::lsecir::TimeExposed<ScalarType>>()[0]                      = 3.2;
        model->parameters.get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[0]           = 2;
        model->parameters.get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[0]             = 5.8;
        model->parameters.get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[0]               = 9.5;
        model->parameters.get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[0]             = 7.1;
        model->parameters.get<mio::lsecir::TransmissionProbabilityOnContact<ScalarType>>()[0] = 0.05;

        mio::ContactMatrixGroup<ScalarType>& contact_matrix =
            model->parameters.get<mio::lsecir::ContactPatterns<ScalarType>>();
        contact_matrix[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(1, 1, 10));
        contact_matrix[0].add_damping(0.7, mio::SimulationTime<ScalarType>(2.));

        model->parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[0] = 0.7;
        model->parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[0] = 0.25;
        model->parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[0] = 0.09;
        model->parameters.get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[0]      = 0.2;
        model->parameters.get<mio::lsecir::CriticalPerSevere<ScalarType>>()[0]              = 0.25;
        model->parameters.get<mio::lsecir::DeathsPerCritical<ScalarType>>()[0]              = 0.3;
    }

    virtual void TearDown()
    {
        delete model;
    }

public:
    Model* model = nullptr;
};

// Test compares a simulation with the result of a previous run stored in a.csv file.
TEST_F(ModelTestLCTSecir, compareWithPreviousRun)
{
    ScalarType tmax                    = 3;
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, ModelTestLCTSecir::Model>(
        0, tmax, 0.5, *model,
        std::make_unique<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>());

    // Compare subcompartments.
    auto compare = load_test_data_csv<ScalarType>("lct-secir-subcompartments-compare.csv");

    ASSERT_EQ(compare.size(), static_cast<size_t>(result.get_num_time_points()));
    for (size_t i = 0; i < compare.size(); i++) {
        ASSERT_EQ(compare[i].size(), static_cast<size_t>(result.get_num_elements()) + 1) << "at row " << i;
        EXPECT_NEAR(result.get_time(i), compare[i][0], 1e-7) << "at row " << i;
        for (size_t j = 1; j < compare[i].size(); j++) {
            EXPECT_NEAR(result.get_value(i)[j - 1], compare[i][j], 1e-7) << " at row " << i;
        }
    }

    // Compare InfectionState compartments.
    mio::TimeSeries<ScalarType> population = model->calculate_compartments(result);
    auto compare_population                = load_test_data_csv<ScalarType>("lct-secir-compartments-compare.csv");

    ASSERT_EQ(compare_population.size(), static_cast<size_t>(population.get_num_time_points()));
    for (size_t i = 0; i < compare_population.size(); i++) {
        ASSERT_EQ(compare_population[i].size(), static_cast<size_t>(population.get_num_elements()) + 1)
            << "at row " << i;
        EXPECT_NEAR(population.get_time(i), compare_population[i][0], 1e-7) << "at row " << i;
        for (size_t j = 1; j < compare_population[i].size(); j++) {
            EXPECT_NEAR(population.get_value(i)[j - 1], compare_population[i][j], 1e-7) << " at row " << i;
        }
    }
}

/* Test compares a simulation with the result of a previous run stored in a .csv file.
 * Here, three groups with the same parametrization as the one used for the previous result is used. 
 * It is tested, if the sum of the results of the groups is equal to the result with one group.
 */
TEST(TestLCTSecir, compareWithPreviousRunThreeGroups)
{
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<ScalarType, InfState, 1, 2, 3, 1, 1, 5, 1, 1>;
    using Model    = mio::lsecir::Model<ScalarType, LctState, LctState, LctState>;

    // Initialize a model.
    Model model;
    size_t num_groups          = Model::num_groups;
    size_t num_subcompartments = (size_t)LctState::Count;

    // Define initial distribution of the population in the subcompartments.
    std::vector<std::vector<ScalarType>> initial_populations = {{750}, {30, 20},          {20, 10, 10}, {50},
                                                                {50},  {10, 10, 5, 3, 2}, {20},         {10}};
    // Transfer the initial values in initial_populations to the model.
    std::vector<ScalarType> flat_initial_populations;
    for (auto&& vec : initial_populations) {
        flat_initial_populations.insert(flat_initial_populations.end(), vec.begin(), vec.end());
    }
    for (size_t group = 0; group < num_groups; group++) {
        for (size_t i = 0; i < LctState::Count; i++) {
            model.populations[num_subcompartments * group + i] = flat_initial_populations[i] / (ScalarType)num_groups;
        }
    }
    // Set parameters.
    for (size_t group = 0; group < num_groups; group++) {
        model.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[group]                      = 3.2;
        model.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[group]           = 2;
        model.parameters.template get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[group]             = 5.8;
        model.parameters.template get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[group]               = 9.5;
        model.parameters.template get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[group]             = 7.1;
        model.parameters.template get<mio::lsecir::TransmissionProbabilityOnContact<ScalarType>>()[group] = 0.05;

        model.parameters.template get<mio::lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[group] = 0.7;
        model.parameters.template get<mio::lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[group] = 0.25;
        model.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[group] = 0.09;
        model.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[group]      = 0.2;
        model.parameters.template get<mio::lsecir::CriticalPerSevere<ScalarType>>()[group]              = 0.25;
        model.parameters.template get<mio::lsecir::DeathsPerCritical<ScalarType>>()[group]              = 0.3;
    }
    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        model.parameters.template get<mio::lsecir::ContactPatterns<ScalarType>>();
    contact_matrix[0] = mio::ContactMatrix<ScalarType>(
        Eigen::MatrixX<ScalarType>::Constant(num_groups, num_groups, 10. / (ScalarType)num_groups));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime<ScalarType>(2.));

    // Compare just one step as the adaptive methods does not necessarily produce the same step sizes.
    auto compare    = load_test_data_csv<ScalarType>("lct-secir-subcompartments-compare.csv");
    ScalarType dt   = compare[1][0] - compare[0][0];
    ScalarType tmax = 2 * dt;

    auto integrator =
        std::make_unique<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // Choose dt_min = dt_max so that we have a fixed time step and can compare to the result with one group.
    integrator->set_dt_min(dt);
    integrator->set_dt_max(dt);
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(0, tmax, dt, model, std::move(integrator));

    // Compare subcompartments.
    EXPECT_NEAR(result.get_time(1), compare[1][0], 1e-7);
    ScalarType sunum_groups = 0;
    for (size_t j = 1; j < compare[1].size(); j++) {
        for (size_t group = 0; group < num_groups; group++) {
            sunum_groups += result.get_value(1)[group * num_subcompartments + j - 1];
        }
        EXPECT_NEAR(sunum_groups, compare[1][j], 1e-4) << "at subcompartment number " << j;
        sunum_groups = 0;
    }

    // Compare InfectionState compartments.
    mio::TimeSeries<ScalarType> population = model.calculate_compartments(result);
    auto compare_population                = load_test_data_csv<ScalarType>("lct-secir-compartments-compare.csv");
    EXPECT_NEAR(population.get_time(1), compare_population[1][0], 1e-7);
    for (size_t j = 1; j < compare_population[1].size(); j++) {
        for (size_t group = 0; group < num_groups; group++) {
            sunum_groups += population.get_value(1)[group * (size_t)InfState::Count + j - 1];
        }
        EXPECT_NEAR(sunum_groups, compare_population[1][j], 1e-4) << "at compartment number " << j;
        sunum_groups = 0;
    }
}

// Test calculate_compartments with a TimeSeries that has an incorrect number of elements.
TEST_F(ModelTestLCTSecir, testCalculatePopWrongSize)
{
    // Deactivate temporarily log output because an error is expected.
    mio::set_log_level(mio::LogLevel::off);
    // TimeSeries has to have LctState::Count elements.
    size_t wrong_size = LctState::Count - 2;
    // Define TimeSeries with wrong_size elements.
    mio::TimeSeries<ScalarType> wrong_num_elements(wrong_size);
    Eigen::VectorX<ScalarType> vec_wrong_size = Eigen::VectorX<ScalarType>::Ones(wrong_size);
    wrong_num_elements.add_time_point(-10, vec_wrong_size);
    wrong_num_elements.add_time_point(-9, vec_wrong_size);
    // Call the calculate_compartments function with the TimeSeries with a wrong number of elements.
    mio::TimeSeries<ScalarType> population = model->calculate_compartments(wrong_num_elements);
    // A TimeSeries of the right size with values -1 is expected.
    ASSERT_EQ(1, population.get_num_time_points());
    for (int i = 0; i < population.get_num_elements(); i++) {
        EXPECT_EQ(-1, population.get_last_value()[i]);
    }
    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}

//Check constraints of Parameters class.
TEST(TestLCTSecir, testConstraintsParameters)
{
    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // Check for exceptions of parameters.
    mio::lsecir::Parameters<ScalarType> parameters_lct(1);
    parameters_lct.get<mio::lsecir::TimeExposed<ScalarType>>()[0]                      = 0;
    parameters_lct.get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[0]           = 3.1;
    parameters_lct.get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[0]             = 6.1;
    parameters_lct.get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[0]               = 11.1;
    parameters_lct.get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[0]             = 17.1;
    parameters_lct.get<mio::lsecir::TransmissionProbabilityOnContact<ScalarType>>()[0] = 0.01;
    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        parameters_lct.get<mio::lsecir::ContactPatterns<ScalarType>>();
    contact_matrix[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(1, 1, 10));

    parameters_lct.get<mio::lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[0] = 1;
    parameters_lct.get<mio::lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[0] = 1;
    parameters_lct.get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[0] = 0.1;
    parameters_lct.get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[0]      = 0.1;
    parameters_lct.get<mio::lsecir::CriticalPerSevere<ScalarType>>()[0]              = 0.1;
    parameters_lct.get<mio::lsecir::DeathsPerCritical<ScalarType>>()[0]              = 0.1;

    // Check improper TimeExposed.
    bool constraint_check = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::TimeExposed<ScalarType>>()[0] = 3.1;

    // Check TimeInfectedNoSymptoms.
    parameters_lct.get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[0] = 0.1;
    constraint_check                                                         = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[0] = 3.1;

    // Check TimeInfectedSymptoms.
    parameters_lct.get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[0] = -0.1;
    constraint_check                                                       = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[0] = 6.1;

    // Check TimeInfectedSevere.
    parameters_lct.get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[0] = 0.5;
    constraint_check                                                     = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[0] = 11.1;

    // Check TimeInfectedCritical.
    parameters_lct.get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[0] = 0.;
    constraint_check                                                       = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[0] = 17.1;

    // Check TransmissionProbabilityOnContact.
    parameters_lct.get<mio::lsecir::TransmissionProbabilityOnContact<ScalarType>>()[0] = -1;
    constraint_check = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::TransmissionProbabilityOnContact<ScalarType>>()[0] = 0.01;

    // Check RelativeTransmissionNoSymptoms.
    parameters_lct.get<mio::lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[0] = 1.5;
    constraint_check = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[0] = 1;

    // Check RiskOfInfectionFromSymptomatic.
    parameters_lct.get<mio::lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[0] = 1.5;
    constraint_check = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[0] = 1;

    // Check RecoveredPerInfectedNoSymptoms.
    parameters_lct.get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[0] = 1.5;
    constraint_check = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[0] = 0.1;

    // Check SeverePerInfectedSymptoms.
    parameters_lct.get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[0] = -1;
    constraint_check                                                            = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[0] = 0.1;

    // Check CriticalPerSevere.
    parameters_lct.get<mio::lsecir::CriticalPerSevere<ScalarType>>()[0] = -1;
    constraint_check                                                    = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::CriticalPerSevere<ScalarType>>()[0] = 0.1;

    // Check DeathsPerCritical.
    parameters_lct.get<mio::lsecir::DeathsPerCritical<ScalarType>>()[0] = -1;
    constraint_check                                                    = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::DeathsPerCritical<ScalarType>>()[0] = 0.1;

    // Check Seasonality.
    parameters_lct.set<mio::lsecir::Seasonality<ScalarType>>(1);
    constraint_check = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.set<mio::lsecir::Seasonality<ScalarType>>(0.1);

    // Check with correct parameters.
    constraint_check = parameters_lct.check_constraints();
    EXPECT_FALSE(constraint_check);

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}

// Check constraints of the Model setup.
TEST(TestLCTSecir, testConstraintsModel)
{
    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    using InfState  = mio::lsecir::InfectionState;
    using LctState1 = mio::LctInfectionState<ScalarType, InfState, 1, 1, 1, 1, 1, 1, 1, 1>;

    // Check for improper number of subcompartments for Susceptible.
    using LctStatewrongSusceptibles = mio::LctInfectionState<ScalarType, InfState, 3, 1, 1, 1, 1, 1, 1, 1>;
    using ModelwrongSusceptibles    = mio::lsecir::Model<ScalarType, LctState1, LctStatewrongSusceptibles>;
    ModelwrongSusceptibles modelwrongSusceptibles;
    bool constraint_check = modelwrongSusceptibles.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check for improper number of subcompartments for Recovered.
    using LctStatewrongRecovered = mio::LctInfectionState<ScalarType, InfState, 1, 1, 1, 1, 1, 1, 10, 1>;
    using ModelwrongRecovered    = mio::lsecir::Model<ScalarType, LctState1, LctStatewrongRecovered>;
    ModelwrongRecovered modelwrongRecovered;
    constraint_check = modelwrongRecovered.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check for improper number of subcompartments for Dead.
    using LctStatewrongDead = mio::LctInfectionState<ScalarType, InfState, 1, 1, 1, 1, 1, 1, 1, 2>;
    using ModelwrongDead    = mio::lsecir::Model<ScalarType, LctState1, LctStatewrongDead>;
    ModelwrongDead modelwrongDead;
    constraint_check = modelwrongDead.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check with a negative number in the initial population distribution.
    using LctStatevalid = mio::LctInfectionState<ScalarType, InfState, 1, 2, 3, 2, 2, 2, 1, 1>;
    using Model         = mio::lsecir::Model<ScalarType, LctState1, LctStatevalid>;
    Model model;
    model.populations[0] = -1000;
    constraint_check     = model.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);

    // Check for valid Setup.
    model.populations[0] = 1000;
    constraint_check     = model.check_constraints();
    EXPECT_FALSE(constraint_check);
}
