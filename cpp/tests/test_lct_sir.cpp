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

#include "lct_sir/model.h"
#include "lct_sir/infection_state.h"
#include "lct_sir/parameters.h"
#include "ode_sir/model.h"
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
TEST(TestLCTSir, simulateDefault)
{
    using InfState  = mio::lsir::InfectionState;
    using LctState  = mio::LctInfectionState<ScalarType, InfState, 1, 1, 1>;
    using Model     = mio::lsir::Model<ScalarType, LctState>;
    ScalarType t0   = 0;
    ScalarType tmax = 1;
    ScalarType dt   = 0.1;

    Eigen::VectorX<ScalarType> init = Eigen::VectorX<ScalarType>::Constant((Eigen::Index)InfState::Count, 15);
    // init[0]                         = 200;
    // init[3]                         = 50;
    // init[5]                         = 30;

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
TEST(TestLCTSir, compareWithOdeSir)
{
    using InfState  = mio::lsir::InfectionState;
    using LctState  = mio::LctInfectionState<ScalarType, InfState, 1, 1, 1>;
    using Model     = mio::lsir::Model<ScalarType, LctState>;
    ScalarType t0   = 0;
    ScalarType tmax = 5;
    ScalarType dt   = 0.1;

    // Initialization vector for both models.
    Eigen::VectorX<ScalarType> init = Eigen::VectorX<ScalarType>::Constant((Eigen::Index)InfState::Count, 15);
    // init[0]                         = 200;
    // init[3]                         = 50;
    // init[5]                         = 30;

    // Define LCT model.
    Model model_lct;
    //Set initial values
    for (size_t i = 0; i < LctState::Count; i++) {
        model_lct.populations[i] = init[i];
    }

    // Set Parameters.
    model_lct.parameters.template get<mio::lsir::TimeInfected<ScalarType>>()[0] = 3.2;

    model_lct.parameters.template get<mio::lsir::TransmissionProbabilityOnContact<ScalarType>>()[0] = 0.05;

    mio::ContactMatrixGroup<ScalarType>& contact_matrix_lct =
        model_lct.parameters.template get<mio::lsir::ContactPatterns<ScalarType>>();
    contact_matrix_lct[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(1, 1, 10));
    contact_matrix_lct[0].add_damping(0.7, mio::SimulationTime<ScalarType>(2.));

    model_lct.parameters.template get<mio::lsir::RiskOfInfectionFromSymptomatic<ScalarType>>()[0] = 1.;

    // Simulate.
    mio::TimeSeries<ScalarType> result_lct = mio::simulate<ScalarType, Model>(
        t0, tmax, dt, model_lct,
        std::make_unique<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>());

    // Initialize ODE model with one age group.
    mio::osir::Model<ScalarType> model_ode(1);
    // Set initial distribution of the population.
    model_ode.populations[{(mio::AgeGroup)0, mio::osir::InfectionState::Infected}] =
        init[(Eigen::Index)InfState::Infected];
    model_ode.populations[{(mio::AgeGroup)0, mio::osir::InfectionState::Recovered}] =
        init[(Eigen::Index)InfState::Recovered];
    model_ode.populations.set_difference_from_total({(mio::AgeGroup)0, mio::osir::InfectionState::Susceptible},
                                                    init.sum());

    // Set parameters according to the parameters of the LCT model.
    model_ode.parameters.template get<mio::osir::TimeInfected<ScalarType>>()[(mio::AgeGroup)0] = 3.2;

    mio::ContactMatrixGroup<ScalarType>& contact_matrix_ode =
        model_ode.parameters.template get<mio::osir::ContactPatterns<ScalarType>>();
    contact_matrix_ode[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(1, 1, 10));
    contact_matrix_ode[0].add_damping(0.7, mio::SimulationTime<ScalarType>(2.));

    model_ode.parameters.template get<mio::osir::TransmissionProbabilityOnContact<ScalarType>>()[(mio::AgeGroup)0] =
        0.05;
    // model_ode.parameters.template get<mio::osir::RiskOfInfectionFromSymptomatic<ScalarType>>()[(mio::AgeGroup)0] = 0.25;

    // Simulate.
    mio::TimeSeries<ScalarType> result_ode = mio::simulate<ScalarType>(
        t0, tmax, dt, model_ode,
        std::make_unique<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>());

    // Simulation results should be equal.
    ASSERT_EQ(result_lct.get_num_time_points(), result_ode.get_num_time_points());
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(result_lct.get_time(i), result_ode.get_time(i), 1e-5);

        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::Susceptible],
                    result_ode[i][(Eigen::Index)mio::osir::InfectionState::Susceptible], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::Infected],
                    result_ode[i][(Eigen::Index)mio::osir::InfectionState::Infected], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::Recovered],
                    result_ode[i][(Eigen::Index)mio::osir::InfectionState::Recovered], 1e-5);
    }
}
