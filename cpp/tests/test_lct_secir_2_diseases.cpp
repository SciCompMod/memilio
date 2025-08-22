/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Annika Jungklaus, Lena Ploetzke
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

#include "lct_secir_2_diseases/model.h"
#include "lct_secir_2_diseases/infection_state.h"
#include "lct_secir_2_diseases/parameters.h"
#include "lct_secir/model.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/data/analyze_result.h"
#include "memilio/compartments/simulation.h"
#include "load_test_data.h"
#include <vector>
#include <gtest/gtest.h>
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"

// Test confirms that default construction of an LCT2D model works.
TEST(TestLCTSecir2d, simulateDefault)
{
    using InfState = mio::lsecir2d::InfectionState;
    using LctState =
        mio::LctInfectionState<InfState, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1>;
    using Model     = mio::lsecir2d::Model<LctState>;
    ScalarType t0   = 0;
    ScalarType tmax = 1;
    ScalarType dt   = 0.1;

    Eigen::VectorX<ScalarType> init = Eigen::VectorX<ScalarType>::Constant((Eigen::Index)InfState::Count, 0);
    init[0]                         = 200;
    init[3]                         = 50; // people infected with disease a
    init[5]                         = 30;
    init[15]                        = 50; // people infected with disease b
    init[17]                        = 30;

    Model model;
    for (size_t i = 0; i < LctState::Count; i++) {
        model.populations[i] = init[i];
    }

    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
    ScalarType sum_pop = init.sum();
    for (Eigen::Index i = 0; i < result.get_num_time_points(); i++) {
        EXPECT_NEAR(sum_pop, result[i].sum(), 1e-5); // check that total pop. is constant
    }
}

/* Tests comparing the result for an LCT SECIR 2 DISEASE model (with transmission prob. 0 for one disease and 2 age groups) 
    with the result of the equivalent LCT SECIR model. */
// 1. Test: First infection for disease a (transmission prob. of disease b = 0)
TEST(TestLCTSecir2d, compareWithLCTSecir1)
{
    using InfState2d = mio::lsecir2d::InfectionState;
    using LctState2d = mio::LctInfectionState<InfState2d, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 1, 1>;
    using Model_2d   = mio::lsecir2d::Model<LctState2d, LctState2d>;

    using InfState_lct = mio::lsecir::InfectionState;
    using LctState_lct = mio::LctInfectionState<InfState_lct, 1, 1, 1, 1, 1, 1, 1, 1>;
    using Model_lct    = mio::lsecir::Model<LctState_lct, LctState_lct>;

    ScalarType t0   = 0;
    ScalarType tmax = 5;
    ScalarType dt   = 0.1;

    // Initialization vector for lct2d model.
    Eigen::VectorX<ScalarType> init_lct2d = Eigen::VectorX<ScalarType>::Constant((Eigen::Index)InfState2d::Count, 0);
    init_lct2d[0]                         = 200; // lct and lct2d use different infection states
    init_lct2d[3]                         = 50; // make sure initial pop. is in the same compartments for lct and lct2d
    init_lct2d[5]                         = 30;

    // Define LCT2D model.
    Model_2d model_lct2d;
    //Set initial values
    for (size_t i = 0; i < LctState2d::Count; i++) {
        model_lct2d.populations[i] = init_lct2d[i];
    }

    // Set Parameters.
    model_lct2d.parameters.get<mio::lsecir2d::TimeExposed_a>()[0]            = 3.2;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedNoSymptoms_a>()[0] = 2;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedSymptoms_a>()[0]   = 5.8;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedSevere_a>()[0]     = 9.5;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedCritical_a>()[0]   = 7.1;
    model_lct2d.parameters.get<mio::lsecir2d::TimeExposed_b>()[0]            = 3.2;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedNoSymptoms_b>()[0] = 2;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedSymptoms_b>()[0]   = 5.8;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedSevere_b>()[0]     = 9.5;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedCritical_b>()[0]   = 7.1;

    model_lct2d.parameters.get<mio::lsecir2d::TransmissionProbabilityOnContact_a>()[0] = 0.05;
    model_lct2d.parameters.get<mio::lsecir2d::TransmissionProbabilityOnContact_b>()[0] = 0.;

    mio::ContactMatrixGroup& contact_matrix_lct2d = model_lct2d.parameters.get<mio::lsecir2d::ContactPatterns>();
    contact_matrix_lct2d[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(2, 2, 10));
    contact_matrix_lct2d[0].add_damping(0.7, mio::SimulationTime(2.));

    model_lct2d.parameters.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_a>()[0] = 0.7;
    model_lct2d.parameters.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_a>()[0] = 0.25;
    model_lct2d.parameters.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_a>()[0] = 0.09;
    model_lct2d.parameters.get<mio::lsecir2d::SeverePerInfectedSymptoms_a>()[0]      = 0.2;
    model_lct2d.parameters.get<mio::lsecir2d::CriticalPerSevere_a>()[0]              = 0.25;
    model_lct2d.parameters.get<mio::lsecir2d::DeathsPerCritical_a>()[0]              = 0.3;
    model_lct2d.parameters.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_b>()[0] = 0.7;
    model_lct2d.parameters.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_b>()[0] = 0.25;
    model_lct2d.parameters.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_b>()[0] = 0.09;
    model_lct2d.parameters.get<mio::lsecir2d::SeverePerInfectedSymptoms_b>()[0]      = 0.2;
    model_lct2d.parameters.get<mio::lsecir2d::CriticalPerSevere_b>()[0]              = 0.25;
    model_lct2d.parameters.get<mio::lsecir2d::DeathsPerCritical_b>()[0]              = 0.3;
    model_lct2d.parameters.get<mio::lsecir2d::StartDay>()                            = 50;
    model_lct2d.parameters.get<mio::lsecir2d::Seasonality>()                         = 0.1;

    // Initialization vector for LCT model.
    Eigen::VectorX<ScalarType> init_lct = Eigen::VectorX<ScalarType>::Constant((Eigen::Index)InfState_lct::Count, 0);
    init_lct[0]                         = 200;
    init_lct[3]                         = 50;
    init_lct[5]                         = 30;

    // Define LCT model.
    Model_lct model_lct;
    //Set initial values
    for (size_t i = 0; i < LctState_lct::Count; i++) {
        model_lct.populations[i] = init_lct[i];
    }

    // Set Parameters.
    model_lct.parameters.get<mio::lsecir::TimeExposed>()[0]            = 3.2;
    model_lct.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[0] = 2;
    model_lct.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[0]   = 5.8;
    model_lct.parameters.get<mio::lsecir::TimeInfectedSevere>()[0]     = 9.5;
    model_lct.parameters.get<mio::lsecir::TimeInfectedCritical>()[0]   = 7.1;

    model_lct.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[0] = 0.05;

    mio::ContactMatrixGroup& contact_matrix_lct = model_lct.parameters.get<mio::lsecir::ContactPatterns>();
    contact_matrix_lct[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(2, 2, 10));
    contact_matrix_lct[0].add_damping(0.7, mio::SimulationTime(2.));

    model_lct.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[0] = 0.7;
    model_lct.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[0] = 0.25;
    model_lct.parameters.get<mio::lsecir::StartDay>()                          = 50;
    model_lct.parameters.get<mio::lsecir::Seasonality>()                       = 0.1;
    model_lct.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[0] = 0.09;
    model_lct.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[0]      = 0.2;
    model_lct.parameters.get<mio::lsecir::CriticalPerSevere>()[0]              = 0.25;
    model_lct.parameters.get<mio::lsecir::DeathsPerCritical>()[0]              = 0.3;

    // Simulate
    mio::TimeSeries<ScalarType> result_lct2d = mio::simulate<ScalarType, Model_2d>(
        t0, tmax, dt, model_lct2d,
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>());

    mio::TimeSeries<ScalarType> result_lct = mio::simulate<ScalarType, Model_lct>(
        t0, tmax, dt, model_lct,
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>());

    // Simulation results should be equal.
    // Compare LCT with infection 1a in LCT2D
    ASSERT_EQ(result_lct.get_num_time_points(), result_lct2d.get_num_time_points());
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(result_lct.get_time(i), result_lct2d.get_time(i), 1e-5);

        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState_lct::Susceptible],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::Susceptible], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState_lct::Exposed],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::Exposed_1a], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState_lct::InfectedNoSymptoms],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::InfectedNoSymptoms_1a], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState_lct::InfectedSymptoms],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::InfectedSymptoms_1a], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState_lct::InfectedCritical],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::InfectedCritical_1a], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState_lct::InfectedSevere],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::InfectedSevere_1a], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState_lct::Recovered],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::Recovered_a], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState_lct::Dead],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::Dead_a], 1e-5);
    }
}

// 2. Test: First infection with disease b (transmission prob. of disease a = 0)
TEST(TestLCTSecir2d, compareWithLCTSecir2)
{
    using InfState2d = mio::lsecir2d::InfectionState;
    using LctState2d = mio::LctInfectionState<InfState2d, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 1, 1>;
    using Model_2d   = mio::lsecir2d::Model<LctState2d, LctState2d>;
    ScalarType t0    = 0;
    ScalarType tmax  = 5;
    ScalarType dt    = 0.1;

    // Initialization vector for lct2d model.
    Eigen::VectorX<ScalarType> init_lct2d = Eigen::VectorX<ScalarType>::Constant((Eigen::Index)InfState2d::Count, 0);
    init_lct2d[0]                         = 200; // lct and lct2d use different infection states
    init_lct2d[15]                        = 50; // make sure initial pop. is in the same compartments for lct and lct2d
    init_lct2d[17]                        = 30;

    // Define LCT2D model.
    Model_2d model_lct2d;
    //Set initial values
    for (size_t i = 0; i < LctState2d::Count; i++) {
        model_lct2d.populations[i] = init_lct2d[i];
    }

    // Set Parameters.
    model_lct2d.parameters.get<mio::lsecir2d::TimeExposed_a>()[0]            = 3.2;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedNoSymptoms_a>()[0] = 2;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedSymptoms_a>()[0]   = 5.8;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedSevere_a>()[0]     = 9.5;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedCritical_a>()[0]   = 7.1;
    model_lct2d.parameters.get<mio::lsecir2d::TimeExposed_b>()[0]            = 3.2;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedNoSymptoms_b>()[0] = 2;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedSymptoms_b>()[0]   = 5.8;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedSevere_b>()[0]     = 9.5;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedCritical_b>()[0]   = 7.1;

    model_lct2d.parameters.get<mio::lsecir2d::TransmissionProbabilityOnContact_a>()[0] = 0.;
    model_lct2d.parameters.get<mio::lsecir2d::TransmissionProbabilityOnContact_b>()[0] = 0.05;

    mio::ContactMatrixGroup& contact_matrix_lct2d = model_lct2d.parameters.get<mio::lsecir2d::ContactPatterns>();
    contact_matrix_lct2d[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(2, 2, 10));
    contact_matrix_lct2d[0].add_damping(0.7, mio::SimulationTime(2.));

    model_lct2d.parameters.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_a>()[0] = 0.7;
    model_lct2d.parameters.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_a>()[0] = 0.25;
    model_lct2d.parameters.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_a>()[0] = 0.09;
    model_lct2d.parameters.get<mio::lsecir2d::SeverePerInfectedSymptoms_a>()[0]      = 0.2;
    model_lct2d.parameters.get<mio::lsecir2d::CriticalPerSevere_a>()[0]              = 0.25;
    model_lct2d.parameters.get<mio::lsecir2d::DeathsPerCritical_a>()[0]              = 0.3;
    model_lct2d.parameters.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_b>()[0] = 0.7;
    model_lct2d.parameters.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_b>()[0] = 0.25;
    model_lct2d.parameters.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_b>()[0] = 0.09;
    model_lct2d.parameters.get<mio::lsecir2d::SeverePerInfectedSymptoms_b>()[0]      = 0.2;
    model_lct2d.parameters.get<mio::lsecir2d::CriticalPerSevere_b>()[0]              = 0.25;
    model_lct2d.parameters.get<mio::lsecir2d::DeathsPerCritical_b>()[0]              = 0.3;
    model_lct2d.parameters.get<mio::lsecir2d::StartDay>()                            = 50;
    model_lct2d.parameters.get<mio::lsecir2d::Seasonality>()                         = 0.1;

    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, 1, 1, 1, 1, 1, 1, 1>;
    using Model    = mio::lsecir::Model<LctState, LctState>;

    // Initialization vector for LCT model.
    Eigen::VectorX<ScalarType> init_lct = Eigen::VectorX<ScalarType>::Constant((Eigen::Index)InfState::Count, 0);
    init_lct[0]                         = 200;
    init_lct[3]                         = 50;
    init_lct[5]                         = 30;

    // Define LCT model.
    Model model_lct;
    //Set initial values
    for (size_t i = 0; i < LctState::Count; i++) {
        model_lct.populations[i] = init_lct[i];
    }

    // Set Parameters.
    model_lct.parameters.get<mio::lsecir::TimeExposed>()[0]            = 3.2;
    model_lct.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[0] = 2;
    model_lct.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[0]   = 5.8;
    model_lct.parameters.get<mio::lsecir::TimeInfectedSevere>()[0]     = 9.5;
    model_lct.parameters.get<mio::lsecir::TimeInfectedCritical>()[0]   = 7.1;

    model_lct.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[0] = 0.05;

    mio::ContactMatrixGroup& contact_matrix_lct = model_lct.parameters.get<mio::lsecir::ContactPatterns>();
    contact_matrix_lct[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(2, 2, 10));
    contact_matrix_lct[0].add_damping(0.7, mio::SimulationTime(2.));

    model_lct.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[0] = 0.7;
    model_lct.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[0] = 0.25;
    model_lct.parameters.get<mio::lsecir::StartDay>()                          = 50;
    model_lct.parameters.get<mio::lsecir::Seasonality>()                       = 0.1;
    model_lct.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[0] = 0.09;
    model_lct.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[0]      = 0.2;
    model_lct.parameters.get<mio::lsecir::CriticalPerSevere>()[0]              = 0.25;
    model_lct.parameters.get<mio::lsecir::DeathsPerCritical>()[0]              = 0.3;

    // Simulate
    mio::TimeSeries<ScalarType> result_lct2d = mio::simulate<ScalarType, Model_2d>(
        t0, tmax, dt, model_lct2d,
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>());

    mio::TimeSeries<ScalarType> result_lct = mio::simulate<ScalarType, Model>(
        t0, tmax, dt, model_lct,
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>());

    // Simulation results should be equal.
    // Compare LCT with Infection 1b in LCT2D
    ASSERT_EQ(result_lct.get_num_time_points(), result_lct2d.get_num_time_points());
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(result_lct.get_time(i), result_lct2d.get_time(i), 1e-5);

        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::Susceptible],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::Susceptible], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::Exposed],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::Exposed_1b], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::InfectedNoSymptoms],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::InfectedNoSymptoms_1b], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::InfectedSymptoms],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::InfectedSymptoms_1b], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::InfectedCritical],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::InfectedCritical_1b], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::InfectedSevere],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::InfectedSevere_1b], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::Recovered],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::Recovered_b], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::Dead],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::Dead_b], 1e-5);
    }
}

// 3. Test: Second infection with disease a (transmission prob. of disease b = 0, start in Recovered_b)*/
TEST(TestLCTSecir2d, compareWithLCTSecir3)
{
    using InfState_2d = mio::lsecir2d::InfectionState;
    using LctState_2d = mio::LctInfectionState<InfState_2d, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                               1, 1, 1, 1, 1, 1>;
    using Model_2d    = mio::lsecir2d::Model<LctState_2d, LctState_2d>;

    using InfState_lct = mio::lsecir::InfectionState;
    using LctState_lct = mio::LctInfectionState<InfState_lct, 1, 1, 1, 1, 1, 1, 1, 1>;
    using Model_lct    = mio::lsecir::Model<LctState_lct, LctState_lct>;

    ScalarType t0   = 0;
    ScalarType tmax = 5;
    ScalarType dt   = 0.1;

    // Initialization vector for lct2d model.
    Eigen::VectorX<ScalarType> init_lct2d = Eigen::VectorX<ScalarType>::Constant((Eigen::Index)InfState_2d::Count, 0);
    init_lct2d[18]                        = 200; // lct and lct2d use different infection states
    init_lct2d[10]                        = 50; // make sure initial pop. is in the same compartments for lct and lct2d
    init_lct2d[12]                        = 30;

    // Define LCT2D model.
    Model_2d model_lct2d;
    //Set initial values
    for (size_t i = 0; i < LctState_2d::Count; i++) {
        model_lct2d.populations[i] = init_lct2d[i];
    }

    // Set Parameters.
    model_lct2d.parameters.get<mio::lsecir2d::TimeExposed_a>()[0]            = 3.2;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedNoSymptoms_a>()[0] = 2;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedSymptoms_a>()[0]   = 5.8;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedSevere_a>()[0]     = 9.5;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedCritical_a>()[0]   = 7.1;
    model_lct2d.parameters.get<mio::lsecir2d::TimeExposed_b>()[0]            = 3.2;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedNoSymptoms_b>()[0] = 2;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedSymptoms_b>()[0]   = 5.8;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedSevere_b>()[0]     = 9.5;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedCritical_b>()[0]   = 7.1;

    model_lct2d.parameters.get<mio::lsecir2d::TransmissionProbabilityOnContact_a>()[0] = 0.05;
    model_lct2d.parameters.get<mio::lsecir2d::TransmissionProbabilityOnContact_b>()[0] = 0.;

    mio::ContactMatrixGroup& contact_matrix_lct2d = model_lct2d.parameters.get<mio::lsecir2d::ContactPatterns>();
    contact_matrix_lct2d[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(2, 2, 10));
    contact_matrix_lct2d[0].add_damping(0.7, mio::SimulationTime(2.));

    model_lct2d.parameters.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_a>()[0] = 0.7;
    model_lct2d.parameters.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_a>()[0] = 0.25;
    model_lct2d.parameters.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_a>()[0] = 0.09;
    model_lct2d.parameters.get<mio::lsecir2d::SeverePerInfectedSymptoms_a>()[0]      = 0.2;
    model_lct2d.parameters.get<mio::lsecir2d::CriticalPerSevere_a>()[0]              = 0.25;
    model_lct2d.parameters.get<mio::lsecir2d::DeathsPerCritical_a>()[0]              = 0.3;
    model_lct2d.parameters.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_b>()[0] = 0.7;
    model_lct2d.parameters.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_b>()[0] = 0.25;
    model_lct2d.parameters.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_b>()[0] = 0.09;
    model_lct2d.parameters.get<mio::lsecir2d::SeverePerInfectedSymptoms_b>()[0]      = 0.2;
    model_lct2d.parameters.get<mio::lsecir2d::CriticalPerSevere_b>()[0]              = 0.25;
    model_lct2d.parameters.get<mio::lsecir2d::DeathsPerCritical_b>()[0]              = 0.3;
    model_lct2d.parameters.get<mio::lsecir2d::StartDay>()                            = 50;
    model_lct2d.parameters.get<mio::lsecir2d::Seasonality>()                         = 0.1;

    // Initialization vector for LCT model.
    Eigen::VectorX<ScalarType> init_lct = Eigen::VectorX<ScalarType>::Constant((Eigen::Index)InfState_lct::Count, 0);
    init_lct[0]                         = 200;
    init_lct[3]                         = 50;
    init_lct[5]                         = 30;

    // Define LCT model.
    Model_lct model_lct;
    //Set initial values
    for (size_t i = 0; i < LctState_lct::Count; i++) {
        model_lct.populations[i] = init_lct[i];
    }

    // Set Parameters.
    model_lct.parameters.get<mio::lsecir::TimeExposed>()[0]            = 3.2;
    model_lct.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[0] = 2;
    model_lct.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[0]   = 5.8;
    model_lct.parameters.get<mio::lsecir::TimeInfectedSevere>()[0]     = 9.5;
    model_lct.parameters.get<mio::lsecir::TimeInfectedCritical>()[0]   = 7.1;

    model_lct.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[0] = 0.05;

    mio::ContactMatrixGroup& contact_matrix_lct = model_lct.parameters.get<mio::lsecir::ContactPatterns>();
    contact_matrix_lct[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(2, 2, 10));
    contact_matrix_lct[0].add_damping(0.7, mio::SimulationTime(2.));

    model_lct.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[0] = 0.7;
    model_lct.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[0] = 0.25;
    model_lct.parameters.get<mio::lsecir::StartDay>()                          = 50;
    model_lct.parameters.get<mio::lsecir::Seasonality>()                       = 0.1;
    model_lct.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[0] = 0.09;
    model_lct.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[0]      = 0.2;
    model_lct.parameters.get<mio::lsecir::CriticalPerSevere>()[0]              = 0.25;
    model_lct.parameters.get<mio::lsecir::DeathsPerCritical>()[0]              = 0.3;

    // Simulate
    mio::TimeSeries<ScalarType> result_lct2d = mio::simulate<ScalarType, Model_2d>(
        t0, tmax, dt, model_lct2d,
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>());

    mio::TimeSeries<ScalarType> result_lct = mio::simulate<ScalarType, Model_lct>(
        t0, tmax, dt, model_lct,
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>());

    // Simulation results should be equal.
    // Compare LCT with Infection 2a in LCT2D
    ASSERT_EQ(result_lct.get_num_time_points(), result_lct2d.get_num_time_points());
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(result_lct.get_time(i), result_lct2d.get_time(i), 1e-5);

        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState_lct::Susceptible],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::Recovered_b], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState_lct::Exposed],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::Exposed_2a], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState_lct::InfectedNoSymptoms],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::InfectedNoSymptoms_2a], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState_lct::InfectedSymptoms],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::InfectedSymptoms_2a], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState_lct::InfectedCritical],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::InfectedCritical_2a], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState_lct::InfectedSevere],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::InfectedSevere_2a], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState_lct::Recovered],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::Recovered_ab], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState_lct::Dead],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::Dead_a], 1e-5);
    }
}

// 4. Test: Second infection with disease b (transmission prob. of disease a = 0, start in Recovered_a)*/
TEST(TestLCTSecir2d, compareWithLCTSecir4)
{
    using InfState2d = mio::lsecir2d::InfectionState;
    using LctState2d = mio::LctInfectionState<InfState2d, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 1, 1>;
    using Model_2d   = mio::lsecir2d::Model<LctState2d, LctState2d>;
    ScalarType t0    = 0;
    ScalarType tmax  = 5;
    ScalarType dt    = 0.1;

    // Initialization vector for lct2d model.
    Eigen::VectorX<ScalarType> init_lct2d = Eigen::VectorX<ScalarType>::Constant((Eigen::Index)InfState2d::Count, 0);
    init_lct2d[6]                         = 200; // lct and lct2d use different infection states
    init_lct2d[22]                        = 50; // make sure initial pop. is in the same compartments for lct and lct2d
    init_lct2d[24]                        = 30;

    // Define LCT2D model.
    Model_2d model_lct2d;
    //Set initial values
    for (size_t i = 0; i < LctState2d::Count; i++) {
        model_lct2d.populations[i] = init_lct2d[i];
    }

    // Set Parameters.
    model_lct2d.parameters.get<mio::lsecir2d::TimeExposed_a>()[0]            = 3.2;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedNoSymptoms_a>()[0] = 2;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedSymptoms_a>()[0]   = 5.8;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedSevere_a>()[0]     = 9.5;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedCritical_a>()[0]   = 7.1;
    model_lct2d.parameters.get<mio::lsecir2d::TimeExposed_b>()[0]            = 3.2;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedNoSymptoms_b>()[0] = 2;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedSymptoms_b>()[0]   = 5.8;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedSevere_b>()[0]     = 9.5;
    model_lct2d.parameters.get<mio::lsecir2d::TimeInfectedCritical_b>()[0]   = 7.1;

    model_lct2d.parameters.get<mio::lsecir2d::TransmissionProbabilityOnContact_a>()[0] = 0.;
    model_lct2d.parameters.get<mio::lsecir2d::TransmissionProbabilityOnContact_b>()[0] = 0.05;

    mio::ContactMatrixGroup& contact_matrix_lct2d = model_lct2d.parameters.get<mio::lsecir2d::ContactPatterns>();
    contact_matrix_lct2d[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(2, 2, 10));
    contact_matrix_lct2d[0].add_damping(0.7, mio::SimulationTime(2.));

    model_lct2d.parameters.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_a>()[0] = 0.7;
    model_lct2d.parameters.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_a>()[0] = 0.25;
    model_lct2d.parameters.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_a>()[0] = 0.09;
    model_lct2d.parameters.get<mio::lsecir2d::SeverePerInfectedSymptoms_a>()[0]      = 0.2;
    model_lct2d.parameters.get<mio::lsecir2d::CriticalPerSevere_a>()[0]              = 0.25;
    model_lct2d.parameters.get<mio::lsecir2d::DeathsPerCritical_a>()[0]              = 0.3;
    model_lct2d.parameters.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_b>()[0] = 0.7;
    model_lct2d.parameters.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_b>()[0] = 0.25;
    model_lct2d.parameters.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_b>()[0] = 0.09;
    model_lct2d.parameters.get<mio::lsecir2d::SeverePerInfectedSymptoms_b>()[0]      = 0.2;
    model_lct2d.parameters.get<mio::lsecir2d::CriticalPerSevere_b>()[0]              = 0.25;
    model_lct2d.parameters.get<mio::lsecir2d::DeathsPerCritical_b>()[0]              = 0.3;
    model_lct2d.parameters.get<mio::lsecir2d::StartDay>()                            = 50;
    model_lct2d.parameters.get<mio::lsecir2d::Seasonality>()                         = 0.1;

    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, 1, 1, 1, 1, 1, 1, 1>;
    using Model    = mio::lsecir::Model<LctState, LctState>;

    // Initialization vector for LCT model.
    Eigen::VectorX<ScalarType> init_lct = Eigen::VectorX<ScalarType>::Constant((Eigen::Index)InfState::Count, 0);
    init_lct[0]                         = 200;
    init_lct[3]                         = 50;
    init_lct[5]                         = 30;

    // Define LCT model.
    Model model_lct;
    //Set initial values
    for (size_t i = 0; i < LctState::Count; i++) {
        model_lct.populations[i] = init_lct[i];
    }

    // Set Parameters.
    model_lct.parameters.get<mio::lsecir::TimeExposed>()[0]            = 3.2;
    model_lct.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[0] = 2;
    model_lct.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[0]   = 5.8;
    model_lct.parameters.get<mio::lsecir::TimeInfectedSevere>()[0]     = 9.5;
    model_lct.parameters.get<mio::lsecir::TimeInfectedCritical>()[0]   = 7.1;

    model_lct.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[0] = 0.05;

    mio::ContactMatrixGroup& contact_matrix_lct = model_lct.parameters.get<mio::lsecir::ContactPatterns>();
    contact_matrix_lct[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(2, 2, 10));
    contact_matrix_lct[0].add_damping(0.7, mio::SimulationTime(2.));

    model_lct.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[0] = 0.7;
    model_lct.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[0] = 0.25;
    model_lct.parameters.get<mio::lsecir::StartDay>()                          = 50;
    model_lct.parameters.get<mio::lsecir::Seasonality>()                       = 0.1;
    model_lct.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[0] = 0.09;
    model_lct.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[0]      = 0.2;
    model_lct.parameters.get<mio::lsecir::CriticalPerSevere>()[0]              = 0.25;
    model_lct.parameters.get<mio::lsecir::DeathsPerCritical>()[0]              = 0.3;

    // Simulate
    mio::TimeSeries<ScalarType> result_lct2d = mio::simulate<ScalarType, Model_2d>(
        t0, tmax, dt, model_lct2d,
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>());

    mio::TimeSeries<ScalarType> result_lct = mio::simulate<ScalarType, Model>(
        t0, tmax, dt, model_lct,
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>());

    // Simulation results should be equal.
    // Compare LCT with Infection 2b in LCT2D
    ASSERT_EQ(result_lct.get_num_time_points(), result_lct2d.get_num_time_points());
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(result_lct.get_time(i), result_lct2d.get_time(i), 1e-5);

        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::Susceptible],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::Recovered_a], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::Exposed],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::Exposed_2b], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::InfectedNoSymptoms],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::InfectedNoSymptoms_2b], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::InfectedSymptoms],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::InfectedSymptoms_2b], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::InfectedCritical],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::InfectedCritical_2b], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::InfectedSevere],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::InfectedSevere_2b], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::Recovered],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::Recovered_ab], 1e-5);
        EXPECT_NEAR(result_lct[i][(Eigen::Index)InfState::Dead],
                    result_lct2d[i][(Eigen::Index)mio::lsecir2d::InfectionState::Dead_b], 1e-5);
    }
}

// Run the model with more than one subcompartment for states E,I,C,H,U
// and calculate the TimeSeries with no subcompartments from the result
TEST(TestLCTSecir2d, testSubcompartments)
{
    using InfState = mio::lsecir2d::InfectionState;
    using LctState =
        mio::LctInfectionState<InfState, 1, 2, 3, 3, 3, 3, 1, 1, 2, 3, 3, 3, 3, 2, 3, 3, 3, 3, 1, 1, 2, 3, 3, 3, 3, 1>;
    using Model     = mio::lsecir2d::Model<LctState>;
    ScalarType t0   = 0;
    ScalarType tmax = 1;
    ScalarType dt   = 0.1;

    // Initial population, split into subcompartments
    std::vector<std::vector<ScalarType>> init = {
        {200},        {0, 0},    {30, 10, 10}, {0, 0, 0}, {10, 10, 10}, {0, 0, 0}, {0},       {0},       {0, 0},
        {30, 10, 10}, {0, 0, 0}, {10, 10, 10}, {0, 0, 0}, {0, 0},       {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0},
        {0},          {0},       {0, 0},       {0, 0, 0}, {0, 0, 0},    {0, 0, 0}, {0, 0, 0}, {0}};

    Model model;

    // Transfer the initial values in initial_populations to the model.
    std::vector<ScalarType> flat_init;
    for (auto&& vec : init) {
        flat_init.insert(flat_init.end(), vec.begin(), vec.end());
    }
    for (size_t i = 0; i < LctState::Count; i++) {
        model.populations[i] = flat_init[i];
    }

    mio::TimeSeries<ScalarType> result                        = mio::simulate<ScalarType, Model>(t0, tmax, dt, model);
    mio::TimeSeries<ScalarType> population_no_subcompartments = model.calculate_compartments(result);
    auto interpolated_results = mio::interpolate_simulation_result(population_no_subcompartments);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
}

// Model setup to compare result with a previous output.
class ModelTestLCTSecir2d : public testing::Test
{
public:
    using InfState = mio::lsecir2d::InfectionState;
    using LctState =
        mio::LctInfectionState<InfState, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1>;
    using Model = mio::lsecir2d::Model<LctState>;

protected:
    virtual void SetUp()
    {
        // Define initial distribution of the population in the subcompartments.
        std::vector<std::vector<ScalarType>> initial_populations = {{200}, {0}, {0}, {30}, {0}, {0}, {0},  {0}, {0},
                                                                    {0},   {0}, {0}, {0},  {0}, {0}, {30}, {0}, {0},
                                                                    {0},   {0}, {0}, {0},  {0}, {0}, {0},  {0}};
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
        model->parameters.get<mio::lsecir2d::TimeExposed_a>()[0]                      = 3.2;
        model->parameters.get<mio::lsecir2d::TimeInfectedNoSymptoms_a>()[0]           = 2;
        model->parameters.get<mio::lsecir2d::TimeInfectedSymptoms_a>()[0]             = 5.8;
        model->parameters.get<mio::lsecir2d::TimeInfectedSevere_a>()[0]               = 9.5;
        model->parameters.get<mio::lsecir2d::TimeInfectedCritical_a>()[0]             = 7.1;
        model->parameters.get<mio::lsecir2d::TransmissionProbabilityOnContact_a>()[0] = 0.05;
        model->parameters.get<mio::lsecir2d::TimeExposed_b>()[0]                      = 3.2;
        model->parameters.get<mio::lsecir2d::TimeInfectedNoSymptoms_b>()[0]           = 2;
        model->parameters.get<mio::lsecir2d::TimeInfectedSymptoms_b>()[0]             = 5.8;
        model->parameters.get<mio::lsecir2d::TimeInfectedSevere_b>()[0]               = 9.5;
        model->parameters.get<mio::lsecir2d::TimeInfectedCritical_b>()[0]             = 7.1;
        model->parameters.get<mio::lsecir2d::TransmissionProbabilityOnContact_b>()[0] = 0.05;

        mio::ContactMatrixGroup& contact_matrix = model->parameters.get<mio::lsecir2d::ContactPatterns>();
        contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
        contact_matrix[0].add_damping(0.7, mio::SimulationTime(2.));

        model->parameters.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_a>()[0] = 0.7;
        model->parameters.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_a>()[0] = 0.25;
        model->parameters.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_a>()[0] = 0.09;
        model->parameters.get<mio::lsecir2d::SeverePerInfectedSymptoms_a>()[0]      = 0.2;
        model->parameters.get<mio::lsecir2d::CriticalPerSevere_a>()[0]              = 0.25;
        model->parameters.get<mio::lsecir2d::DeathsPerCritical_a>()[0]              = 0.3;
        model->parameters.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_b>()[0] = 0.7;
        model->parameters.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_b>()[0] = 0.25;
        model->parameters.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_b>()[0] = 0.09;
        model->parameters.get<mio::lsecir2d::SeverePerInfectedSymptoms_b>()[0]      = 0.2;
        model->parameters.get<mio::lsecir2d::CriticalPerSevere_b>()[0]              = 0.25;
        model->parameters.get<mio::lsecir2d::DeathsPerCritical_b>()[0]              = 0.3;
    }

    virtual void TearDown()
    {
        delete model;
    }

public:
    Model* model = nullptr;
};

// Test calculate_compartments with a TimeSeries that has an incorrect number of elements.
TEST_F(ModelTestLCTSecir2d, testCalculatePopWrongSize)
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
TEST(TestLCTSecir2d, testConstraintsParameters)
{
    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // Check for exceptions of parameters.
    mio::lsecir2d::Parameters parameters_lct2d(1);
    parameters_lct2d.get<mio::lsecir2d::TimeExposed_a>()[0]                      = 0;
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedNoSymptoms_a>()[0]           = 3.1;
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedSymptoms_a>()[0]             = 6.1;
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedSevere_a>()[0]               = 11.1;
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedCritical_a>()[0]             = 17.1;
    parameters_lct2d.get<mio::lsecir2d::TransmissionProbabilityOnContact_a>()[0] = 0.01;
    parameters_lct2d.get<mio::lsecir2d::TimeExposed_b>()[0]                      = 3.1;
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedNoSymptoms_b>()[0]           = 3.1;
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedSymptoms_b>()[0]             = 6.1;
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedSevere_b>()[0]               = 11.1;
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedCritical_b>()[0]             = 17.1;
    parameters_lct2d.get<mio::lsecir2d::TransmissionProbabilityOnContact_b>()[0] = 0.01;
    mio::ContactMatrixGroup& contact_matrix = parameters_lct2d.get<mio::lsecir2d::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));

    parameters_lct2d.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_a>()[0] = 1;
    parameters_lct2d.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_a>()[0] = 1;
    parameters_lct2d.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_a>()[0] = 0.1;
    parameters_lct2d.get<mio::lsecir2d::SeverePerInfectedSymptoms_a>()[0]      = 0.1;
    parameters_lct2d.get<mio::lsecir2d::CriticalPerSevere_a>()[0]              = 0.1;
    parameters_lct2d.get<mio::lsecir2d::DeathsPerCritical_a>()[0]              = 0.1;
    parameters_lct2d.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_b>()[0] = 1;
    parameters_lct2d.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_b>()[0] = 1;
    parameters_lct2d.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_b>()[0] = 0.1;
    parameters_lct2d.get<mio::lsecir2d::SeverePerInfectedSymptoms_b>()[0]      = 0.1;
    parameters_lct2d.get<mio::lsecir2d::CriticalPerSevere_b>()[0]              = 0.1;
    parameters_lct2d.get<mio::lsecir2d::DeathsPerCritical_b>()[0]              = 0.1;

    // Check improper TimeExposed.
    bool constraint_check = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::TimeExposed_a>()[0] = 3.1;

    parameters_lct2d.get<mio::lsecir2d::TimeExposed_b>()[0] = 0.1;
    constraint_check                                        = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::TimeExposed_b>()[0] = 3.1;

    // Check TimeInfectedNoSymptoms.
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedNoSymptoms_a>()[0] = 0.1;
    constraint_check                                                   = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedNoSymptoms_a>()[0] = 3.1;

    parameters_lct2d.get<mio::lsecir2d::TimeInfectedNoSymptoms_b>()[0] = 0.1;
    constraint_check                                                   = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedNoSymptoms_b>()[0] = 3.1;

    // Check TimeInfectedSymptoms.
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedSymptoms_a>()[0] = -0.1;
    constraint_check                                                 = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedSymptoms_a>()[0] = 6.1;

    parameters_lct2d.get<mio::lsecir2d::TimeInfectedSymptoms_b>()[0] = -0.1;
    constraint_check                                                 = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedSymptoms_b>()[0] = 6.1;

    // Check TimeInfectedSevere.
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedSevere_a>()[0] = 0.5;
    constraint_check                                               = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedSevere_a>()[0] = 11.1;

    parameters_lct2d.get<mio::lsecir2d::TimeInfectedSevere_b>()[0] = 0.5;
    constraint_check                                               = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedSevere_b>()[0] = 11.1;

    // Check TimeInfectedCritical.
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedCritical_a>()[0] = 0.;
    constraint_check                                                 = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedCritical_a>()[0] = 17.1;

    parameters_lct2d.get<mio::lsecir2d::TimeInfectedCritical_b>()[0] = 0.;
    constraint_check                                                 = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::TimeInfectedCritical_b>()[0] = 17.1;

    // Check TransmissionProbabilityOnContact.
    parameters_lct2d.get<mio::lsecir2d::TransmissionProbabilityOnContact_a>()[0] = -1;
    constraint_check                                                             = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::TransmissionProbabilityOnContact_a>()[0] = 0.01;

    parameters_lct2d.get<mio::lsecir2d::TransmissionProbabilityOnContact_b>()[0] = -1;
    constraint_check                                                             = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::TransmissionProbabilityOnContact_b>()[0] = 0.01;

    // Check RelativeTransmissionNoSymptoms.
    parameters_lct2d.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_a>()[0] = 1.5;
    constraint_check                                                           = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_a>()[0] = 1;

    parameters_lct2d.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_b>()[0] = 1.5;
    constraint_check                                                           = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_b>()[0] = 1;

    // Check RiskOfInfectionFromSymptomatic.
    parameters_lct2d.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_a>()[0] = 1.5;
    constraint_check                                                           = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_a>()[0] = 1;

    parameters_lct2d.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_b>()[0] = 1.5;
    constraint_check                                                           = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_b>()[0] = 1;

    // Check RecoveredPerInfectedNoSymptoms.
    parameters_lct2d.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_a>()[0] = 1.5;
    constraint_check                                                           = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_a>()[0] = 0.1;

    parameters_lct2d.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_b>()[0] = 1.5;
    constraint_check                                                           = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_b>()[0] = 0.1;

    // Check SeverePerInfectedSymptoms.
    parameters_lct2d.get<mio::lsecir2d::SeverePerInfectedSymptoms_a>()[0] = -1;
    constraint_check                                                      = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::SeverePerInfectedSymptoms_a>()[0] = 0.1;

    parameters_lct2d.get<mio::lsecir2d::SeverePerInfectedSymptoms_b>()[0] = -1;
    constraint_check                                                      = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::SeverePerInfectedSymptoms_b>()[0] = 0.1;

    // Check CriticalPerSevere.
    parameters_lct2d.get<mio::lsecir2d::CriticalPerSevere_a>()[0] = -1;
    constraint_check                                              = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::CriticalPerSevere_a>()[0] = 0.1;

    parameters_lct2d.get<mio::lsecir2d::CriticalPerSevere_b>()[0] = -1;
    constraint_check                                              = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::CriticalPerSevere_b>()[0] = 0.1;

    // Check DeathsPerCritical.
    parameters_lct2d.get<mio::lsecir2d::DeathsPerCritical_a>()[0] = -1;
    constraint_check                                              = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::DeathsPerCritical_a>()[0] = 0.1;

    parameters_lct2d.get<mio::lsecir2d::DeathsPerCritical_b>()[0] = -1;
    constraint_check                                              = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.get<mio::lsecir2d::DeathsPerCritical_b>()[0] = 0.1;

    // Check Seasonality.
    parameters_lct2d.set<mio::lsecir2d::Seasonality>(1);
    constraint_check = parameters_lct2d.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct2d.set<mio::lsecir2d::Seasonality>(0.1);

    // Check with correct parameters.
    constraint_check = parameters_lct2d.check_constraints();
    EXPECT_FALSE(constraint_check);

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}

// Check constraints of the Model setup.
TEST(TestLCTSecir2d, testConstraintsModel)
{
    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    using InfState = mio::lsecir2d::InfectionState;

    // Check for improper number of subcompartments for Susceptible.
    using LctStatewrongSusceptibles =
        mio::LctInfectionState<InfState, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1>;
    using ModelwrongSusceptibles = mio::lsecir2d::Model<LctStatewrongSusceptibles>;
    ModelwrongSusceptibles modelwrongSusceptibles;
    bool constraint_check = modelwrongSusceptibles.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check for improper number of subcompartments for Recovered.
    using LctStatewrongRecovered_a =
        mio::LctInfectionState<InfState, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1>;
    using ModelwrongRecovered_a = mio::lsecir2d::Model<LctStatewrongRecovered_a>;
    ModelwrongRecovered_a modelwrongRecovered_a;
    constraint_check = modelwrongRecovered_a.check_constraints();
    EXPECT_TRUE(constraint_check);

    using LctStatewrongRecovered_b =
        mio::LctInfectionState<InfState, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1>;
    using ModelwrongRecovered_b = mio::lsecir2d::Model<LctStatewrongRecovered_b>;
    ModelwrongRecovered_b modelwrongRecovered_b;
    constraint_check = modelwrongRecovered_b.check_constraints();
    EXPECT_TRUE(constraint_check);

    using LctStatewrongRecovered_ab =
        mio::LctInfectionState<InfState, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3>;
    using ModelwrongRecovered_ab = mio::lsecir2d::Model<LctStatewrongRecovered_ab>;
    ModelwrongRecovered_ab modelwrongRecovered_ab;
    constraint_check = modelwrongRecovered_ab.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check for improper number of subcompartments for Dead.
    using LctStatewrongDead_a =
        mio::LctInfectionState<InfState, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1>;
    using ModelwrongDead_a = mio::lsecir2d::Model<LctStatewrongDead_a>;
    ModelwrongDead_a modelwrongDead_a;
    constraint_check = modelwrongDead_a.check_constraints();
    EXPECT_TRUE(constraint_check);

    using LctStatewrongDead_b =
        mio::LctInfectionState<InfState, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1>;
    using ModelwrongDead_b = mio::lsecir2d::Model<LctStatewrongDead_b>;
    ModelwrongDead_b modelwrongDead_b;
    constraint_check = modelwrongDead_b.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check with a negative number in the initial population distribution.
    using LctStatevalid =
        mio::LctInfectionState<InfState, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1>;
    using Model = mio::lsecir2d::Model<LctStatevalid>;
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
