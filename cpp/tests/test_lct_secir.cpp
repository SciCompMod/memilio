/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

#include "load_test_data.h"
#include "lct_secir/model.h"
#include "lct_secir/infection_state.h"
#include "lct_secir/simulation.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/eigen.h"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"

#include "ode_secir/model.h"
#include "ode_secir/parameter_space.h"
#include "ode_secir/analyze_result.h"
#include "ode_secir/parameters.h"

#include <gtest/gtest.h>

#include <iostream>

TEST(TestLCTSecir, simulateDefault)
{
    ScalarType t0   = 0;
    ScalarType tmax = 1;
    ScalarType dt   = 0.1;

    mio::lsecir::Model model(Eigen::VectorXd::Constant((int)mio::lsecir::InfectionStateBase::Count, 15));
    mio::TimeSeries<ScalarType> result = mio::lsecir::simulate(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
}

TEST(TestLCTSecir, compareWithOdeSecir)
{
    ScalarType t0   = 0;
    ScalarType tmax = 5;
    ScalarType dt   = 0.1;

    Eigen::VectorXd init = Eigen::VectorXd::Constant((int)mio::lsecir::InfectionStateBase::Count, 15);
    init[0]              = 200;
    init[3]              = 50;
    init[5]              = 30;

    mio::lsecir::Model model_lct(init);
    // Set Parameters of the model
    model_lct.parameters.get<mio::lsecir::TimeExposed>()            = 2 * 4.2 - 5.2;
    model_lct.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>() = 2 * (5.2 - 4.2);
    model_lct.parameters.get<mio::lsecir::TimeInfectedSymptoms>()   = 5.8;
    model_lct.parameters.get<mio::lsecir::TimeInfectedSevere>()     = 9.5;
    model_lct.parameters.get<mio::lsecir::TimeInfectedCritical>()   = 7.1;

    model_lct.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.05;

    mio::ContactMatrixGroup& contact_matrix_lct = model_lct.parameters.get<mio::lsecir::ContactPatterns>();
    contact_matrix_lct[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
    contact_matrix_lct[0].add_damping(0.7, mio::SimulationTime(5.));

    model_lct.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 0.7;
    model_lct.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 0.25;
    model_lct.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.09;
    model_lct.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.2;
    model_lct.parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.25;
    model_lct.parameters.get<mio::lsecir::DeathsPerCritical>()              = 0.3;

    mio::TimeSeries<ScalarType> result_lct = mio::lsecir::simulate(
        t0, tmax, dt, model_lct,
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>());
    mio::lsecir::print_TimeSeries(result_lct);

    // Initialize ODE model with one single age group
    mio::osecir::Model model_ode(1);

    //Set population
    model_ode.populations.set_total(init.sum());
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}] =
        init[Eigen::Index(mio::lsecir::InfectionStateBase::Exposed)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] =
        init[Eigen::Index(mio::lsecir::InfectionStateBase::InfectedNoSymptoms)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}] =
        init[Eigen::Index(mio::lsecir::InfectionStateBase::InfectedSymptoms)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}] =
        init[Eigen::Index(mio::lsecir::InfectionStateBase::InfectedSevere)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}] =
        init[Eigen::Index(mio::lsecir::InfectionStateBase::InfectedCritical)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}] =
        init[Eigen::Index(mio::lsecir::InfectionStateBase::Recovered)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}] =
        init[Eigen::Index(mio::lsecir::InfectionStateBase::Dead)];
    model_ode.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                    init.sum());

    // set parameters fitting to these of the lct model
    model_ode.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0] =
        5.2; // TimeExposed = 2 * SerialInterval - IncubationTime
    model_ode.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0] =
        4.2; // TimeInfectedNoSymptoms = 2* (IncubationTime - SerialInterval)
    model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0] = 5.8;
    model_ode.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]   = 9.5;
    model_ode.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0] = 7.1;

    mio::ContactMatrixGroup& contact_matrix_ode = model_ode.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix_ode[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
    contact_matrix_ode[0].add_damping(0.7, mio::SimulationTime(5.));

    model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0] = 0.05;
    model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]   = 0.7;
    model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]   = 0.09;
    model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]   = 0.25;
    model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]        = 0.2;
    model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]                = 0.25;
    model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]                = 0.3;
    // zu R großer unterschied sowie S, E,C
    // evtl erstmal für LCT testen ob immer n0 ergibt.
    mio::TimeSeries<double> result_ode = mio::osecir::simulate(
        t0, tmax, dt, model_ode,
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>());
    mio::lsecir::print_TimeSeries(result_ode);
    ASSERT_EQ(result_lct.get_num_time_points(), result_ode.get_num_time_points());
}