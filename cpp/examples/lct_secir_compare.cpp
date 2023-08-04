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

#include "boost/filesystem/operations.hpp"
#include "lct_secir/model.h"
#include "lct_secir/infection_state.h"
#include "lct_secir/simulation.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/eigen.h"
#include "memilio/io/result_io.h"

#include "ode_secir/model.h"
#include <string>
#include <vector>
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "memilio/io/io.h"
#include "boost/filesystem.hpp"

int main()
{ /* Runs a simulation with comparable parameters and initial data for a LCT SECIR model and an ODE model. 
    Results are automatically stored in hdf5 files.*/
    bool print = true;
    bool save  = true;

    ScalarType t0   = 0;
    ScalarType tmax = 20;
    ScalarType dt   = 0.1;
    Eigen::VectorXd init((int)mio::osecir::InfectionState::Count);
    init << 7500, 90, 50, 70, 18, 8, 0, 0;
    init = init * (10000 / init.sum());

    // ---Perform simulation for the LCT SECIR model.---
    // Set vector that specifies the number of subcompartments
    std::vector<int> num_subcompartments((int)mio::lsecir::InfectionStateBase::Count, 1);
    num_subcompartments[(int)mio::lsecir::InfectionStateBase::Exposed]            = 20;
    num_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] = 20;
    num_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms]   = 20;
    num_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSevere]     = 20;
    num_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedCritical]   = 20;
    mio::lsecir::InfectionState infectionStates(num_subcompartments);

    // define initial population distribution in infection states, one entry per Subcompartment
    Eigen::VectorXd init_lct(infectionStates.get_count());
    // if jumpstart is true, lct model will be initialized just with individuals in the first subcompartments
    bool jumpstart = false;
    if (jumpstart) {
        for (int i = 0; i < (int)mio::lsecir::InfectionStateBase::Count; i++) {
            init_lct[infectionStates.get_firstindex(i)] = init[i];
            for (int j = infectionStates.get_firstindex(i) + 1;
                 j < infectionStates.get_firstindex(i) + infectionStates.get_number(i); j++) {
                init_lct[j] = 0;
            }
        }
    }
    else {
        for (int i = 0; i < (int)mio::lsecir::InfectionStateBase::Count; i++) {
            for (int j = infectionStates.get_firstindex(i);
                 j < infectionStates.get_firstindex(i) + infectionStates.get_number(i); j++) {
                init_lct[j] = init[i] / infectionStates.get_number(i);
            }
        }
    }

    // initialize model
    mio::lsecir::Model model_lct(std::move(init_lct), infectionStates);

    // Set Parameters of the model
    model_lct.parameters.get<mio::lsecir::TimeExposed>()            = 2 * 4.2 - 5.2;
    model_lct.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>() = 2 * (5.2 - 4.2);
    model_lct.parameters.get<mio::lsecir::TimeInfectedSymptoms>()   = 5.8;
    model_lct.parameters.get<mio::lsecir::TimeInfectedSevere>()     = 9.5;
    model_lct.parameters.get<mio::lsecir::TimeInfectedCritical>()   = 7.1;

    model_lct.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.05;

    mio::ContactMatrixGroup& contact_matrix = model_lct.parameters.get<mio::lsecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 9));

    model_lct.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 0.7;
    model_lct.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 0.25;
    model_lct.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.09;
    model_lct.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.2;
    model_lct.parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.25;
    model_lct.parameters.get<mio::lsecir::DeathsPerCritical>()              = 0.3;

    // Perform simulation.
    mio::TimeSeries<ScalarType> result_lct = mio::lsecir::simulate(
        t0, tmax, dt, model_lct,
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>());
    // Calculate the distribution in Infectionstates without subcompartments of the result.
    mio::TimeSeries<ScalarType> populations_lct = model_lct.calculate_populations(result_lct);

    // Save result in HDF5 file
    if (save) {
        auto save_result_status_subcompartments =
            mio::save_result({result_lct}, {0}, 1, "result_lct_subcompartments.h5");
        auto save_result_status_lct = mio::save_result({populations_lct}, {0}, 1, "result_lct.h5");
    }
    if (print) {
        mio::lsecir::print_TimeSeries(populations_lct, model_lct.get_heading_CompartmentsBase());
    }

    // ---Perform simulation for the ODE SECIR model.---
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
    // no restrictions by additional parameters
    model_ode.parameters.set<mio::osecir::StartDay>(0);
    model_ode.parameters.set<mio::osecir::Seasonality>(0);
    model_ode.parameters.get<mio::osecir::TestAndTraceCapacity>() = std::numeric_limits<double>::max();
    model_ode.parameters.get<mio::osecir::ICUCapacity>()          = std::numeric_limits<double>::max();

    model_ode.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0] =
        5.2; // TimeExposed = 2 * SerialInterval - IncubationTime
    model_ode.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0] =
        4.2; // TimeInfectedNoSymptoms = 2* (IncubationTime - SerialInterval)
    model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0] = 5.8;
    model_ode.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]   = 9.5;
    model_ode.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0] = 7.1;

    mio::ContactMatrixGroup& contact_matrix_ode = model_ode.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix_ode[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 9));

    model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0] = 0.05;
    model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]   = 0.7;
    model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]   = 0.09;
    model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]   = 0.25;
    model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]        = 0.2;
    model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]                = 0.25;
    model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]                = 0.3;

    mio::TimeSeries<double> result_ode =
        simulate(t0, tmax, dt, model_ode,
                 std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>());

    // Save result in HDF5 file
    if (save) {
        auto save_result_status_ode = mio::save_result({result_ode}, {0}, 1, "result_ode.h5");
    }
    if (print) {
        mio::lsecir::print_TimeSeries(result_ode, "S | E | C | I | H | U | R | D");
    }
}