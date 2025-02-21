/* 
* Copyright (C) 2020-2025 MEmilio
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

#include "memilio/config.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"
#include "memilio/compartments/simulation.h"

#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <iostream>
#include <vector>
#include <omp.h>

namespace params
{
// num_subcompartments is used as a template argument and has to be a constexpr.
constexpr size_t num_subcompartments = NUM_SUBCOMPARTMENTS;

// Define (non-age-resolved) parameters.
const ScalarType dt                             = 0.01;
const ScalarType seasonality                    = 0.;
const ScalarType relativeTransmissionNoSymptoms = 1.;
const ScalarType riskOfInfectionFromSymptomatic = 0.3;
const ScalarType total_population               = 83155031.0;

const ScalarType timeExposed                      = 3.335;
const ScalarType timeInfectedNoSymptoms           = 2.58916;
const ScalarType timeInfectedSymptoms             = 6.94547;
const ScalarType timeInfectedSevere               = 7.28196;
const ScalarType timeInfectedCritical             = 13.066;
const ScalarType transmissionProbabilityOnContact = 0.07333;
const ScalarType recoveredPerInfectedNoSymptoms   = 0.206901;
const ScalarType severePerInfectedSymptoms        = 0.07864;
const ScalarType criticalPerSevere                = 0.17318;
const ScalarType deathsPerCritical                = 0.21718;
} // namespace params

/** 
* @brief Performs multiple simulations with one model to get an average run time.
*
*   This function measures the run time taken for the simulation execution.
*   The model setup is not included in the run time. The run time is averaged over several runs.
*   The simulation uses (non-age-resolved) LCT models with Covid-19 inspired parameters and a contact rate for Germany.
*   The initial values are set to some realistic values. 
*
* @param[in] num_runs Number of runs with run time measurement. 
* @param[in] num_warm_up_runs Number of warm-up runs before actual measurements begin. 
*       These runs are used to allow the system to stabilize.
* @param[in] tmax Time horizon of the simulation.
* @param[in] use_adaptive_solver Determines whether to use an adaptive solver. If false, a fixed step size is used.
*            Default is false. 
*/
void simulate(size_t num_runs, size_t num_warm_up_runs, ScalarType tmax, bool use_adaptive_solver = false)
{
    using namespace params;
    std::cout << "{ \"Subcompartments\": " << num_subcompartments << ", " << std::endl;
    // Initialize (non-age-resolved) LCT model.
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, num_subcompartments, num_subcompartments, num_subcompartments,
                                            num_subcompartments, num_subcompartments, 1, 1>;
    using Model    = mio::lsecir::Model<LctState>;
    Model model;

    // Set parameters.
    model.parameters.get<mio::lsecir::TimeExposed>()[0]                      = timeExposed;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[0]           = timeInfectedNoSymptoms;
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[0]             = timeInfectedSymptoms;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()[0]               = timeInfectedSevere;
    model.parameters.get<mio::lsecir::TimeInfectedCritical>()[0]             = timeInfectedCritical;
    model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[0] = transmissionProbabilityOnContact;

    model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[0] = relativeTransmissionNoSymptoms;
    model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[0] = riskOfInfectionFromSymptomatic;
    model.parameters.get<mio::lsecir::Seasonality>()                       = seasonality;

    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[0] = recoveredPerInfectedNoSymptoms;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[0]      = severePerInfectedSymptoms;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()[0]              = criticalPerSevere;
    model.parameters.get<mio::lsecir::DeathsPerCritical>()[0]              = deathsPerCritical;
    // Realistic average number of contacts.
    mio::ContactMatrixGroup contact_matrix               = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                                    = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 7.69129));
    model.parameters.get<mio::lsecir::ContactPatterns>() = mio::UncertainContactMatrix<ScalarType>(contact_matrix);

    // Set initial values.
    // Vector is a "random vector" taken from another example. Just need some realistic values.
    std::vector<ScalarType> init = {3966564.2110, 664.2367, 545.5523, 1050.3946, 5.6045, 0.5844, 307.4165, 0.};
    // Use init as a basis to define appropriate initial values.
    // Compartment values are distributed equally to subcompartments.
    model.populations[0]                   = init[0]; // Susceptible.
    model.populations[LctState::Count - 2] = init[6]; // Recovered.
    model.populations[LctState::Count - 1] = init[7]; // Dead.
    for (size_t i = 1; i < (size_t)InfState::Count - 2; i++) {
        for (size_t subcomp = 0; subcomp < num_subcompartments; subcomp++) {
            model.populations[(i - 1) * num_subcompartments + 1 + subcomp] = init[i] / (ScalarType)num_subcompartments;
        }
    }

    // Set integrator of fifth order.
    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    if (!use_adaptive_solver) {
        // Choose dt_min = dt_max to get a fixed step size.
        integrator->set_dt_min(dt);
        integrator->set_dt_max(dt);
    }

    // Perform simulation several times and measure run times.
    // Warm up runs.
    mio::set_log_level(mio::LogLevel::off);
    for (size_t i = 0; i < num_warm_up_runs; i++) {
        mio::simulate<ScalarType, Model>(0, tmax, dt, model, integrator);
    }
    // Simulate one time to track the number of steps.
    auto result = mio::simulate<ScalarType, Model>(0, tmax, dt, model, integrator);
    std::cout << "\"Steps\": " << result.get_num_time_points() << "," << std::endl;

    // Runs with timing.
    ScalarType total = 0;
    for (size_t i = 0; i < num_runs; i++) {
        total -= omp_get_wtime();
        mio::simulate<ScalarType, Model>(0, tmax, dt, model, integrator);
        total += omp_get_wtime();
    }
    std::cout << "\"Time\": " << total / num_runs << "\n}," << std::endl;
    mio::set_log_level(mio::LogLevel::warn);
}

/**
* Usage: lct_timing <num_runs> <num_warm_up_runs> <use_adaptive_solver> 
*   All command line arguments are optional. Simple default values are provided if not specified.
*   All parameters are passed to the simulation() function. See the documentation for a description of the parameters.
*
*  The numbers of subcompartments used in the LCT model is determined by the preprocessor macro NUM_SUBCOMPARTMENTS.
*   You can set the number via the flag -DNUM_SUBCOMPARTMENTS=... . 
*/
int main(int argc, char** argv)
{
    const ScalarType tmax    = 20;
    size_t num_runs          = 100;
    size_t num_warm_up_runs  = 10;
    bool use_adaptive_solver = false;

    switch (argc) {
    case 4:
        use_adaptive_solver = std::stoi(argv[3]);
        [[fallthrough]];
    case 3:
        num_warm_up_runs = std::stod(argv[2]);
        [[fallthrough]];
    case 2:
        num_runs = std::stod(argv[1]);
    }
    simulate(num_runs, num_warm_up_runs, tmax, use_adaptive_solver);
    return 0;
}
