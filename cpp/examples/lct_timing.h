/* 
* Copyright (C) 2020-2024 MEmilio
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
#include "lct_secir/initializer_flows.h"

#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"
#include "memilio/compartments/simulation.h"
#include "memilio/data/analyze_result.h"

#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <string>
#include <iostream>
#include <vector>
#include <omp.h>

namespace params
{
constexpr size_t num_groups = 6;

// Parameters
const ScalarType seasonality                    = 0.;
const ScalarType RelativeTransmissionNoSymptoms = 1.;
const ScalarType RiskOfInfectionFromSymptomatic = 0.3;
const ScalarType tmax                           = 100;

const ScalarType dt                = 0.01;
const ScalarType age_group_sizes[] = {3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434};
const ScalarType total_population  = 83155031.0;

const ScalarType TransmissionProbabilityOnContact[] = {0.03, 0.06, 0.06, 0.06, 0.09, 0.175};

const ScalarType TimeExposed[]            = {3.335, 3.335, 3.335, 3.335, 3.335, 3.335};
const ScalarType TimeInfectedNoSymptoms[] = {2.74, 2.74, 2.565, 2.565, 2.565, 2.565};
const ScalarType TimeInfectedSymptoms[]   = {7.02625, 7.02625, 7.0665, 6.9385, 6.835, 6.775};
const ScalarType TimeInfectedSevere[]     = {5., 5., 5.925, 7.55, 8.5, 11.};
const ScalarType TimeInfectedCritical[]   = {6.95, 6.95, 6.86, 17.36, 17.1, 11.6};

const ScalarType RecoveredPerInfectedNoSymptoms[] = {1 - 0.75, 1 - 0.75, 1 - 0.8, 1 - 0.8, 1 - 0.8, 1 - 0.8};
const ScalarType SeverePerInfectedSymptoms[]      = {0.0075, 0.0075, 0.019, 0.0615, 0.165, 0.225};
const ScalarType CriticalPerSevere[]              = {0.075, 0.075, 0.075, 0.15, 0.3, 0.4};
const ScalarType DeathsPerCritical[]              = {0.05, 0.05, 0.14, 0.14, 0.4, 0.6};

} // namespace params

/** @brief Returns transitions that can be used to initialize an IDE model or 
*    to calculate initial values for a LCT model.
*/
mio::TimeSeries<ScalarType> get_initial_flows()
{
    using namespace params;
    // The initialization vector for the LCT model is calculated by defining transitions.
    // Create TimeSeries with num_transitions elements for each group.
    using InfTransition = mio::lsecir::InfectionTransition;
    int num_transitions = (int)InfTransition::Count;
    mio::TimeSeries<ScalarType> init(num_transitions * num_groups);

    // Add time points for initialization of transitions.
    /* For this example, the intention is to create nearly constant values for SusceptiblesToExposed flow
    at the beginning of the simulation. Therefore we initialize the flows accordingly constant for
    SusceptiblesToExposed and derive matching values for the other flows.*/
    const ScalarType SusceptibleToExposed_dayinit[] = {199.471, 633.29, 3779.09, 4140.7, 1385.04, 479.161};
    Eigen::VectorXd init_transitions(num_transitions * num_groups);
    for (size_t group = 0; group < num_groups; group++) {
        init_transitions[(int)InfTransition::SusceptibleToExposed + group * num_transitions] =
            SusceptibleToExposed_dayinit[group];
        init_transitions[(int)InfTransition::ExposedToInfectedNoSymptoms + group * num_transitions] =
            SusceptibleToExposed_dayinit[group];
        init_transitions[(int)InfTransition::InfectedNoSymptomsToInfectedSymptoms + group * num_transitions] =
            SusceptibleToExposed_dayinit[group] * (1 - RecoveredPerInfectedNoSymptoms[group]);
        init_transitions[(int)InfTransition::InfectedNoSymptomsToRecovered + group * num_transitions] =
            SusceptibleToExposed_dayinit[group] * RecoveredPerInfectedNoSymptoms[group];
        init_transitions[(int)InfTransition::InfectedSymptomsToInfectedSevere + group * num_transitions] =
            init_transitions[(int)InfTransition::InfectedNoSymptomsToInfectedSymptoms + group * num_transitions] *
            SeverePerInfectedSymptoms[group];
        init_transitions[(int)InfTransition::InfectedSymptomsToRecovered + group * num_transitions] =
            init_transitions[(int)InfTransition::InfectedNoSymptomsToInfectedSymptoms + group * num_transitions] *
            (1 - SeverePerInfectedSymptoms[group]);
        init_transitions[(int)InfTransition::InfectedSevereToInfectedCritical + group * num_transitions] =
            init_transitions[(int)InfTransition::InfectedSymptomsToInfectedSevere + group * num_transitions] *
            CriticalPerSevere[group];
        init_transitions[(int)InfTransition::InfectedSevereToRecovered + group * num_transitions] =
            init_transitions[(int)InfTransition::InfectedSymptomsToInfectedSevere + group * num_transitions] *
            (1 - CriticalPerSevere[group]);
        init_transitions[(int)InfTransition::InfectedCriticalToDead + group * num_transitions] =
            init_transitions[(int)InfTransition::InfectedSevereToInfectedCritical + group * num_transitions] *
            DeathsPerCritical[group];
        init_transitions[(int)InfTransition::InfectedCriticalToRecovered + group * num_transitions] =
            init_transitions[(int)InfTransition::InfectedSevereToInfectedCritical + group * num_transitions] *
            (1 - DeathsPerCritical[group]);
    }
    init_transitions = init_transitions * dt;
    // Add initial time point to time series.
    init.add_time_point(-350, init_transitions);
    // Add further time points until time 0 with constant values.
    while (init.get_last_time() < -dt + 1e-10) {
        init.add_time_point(init.get_last_time() + dt, init_transitions);
    }
    return init;
}

/**
 * @brief Performs a simulation of a real scenario with an LCT and an ODE model.
 *
 */
template <size_t num_subcompartments>
void simulate()
{
    using namespace params;
    // ----- Initialize age resolved model. -----
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, num_subcompartments, num_subcompartments, num_subcompartments,
                                            num_subcompartments, num_subcompartments, 1, 1>;
    using Model    = mio::lsecir::Model<LctState, LctState, LctState, LctState, LctState, LctState>;
    Model model;
    std::cout << num_subcompartments << std::endl;

    // Define parameters used for simulation and initialization.
    for (size_t group = 0; group < num_groups; group++) {
        model.parameters.template get<mio::lsecir::TimeExposed>()[group]            = TimeExposed[group];
        model.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms>()[group] = TimeInfectedNoSymptoms[group];
        model.parameters.template get<mio::lsecir::TimeInfectedSymptoms>()[group]   = TimeInfectedSymptoms[group];
        model.parameters.template get<mio::lsecir::TimeInfectedSevere>()[group]     = TimeInfectedSevere[group];
        model.parameters.template get<mio::lsecir::TimeInfectedCritical>()[group]   = TimeInfectedCritical[group];
        model.parameters.template get<mio::lsecir::TransmissionProbabilityOnContact>()[group] =
            TransmissionProbabilityOnContact[group];

        model.parameters.template get<mio::lsecir::RelativeTransmissionNoSymptoms>()[group] =
            RelativeTransmissionNoSymptoms;
        model.parameters.template get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[group] =
            RiskOfInfectionFromSymptomatic;

        model.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[group] =
            RecoveredPerInfectedNoSymptoms[group];
        model.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms>()[group] =
            SeverePerInfectedSymptoms[group];
        model.parameters.template get<mio::lsecir::CriticalPerSevere>()[group] = CriticalPerSevere[group];
        model.parameters.template get<mio::lsecir::DeathsPerCritical>()[group] = DeathsPerCritical[group];
    }
    mio::ContactMatrixGroup& contact_matrix = model.parameters.template get<mio::lsecir::ContactPatterns>();
    /**3.9547 1.1002 2.9472   2.05 0.3733 0.0445
    0.3327 3.5892  1.236 1.9208 0.2681 0.0161
    0.246 0.7124 5.6518 3.2939 0.2043 0.0109
    0.1742 0.8897 3.3124 4.5406 0.4262 0.0214
    0.0458 0.1939 0.5782 1.3825  1.473 0.0704
    0.1083 0.1448 0.4728 0.9767 0.6266 0.1724 */
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, 3.5));

    model.parameters.template get<mio::lsecir::ContactPatterns>() = contact_matrix;
    model.parameters.template get<mio::lsecir::Seasonality>()     = seasonality;

    mio::Vector<ScalarType> total_confirmed_cases(num_groups);
    total_confirmed_cases << 6820., 19164., 122877., 145125., 53235., 26468.;
    total_confirmed_cases = total_confirmed_cases * 0.2;
    mio::Vector<ScalarType> deaths(num_groups);
    deaths << 0., 0., 0., 0., 0., 0.;
    Eigen::VectorXd age_group_sizes_eig(num_groups);
    age_group_sizes_eig << 3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434;
    mio::lsecir::Initializer<Model> initializer(std::move(get_initial_flows()), model);
    initializer.set_tol_for_support_max(1e-6);
    initializer.compute_initialization_vector(age_group_sizes_eig, deaths, total_confirmed_cases);

    // Perform simulation.
    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // Choose dt_min = dt_max so that we have a fixed time step and can compare to the result with one group.
    integrator->set_dt_min(dt);
    integrator->set_dt_max(dt);
    double total = 0;

    total -= omp_get_wtime();
    mio::simulate<ScalarType, Model>(0, tmax, dt, model, integrator);
    total += omp_get_wtime();

    std::cout << "Simulation took " << total << " seconds in average!" << std::endl;
    // TODO: Add initialization with vector !!
    // TODO: Add realistic contacts.
}
