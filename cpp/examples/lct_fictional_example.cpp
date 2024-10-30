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
#include "lct_secir/initializer_flows.h"

#include "memilio/config.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"
#include "memilio/compartments/simulation.h"
#include "memilio/data/analyze_result.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"

#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <string>
#include <iostream>
#include <vector>
#include <omp.h>

namespace params
{
// Necessary because num_subcompartments is used as a template argument and has to be a constexpr.
constexpr int num_subcompartments = 3;
constexpr size_t num_groups       = 6;

// Parameters
const ScalarType dt                             = 0.01;
const ScalarType seasonality                    = 0.;
const ScalarType RelativeTransmissionNoSymptoms = 1.;
const ScalarType RiskOfInfectionFromSymptomatic = 0.3;

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

mio::TimeSeries<ScalarType> add_age_groups(mio::TimeSeries<ScalarType> ageres)
{
    using namespace params;
    size_t infstatecount = (size_t)mio::lsecir::InfectionState::Count;
    mio::TimeSeries<ScalarType> noage(infstatecount);

    for (Eigen::Index timepoint = 0; timepoint < ageres.get_num_time_points(); ++timepoint) {
        Eigen::VectorXd result = Eigen::VectorXd::Zero(infstatecount);
        for (size_t infstate = 0; infstate < infstatecount; infstate++) {
            for (size_t group = 0; group < num_groups; group++) {
                result[infstate] += ageres.get_value(timepoint)[group * infstatecount + infstate];
            }
        }
        noage.add_time_point(ageres.get_time(timepoint), result);
    }

    return noage;
}
/** @brief Returns contact matrix in relation to defined R0.
* Contacts are defined such that R0 equals 1 at the beginning of the simulation and jumps to R0 in 
* the time interval [1.9,2.0].
*/
mio::UncertainContactMatrix<ScalarType> get_contact_matrix(ScalarType R0)
{
    using namespace params; // Contacts for R=1 for every age group:
    ScalarType contacts_R1[6];
    for (size_t group = 0; group < num_groups; group++) {
        contacts_R1[group] = 1. / (num_groups * TransmissionProbabilityOnContact[group] *
                                   (TimeInfectedNoSymptoms[group] * RelativeTransmissionNoSymptoms +
                                    (1 - RecoveredPerInfectedNoSymptoms[group]) * TimeInfectedSymptoms[group] *
                                        RiskOfInfectionFromSymptomatic));
    }
    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, num_groups);
    ScalarType scale                       = 1.;
    if (R0 > 1.) {
        scale = R0;
    }
    Eigen::MatrixXd a = Eigen::MatrixXd::Zero(num_groups, num_groups);
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_groups, num_groups, scale * contacts_R1[0]));
    for (size_t group = 0; group < num_groups; group++) {
        //a.row(group) = Eigen::MatrixXd::Constant(1, num_groups, scale * contacts_R1[group]);
        a(group, group) = scale * contacts_R1[group] * num_groups;
    }
    contact_matrix[0] = mio::ContactMatrix(a);

    if (R0 <= 1.) {
        // Perform simulation with dropping R0.
        contact_matrix[0].add_damping(0., mio::SimulationTime(1.9));
        contact_matrix[0].add_damping(R0, mio::SimulationTime(2.));
    }
    else {
        // Perform simulation with rising R0.
        contact_matrix[0].add_damping(1 - 1. / R0, mio::SimulationTime(-1.));
        contact_matrix[0].add_damping(1 - 1. / R0, mio::SimulationTime(1.9));
        contact_matrix[0].add_damping(0., mio::SimulationTime(2.));
    }
    return mio::UncertainContactMatrix<ScalarType>(contact_matrix);
}

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
* @brief Perform a fictive simulation with realistic parameters and contacts, such that the reproduction number 
*   is approximately 1 at the beginning and rising or dropping at simulation time 2.
*   
*   This scenario should enable a comparison of the qualitative behavior of different LCT models.
*   
* @param[in] R0 Define R0 from simulation time 2 on. Please use a number > 0.
* @param[in] tmax End time of the simulation.
* @param[in] save_dir Specifies the directory where the results should be stored. Provide an empty string if results should not be saved.
* @returns Any io errors that happen during saving the results.
*/
mio::IOResult<void> simulate_lct_model(ScalarType R0, ScalarType tmax, std::string save_dir = "")
{
    using namespace params;
    Eigen::VectorXd age_group_sizes(num_groups);
    age_group_sizes << 3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434;
    Eigen::VectorXd total_confirmed_cases(num_groups);
    total_confirmed_cases << 6820., 19164., 122877., 145125., 53235., 26468.;
    total_confirmed_cases = total_confirmed_cases * 0.2;
    Eigen::VectorXd deaths(num_groups);
    deaths << 0., 0., 0., 0., 0., 0.;
    // deaths << 2., 1., 25., 528., 3459., 6625.;

    std::cout << "Simulation with LCT model and " << num_subcompartments << " subcompartments." << std::endl;
    // Initialize model.
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, num_subcompartments, num_subcompartments, num_subcompartments,
                                            num_subcompartments, num_subcompartments, 1, 1>;
    using Model    = mio::lsecir::Model<LctState, LctState, LctState, LctState, LctState, LctState>;
    Model model;

    // Define parameters used for simulation and initialization.
    for (size_t group = 0; group < num_groups; group++) {
        model.parameters.get<mio::lsecir::TimeExposed>()[group]            = TimeExposed[group];
        model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[group] = TimeInfectedNoSymptoms[group];
        model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[group]   = TimeInfectedSymptoms[group];
        model.parameters.get<mio::lsecir::TimeInfectedSevere>()[group]     = TimeInfectedSevere[group];
        model.parameters.get<mio::lsecir::TimeInfectedCritical>()[group]   = TimeInfectedCritical[group];
        model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[group] =
            TransmissionProbabilityOnContact[group];

        model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[group] = RelativeTransmissionNoSymptoms;
        model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[group] = RiskOfInfectionFromSymptomatic;

        model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[group] =
            RecoveredPerInfectedNoSymptoms[group];
        model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[group] = SeverePerInfectedSymptoms[group];
        model.parameters.get<mio::lsecir::CriticalPerSevere>()[group]         = CriticalPerSevere[group];
        model.parameters.get<mio::lsecir::DeathsPerCritical>()[group]         = DeathsPerCritical[group];
    }

    model.parameters.get<mio::lsecir::ContactPatterns>() = get_contact_matrix(R0);
    model.parameters.get<mio::lsecir::Seasonality>()     = seasonality;

    // Get initialization vector for LCT model with num_subcompartments subcompartments.
    mio::lsecir::Initializer<Model> initializer(std::move(get_initial_flows()), model);
    initializer.set_tol_for_support_max(1e-6);

    auto status = initializer.compute_initialization_vector(age_group_sizes, deaths, total_confirmed_cases);
    if (status) {
        return mio::failure(mio::StatusCode::InvalidValue,
                            "One of the model constraints are not fulfilled using the initialization method.");
    }

    // Perform simulation.
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(
        0, tmax, dt, model,
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>(
            1e-10, 1e-5, 0, dt));
    // Calculate result without division in subcompartments.
    mio::TimeSeries<ScalarType> populations = add_age_groups(model.calculate_compartments(result));

    if (!save_dir.empty()) {
        std::string R0string = std::to_string(R0);
        std::string filename = save_dir + "fictional_lct_" + R0string.substr(0, R0string.find(".") + 2) + "_" +
                               std::to_string(num_subcompartments);
        if (tmax > 50) {
            filename = filename + "_long";
        }
        filename                               = filename + ".h5";
        mio::IOResult<void> save_result_status = mio::save_result({populations}, {0}, 1, filename);
    }

    return mio::success();
}

int main()
{
    std::string save_dir = "../../data/simulation_lct/dropR0/";
    auto result          = simulate_lct_model(0.5, 12, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}
