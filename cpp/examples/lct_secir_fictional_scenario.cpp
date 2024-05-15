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
#include "lct_secir/parameters.h"
#include "lct_secir/simulation.h"

#include "ide_secir/model.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/parameters.h"
#include "ide_secir/simulation.h"

#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"
#include "memilio/utils/time_series.h"
#include "memilio/math/eigen.h"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <string>
#include <map>
#include <iostream>

// necessary because num_subcompartments is used as a template argument and has ti be
constexpr int num_subcompartments = 20;

// Parameters are calculated via examples/compute_parameters.cpp.
std::map<std::string, ScalarType> simulation_parameter = {{"dt_flows", 0.1},
                                                          {"total_population", 83155031.},
                                                          {"total_confirmed_cases", 341223.},
                                                          {"deaths", 9710.},
                                                          {"TimeExposed", 3.335},
                                                          {"TimeInfectedNoSymptoms", 3.31331},
                                                          {"TimeInfectedSymptoms", 6.94547},
                                                          {"TimeInfectedSevere", 11.634346},
                                                          {"TimeInfectedCritical", 17.476959},
                                                          {"TransmissionProbabilityOnContact", 0.0733271},
                                                          {"RelativeTransmissionNoSymptoms", 1},
                                                          {"RiskOfInfectionFromSymptomatic", 0.3},
                                                          {"Seasonality", 0.},
                                                          {"RecoveredPerInfectedNoSymptoms", 0.206901},
                                                          {"SeverePerInfectedSymptoms", 0.0786429},
                                                          {"CriticalPerSevere", 0.173176},
                                                          {"DeathsPerCritical", 0.217177}};

mio::UncertainContactMatrix get_contact_matrix(ScalarType R0)
{
    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    if (R0 <= 1.) {
        // Perform simulation with dropping R0.
        contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 2.7463));
        contact_matrix[0].add_damping(0., mio::SimulationTime(1.9));
        contact_matrix[0].add_damping(R0, mio::SimulationTime(2.));
    }
    else {
        // Perform simulation with rising R0.
        contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, R0 * 2.7463));
        contact_matrix[0].add_damping(1 - 1. / R0, mio::SimulationTime(-1.));
        contact_matrix[0].add_damping(1 - 1. / R0, mio::SimulationTime(1.9));
        contact_matrix[0].add_damping(0., mio::SimulationTime(2.));
    }

    return mio::UncertainContactMatrix(contact_matrix);
}

mio::TimeSeries<ScalarType> get_initial_flows()
{
    // The initialization vector for the LCT model is calculated by defining transitions.
    // Create TimeSeries with num_transitions elements.
    int num_transitions = (int)mio::lsecir::InfectionTransition::Count;
    mio::TimeSeries<ScalarType> init(num_transitions);

    // Add time points for initialization of transitions.
    /* For this example, the intention is to create nearly constant values for SusceptiblesToExposed flow 
    at the beginning of the simulation. Therefore we initalize the flows accordingly constant for 
    SusceptiblesToExposed and derive matching values for the other flows.*/
    // 7-Tage-Inzidenz at 15.10.2020 was 34.1, see https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/Okt_2020/2020-10-15-de.pdf?__blob=publicationFile.
    ScalarType SusceptibleToExposed_const = (34.1 / 7.) * simulation_parameter["total_population"] / 100000.;
    Eigen::VectorXd init_transitions(num_transitions);
    init_transitions[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]        = SusceptibleToExposed_const;
    init_transitions[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms] = SusceptibleToExposed_const;
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] =
        SusceptibleToExposed_const * (1 - simulation_parameter["RecoveredPerInfectedNoSymptoms"]);
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered] =
        SusceptibleToExposed_const * simulation_parameter["RecoveredPerInfectedNoSymptoms"];
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
        simulation_parameter["SeverePerInfectedSymptoms"];
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
        (1 - simulation_parameter["SeverePerInfectedSymptoms"]);
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] *
        simulation_parameter["CriticalPerSevere"];
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] *
        (1 - simulation_parameter["CriticalPerSevere"]);
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] *
        simulation_parameter["DeathsPerCritical"];
    init_transitions[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered] =
        init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] *
        (1 - simulation_parameter["DeathsPerCritical"]);
    init_transitions = init_transitions * simulation_parameter["dt_flows"];

    // Add initial time point to time series.
    init.add_time_point(-350, init_transitions);
    // Add further time points until time 0 with constant values.
    while (init.get_last_time() < -1e-10) {
        init.add_time_point(init.get_last_time() + simulation_parameter["dt_flows"], init_transitions);
    }
    return init;
}

mio::IOResult<void> simulate_ide_model(ScalarType R0, ScalarType tmax, std::string save_dir = "")
{
    // Initialize model.
    mio::isecir::Model model_ide(std::move(get_initial_flows()), simulation_parameter["total_population"],
                                 simulation_parameter["deaths"], simulation_parameter["total_confirmed_cases"]);

    // Set working parameters.
    // Set TransitionDistributions.
    mio::ConstantFunction initialfunc(0);
    mio::StateAgeFunctionWrapper delaydistributioninit(initialfunc);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib((int)mio::isecir::InfectionTransition::Count,
                                                               delaydistributioninit);

    mio::ExponentialSurvivalFunction expInfectedSevereToInfectedCritical(1. / 9.36);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical].set_state_age_function(
        expInfectedSevereToInfectedCritical);

    mio::LognormSurvivalFunction lognInfectedSevereToRecovered(0.76, -0.45, 9.41);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered].set_state_age_function(
        lognInfectedSevereToRecovered);

    mio::ExponentialSurvivalFunction expInfectedCriticalToDeath(1. / 14.88, 1);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead].set_state_age_function(
        expInfectedCriticalToDeath);
    expInfectedCriticalToDeath.set_distribution_parameter(1 / 16.92);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered].set_state_age_function(
        expInfectedCriticalToDeath);

    mio::GammaSurvivalFunction erlangExposedToInfectedNoSymptoms(num_subcompartments, 0, 3.335 / num_subcompartments);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms].set_state_age_function(
        erlangExposedToInfectedNoSymptoms);

    mio::GammaSurvivalFunction erlangInfectedNoSymptomsToInfectedSymptoms(num_subcompartments, 0,
                                                                          1.865 / num_subcompartments);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
        .set_state_age_function(erlangInfectedNoSymptomsToInfectedSymptoms);
    erlangInfectedNoSymptomsToInfectedSymptoms.set_scale(8.865 / num_subcompartments);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered].set_state_age_function(
        erlangInfectedNoSymptomsToInfectedSymptoms);

    mio::GammaSurvivalFunction erlangInfectedSymptomsToInfectedSevere(num_subcompartments, 0,
                                                                      6.30662 / num_subcompartments);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere].set_state_age_function(
        erlangInfectedSymptomsToInfectedSevere);
    erlangInfectedSymptomsToInfectedSevere.set_scale(7. / num_subcompartments);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered].set_state_age_function(
        erlangInfectedSymptomsToInfectedSevere);

    model_ide.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    // Set other parameters.
    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 1.);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
        1 - simulation_parameter["RecoveredPerInfectedNoSymptoms"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] =
        simulation_parameter["RecoveredPerInfectedNoSymptoms"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)] =
        simulation_parameter["SeverePerInfectedSymptoms"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)] =
        1 - simulation_parameter["SeverePerInfectedSymptoms"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)] =
        simulation_parameter["CriticalPerSevere"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)] =
        1 - simulation_parameter["CriticalPerSevere"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)] =
        simulation_parameter["DeathsPerCritical"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)] =
        1 - simulation_parameter["DeathsPerCritical"];

    model_ide.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    model_ide.parameters.get<mio::isecir::ContactPatterns>() = get_contact_matrix(R0);

    mio::ConstantFunction constfunc(simulation_parameter["TransmissionProbabilityOnContact"]);
    mio::StateAgeFunctionWrapper StateAgeFunctionWrapperide(constfunc);
    model_ide.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(StateAgeFunctionWrapperide);
    StateAgeFunctionWrapperide.set_distribution_parameter(simulation_parameter["RelativeTransmissionNoSymptoms"]);
    model_ide.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(StateAgeFunctionWrapperide);
    StateAgeFunctionWrapperide.set_distribution_parameter(simulation_parameter["RiskOfInfectionFromSymptomatic"]);
    model_ide.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(StateAgeFunctionWrapperide);

    model_ide.set_tol_for_support_max(1e-6);
    model_ide.check_constraints(simulation_parameter["dt_flows"]);

    // Simulate.
    mio::isecir::Simulation sim(model_ide, simulation_parameter["dt_flows"]);
    sim.advance(tmax);

    if (!save_dir.empty()) {
        std::string R0string     = std::to_string(R0);
        std::string filename_ide = save_dir + "fictional_ide_" + R0string.substr(0, R0string.find(".") + 2) + "_" +
                                   std::to_string(num_subcompartments);
        if (tmax > 50) {
            filename_ide = filename_ide + "_long";
        }
        filename_ide                           = filename_ide + ".h5";
        mio::IOResult<void> save_result_status = mio::save_result({sim.get_result()}, {0}, 1, filename_ide);
    }
    return mio::success();
}

/** 
* @brief Perform a fictive simulation with realistic parameters and contacts, such that the reproduction number 
*   is approximately 1 at the beginning and rising or dropping at simulationtime 2.
*   
*   This scenario should enable a comparison of the qualitative behavior of different models.
*   
* @param[in] R0 Define R0 from simulationtime 2 on. Please use a number > 0.
* @param[in] num_subcompartments Number of subcompartments for each compartment where subcompartments make sense. 
*        Number is also used for some distributions of the IDE model.
*        Set num_subcompartments = 0 to use subcompartments with an expected sojourn time of approximately 1.
* @param[in] save_dir Specifies the directory where the results should be stored. Provide an empty string if results should not be saved.
* @param[in] tmax End time of the simulation.
* @returns Any io errors that happen during saving the results.
*/
mio::IOResult<void> simulate_lct_model(ScalarType R0, ScalarType tmax, std::string save_dir = "")
{
    using Model = mio::lsecir::Model<num_subcompartments, num_subcompartments, num_subcompartments, num_subcompartments,
                                     num_subcompartments>;
    using LctState = Model::LctState;
    // Initialize model and perform simulation.
    Model model(std::move(Eigen::VectorXd::Zero(LctState::Count)));

    // Define parameters used for simulation and initialization.

    model.parameters.get<mio::lsecir::TimeExposed>()            = simulation_parameter["TimeExposed"];
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>() = simulation_parameter["TimeInfectedNoSymptoms"];
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()   = simulation_parameter["TimeInfectedSymptoms"];
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()     = simulation_parameter["TimeInfectedSevere"];
    model.parameters.get<mio::lsecir::TimeInfectedCritical>()   = simulation_parameter["TimeInfectedCritical"];
    model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>() =
        simulation_parameter["TransmissionProbabilityOnContact"];

    model.parameters.get<mio::lsecir::ContactPatterns>() = get_contact_matrix(R0);

    model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>() =
        simulation_parameter["RelativeTransmissionNoSymptoms"];
    model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() =
        simulation_parameter["RiskOfInfectionFromSymptomatic"];
    model.parameters.get<mio::lsecir::Seasonality>() = simulation_parameter["Seasonality"];
    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() =
        simulation_parameter["RecoveredPerInfectedNoSymptoms"];
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>() = simulation_parameter["SeverePerInfectedSymptoms"];
    model.parameters.get<mio::lsecir::CriticalPerSevere>()         = simulation_parameter["CriticalPerSevere"];
    model.parameters.get<mio::lsecir::DeathsPerCritical>()         = simulation_parameter["DeathsPerCritical"];

    // Get initialization vector for LCT model with num_subcompartments subcompartments.
    mio::lsecir::Initializer<Model> initializer(std::move(get_initial_flows()), model);
    initializer.set_tol_for_support_max(1e-6);
    auto status = initializer.compute_initialization_vector(simulation_parameter["total_population"],
                                                            simulation_parameter["deaths"],
                                                            simulation_parameter["total_confirmed_cases"]);
    if (status) {
        return mio::failure(mio::StatusCode::InvalidValue,
                            "One ofthe model constraints are not fulfilled using the initialization method.");
    }

    // Perform simulation.
    mio::TimeSeries<ScalarType> result = mio::lsecir::simulate(
        0, tmax, 0.1, model,
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>(1e-10, 1e-5, 0,
                                                                                                         0.1));
    // Calculate result without division in subcompartments.
    mio::TimeSeries<ScalarType> populations = model.calculate_populations(result);

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
    // Options used: For R0=2 epidemic peak use tmax=150,
    // for R0=4 epidemic peak use tmax = 75.
    // For short things: 10 days and R0=0.5 or 2
    ScalarType R0 = 0.5;
    // Paths are valid if file is executed eg in memilio/build/bin.
    // Folders have to exist beforehand.
    std::string save_dir = "../../data/simulation_lct/dropR0short/";

    auto result = simulate_lct_model(R0, 10, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }

    /*result = simulate_ide_model(R0, 10, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }

    R0       = 4.;
    save_dir = "../../data/simulation_lct/riseR04long/";

    result = simulate_lct_model(R0, 75,save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }

    result = simulate_ide_model(R0, 75,save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    } */

    return 0;
}