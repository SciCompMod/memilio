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

#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/math/eigen.h"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <string>
#include <map>
#include <iostream>

// Necessary because num_subcompartments is used as a template argument and has to be a constexpr.
constexpr int num_subcompartments = 3;

// Parameters are calculated via examples/compute_parameters.cpp.
std::map<std::string, ScalarType> simulation_parameter = {{"dt_flows", 0.1},
                                                          {"total_population", 83155031.},
                                                          {"total_confirmed_cases", 341223.},
                                                          {"deaths", 9710.},
                                                          {"TimeExposed", 3.335},
                                                          {"TimeInfectedNoSymptomsToInfectedSymptoms", 1.865},
                                                          {"TimeInfectedNoSymptomsToRecovered", 8.865},
                                                          {"TimeInfectedSymptomsToInfectedSevere", 6.30663},
                                                          {"TimeInfectedSymptomsToRecovered", 7.},
                                                          {"TimeInfectedSevereToInfectedCritical", 9.36},
                                                          {"TimeInfectedSevereToRecovered", 12.110701},
                                                          {"TimeInfectedCriticalToDead", 15.88},
                                                          {"TimeInfectedCriticalToRecovered", 17.92},
                                                          {"TransmissionProbabilityOnContact", 0.0733271},
                                                          {"RelativeTransmissionNoSymptoms", 1},
                                                          {"RiskOfInfectionFromSymptomatic", 0.3},
                                                          {"Seasonality", 0.},
                                                          {"RecoveredPerInfectedNoSymptoms", 0.206901},
                                                          {"SeverePerInfectedSymptoms", 0.0786429},
                                                          {"CriticalPerSevere", 0.173176},
                                                          {"DeathsPerCritical", 0.217177}};

/** @brief Returns transitions that can be used to inizialize an IDE model or 
*    to calcuate initial values for a LCT model.
*/
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
    init_transitions[(int)mio::lsecir::InfectionTransition::SusceptibleToExposed]        = SusceptibleToExposed_const;
    init_transitions[(int)mio::lsecir::InfectionTransition::ExposedToInfectedNoSymptoms] = SusceptibleToExposed_const;
    init_transitions[(int)mio::lsecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] =
        SusceptibleToExposed_const * (1 - simulation_parameter["RecoveredPerInfectedNoSymptoms"]);
    init_transitions[(int)mio::lsecir::InfectionTransition::InfectedNoSymptomsToRecovered] =
        SusceptibleToExposed_const * simulation_parameter["RecoveredPerInfectedNoSymptoms"];
    init_transitions[(int)mio::lsecir::InfectionTransition::InfectedSymptomsToInfectedSevere] =
        init_transitions[(int)mio::lsecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
        simulation_parameter["SeverePerInfectedSymptoms"];
    init_transitions[(int)mio::lsecir::InfectionTransition::InfectedSymptomsToRecovered] =
        init_transitions[(int)mio::lsecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
        (1 - simulation_parameter["SeverePerInfectedSymptoms"]);
    init_transitions[(int)mio::lsecir::InfectionTransition::InfectedSevereToInfectedCritical] =
        init_transitions[(int)mio::lsecir::InfectionTransition::InfectedSymptomsToInfectedSevere] *
        simulation_parameter["CriticalPerSevere"];
    init_transitions[(int)mio::lsecir::InfectionTransition::InfectedSevereToRecovered] =
        init_transitions[(int)mio::lsecir::InfectionTransition::InfectedSymptomsToInfectedSevere] *
        (1 - simulation_parameter["CriticalPerSevere"]);
    init_transitions[(int)mio::lsecir::InfectionTransition::InfectedCriticalToDead] =
        init_transitions[(int)mio::lsecir::InfectionTransition::InfectedSevereToInfectedCritical] *
        simulation_parameter["DeathsPerCritical"];
    init_transitions[(int)mio::lsecir::InfectionTransition::InfectedCriticalToRecovered] =
        init_transitions[(int)mio::lsecir::InfectionTransition::InfectedSevereToInfectedCritical] *
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

/** 
* @brief 
*   
* @param[in] R0 Define R0 from simulationtime 2 on. Please use a number > 0.
* @returns Iniutal value vecotr.
*/
Eigen::VectorXd calculate_inital_values(ScalarType R0)
{
    // Initialize model.
    using Model = mio::lsecir::Model<num_subcompartments, num_subcompartments, num_subcompartments, num_subcompartments,
                                     num_subcompartments>;
    using LctState = Model::LctState;
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
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>(
            1e-10, 1e-5, 0, 0.1));
    // Calculate result without division in subcompartments.
    mio::TimeSeries<ScalarType> populations = model.calculate_populations(result);

    if (!save_dir.empty()) {
        auto interpolated_result = mio::interpolate_simulation_result(populations, 0.1);
        interpolated_result.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 16, 8);
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

    /*auto result = simulate_lct_model(R0, 10, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }*/

    auto result = simulate_ide_model(R0, 10, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }
    /*
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