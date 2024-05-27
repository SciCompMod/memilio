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
#include "memilio/utils/logging.h"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <string>
#include <map>
#include <iostream>

// Necessary because num_subcompartments is used as a template argument and has to be a constexpr.
constexpr int num_subcompartments = 10;

// Parameters are calculated via examples/compute_parameters.cpp.
std::map<std::string, ScalarType> simulation_parameter = {{"dt_flows", 0.001},
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
    while (init.get_last_time() < -simulation_parameter["dt_flows"] + 1e-10) {
        init.add_time_point(init.get_last_time() + simulation_parameter["dt_flows"], init_transitions);
    }
    return init;
}

/** 
* @brief 
*   
* @returns Iniutal value vecotr.
*/
void calculate_inital_values()
{
    // Initialize model.
    using Model    = mio::lsecir::Model<num_subcompartments, num_subcompartments, num_subcompartments, 1, 1>;
    using LctState = Model::LctState;
    mio::lsecir::Parameters parameters_lct;

    // Define parameters used for simulation and initialization.
    parameters_lct.get<mio::lsecir::TimeExposed>() = simulation_parameter["TimeExposed"];
    parameters_lct.get<mio::lsecir::TimeInfectedNoSymptoms>() =
        simulation_parameter["TimeInfectedNoSymptomsToInfectedSymptoms"];
    parameters_lct.get<mio::lsecir::TimeInfectedSymptoms>() =
        simulation_parameter["TimeInfectedSymptomsToInfectedSevere"];
    parameters_lct.get<mio::lsecir::TimeInfectedSevere>() =
        simulation_parameter["TimeInfectedSevereToInfectedCritical"];
    parameters_lct.get<mio::lsecir::TimeInfectedCritical>() = simulation_parameter["TimeInfectedCriticalToDead"];
    parameters_lct.get<mio::lsecir::TransmissionProbabilityOnContact>() =
        simulation_parameter["TransmissionProbabilityOnContact"];

    mio::ContactMatrixGroup& contact_matrix = parameters_lct.get<mio::lsecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 2.7463));

    parameters_lct.get<mio::lsecir::RelativeTransmissionNoSymptoms>() =
        simulation_parameter["RelativeTransmissionNoSymptoms"];
    parameters_lct.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() =
        simulation_parameter["RiskOfInfectionFromSymptomatic"];
    parameters_lct.get<mio::lsecir::Seasonality>() = simulation_parameter["Seasonality"];
    parameters_lct.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() =
        simulation_parameter["RecoveredPerInfectedNoSymptoms"];
    parameters_lct.get<mio::lsecir::SeverePerInfectedSymptoms>() = simulation_parameter["SeverePerInfectedSymptoms"];
    parameters_lct.get<mio::lsecir::CriticalPerSevere>()         = simulation_parameter["CriticalPerSevere"];
    parameters_lct.get<mio::lsecir::DeathsPerCritical>()         = simulation_parameter["DeathsPerCritical"];

    Model model(std::move(Eigen::VectorXd::Zero(LctState::Count)), std::move(parameters_lct));

    // Get initialization vector for LCT model with num_subcompartments subcompartments.
    mio::lsecir::Initializer<Model> initializer(std::move(get_initial_flows()), model);
    initializer.set_tol_for_support_max(1e-6);
    auto status = initializer.compute_initialization_vector(simulation_parameter["total_population"],
                                                            simulation_parameter["deaths"],
                                                            simulation_parameter["total_confirmed_cases"]);
    if (status) {
        mio::log_error("One of the model constraints are not fulfilled using the initialization method.");
    }
    auto result_not_Recovered = model.get_initial_values();

    // --- Other initial values. ---
    // As initializer only uses the incoming flows, we can calculate remaining results at one time.
    parameters_lct.get<mio::lsecir::TimeInfectedNoSymptoms>() =
        simulation_parameter["TimeInfectedNoSymptomsToRecovered"];
    parameters_lct.get<mio::lsecir::TimeInfectedSymptoms>() = simulation_parameter["TimeInfectedSymptomsToRecovered"];
    parameters_lct.get<mio::lsecir::TimeInfectedSevere>()   = simulation_parameter["TimeInfectedSevereToRecovered"];
    parameters_lct.get<mio::lsecir::TimeInfectedCritical>() = simulation_parameter["TimeInfectedCriticalToRecovered"];
    Model modelToRecovered(std::move(Eigen::VectorXd::Zero(LctState::Count)), std::move(parameters_lct));

    // Get initialization vector.
    mio::lsecir::Initializer<Model> initializerToRecovered(std::move(get_initial_flows()), modelToRecovered);
    initializerToRecovered.set_tol_for_support_max(1e-6);
    status = initializerToRecovered.compute_initialization_vector(simulation_parameter["total_population"],
                                                                  simulation_parameter["deaths"],
                                                                  simulation_parameter["total_confirmed_cases"]);
    auto result_ToRecovered = modelToRecovered.get_initial_values();

    if (status) {
        mio::log_error("One of the model constraints are not fulfilled using the initialization method.");
    }
    // All of the components are calculated, now print result.
    ScalarType recovered    = simulation_parameter["total_confirmed_cases"] - simulation_parameter["deaths"];
    ScalarType Susceptibles = simulation_parameter["total_population"] - simulation_parameter["deaths"];
    // Exposed.
    std::cout << "  Exposed" << std::endl;
    for (int i = LctState::get_first_index<LctState::InfectionState::Exposed>();
         i < LctState::get_first_index<LctState::InfectionState::InfectedNoSymptoms>(); i++) {
        std::cout << std::fixed << std::setprecision(8) << result_not_Recovered[i] << ", ";
    }
    Susceptibles -= result_not_Recovered
                        .segment(LctState::get_first_index<LctState::InfectionState::Exposed>(), num_subcompartments)
                        .sum();
    // InfectedNoSymptoms.
    std::cout << "\n  InfectedNoSymptomsToInfectedSymptoms" << std::endl;
    for (int i = LctState::get_first_index<LctState::InfectionState::InfectedNoSymptoms>();
         i < LctState::get_first_index<LctState::InfectionState::InfectedSymptoms>(); i++) {
        std::cout << std::fixed << std::setprecision(8)
                  << (1 - simulation_parameter["RecoveredPerInfectedNoSymptoms"]) * result_not_Recovered[i] << ", ";
    }
    Susceptibles -=
        (1 - simulation_parameter["RecoveredPerInfectedNoSymptoms"]) *
        result_not_Recovered
            .segment(LctState::get_first_index<LctState::InfectionState::InfectedNoSymptoms>(), num_subcompartments)
            .sum();
    std::cout << "\n  InfectedNoSymptomsToRecovered" << std::endl;
    for (int i = LctState::get_first_index<LctState::InfectionState::InfectedNoSymptoms>();
         i < LctState::get_first_index<LctState::InfectionState::InfectedSymptoms>(); i++) {
        std::cout << std::fixed << std::setprecision(8)
                  << simulation_parameter["RecoveredPerInfectedNoSymptoms"] * result_ToRecovered[i] << ", ";
    }
    Susceptibles -=
        simulation_parameter["RecoveredPerInfectedNoSymptoms"] *
        result_ToRecovered
            .segment(LctState::get_first_index<LctState::InfectionState::InfectedNoSymptoms>(), num_subcompartments)
            .sum();
    // InfectedSymptoms.
    std::cout << "\n  InfectedSymptomsToInfectedSevere" << std::endl;
    for (int i = LctState::get_first_index<LctState::InfectionState::InfectedSymptoms>();
         i < LctState::get_first_index<LctState::InfectionState::InfectedSevere>(); i++) {
        std::cout << std::fixed << std::setprecision(8)
                  << simulation_parameter["SeverePerInfectedSymptoms"] * result_not_Recovered[i] << ", ";
    }
    ScalarType dummy =
        simulation_parameter["SeverePerInfectedSymptoms"] *
        result_not_Recovered
            .segment(LctState::get_first_index<LctState::InfectionState::InfectedSymptoms>(), num_subcompartments)
            .sum();
    Susceptibles -= dummy;
    recovered -= dummy;
    std::cout << "\n  InfectedSymptomsToRecovered" << std::endl;
    for (int i = LctState::get_first_index<LctState::InfectionState::InfectedSymptoms>();
         i < LctState::get_first_index<LctState::InfectionState::InfectedSevere>(); i++) {
        std::cout << std::fixed << std::setprecision(8)
                  << (1 - simulation_parameter["SeverePerInfectedSymptoms"]) * result_ToRecovered[i] << ", ";
    }
    dummy = (1 - simulation_parameter["SeverePerInfectedSymptoms"]) *
            result_ToRecovered
                .segment(LctState::get_first_index<LctState::InfectionState::InfectedSymptoms>(), num_subcompartments)
                .sum();
    Susceptibles -= dummy;
    recovered -= dummy;
    // InfectedSevere.
    std::cout << "\n  InfectedSevereToInfectedCritical" << std::endl;
    dummy = simulation_parameter["CriticalPerSevere"] *
            result_not_Recovered[LctState::get_first_index<LctState::InfectionState::InfectedSevere>()];
    std::cout << std::fixed << std::setprecision(8) << dummy << ", ";
    Susceptibles -= dummy;
    recovered -= dummy;
    std::cout << "\n  InfectedSevereToRecovered" << std::endl;
    dummy = (1 - simulation_parameter["CriticalPerSevere"]) *
            result_ToRecovered[LctState::get_first_index<LctState::InfectionState::InfectedSevere>()];
    std::cout << std::fixed << std::setprecision(8) << dummy << ", ";
    Susceptibles -= dummy;
    recovered -= dummy;
    // InfectedCritical.
    std::cout << "\n  InfectedCriticalToDead" << std::endl;
    dummy = simulation_parameter["DeathsPerCritical"] *
            result_not_Recovered[LctState::get_first_index<LctState::InfectionState::InfectedCritical>()];
    std::cout << std::fixed << std::setprecision(8) << dummy << ", ";
    Susceptibles -= dummy;
    recovered -= dummy;
    std::cout << "\n  InfectedCriticalToRecovered" << std::endl;
    dummy = (1 - simulation_parameter["DeathsPerCritical"]) *
            result_ToRecovered[LctState::get_first_index<LctState::InfectionState::InfectedCritical>()];

    std::cout << std::fixed << std::setprecision(8) << dummy << ", ";
    Susceptibles -= dummy;
    recovered -= dummy;
    Susceptibles -= recovered;
    // Recovered, Dead, Susceptibles.
    std::cout << "\n  Recovered " << std::fixed << std::setprecision(8) << recovered << ", " << std::endl;
    std::cout << "  Dead " << std::fixed << std::setprecision(8) << simulation_parameter["deaths"] << "; " << std::endl;
    std::cout << "  Susceptibles " << std::fixed << std::setprecision(8) << Susceptibles << ", " << std::endl;
}

int main()
{
    calculate_inital_values();
}