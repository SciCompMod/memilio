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

#include "glct_secir/model.h"
#include "glct_secir/infection_state.h"
#include "glct_secir/simulation.h"

#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"
#include "memilio/data/analyze_result.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"

#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <string>
#include <map>

// Necessary because num_subcompartments is used as a template argument and has to be a constexpr.
constexpr int num_subcompartments       = 3;
constexpr int twice_num_subcompartments = 6;
Eigen::VectorXd init_vector()
{
    Eigen::VectorXd vec(8);
    vec << 82787281.86349970, 13306.99948374, 13219.13701656, 22153.17301890, 2926.86690340, 762.50129842,
        305670.45877928, 9710.;

    return vec;
}

// Parameters are calculated via examples/compute_parameters.cpp.
std::map<std::string, ScalarType> simulation_parameter = {{"TimeExposed", 3.335},
                                                          {"TransmissionProbabilityOnContact", 0.0733271},
                                                          {"RelativeTransmissionNoSymptoms", 1},
                                                          {"RiskOfInfectionFromSymptomatic", 0.3},
                                                          {"Seasonality", 0.},
                                                          {"RecoveredPerInfectedNoSymptoms", 0.206901},
                                                          {"SeverePerInfectedSymptoms", 0.0786429},
                                                          {"CriticalPerSevere", 0.173176},
                                                          {"DeathsPerCritical", 0.217177}};

/** @brief Returns contact matrix in relation to defined R0.
* Contacts are defined such that R0 equals 1 at the beginning of the simulation and jumps to R0 in 
* the time interval [1.9,2.0].
*/
mio::UncertainContactMatrix<ScalarType> get_contact_matrix(ScalarType R0)
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

    return mio::UncertainContactMatrix<ScalarType>(contact_matrix);
}

/** 
* @brief Perform a fictive simulation with realistic parameters and contacts, such that the reproduction number 
*   is approximately 1 at the beginning and rising or dropping at simulationtime 2.
*   
*   This scenario should enable a comparison of the qualitative behavior of different LCT models.
*   
* @param[in] R0 Define R0 from simulationtime 2 on. Please use a number > 0.
* @param[in] tmax End time of the simulation.
* @param[in] save_dir Specifies the directory where the results should be stored. Provide an empty string if results should not be saved.
* @returns Any io errors that happen during saving the results.
*/
mio::IOResult<void> simulate_glct_model(ScalarType R0, ScalarType tmax, std::string save_dir = "")
{
    // Initialize model.
    using Model = mio::glsecir::Model<num_subcompartments, twice_num_subcompartments, twice_num_subcompartments, 2, 2>;
    using LctState = Model::LctState;
    Model model(std::move(Eigen::VectorXd::Zero(LctState::Count)));

    // Set Parameters.
    // Exposed.
    model.parameters.get<mio::glsecir::StartingProbabilitiesExposed>() =
        mio::glsecir::StartingProbabilitiesExposed().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::Exposed>());
    model.parameters.get<mio::glsecir::TransitionMatrixExposed>() = mio::glsecir::TransitionMatrixExposed().get_default(
        LctState::get_num_subcompartments<LctState::InfectionState::Exposed>(), 3.335);

    // InfectedNoSymptoms.
    Eigen::VectorXd StartingProbabilitiesInfectedNoSymptoms = Eigen::VectorXd::Zero(twice_num_subcompartments);
    StartingProbabilitiesInfectedNoSymptoms[0] = 1 - simulation_parameter["RecoveredPerInfectedNoSymptoms"];
    StartingProbabilitiesInfectedNoSymptoms[num_subcompartments] =
        simulation_parameter["RecoveredPerInfectedNoSymptoms"];
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedNoSymptoms>() =
        StartingProbabilitiesInfectedNoSymptoms;

    Eigen::MatrixXd TransitionMatrixInfectedNoSymptoms =
        Eigen::VectorXd::Constant(twice_num_subcompartments, -num_subcompartments / 1.865).asDiagonal();
    TransitionMatrixInfectedNoSymptoms.diagonal(1).setConstant(num_subcompartments / 1.865);

    for (int i = num_subcompartments; i < twice_num_subcompartments - 1; i++) {
        TransitionMatrixInfectedNoSymptoms(i, i)     = -num_subcompartments / 8.865;
        TransitionMatrixInfectedNoSymptoms(i, i + 1) = num_subcompartments / 8.865;
    }
    TransitionMatrixInfectedNoSymptoms(num_subcompartments - 1, num_subcompartments) = 0;
    TransitionMatrixInfectedNoSymptoms(twice_num_subcompartments - 1, twice_num_subcompartments - 1) =
        -num_subcompartments / 8.865;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptoms>() = TransitionMatrixInfectedNoSymptoms;

    // InfectedSymptoms.
    Eigen::VectorXd StartingProbabilitiesInfectedSymptoms      = Eigen::VectorXd::Zero(twice_num_subcompartments);
    StartingProbabilitiesInfectedSymptoms[0]                   = simulation_parameter["SeverePerInfectedSymptoms"];
    StartingProbabilitiesInfectedSymptoms[num_subcompartments] = 1 - simulation_parameter["SeverePerInfectedSymptoms"];
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedSymptoms>() = StartingProbabilitiesInfectedSymptoms;

    Eigen::MatrixXd TransitionMatrixInfectedSymptoms =
        Eigen::VectorXd::Constant(twice_num_subcompartments, -num_subcompartments / 6.30662).asDiagonal();
    TransitionMatrixInfectedSymptoms.diagonal(1).setConstant(num_subcompartments / 6.30662);

    for (int i = num_subcompartments; i < twice_num_subcompartments - 1; i++) {
        TransitionMatrixInfectedSymptoms(i, i)     = -num_subcompartments / 7.;
        TransitionMatrixInfectedSymptoms(i, i + 1) = num_subcompartments / 7.;
    }
    TransitionMatrixInfectedSymptoms(num_subcompartments - 1, num_subcompartments) = 0;
    TransitionMatrixInfectedSymptoms(twice_num_subcompartments - 1, twice_num_subcompartments - 1) =
        -num_subcompartments / 7.;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSymptoms>() = TransitionMatrixInfectedSymptoms;

    // InfectedSevere.
    Eigen::VectorXd StartingProbabilitiesInfectedSevere = Eigen::VectorXd::Zero(2);
    StartingProbabilitiesInfectedSevere[0]              = simulation_parameter["CriticalPerSevere"];
    StartingProbabilitiesInfectedSevere[1]              = 1 - simulation_parameter["CriticalPerSevere"];
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>() = StartingProbabilitiesInfectedSevere;

    Eigen::MatrixXd TransitionMatrixInfectedSevere                       = Eigen::MatrixXd::Zero(2, 2);
    TransitionMatrixInfectedSevere(0, 0)                                 = 1. / 9.36;
    TransitionMatrixInfectedSevere(1, 1)                                 = 1. / 12.110701;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSevere>() = TransitionMatrixInfectedSevere;

    // InfectedCritical.
    Eigen::VectorXd StartingProbabilitiesInfectedCritical = Eigen::VectorXd::Zero(2);
    StartingProbabilitiesInfectedCritical[0]              = simulation_parameter["DeathsPerCritical"];
    StartingProbabilitiesInfectedCritical[1]              = 1 - simulation_parameter["DeathsPerCritical"];
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedCritical>() = StartingProbabilitiesInfectedCritical;

    Eigen::MatrixXd TransitionMatrixInfectedCritical                       = Eigen::MatrixXd::Zero(2, 2);
    TransitionMatrixInfectedCritical(0, 0)                                 = 1. / 15.88;
    TransitionMatrixInfectedCritical(1, 1)                                 = 1. / 17.92;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedCritical>() = TransitionMatrixInfectedCritical;

    // Remaining parameters.
    model.parameters.get<mio::glsecir::ContactPatterns>() = get_contact_matrix(R0);

    model.parameters.get<mio::glsecir::TransmissionProbabilityOnContact>() =
        simulation_parameter["TransmissionProbabilityOnContact"];
    model.parameters.get<mio::glsecir::RelativeTransmissionNoSymptoms>() =
        simulation_parameter["RelativeTransmissionNoSymptoms"];
    model.parameters.get<mio::glsecir::RiskOfInfectionFromSymptomatic>() =
        simulation_parameter["RiskOfInfectionFromSymptomatic"];
    model.parameters.get<mio::glsecir::Seasonality>() = simulation_parameter["Seasonality"];
    model.parameters.get<mio::glsecir::RecoveredPerInfectedNoSymptoms>() =
        simulation_parameter["RecoveredPerInfectedNoSymptoms"];
    model.parameters.get<mio::glsecir::SeverePerInfectedSymptoms>() = simulation_parameter["SeverePerInfectedSymptoms"];
    model.parameters.get<mio::glsecir::CriticalPerSevere>()         = simulation_parameter["CriticalPerSevere"];
    model.parameters.get<mio::glsecir::DeathsPerCritical>()         = simulation_parameter["DeathsPerCritical"];

    // Set initialization vector for GLCT model.
    Eigen::VectorXd vec_comp = init_vector();
    Eigen::VectorXd init_subcomp(LctState::Count);
    init_subcomp[0] = vec_comp[0];
    for (int i = LctState::get_first_index<LctState::InfectionState::Exposed>();
         i < LctState::get_first_index<LctState::InfectionState::InfectedNoSymptoms>(); i++) {
        init_subcomp[i] = vec_comp[1] / LctState::get_num_subcompartments<LctState::InfectionState::Exposed>();
    }
    for (int i = LctState::get_first_index<LctState::InfectionState::InfectedNoSymptoms>();
         i < LctState::get_first_index<LctState::InfectionState::InfectedNoSymptoms>() + num_subcompartments; i++) {
        init_subcomp[i] =
            (1 - simulation_parameter["RecoveredPerInfectedNoSymptoms"]) * vec_comp[2] / num_subcompartments;
    }
    for (int i = LctState::get_first_index<LctState::InfectionState::InfectedNoSymptoms>() + num_subcompartments;
         i < LctState::get_first_index<LctState::InfectionState::InfectedSymptoms>(); i++) {
        init_subcomp[i] = simulation_parameter["RecoveredPerInfectedNoSymptoms"] * vec_comp[2] / num_subcompartments;
    }
    for (int i = LctState::get_first_index<LctState::InfectionState::InfectedSymptoms>();
         i < LctState::get_first_index<LctState::InfectionState::InfectedSymptoms>() + num_subcompartments; i++) {
        init_subcomp[i] = simulation_parameter["SeverePerInfectedSymptoms"] * vec_comp[3] / num_subcompartments;
    }
    for (int i = LctState::get_first_index<LctState::InfectionState::InfectedSymptoms>() + num_subcompartments;
         i < LctState::get_first_index<LctState::InfectionState::InfectedSevere>(); i++) {
        init_subcomp[i] = (1 - simulation_parameter["SeverePerInfectedSymptoms"]) * vec_comp[3] / num_subcompartments;
    }
    init_subcomp[LctState::get_first_index<LctState::InfectionState::InfectedSevere>()] =
        simulation_parameter["CriticalPerSevere"] * vec_comp[4];
    init_subcomp[LctState::get_first_index<LctState::InfectionState::InfectedSevere>() + 1] =
        (1 - simulation_parameter["CriticalPerSevere"]) * vec_comp[4];
    init_subcomp[LctState::get_first_index<LctState::InfectionState::InfectedCritical>()] =
        simulation_parameter["DeathsPerCritical"] * vec_comp[5];
    init_subcomp[LctState::get_first_index<LctState::InfectionState::InfectedCritical>() + 1] =
        (1 - simulation_parameter["DeathsPerCritical"]) * vec_comp[5];
    init_subcomp[LctState::get_first_index<LctState::InfectionState::Recovered>()] = vec_comp[6];
    init_subcomp[LctState::get_first_index<LctState::InfectionState::Dead>()]      = vec_comp[7];

    model.set_initial_values(init_subcomp);

    // Perform simulation.
    mio::TimeSeries<ScalarType> result = mio::glsecir::simulate(
        0, tmax, 0.1, model,
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>(
            1e-10, 1e-5, 1e-10, 0.1));
    // Calculate result without division in subcompartments.
    mio::TimeSeries<ScalarType> populations = model.calculate_populations(result);

    if (!save_dir.empty()) {
        std::string R0string = std::to_string(R0);
        std::string filename = save_dir + "fictional_glct_" + R0string.substr(0, R0string.find(".") + 2) + "_" +
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
    ScalarType R0 = 0.5;
    // Path is valid if file is executed eg in memilio/build/bin.
    // Folder has to exist beforehand.
    std::string save_dir = "../../data/simulation_lct/dropR0short/";
    //std::string save_dir = "../../data/simulation_lct/riseR02long/";

    auto result = simulate_glct_model(R0, 10, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}