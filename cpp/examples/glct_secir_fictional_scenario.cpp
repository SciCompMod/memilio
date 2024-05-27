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
#include "glct_secir/parameters.h"
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
// Valid initialization vectors only available for 3 and 10 subcompartments.
// For other numbers of subcompartments, an initialization vector is used but the result is maybe distorted.
constexpr int num_subcompartments       = 3;
constexpr int twice_num_subcompartments = 6;

// Parameters are calculated via examples/compute_parameters.cpp.
std::map<std::string, ScalarType> simulation_parameter = {{"TimeExposed", 3.335},
                                                          {"TransmissionProbabilityOnContact", 0.0733271},
                                                          {"RelativeTransmissionNoSymptoms", 1.},
                                                          {"RiskOfInfectionFromSymptomatic", 0.3},
                                                          {"Seasonality", 0.},
                                                          {"RecoveredPerInfectedNoSymptoms", 0.206901},
                                                          {"SeverePerInfectedSymptoms", 0.0786429},
                                                          {"CriticalPerSevere", 0.173176},
                                                          {"DeathsPerCritical", 0.217177}};

Eigen::VectorXd get_initialization_vector()
{
    using Model = mio::glsecir::Model<num_subcompartments, twice_num_subcompartments, twice_num_subcompartments, 2, 2>;
    using LctState = Model::LctState;
    Eigen::VectorXd init_subcomp(LctState::Count);
    if (num_subcompartments == 3) {
        init_subcomp << 82786880.87971967, //Susceptibles
            4501.15138702, 4503.17589402, 4503.17594660, // Exposed
            1995.63097927, 1997.23640301, 1997.23677271, // InfectedNoSymptomsToInfectedSymptoms
            2476.22539714, 2476.64393643, 2476.64356410, //  InfectedNoSymptomsToRecovered
            531.01119455, 531.13743225, 531.13738356, //  InfectedSymptomsToInfectedSevere
            6905.30655382, 6906.78534187, 6906.78459191, // InfectedSymptomsToRecovered
            409.51336414, //InfectedSevereToInfectedCritical
            2529.82777984, // InfectedSevereToRecovered
            150.89097051, //InfectedCriticalToDead
            613.76366884, //InfectedCriticalToRecovered
            305496.84171871, //Recovered
            9710.00000000; //  Dead
    }
    else if (num_subcompartments == 10) {
        init_subcomp << 82786880.85122260, //Susceptibles
            1348.92959502, 1350.95296417, 1350.95395534, 1350.95393656, 1350.95391844, 1350.95390254, 1350.95388719,
            1350.95387321, 1350.95385957, 1350.95384662, // Exposed
            597.56641004, 599.16989113, 599.17132094, 599.17131635, 599.17131265, 599.17130876, 599.17130487,
            599.17130154, 599.17129810, 599.17129482, //  InfectedNoSymptomsToInfectedSymptoms
            742.57588476, 742.99474927, 742.99479428, 742.99476528, 742.99473822, 742.99471383, 742.99469033,
            742.99466849, 742.99464776, 742.99462731, //  InfectedNoSymptomsToRecovered
            159.21519311, 159.34144908, 159.34147749, 159.34147315, 159.34146915, 159.34146538, 159.34146192,
            159.34145862, 159.34145549, 159.34145249, //  InfectedSymptomsToInfectedSevere
            2070.55965000, 2072.03889013, 2072.03916970, 2072.03910523, 2072.03904662, 2072.03899265, 2072.03894331,
            2072.03889426, 2072.03885019, 2072.03880558, //  InfectedSymptomsToRecovered
            409.51336414, //  InfectedSevereToInfectedCritical
            2529.82777984, //  InfectedSevereToRecovered
            150.89097051, //  InfectedCriticalToDead
            613.76366884, //  InfectedCriticalToRecovered
            305496.80551313, //  Recovered
            9710.00000000; //  Dead
    }
    else {
        Eigen::VectorXd vec_comp(8);
        vec_comp << 82787281.86349970, 13306.99948374, 13219.13701656, 22153.17301890, 2926.86690340, 762.50129842,
            305670.45877928, 9710.;
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
            init_subcomp[i] =
                simulation_parameter["RecoveredPerInfectedNoSymptoms"] * vec_comp[2] / num_subcompartments;
        }
        for (int i = LctState::get_first_index<LctState::InfectionState::InfectedSymptoms>();
             i < LctState::get_first_index<LctState::InfectionState::InfectedSymptoms>() + num_subcompartments; i++) {
            init_subcomp[i] = simulation_parameter["SeverePerInfectedSymptoms"] * vec_comp[3] / num_subcompartments;
        }
        for (int i = LctState::get_first_index<LctState::InfectionState::InfectedSymptoms>() + num_subcompartments;
             i < LctState::get_first_index<LctState::InfectionState::InfectedSevere>(); i++) {
            init_subcomp[i] =
                (1 - simulation_parameter["SeverePerInfectedSymptoms"]) * vec_comp[3] / num_subcompartments;
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
    }
    return init_subcomp;
}

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
    std::cout << "Simulation with GLCT model and " << num_subcompartments << " subcompartments." << std::endl;
    // Initialize model.
    using Model = mio::glsecir::Model<num_subcompartments, twice_num_subcompartments, twice_num_subcompartments, 2, 2>;
    using LctState = Model::LctState;
    Model model(std::move(Eigen::VectorXd::Zero(LctState::Count)));

    // Set Parameters.
    // Exposed.
    model.parameters.get<mio::glsecir::StartingProbabilitiesExposed>() =
        mio::glsecir::StartingProbabilitiesExposed().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::Exposed>());
    model.parameters.get<mio::glsecir::TransitionMatrixExposedToInfectedNoSymptoms>() =
        mio::glsecir::TransitionMatrixExposedToInfectedNoSymptoms().get_default(num_subcompartments, 3.335);

    // InfectedNoSymptoms.
    Eigen::VectorXd StartingProbabilitiesInfectedNoSymptoms = Eigen::VectorXd::Zero(twice_num_subcompartments);
    StartingProbabilitiesInfectedNoSymptoms[0] = 1 - simulation_parameter["RecoveredPerInfectedNoSymptoms"];
    StartingProbabilitiesInfectedNoSymptoms[num_subcompartments] =
        simulation_parameter["RecoveredPerInfectedNoSymptoms"];
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedNoSymptoms>() =
        StartingProbabilitiesInfectedNoSymptoms;

    model.parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>() =
        mio::glsecir::TransitionMatrixInfectedNoSymptomsToInfectedSymptoms().get_default(num_subcompartments, 1.865);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedNoSymptomsToRecovered().get_default(num_subcompartments, 8.865);

    // InfectedSymptoms.
    Eigen::VectorXd StartingProbabilitiesInfectedSymptoms      = Eigen::VectorXd::Zero(twice_num_subcompartments);
    StartingProbabilitiesInfectedSymptoms[0]                   = simulation_parameter["SeverePerInfectedSymptoms"];
    StartingProbabilitiesInfectedSymptoms[num_subcompartments] = 1 - simulation_parameter["SeverePerInfectedSymptoms"];
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedSymptoms>() = StartingProbabilitiesInfectedSymptoms;

    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToInfectedSevere>() =
        mio::glsecir::TransitionMatrixInfectedSymptomsToInfectedSevere().get_default(num_subcompartments, 6.30662);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered().get_default(num_subcompartments, 7.);

    // InfectedSevere.
    Eigen::VectorXd StartingProbabilitiesInfectedSevere = Eigen::VectorXd::Zero(2);
    StartingProbabilitiesInfectedSevere[0]              = simulation_parameter["CriticalPerSevere"];
    StartingProbabilitiesInfectedSevere[1]              = 1 - simulation_parameter["CriticalPerSevere"];
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>() = StartingProbabilitiesInfectedSevere;

    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToInfectedCritical>() =
        mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered().get_default(1, 9.36);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered().get_default(1, 12.110701);

    // InfectedCritical.
    Eigen::VectorXd StartingProbabilitiesInfectedCritical = Eigen::VectorXd::Zero(2);
    StartingProbabilitiesInfectedCritical[0]              = simulation_parameter["DeathsPerCritical"];
    StartingProbabilitiesInfectedCritical[1]              = 1 - simulation_parameter["DeathsPerCritical"];
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedCritical>() = StartingProbabilitiesInfectedCritical;

    model.parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToDead>() =
        mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered().get_default(1, 15.88);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered().get_default(1, 17.92);

    // Remaining parameters.
    model.parameters.get<mio::glsecir::ContactPatterns>() = get_contact_matrix(R0);

    model.parameters.get<mio::glsecir::TransmissionProbabilityOnContact>() =
        simulation_parameter["TransmissionProbabilityOnContact"];
    model.parameters.get<mio::glsecir::RelativeTransmissionNoSymptoms>() =
        simulation_parameter["RelativeTransmissionNoSymptoms"];
    model.parameters.get<mio::glsecir::RiskOfInfectionFromSymptomatic>() =
        simulation_parameter["RiskOfInfectionFromSymptomatic"];
    model.parameters.get<mio::glsecir::Seasonality>() = simulation_parameter["Seasonality"];

    model.set_initial_values(get_initialization_vector());

    // Perform simulation.
    mio::TimeSeries<ScalarType> result = mio::glsecir::simulate(
        0, tmax, 0.1, model,
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>(
            1e-10, 1e-5, 0.001, 0.001));
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
    ScalarType R0 = 2.;
    // Path is valid if file is executed eg in memilio/build/bin.
    // Folder has to exist beforehand.
    //std::string save_dir = "../../data/simulation_lct/dropR0short/";
    std::string save_dir = "../../data/simulation_lct/riseR02short/";

    auto result = simulate_glct_model(R0, 10, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}