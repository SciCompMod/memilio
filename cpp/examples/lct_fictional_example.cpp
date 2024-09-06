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

// Necessary because num_subcompartments is used as a template argument and has to be a constexpr.
constexpr int num_subcompartments = 1;
constexpr size_t num_groups       = 6;

// Parameters
const ScalarType dt                             = 0.01;
const ScalarType seasonality                    = 0.;
const ScalarType RelativeTransmissionNoSymptoms = 1.;
const ScalarType RiskOfInfectionFromSymptomatic = 0.3;

const ScalarType age_group_sizes[] = {3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434};

const ScalarType total_confirmed_cases              = 341223.;
const ScalarType deaths                             = 9710.;
const ScalarType TransmissionProbabilityOnContact[] = {0.03, 0.06, 0.06, 0.06, 0.09, 0.175};

const ScalarType TimeExposed[]            = {3.335, 3.335, 3.335, 3.335, 3.335, 3.335};
const ScalarType TimeInfectedNoSymptoms[] = {2.74, 2.74, 2.565, 2.565, 2.565, 2.565};
const ScalarType TimeInfectedSymptoms[]   = {7.02625, 7.02625, 7.0665, 6.9385, 6.835, 6.775};
const ScalarType TimeInfectedSevere[]     = {5, 5, 5.925, 7.55, 8.5, 11};
const ScalarType TimeInfectedCritical[]   = {6.95, 6.95, 6.86, 17.36, 17.1, 11.6};

const ScalarType RecoveredPerInfectedNoSymptoms[] = {1 - 0.75, 1 - 0.75, 1 - 0.8, 1 - 0.8, 1 - 0.8, 1 - 0.8};
const ScalarType SeverePerInfectedSymptoms[]      = {0.0075, 0.0075, 0.019, 0.0615, 0.165, 0.225};
const ScalarType CriticalPerSevere[]              = {0.075, 0.075, 0.075, 0.15, 0.3, 0.4};
const ScalarType DeathsPerCritical[]              = {0.05, 0.05, 0.14, 0.14, 0.4, 0.6};

const std::string contact_locations[] = {"home", "school_pf_eig", "work", "other"};

/** @brief Returns contact matrix in relation to defined R0.
* Contacts are defined such that R0 equals 1 at the beginning of the simulation and jumps to R0 in 
* the time interval [1.9,2.0].
*/
mio::UncertainContactMatrix<ScalarType> get_contact_matrix(ScalarType R0)
{
    std::string data_dir                   = "../../data";
    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, num_groups);
    contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_groups, num_groups, 0.));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(baseline, mio::read_mobility_plain(
                                        (data_dir / "contacts" / ("baseline_" + contact_location + ".txt")).string()));
        contact_matrix[0] += baseline;
    }
    if (R0 <= 1.) {
        // Perform simulation with dropping R0.
        //contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_groups, num_groups, 2.7463));
        contact_matrix[0].add_damping(0., mio::SimulationTime(1.9));
        contact_matrix[0].add_damping(R0, mio::SimulationTime(2.));
    }
    else {
        // Perform simulation with rising R0.
        contact_matrix[0] = R0 * contact_matrix[0];
        // contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, R0 * 2.7463));
        contact_matrix[0].add_damping(1 - 1. / R0, mio::SimulationTime(-1.));
        contact_matrix[0].add_damping(1 - 1. / R0, mio::SimulationTime(1.9));
        contact_matrix[0].add_damping(0., mio::SimulationTime(2.));
    }

    return mio::UncertainContactMatrix<ScalarType>(contact_matrix);
}

/** @brief Returns transitions that can be used to inizialize an IDE model or 
*    to calcuate initial values for a LCT model.
*/
// mio::TimeSeries<ScalarType> get_initial_flows()
// {
//     // The initialization vector for the LCT model is calculated by defining transitions.
//     // Create TimeSeries with num_transitions elements.
//     int num_transitions = (int)mio::isecir::InfectionTransition::Count;
//     mio::TimeSeries<ScalarType> init(num_transitions);

//     // Add time points for initialization of transitions.
//     /* For this example, the intention is to create nearly constant values for SusceptiblesToExposed flow
//     at the beginning of the simulation. Therefore we initalize the flows accordingly constant for
//     SusceptiblesToExposed and derive matching values for the other flows.*/
//     // 7-Tage-Inzidenz at 15.10.2020 was 34.1, see https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/Okt_2020/2020-10-15-de.pdf?__blob=publicationFile.
//     ScalarType SusceptibleToExposed_const = (34.1 / 7.) * simulation_parameter["total_population"] / 100000.;
//     Eigen::VectorXd init_transitions(num_transitions);
//     init_transitions[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]        = SusceptibleToExposed_const;
//     init_transitions[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms] = SusceptibleToExposed_const;
//     init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] =
//         SusceptibleToExposed_const * (1 - simulation_parameter["RecoveredPerInfectedNoSymptoms"]);
//     init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered] =
//         SusceptibleToExposed_const * simulation_parameter["RecoveredPerInfectedNoSymptoms"];
//     init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] =
//         init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
//         simulation_parameter["SeverePerInfectedSymptoms"];
//     init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered] =
//         init_transitions[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
//         (1 - simulation_parameter["SeverePerInfectedSymptoms"]);
//     init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] =
//         init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] *
//         simulation_parameter["CriticalPerSevere"];
//     init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered] =
//         init_transitions[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] *
//         (1 - simulation_parameter["CriticalPerSevere"]);
//     init_transitions[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead] =
//         init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] *
//         simulation_parameter["DeathsPerCritical"];
//     init_transitions[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered] =
//         init_transitions[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] *
//         (1 - simulation_parameter["DeathsPerCritical"]);
//     init_transitions = init_transitions * simulation_parameter["dt_flows"];

//     // Add initial time point to time series.
//     init.add_time_point(-350, init_transitions);
//     // Add further time points until time 0 with constant values.
//     while (init.get_last_time() < -simulation_parameter["dt_flows"] + 1e-10) {
//         init.add_time_point(init.get_last_time() + simulation_parameter["dt_flows"], init_transitions);
//     }
//     return init;
// }

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
mio::IOResult<void> simulate_lct_model(ScalarType R0, ScalarType tmax, std::string save_dir = "")
{
    mio::unused(R0);
    mio::unused(tmax);
    mio::unused(save_dir);
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

    //model.parameters.get<mio::lsecir::ContactPatterns>() = get_contact_matrix(R0);
    model.parameters.get<mio::lsecir::Seasonality>() = seasonality;

    // // Get initialization vector for LCT model with num_subcompartments subcompartments.
    // mio::lsecir::Initializer<Model> initializer(std::move(get_initial_flows()), model);
    // initializer.set_tol_for_support_max(1e-6);
    // auto status = initializer.compute_initialization_vector(simulation_parameter["total_population"],
    //                                                         simulation_parameter["deaths"],
    //                                                         simulation_parameter["total_confirmed_cases"]);
    // if (status) {
    //     return mio::failure(mio::StatusCode::InvalidValue,
    //                         "One ofthe model constraints are not fulfilled using the initialization method.");
    // }

    // // Perform simulation.
    // mio::TimeSeries<ScalarType> result = mio::lsecir::simulate(
    //     0, tmax, 0.1, model,
    //     std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>(
    //         1e-10, 1e-5, 0, 0.1));
    // // Calculate result without division in subcompartments.
    // mio::TimeSeries<ScalarType> populations = model.calculate_populations(result);

    // if (!save_dir.empty()) {
    //     //auto interpolated_result = mio::interpolate_simulation_result(populations, 0.1);
    //     //interpolated_result.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 16, 8);
    //     std::string R0string = std::to_string(R0);
    //     std::string filename = save_dir + "fictional_lct_" + R0string.substr(0, R0string.find(".") + 2) + "_" +
    //                            std::to_string(num_subcompartments);
    //     if (tmax > 50) {
    //         filename = filename + "_long";
    //     }
    //     filename                               = filename + ".h5";
    //     mio::IOResult<void> save_result_status = mio::save_result({populations}, {0}, 1, filename);
    // }

    return mio::success();
}

int main()
{
    auto result = simulate_lct_model(2., 10);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}
