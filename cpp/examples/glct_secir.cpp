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

#include <vector>

int main()
{
    // Simple example to demonstrate how to run a simulation using an GLCT SECIR model.
    // Parameters, initial values and the number of subcompartments are not meant to represent a realistic scenario.

    using Model    = mio::glsecir::Model<2, 6, 2, 2, 10>;
    using LctState = Model::LctState;

    ScalarType tmax = 10;

    // Define the initial value vector init with the distribution of the population into subcompartments.
    // This method of defining the vector using a vector of vectors is a bit of overhead, but should remind you how
    // the entries of the initial value vector relate to the defined template parameters of the model or the number of subcompartments.
    // It is also possible to define the initial value vector directly.
    std::vector<std::vector<ScalarType>> initial_populations = {
        {750},
        {30, 20},
        {20 * 0.91, 10 * 0.91, 10 * 0.91, 20 * 0.09, 10 * 0.09, 10 * 0.09},
        {50 * 0.2, 50 * 0.8},
        {50 * 0.25, 50 * 0.75},
        {10 * 0.3, 10 * 0.3, 5 * 0.3, 3 * 0.3, 2 * 0.3, 10 * 0.7, 10 * 0.7, 5 * 0.7, 3 * 0.7, 2 * 0.7},
        {20},
        {10}};

    // Assert that initial_populations has the right shape.
    if (initial_populations.size() != (size_t)LctState::InfectionState::Count) {
        mio::log_error("The number of vectors in initial_populations does not match the number of InfectionStates.");
        return 1;
    }
    if ((initial_populations[(int)LctState::InfectionState::Susceptible].size() !=
         (size_t)LctState::get_num_subcompartments<LctState::InfectionState::Susceptible>()) ||
        (initial_populations[(int)LctState::InfectionState::Exposed].size() !=
         (size_t)LctState::get_num_subcompartments<LctState::InfectionState::Exposed>()) ||
        (initial_populations[(int)LctState::InfectionState::InfectedNoSymptoms].size() !=
         (size_t)LctState::get_num_subcompartments<LctState::InfectionState::InfectedNoSymptoms>()) ||
        (initial_populations[(int)LctState::InfectionState::InfectedSymptoms].size() !=
         (size_t)LctState::get_num_subcompartments<LctState::InfectionState::InfectedSymptoms>()) ||
        (initial_populations[(int)LctState::InfectionState::InfectedSevere].size() !=
         (size_t)LctState::get_num_subcompartments<LctState::InfectionState::InfectedSevere>()) ||
        (initial_populations[(int)LctState::InfectionState::InfectedCritical].size() !=
         (size_t)LctState::get_num_subcompartments<LctState::InfectionState::InfectedCritical>()) ||
        (initial_populations[(int)LctState::InfectionState::Recovered].size() !=
         (size_t)LctState::get_num_subcompartments<LctState::InfectionState::Recovered>()) ||
        (initial_populations[(int)LctState::InfectionState::Dead].size() !=
         (size_t)LctState::get_num_subcompartments<LctState::InfectionState::Dead>())) {
        mio::log_error("The length of at least one vector in initial_populations does not match the related number of "
                       "subcompartments.");
        return 1;
    }

    // Transfer the initial values in initial_populations to the vector init.
    Eigen::VectorXd init = Eigen::VectorXd::Zero(LctState::Count);
    init[LctState::get_first_index<LctState::InfectionState::Susceptible>()] =
        initial_populations[(int)LctState::InfectionState::Susceptible][0];
    for (int i = 0; i < LctState::get_num_subcompartments<LctState::InfectionState::Exposed>(); i++) {
        init[LctState::get_first_index<LctState::InfectionState::Exposed>() + i] =
            initial_populations[(int)LctState::InfectionState::Exposed][i];
    }
    for (int i = 0; i < LctState::get_num_subcompartments<LctState::InfectionState::InfectedNoSymptoms>(); i++) {
        init[LctState::get_first_index<LctState::InfectionState::InfectedNoSymptoms>() + i] =
            initial_populations[(int)LctState::InfectionState::InfectedNoSymptoms][i];
    }
    for (int i = 0; i < LctState::get_num_subcompartments<LctState::InfectionState::InfectedSymptoms>(); i++) {
        init[LctState::get_first_index<LctState::InfectionState::InfectedSymptoms>() + i] =
            initial_populations[(int)LctState::InfectionState::InfectedSymptoms][i];
    }
    for (int i = 0; i < LctState::get_num_subcompartments<LctState::InfectionState::InfectedSevere>(); i++) {
        init[LctState::get_first_index<LctState::InfectionState::InfectedSevere>() + i] =
            initial_populations[(int)LctState::InfectionState::InfectedSevere][i];
    }
    for (int i = 0; i < LctState::get_num_subcompartments<LctState::InfectionState::InfectedCritical>(); i++) {
        init[LctState::get_first_index<LctState::InfectionState::InfectedCritical>() + i] =
            initial_populations[(int)LctState::InfectionState::InfectedCritical][i];
    }
    init[LctState::get_first_index<LctState::InfectionState::Recovered>()] =
        initial_populations[(int)LctState::InfectionState::Recovered][0];
    init[LctState::get_first_index<LctState::InfectionState::Dead>()] =
        initial_populations[(int)LctState::InfectionState::Dead][0];

    // Initialize model.
    Model model(std::move(init));

    // Set Parameters.
    // Exposed.
    model.parameters.get<mio::glsecir::StartingProbabilitiesExposed>() =
        mio::glsecir::StartingProbabilitiesExposed().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::Exposed>());
    model.parameters.get<mio::glsecir::TransitionMatrixExposedToInfectedNoSymptoms>() =
        mio::glsecir::TransitionMatrixExposedToInfectedNoSymptoms().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::Exposed>(), 3.2);
    // InfectedNoSymptoms.
    Eigen::VectorXd StartingProbabilitiesInfectedNoSymptoms =
        Eigen::VectorXd::Zero(LctState::get_num_subcompartments<LctState::InfectionState::InfectedNoSymptoms>());
    StartingProbabilitiesInfectedNoSymptoms[0] = 1 - 0.09;
    StartingProbabilitiesInfectedNoSymptoms
        [LctState::get_num_subcompartments<LctState::InfectionState::InfectedNoSymptoms>() / 2.] = 0.09;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedNoSymptoms>() =
        StartingProbabilitiesInfectedNoSymptoms;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>() =
        mio::glsecir::TransitionMatrixInfectedNoSymptomsToInfectedSymptoms().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::InfectedNoSymptoms>() / 2., 2.);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedNoSymptomsToRecovered().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::InfectedNoSymptoms>() / 2., 2.);
    // InfectedSymptoms.
    Eigen::VectorXd StartingProbabilitiesInfectedSymptoms =
        Eigen::VectorXd::Zero(LctState::get_num_subcompartments<LctState::InfectionState::InfectedSymptoms>());
    StartingProbabilitiesInfectedSymptoms[0] = 0.2;
    StartingProbabilitiesInfectedSymptoms
        [LctState::get_num_subcompartments<LctState::InfectionState::InfectedSymptoms>() / 2.] = 1 - 0.2;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedSymptoms>() = StartingProbabilitiesInfectedSymptoms;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToInfectedSevere>() =
        mio::glsecir::TransitionMatrixInfectedSymptomsToInfectedSevere().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::InfectedSymptoms>() / 2., 5.8);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::InfectedSymptoms>() / 2., 5.8);
    // InfectedSevere.
    Eigen::VectorXd StartingProbabilitiesInfectedSevere =
        Eigen::VectorXd::Zero(LctState::get_num_subcompartments<LctState::InfectionState::InfectedSevere>());
    StartingProbabilitiesInfectedSevere[0]                                    = 0.25;
    StartingProbabilitiesInfectedSevere[LctState::get_num_subcompartments<LctState::InfectionState::InfectedSevere>() /
                                        2.]                                   = 1 - 0.25;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>() = StartingProbabilitiesInfectedSevere;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToInfectedCritical>() =
        mio::glsecir::TransitionMatrixInfectedSevereToInfectedCritical().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::InfectedSevere>() / 2., 9.5);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedSevereToRecovered().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::InfectedSevere>() / 2., 9.5);
    // InfectedCritical.
    Eigen::VectorXd StartingProbabilitiesInfectedCritical =
        Eigen::VectorXd::Zero(LctState::get_num_subcompartments<LctState::InfectionState::InfectedCritical>());
    StartingProbabilitiesInfectedCritical[0] = 0.3;
    StartingProbabilitiesInfectedCritical
        [LctState::get_num_subcompartments<LctState::InfectionState::InfectedCritical>() / 2.] = 1 - 0.3;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedCritical>() = StartingProbabilitiesInfectedCritical;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToDead>() =
        mio::glsecir::TransitionMatrixInfectedCriticalToDead().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::InfectedCritical>() / 2., 7.1);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedCriticalToRecovered().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::InfectedCritical>() / 2., 7.1);

    model.parameters.get<mio::glsecir::TransmissionProbabilityOnContact>() = 0.05;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::glsecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
    // From SimulationTime 5, the contact pattern is reduced to 30% of the initial value.
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(5.));

    model.parameters.get<mio::glsecir::RelativeTransmissionNoSymptoms>() = 0.7;
    model.parameters.get<mio::glsecir::RiskOfInfectionFromSymptomatic>() = 0.25;

    // Perform a simulation.
    mio::TimeSeries<ScalarType> result = mio::glsecir::simulate(0, tmax, 0.5, model);
    // Calculate the distribution in the InfectionState%s without subcompartments of the result and print it.
    mio::TimeSeries<ScalarType> population_no_subcompartments = model.calculate_populations(result);
    auto interpolated_result = mio::interpolate_simulation_result(population_no_subcompartments, 0.1);
    interpolated_result.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 16, 8);
}