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
#include "glct_secir/parameters.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"
#include "memilio/compartments/simulation.h"
#include "memilio/data/analyze_result.h"

#include <vector>

int main()
{
    // Simple example to demonstrate how to run a simulation using an GLCT SECIR model.
    // This example recreates the lct_secir.cpp example using the GLCT model.
    // This means, that we use the corresponding initial numbers, parameters and numbers of subcompartments that are
    // equivalent to the choices in the LCT example.
    // Parameters, initial values and the number of subcompartments are not meant to represent a realistic scenario.
    // We need to double the choices of the number of subcompartments for some compartments,
    // as we define different strains for different transition possibilities. For both strains, the same number of
    // subcompartments is chosen. The transition probabilities are defined in the StartingProbabilities.
    constexpr size_t NumExposed = 2, NumInfectedNoSymptoms = 6, NumInfectedSymptoms = 2, NumInfectedSevere = 2,
                     NumInfectedCritical = 10;
    using Model    = mio::glsecir::Model<NumExposed, NumInfectedNoSymptoms, NumInfectedSymptoms, NumInfectedSevere,
                                      NumInfectedCritical>;
    using LctState = Model::LctState;
    using InfectionState = LctState::InfectionState;

    Model model;

    const ScalarType tmax    = 10;
    const ScalarType t0      = 0;
    const ScalarType dt_init = 10; ///< Initially used step size for adaptive method.

    // Define some epidemiological parameters needed throughput the model definition and initialization.
    const ScalarType TimeExposed                    = 3.2;
    const ScalarType TimeInfectedNoSymptoms         = 2.;
    const ScalarType TimeInfectedSymptoms           = 5.8;
    const ScalarType TimeInfectedSevere             = 9.5;
    const ScalarType TimeInfectedCritical           = 7.1;
    const ScalarType RecoveredPerInfectedNoSymptoms = 0.09;
    const ScalarType SeverePerInfectedSymptoms      = 0.2;
    const ScalarType CriticalPerSevere              = 0.25;
    const ScalarType DeathsPerCritical              = 0.3;

    // Define the initial values with the distribution of the population into subcompartments.
    // This method of defining the initial values using a vector of vectors is not necessary, but should remind you
    // how the entries of the initial value vector relate to the defined template parameters of the model or the number
    // of subcompartments. It is also possible to define the initial values directly.
    // Split the initial population defined in the lct_secir.cpp example according to the transition probabilities in
    // the two strains per compartment.
    std::vector<std::vector<ScalarType>> initial_populations = {
        {750},
        {30, 20},
        {20 * (1 - RecoveredPerInfectedNoSymptoms), 10 * (1 - RecoveredPerInfectedNoSymptoms),
         10 * (1 - RecoveredPerInfectedNoSymptoms), 20 * RecoveredPerInfectedNoSymptoms,
         10 * RecoveredPerInfectedNoSymptoms, 10 * RecoveredPerInfectedNoSymptoms},
        {50 * SeverePerInfectedSymptoms, 50 * (1 - SeverePerInfectedSymptoms)},
        {50 * CriticalPerSevere, 50 * (1 - CriticalPerSevere)},
        {10 * DeathsPerCritical, 10 * DeathsPerCritical, 5 * DeathsPerCritical, 3 * DeathsPerCritical,
         2 * DeathsPerCritical, 10 * (1 - DeathsPerCritical), 10 * (1 - DeathsPerCritical), 5 * (1 - DeathsPerCritical),
         3 * (1 - DeathsPerCritical), 2 * (1 - DeathsPerCritical)},
        {20},
        {10}};

    // Assert that initial_populations has the right shape.
    if (initial_populations.size() != (size_t)InfectionState::Count) {
        mio::log_error("The number of vectors in initial_populations does not match the number of InfectionStates.");
        return 1;
    }
    if ((initial_populations[(size_t)InfectionState::Susceptible].size() !=
         LctState::get_num_subcompartments<InfectionState::Susceptible>()) ||
        (initial_populations[(size_t)InfectionState::Exposed].size() != NumExposed) ||
        (initial_populations[(size_t)InfectionState::InfectedNoSymptoms].size() != NumInfectedNoSymptoms) ||
        (initial_populations[(size_t)InfectionState::InfectedSymptoms].size() != NumInfectedSymptoms) ||
        (initial_populations[(size_t)InfectionState::InfectedSevere].size() != NumInfectedSevere) ||
        (initial_populations[(size_t)InfectionState::InfectedCritical].size() != NumInfectedCritical) ||
        (initial_populations[(size_t)InfectionState::Recovered].size() !=
         LctState::get_num_subcompartments<InfectionState::Recovered>()) ||
        (initial_populations[(size_t)InfectionState::Dead].size() !=
         LctState::get_num_subcompartments<InfectionState::Dead>())) {
        mio::log_error("The length of at least one vector in initial_populations does not match the related number of "
                       "subcompartments.");
        return 1;
    }

    // Transfer the initial values in initial_populations to the model.
    std::vector<ScalarType> flat_initial_populations;
    for (auto&& vec : initial_populations) {
        flat_initial_populations.insert(flat_initial_populations.end(), vec.begin(), vec.end());
    }
    for (size_t i = 0; i < LctState::Count; i++) {
        model.populations[mio::Index<LctState>(i)] = flat_initial_populations[i];
    }

    // Set Parameters according to the LCT model definitions, e.g. use Erlang-distributions.
    // Exposed.
    // The get_default of the StartingProbabilities returns the first unit vector of the defined size.
    // It is necessary to set it although the default method is used to define the length of the vector.
    model.parameters.get<mio::glsecir::StartingProbabilitiesExposed>() =
        mio::glsecir::StartingProbabilitiesExposed().get_default(
            LctState::get_num_subcompartments<InfectionState::Exposed>());
    // The get_default function returns the TransitionMatrix that is required to have an Erlang-distributed
    // stay time with an average of TimeExposed.
    model.parameters.get<mio::glsecir::TransitionMatrixExposedToInfectedNoSymptoms>() =
        mio::glsecir::TransitionMatrixExposedToInfectedNoSymptoms().get_default(
            LctState::get_num_subcompartments<InfectionState::Exposed>(), TimeExposed);
    // This definition of the StartingProbability and the TransitionMatrix lead to an Erlang-distributed latent stage.

    // InfectedNoSymptoms.
    // For InfectedNoSymptoms, two strains has to be defined, one for the Transition
    // InfectedNoSymptomsToInfectedSymptoms and one for InfectedNoSymptomsToRecovered.
    // The strains have a length of NumInfectedNoSymptoms/2. each.
    // The transition probability is included in the StartingProbability vector.
    mio::Vector<ScalarType> StartingProbabilitiesInfectedNoSymptoms =
        mio::Vector<ScalarType>::Zero(LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>());
    StartingProbabilitiesInfectedNoSymptoms[0] = 1 - RecoveredPerInfectedNoSymptoms;
    StartingProbabilitiesInfectedNoSymptoms[(Eigen::Index)(
        LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>() / 2.)] = RecoveredPerInfectedNoSymptoms;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedNoSymptoms>() =
        StartingProbabilitiesInfectedNoSymptoms;
    // Define equal TransitionMatrices for the strains.
    // They follow the same Erlang-distribution such that we get the same result as with one strain in the LCT model.
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>() =
        mio::glsecir::TransitionMatrixInfectedNoSymptomsToInfectedSymptoms().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>() / 2.),
            TimeInfectedNoSymptoms);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedNoSymptomsToRecovered().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>() / 2.),
            TimeInfectedNoSymptoms);
    // Do the same for all compartments.
    // InfectedSymptoms.
    mio::Vector<ScalarType> StartingProbabilitiesInfectedSymptoms =
        mio::Vector<ScalarType>::Zero(LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>());
    StartingProbabilitiesInfectedSymptoms[0]                                         = SeverePerInfectedSymptoms;
    StartingProbabilitiesInfectedSymptoms[(Eigen::Index)(
        LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>() / 2.)] = 1 - SeverePerInfectedSymptoms;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedSymptoms>() = StartingProbabilitiesInfectedSymptoms;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToInfectedSevere>() =
        mio::glsecir::TransitionMatrixInfectedSymptomsToInfectedSevere().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>() / 2.), TimeInfectedSymptoms);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>() / 2.), TimeInfectedSymptoms);
    // InfectedSevere.
    mio::Vector<ScalarType> StartingProbabilitiesInfectedSevere =
        mio::Vector<ScalarType>::Zero(LctState::get_num_subcompartments<InfectionState::InfectedSevere>());
    StartingProbabilitiesInfectedSevere[0]                                         = CriticalPerSevere;
    StartingProbabilitiesInfectedSevere[(Eigen::Index)(
        LctState::get_num_subcompartments<InfectionState::InfectedSevere>() / 2.)] = 1 - CriticalPerSevere;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>() = StartingProbabilitiesInfectedSevere;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToInfectedCritical>() =
        mio::glsecir::TransitionMatrixInfectedSevereToInfectedCritical().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedSevere>() / 2.), TimeInfectedSevere);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedSevereToRecovered().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedSevere>() / 2.), TimeInfectedSevere);
    // InfectedCritical.
    mio::Vector<ScalarType> StartingProbabilitiesInfectedCritical =
        mio::Vector<ScalarType>::Zero(LctState::get_num_subcompartments<InfectionState::InfectedCritical>());
    StartingProbabilitiesInfectedCritical[0]                                         = DeathsPerCritical;
    StartingProbabilitiesInfectedCritical[(Eigen::Index)(
        LctState::get_num_subcompartments<InfectionState::InfectedCritical>() / 2.)] = 1 - DeathsPerCritical;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedCritical>() = StartingProbabilitiesInfectedCritical;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToDead>() =
        mio::glsecir::TransitionMatrixInfectedCriticalToDead().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedCritical>() / 2.), TimeInfectedCritical);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedCriticalToRecovered().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedCritical>() / 2.), TimeInfectedCritical);

    model.parameters.get<mio::glsecir::TransmissionProbabilityOnContact>() = 0.05;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::glsecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
    // From SimulationTime 5, the contact pattern is reduced to 30% of the initial value.
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(5.));

    model.parameters.get<mio::glsecir::RelativeTransmissionNoSymptoms>() = 0.7;
    model.parameters.get<mio::glsecir::RiskOfInfectionFromSymptomatic>() = 0.25;

    // Perform a simulation.
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(t0, tmax, dt_init, model);
    // The simulation result is divided by subcompartments as in the LCT model.
    // We call the function calculate_compartments to get a result according to the InfectionStates.
    mio::TimeSeries<ScalarType> population_no_subcompartments = model.calculate_compartments(result);
    auto interpolated_result = mio::interpolate_simulation_result(population_no_subcompartments, 0.1);
    interpolated_result.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 10, 4);
}
