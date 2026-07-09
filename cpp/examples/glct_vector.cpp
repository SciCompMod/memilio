/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Annika Jungklaus
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

#include "glct_vector/model.h"
#include "glct_vector/parameters.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/compartments/simulation.h"
#include "memilio/data/analyze_result.h"

#include <vector>

int main()
{
    // Simple example to demonstrate how to run a simulation using an GLCT VECTOR model.
    // Parameters, initial values and the number of subcompartments are not meant to represent a realistic scenario.
    // We need to double the choices of the number of subcompartments for some compartments,
    // as we define different strains for different transition possibilities. For both strains, the same number of
    // subcompartments is chosen. The transition probabilities are defined in the StartingProbabilities.
    constexpr size_t NumExposedHuman = 2, NumInfectedHuman = 4, NumExposedVector = 3;

    using Model          = mio::glvector::Model<ScalarType, NumExposedHuman, NumInfectedHuman, NumExposedVector>;
    using LctState       = Model::LctState;
    using InfectionState = LctState::InfectionState;

    Model model;

    const ScalarType tmax    = 10;
    const ScalarType t0      = 0;
    const ScalarType dt_init = 10; ///< Initially used step size for adaptive method.

    // Define some epidemiological parameters needed throughput the model definition and initialization.
    const ScalarType timeExposedHuman      = 3.2;
    const ScalarType timeInfectedHuman     = 2.;
    const ScalarType deathPerInfectedHuman = 0.0;
    const ScalarType timeExposedVector      = NumExposedVector;

    // Define the initial values with the distribution of the population into subcompartments.
    // This method of defining the initial values using a vector of vectors is not necessary, but should remind you
    // how the entries of the initial value vector relate to the defined template parameters of the model or the number
    // of subcompartments. It is also possible to define the initial values directly.
    std::vector<std::vector<ScalarType>> initial_populations = {
        {0},
        {100},
        {10,0},
        {0,0,0,0},
        {0},
        {0},
        {0},
        {0},
        {30},
        {500},
        {20,0,0},
        {0},
        {0}};


    // Assert that initial_populations has the right shape.
    if (initial_populations.size() != (size_t)InfectionState::Count) {
        mio::log_error("The number of vectors in initial_populations does not match the number of InfectionStates.");
        return 1;
    }
    if ((initial_populations[(size_t)InfectionState::UnbornHuman].size() !=
         LctState::get_num_subcompartments<InfectionState::UnbornHuman>()) ||
        (initial_populations[(size_t)InfectionState::SusceptibleHuman].size() !=
         LctState::get_num_subcompartments<InfectionState::SusceptibleHuman>()) ||
        (initial_populations[(size_t)InfectionState::ExposedHuman].size() != NumExposedHuman) ||
        (initial_populations[(size_t)InfectionState::InfectedHuman].size() != NumInfectedHuman) ||
        (initial_populations[(size_t)InfectionState::RecoveredHuman].size() !=
         LctState::get_num_subcompartments<InfectionState::RecoveredHuman>()) ||
        (initial_populations[(size_t)InfectionState::DeadHuman].size() !=
         LctState::get_num_subcompartments<InfectionState::DeadHuman>()) ||
         (initial_populations[(size_t)InfectionState::NaturalDeathHuman].size() !=
         LctState::get_num_subcompartments<InfectionState::NaturalDeathHuman>()) ||
        (initial_populations[(size_t)InfectionState::UnbornVector].size() !=
         LctState::get_num_subcompartments<InfectionState::UnbornVector>()) ||
         (initial_populations[(size_t)InfectionState::AquaticVector].size() !=
         LctState::get_num_subcompartments<InfectionState::AquaticVector>()) ||
        (initial_populations[(size_t)InfectionState::SusceptibleVector].size() !=
         LctState::get_num_subcompartments<InfectionState::SusceptibleVector>()) ||
        (initial_populations[(size_t)InfectionState::ExposedVector].size() != NumExposedVector) ||
        (initial_populations[(size_t)InfectionState::TransmitterVector].size() !=  
        LctState::get_num_subcompartments<InfectionState::TransmitterVector>()) ||
        (initial_populations[(size_t)InfectionState::DeadVector].size() !=
         LctState::get_num_subcompartments<InfectionState::DeadVector>())) {
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



    // Set Parameters for the transitions between compartments.
    // ExposedHuman.
    // The get_default of the StartingProbabilities returns the first unit vector of the defined size.
    // It is necessary to set it although the default method is used to define the length of the vector.
    model.parameters.get<mio::glvector::StartingProbabilitiesExposedHuman<ScalarType>>() =
        mio::glvector::StartingProbabilitiesExposedHuman<ScalarType>().get_default(
            LctState::get_num_subcompartments<InfectionState::ExposedHuman>());
    // The get_default function returns the TransitionMatrix that is required to have an Erlang-distributed
    // stay time with an average of timeExposed.
    model.parameters.get<mio::glvector::TransitionMatrixExposedToInfectedHuman<ScalarType>>() =
        mio::glvector::TransitionMatrixExposedToInfectedHuman<ScalarType>().get_default(
            LctState::get_num_subcompartments<InfectionState::ExposedHuman>(), timeExposedHuman);
  
    // InfectedHuman.
    // For InfectedHuman, two strains has to be defined, one for the Transition
    // InfectedHumanToRecovered and one for InfectedHumanToDead.
    // The strains have a length of NumInfectedHuman/2. each.
    // The transition probability is included in the StartingProbability vector.
    Eigen::VectorX<ScalarType> StartingProbabilitiesInfectedHuman =
        Eigen::VectorX<ScalarType>::Zero(LctState::get_num_subcompartments<InfectionState::InfectedHuman>());
    StartingProbabilitiesInfectedHuman[0] = 1 - deathPerInfectedHuman;
    StartingProbabilitiesInfectedHuman[(Eigen::Index)(
        LctState::get_num_subcompartments<InfectionState::InfectedHuman>() / 2.)] = deathPerInfectedHuman;
    model.parameters.get<mio::glvector::StartingProbabilitiesInfectedHuman<ScalarType>>() =
        StartingProbabilitiesInfectedHuman;
    // Define equal TransitionMatrices for the strains.
    model.parameters.get<mio::glvector::TransitionMatrixInfectedToRecoveredHuman<ScalarType>>() =
        mio::glvector::TransitionMatrixInfectedToRecoveredHuman<ScalarType>().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedHuman>() / 2.),
            timeInfectedHuman);
    model.parameters.get<mio::glvector::TransitionMatrixInfectedToDeadHuman<ScalarType>>() =
        mio::glvector::TransitionMatrixInfectedToDeadHuman<ScalarType>().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedHuman>() / 2.),
            timeInfectedHuman);
    
    // ExposedVector.
    // Set the StartingProbabilities.
    Eigen::VectorX<ScalarType> StartingProbabilitiesExposedVector =
        Eigen::VectorX<ScalarType>::Zero(LctState::get_num_subcompartments<InfectionState::ExposedVector>());
    StartingProbabilitiesExposedVector[0]                                         = 1;
    model.parameters.get<mio::glvector::StartingProbabilitiesExposedVector<ScalarType>>() =
        StartingProbabilitiesExposedVector;
    // Set the TransitionMatrixExposedToTransmitterVector.
    model.parameters.get<mio::glvector::TransitionMatrixExposedToTransmitterVector<ScalarType>>() =
        mio::glvector::TransitionMatrixExposedToTransmitterVector<ScalarType>().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::ExposedVector>()),
             timeExposedVector);
    
    // Setting other parameters.
    ScalarType probMosquitoInfected = 0.5;
    ScalarType probMosquitoDisseminated = 1;
    ScalarType probMosquitoTransmitter = 1;
    model.parameters.get<mio::glvector::TransmissionProbabilityOnContactHumanToVector<ScalarType>>() = 
        probMosquitoInfected * probMosquitoDisseminated * probMosquitoTransmitter;

    model.parameters.get<mio::glvector::TransmissionProbabilityOnContactVectorToHuman<ScalarType>>() = 0.2;
    model.parameters.get<mio::glvector::BitingRateVector<ScalarType>>() = 0.5;

    model.parameters.get<mio::glvector::BirthRateHuman<ScalarType>>() = 0.0000046;
    model.parameters.get<mio::glvector::MortalityRateHuman<ScalarType>>() = 0.000046;

    model.parameters.get<mio::glvector::OvipositionRateVector<ScalarType>>() = 5.;
    model.parameters.get<mio::glvector::MortalityRateAquaticVector<ScalarType>>() = 0.1;
    model.parameters.get<mio::glvector::MortalityRateAdultVector<ScalarType>>() = 0.1;
    model.parameters.get<mio::glvector::TransitionRateAquaticToAdultVector<ScalarType>>() = 0.15;
    model.parameters.get<mio::glvector::SexRatioVector<ScalarType>>() = 0.5;

    model.parameters.get<mio::glvector::SeasonalVariationImpact<ScalarType>>() = 0.;
    model.parameters.get<mio::glvector::SeasonalCycleDuration<ScalarType>>() = 20;
    model.parameters.get<mio::glvector::RatioVectorToHuman<ScalarType>>() = 5.;

    // Perform a simulation.
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(t0, tmax, dt_init, model);
    // The simulation result is divided by subcompartments.
    // We call the function calculate_compartments to get a result according to the InfectionStates.
    mio::TimeSeries<ScalarType> population_no_subcompartments = model.calculate_compartments(result);
    auto interpolated_result = mio::interpolate_simulation_result(population_no_subcompartments, 0.1);
    interpolated_result.print_table({"Uh", "Sh", "Eh", "Ih", "Rh", "Dh", "NDh", "Uv", "Av", "Sv", "Ev", "Tv", "Dv"}, 10, 4);
}
