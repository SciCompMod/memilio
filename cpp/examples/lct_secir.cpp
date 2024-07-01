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
    // Simple example to demonstrate how to run a simulation using an LCT SECIR model.
    // Parameters, initial values and the number of subcompartments are not meant to represent a realistic scenario.
    constexpr size_t NumExposed = 2, NumInfectedNoSymptoms = 3, NumInfectedSymptoms = 1, NumInfectedSevere = 1,
                     NumInfectedCritical = 5;
    using Model          = mio::lsecir::Model<NumExposed, NumInfectedNoSymptoms, NumInfectedSymptoms, NumInfectedSevere,
                                     NumInfectedCritical>;
    using LctState       = Model::LctState;
    using InfectionState = LctState::InfectionState;

    Model model;

    ScalarType tmax = 20;

    // Define the initial values with the distribution of the population into subcompartments.
    // This method of defining the initial values using a vector of vectors is not necessary, but should remind you
    // how the entries of the initial value vector relate to the defined template parameters of the model or the number
    // of subcompartments. It is also possible to define the initial values directly.
    std::vector<std::vector<ScalarType>> initial_populations = {{750}, {30, 20},          {20, 10, 10}, {50},
                                                                {50},  {10, 10, 5, 3, 2}, {20},         {10}};

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

    // Set Parameters.
    model.parameters.get<mio::lsecir::TimeExposed>()            = 3.2;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>() = 2.;
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()   = 5.8;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()     = 9.5;
    // It is also possible to change values with the set function.
    model.parameters.set<mio::lsecir::TimeInfectedCritical>(7.1);

    model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.05;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::lsecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
    // From SimulationTime 5, the contact pattern is reduced to 30% of the initial value.
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(5.));

    model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 0.7;
    model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 0.25;
    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.09;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.2;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.25;
    model.parameters.set<mio::lsecir::DeathsPerCritical>(0.3);

    // Perform a simulation.
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(0, tmax, 0.5, model);
    // The simulation result is divided by subcompartments.
    // We call the function calculate_comparttments to get a result according to the InfectionStates.
    mio::TimeSeries<ScalarType> population_no_subcompartments = model.calculate_compartments(result);
    auto interpolated_results = mio::interpolate_simulation_result(population_no_subcompartments);
    interpolated_results.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 16, 8);
}