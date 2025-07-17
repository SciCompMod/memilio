/* 
* Copyright (C) 2020-2025 MEmilio
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

#include "lct_secir_2_disease/model.h"
#include "lct_secir_2_disease/infection_state.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"
#include "memilio/compartments/simulation.h"
#include "memilio/data/analyze_result.h"

#include <vector>

int main()
{
    // Simple example to demonstrate how to run a simulation using an LCT-SECIR model.
    // One single AgeGroup/Category member is used here.
    // Parameters, initial values and the number of subcompartments are not meant to represent a realistic scenario.
    constexpr size_t NumExposed = 2, NumInfectedNoSymptoms = 3, NumInfectedSymptoms = 1, NumInfectedSevere = 1,
                     NumInfectedCritical = 5;
    using InfState                       = mio::lsecir2d::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, NumExposed, NumInfectedNoSymptoms, NumInfectedSymptoms,
                                            NumInfectedSevere, NumInfectedCritical, 1, 1>;
    using Model    = mio::lsecir2d::Model<LctState>;
    Model model;

    // Variable defines whether the class Initializer is used to define an initial vector from flows or whether a manually
    // defined initial vector is used to initialize the LCT model.

    ScalarType tmax = 10;

    // Set Parameters.
    model.parameters.get<mio::lsecir2d::TimeExposed>()[0]            = 3.2;
    model.parameters.get<mio::lsecir2d::TimeInfectedNoSymptoms>()[0] = 2.;
    model.parameters.get<mio::lsecir2d::TimeInfectedSymptoms>()[0]   = 5.8;
    model.parameters.get<mio::lsecir2d::TimeInfectedSevere>()[0]     = 9.5;
    model.parameters.get<mio::lsecir2d::TimeInfectedCritical>()[0]   = 7.1;

    model.parameters.get<mio::lsecir2d::TransmissionProbabilityOnContact>()[0] = 0.05;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::lsecir2d::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
    // From SimulationTime 5, the contact pattern is reduced to 30% of the initial value.
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(5.));

    model.parameters.get<mio::lsecir2d::RelativeTransmissionNoSymptoms>()[0] = 0.7;
    model.parameters.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic>()[0] = 0.25;
    model.parameters.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms>()[0] = 0.09;
    model.parameters.get<mio::lsecir2d::SeverePerInfectedSymptoms>()[0]      = 0.2;
    model.parameters.get<mio::lsecir2d::CriticalPerSevere>()[0]              = 0.25;
    model.parameters.get<mio::lsecir2d::DeathsPerCritical>()[0]              = 0.3;

    // Simple example how to initialize model without flows.
    // Define the initial values with the distribution of the population into subcompartments.
    // This method of defining the initial values using a vector of vectors is not necessary, but should remind you
    // how the entries of the initial value vector relate to the defined template parameters of the model or the number
    // of subcompartments. It is also possible to define the initial values directly.
    std::vector<std::vector<ScalarType>> initial_populations = {{750}, {30, 20},          {20, 10, 10}, {50},
                                                                {50},  {10, 10, 5, 3, 2}, {20},         {10}};

    // Assert that initial_populations has the right shape.
    if (initial_populations.size() != (size_t)InfState::Count) {
        mio::log_error("The number of vectors in initial_populations does not match the number of InfectionStates.");
        return 1;
    }
    if ((initial_populations[(size_t)InfState::Susceptible].size() !=
         LctState::get_num_subcompartments<InfState::Susceptible>()) ||
        (initial_populations[(size_t)InfState::Exposed].size() != NumExposed) ||
        (initial_populations[(size_t)InfState::InfectedNoSymptoms].size() != NumInfectedNoSymptoms) ||
        (initial_populations[(size_t)InfState::InfectedSymptoms].size() != NumInfectedSymptoms) ||
        (initial_populations[(size_t)InfState::InfectedSevere].size() != NumInfectedSevere) ||
        (initial_populations[(size_t)InfState::InfectedCritical].size() != NumInfectedCritical) ||
        (initial_populations[(size_t)InfState::Recovered].size() !=
         LctState::get_num_subcompartments<InfState::Recovered>()) ||
        (initial_populations[(size_t)InfState::Dead].size() != LctState::get_num_subcompartments<InfState::Dead>())) {
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
        model.populations[i] = flat_initial_populations[i];
    }

    // Perform a simulation.
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(0, tmax, 0.5, model);
    // The simulation result is divided by subcompartments.
    // We call the function calculate_compartments to get a result according to the InfectionStates.
    mio::TimeSeries<ScalarType> population_no_subcompartments = model.calculate_compartments(result);
    auto interpolated_results = mio::interpolate_simulation_result(population_no_subcompartments);
    interpolated_results.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 12, 4);
}
