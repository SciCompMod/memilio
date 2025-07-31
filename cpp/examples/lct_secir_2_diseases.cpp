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

#include "lct_secir_2_diseases/model.h"
#include "lct_secir_2_diseases/infection_state.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
//#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/lct2d_infection_state.h"
//#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"
#include "memilio/compartments/simulation.h"
#include "memilio/data/analyze_result.h"

#include <vector>

int main()
{
    // Simple example to demonstrate how to run a simulation using an LCT-SECIR model.
    // One single AgeGroup/Category member is used here.
    // Parameters, initial values and the number of subcompartments are not meant to represent a realistic scenario.
    constexpr size_t NumExposed_1a = 1, NumInfectedNoSymptoms_1a = 1, NumInfectedSymptoms_1a = 1,
                     NumInfectedSevere_1a = 1, NumInfectedCritical_1a = 1, NumExposed_2a = 1,
                     NumInfectedNoSymptoms_2a = 1, NumInfectedSymptoms_2a = 1, NumInfectedSevere_2a = 1,
                     NumInfectedCritical_2a = 1, NumExposed_1b = 1, NumInfectedNoSymptoms_1b = 1,
                     NumInfectedSymptoms_1b = 1, NumInfectedSevere_1b = 1, NumInfectedCritical_1b = 1,
                     NumExposed_2b = 1, NumInfectedNoSymptoms_2b = 1, NumInfectedSymptoms_2b = 1,
                     NumInfectedSevere_2b = 1, NumInfectedCritical_2b = 1;
    using InfState = mio::lsecir2d::InfectionState;
    using LctState = mio::LctInfectionState<
        InfState, 3, NumExposed_1a, NumInfectedNoSymptoms_1a, NumInfectedSymptoms_1a, NumInfectedSevere_1a,
        NumInfectedCritical_1a, 1, 1, NumExposed_2a, NumInfectedNoSymptoms_2a, NumInfectedSymptoms_2a,
        NumInfectedSevere_2a, NumInfectedCritical_2a, NumExposed_1b, NumInfectedNoSymptoms_1b, NumInfectedSymptoms_1b,
        NumInfectedSevere_1b, NumInfectedCritical_1b, 1, 1, NumExposed_2b, NumInfectedNoSymptoms_2b,
        NumInfectedSymptoms_2b, NumInfectedSevere_2b, NumInfectedCritical_2b, 1>;
    using Model = mio::lsecir2d::Model<LctState>;
    Model model;

    // Variable defines whether the class Initializer is used to define an initial vector from flows or whether a manually
    // defined initial vector is used to initialize the LCT model.

    ScalarType tmax = 5;

    // Set Parameters.
    model.parameters.get<mio::lsecir2d::TimeExposed_a>()[0]            = 3.;
    model.parameters.get<mio::lsecir2d::TimeExposed_b>()[0]            = 3.;
    model.parameters.get<mio::lsecir2d::TimeInfectedNoSymptoms_a>()[0] = 3.;
    model.parameters.get<mio::lsecir2d::TimeInfectedNoSymptoms_b>()[0] = 3.;
    model.parameters.get<mio::lsecir2d::TimeInfectedSymptoms_a>()[0]   = 3.;
    model.parameters.get<mio::lsecir2d::TimeInfectedSymptoms_b>()[0]   = 3.;
    model.parameters.get<mio::lsecir2d::TimeInfectedSevere_a>()[0]     = 3.;
    model.parameters.get<mio::lsecir2d::TimeInfectedSevere_b>()[0]     = 3.;
    model.parameters.get<mio::lsecir2d::TimeInfectedCritical_a>()[0]   = 3.;
    model.parameters.get<mio::lsecir2d::TimeInfectedCritical_b>()[0]   = 3.;

    model.parameters.get<mio::lsecir2d::TransmissionProbabilityOnContact_a>()[0] = 0.1;
    model.parameters.get<mio::lsecir2d::TransmissionProbabilityOnContact_b>()[0] = 0.1;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::lsecir2d::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
    // From SimulationTime 5, the contact pattern is reduced to 30% of the initial value.
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(5.));

    model.parameters.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_a>()[0] = 0.7;
    model.parameters.get<mio::lsecir2d::RelativeTransmissionNoSymptoms_b>()[0] = 0.7;
    model.parameters.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_a>()[0] = 0.25;
    model.parameters.get<mio::lsecir2d::RiskOfInfectionFromSymptomatic_b>()[0] = 0.25;
    model.parameters.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_a>()[0] = 0.09;
    model.parameters.get<mio::lsecir2d::RecoveredPerInfectedNoSymptoms_b>()[0] = 0.09;
    model.parameters.get<mio::lsecir2d::SeverePerInfectedSymptoms_a>()[0]      = 0.2;
    model.parameters.get<mio::lsecir2d::SeverePerInfectedSymptoms_b>()[0]      = 0.2;
    model.parameters.get<mio::lsecir2d::CriticalPerSevere_a>()[0]              = 0.25;
    model.parameters.get<mio::lsecir2d::CriticalPerSevere_b>()[0]              = 0.25;
    model.parameters.get<mio::lsecir2d::DeathsPerCritical_a>()[0]              = 0.3;
    model.parameters.get<mio::lsecir2d::DeathsPerCritical_b>()[0]              = 0.3;

    // Simple example how to initialize model without flows.
    // Define the initial values with the distribution of the population into subcompartments.
    // This method of defining the initial values using a vector of vectors is not necessary, but should remind you
    // how the entries of the initial value vector relate to the defined template parameters of the model or the number
    // of subcompartments. It is also possible to define the initial values directly.
    std::vector<std::vector<ScalarType>> initial_populations = {{0, 0, 0}, {0}, {0},    {0}, {0}, {0}, {9000}, {0}, {0},
                                                                {0},       {0}, {0},    {0}, {0}, {0}, {0},    {0}, {0},
                                                                {0},       {0}, {1000}, {0}, {0}, {0}, {0},    {0}};

    // Assert that initial_populations has the right shape.
    if (initial_populations.size() != (size_t)InfState::Count) {
        mio::log_error("The number of vectors in initial_populations does not match the number of InfectionStates.");
        return 1;
    }
    if ((initial_populations[(size_t)InfState::Susceptible].size() !=
         LctState::get_num_subcompartments<InfState::Susceptible>()) ||
        (initial_populations[(size_t)InfState::Exposed_1a].size() != NumExposed_1a) ||
        (initial_populations[(size_t)InfState::InfectedNoSymptoms_1a].size() != NumInfectedNoSymptoms_1a) ||
        (initial_populations[(size_t)InfState::InfectedSymptoms_1a].size() != NumInfectedSymptoms_1a) ||
        (initial_populations[(size_t)InfState::InfectedSevere_1a].size() != NumInfectedSevere_1a) ||
        (initial_populations[(size_t)InfState::InfectedCritical_1a].size() != NumInfectedCritical_1a) ||
        (initial_populations[(size_t)InfState::Exposed_2a].size() != NumExposed_2a) ||
        (initial_populations[(size_t)InfState::InfectedNoSymptoms_2a].size() != NumInfectedNoSymptoms_2a) ||
        (initial_populations[(size_t)InfState::InfectedSymptoms_2a].size() != NumInfectedSymptoms_2a) ||
        (initial_populations[(size_t)InfState::InfectedSevere_2a].size() != NumInfectedSevere_2a) ||
        (initial_populations[(size_t)InfState::InfectedCritical_2a].size() != NumInfectedCritical_2a) ||
        (initial_populations[(size_t)InfState::Recovered_a].size() !=
         LctState::get_num_subcompartments<InfState::Recovered_a>()) ||
        (initial_populations[(size_t)InfState::Dead_a].size() !=
         LctState::get_num_subcompartments<InfState::Dead_a>()) ||
        (initial_populations[(size_t)InfState::Exposed_1b].size() != NumExposed_1b) ||
        (initial_populations[(size_t)InfState::InfectedNoSymptoms_1b].size() != NumInfectedNoSymptoms_1b) ||
        (initial_populations[(size_t)InfState::InfectedSymptoms_1b].size() != NumInfectedSymptoms_1b) ||
        (initial_populations[(size_t)InfState::InfectedSevere_1b].size() != NumInfectedSevere_1b) ||
        (initial_populations[(size_t)InfState::InfectedCritical_1b].size() != NumInfectedCritical_1b) ||
        (initial_populations[(size_t)InfState::Exposed_2b].size() != NumExposed_2b) ||
        (initial_populations[(size_t)InfState::InfectedNoSymptoms_2b].size() != NumInfectedNoSymptoms_2b) ||
        (initial_populations[(size_t)InfState::InfectedSymptoms_2b].size() != NumInfectedSymptoms_2b) ||
        (initial_populations[(size_t)InfState::InfectedSevere_2b].size() != NumInfectedSevere_2b) ||
        (initial_populations[(size_t)InfState::InfectedCritical_2b].size() != NumInfectedCritical_2b) ||
        (initial_populations[(size_t)InfState::Recovered_ab].size() !=
         LctState::get_num_subcompartments<InfState::Recovered_ab>())) {
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
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(0, tmax, 0.1, model);
    // The simulation result is divided by subcompartments.
    // We call the function calculate_compartments to get a result according to the InfectionStates.
    mio::TimeSeries<ScalarType> population_no_subcompartments = model.calculate_compartments(result);
    auto interpolated_results = mio::interpolate_simulation_result(population_no_subcompartments);

    interpolated_results.print_table({"   S",   "   E1a", "   C1a", "   I1a", "   H1a", "   U1a", "   Ra",
                                      "   Da",  "   E2a", "   C2a", "   I2a", "   H2a", "   U2a", "   E1b",
                                      "   C1b", "   I1b", "   H1b", "   U1b", "   Rb",  "   Db",  "   E2b",
                                      "   C2b", "   I2b", "   H2b", "   U2b", "   Rab"},
                                     7, 3);
}
