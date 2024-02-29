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
#include "lct_secir/simulation.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"

#include <vector>

/** Simple example to demonstrate how to run a simulation using an LCT SECIR model. 
    Parameters, initial values and the number of subcompartments are not meant to represent a realistic scenario. */
int main()
{
    using Model    = mio::lsecir::Model<2, 3, 1, 1, 5>;
    using LctState = Model::LctState;

    ScalarType tmax = 20;

    // Define the initial value vector init with the distribution of the population into subcompartments.
    // This method of defining the vector using a vector of vectors is a bit of overhead, but should remind you how
    // the entries of the initial value vector relate to the defined template parameters of the model or the number of subcompartments.
    // It is also possible to define the initial value vector directly.
    std::vector<std::vector<ScalarType>> initial_populations = {{750}, {30, 20},          {20, 10, 10}, {50},
                                                                {50},  {10, 10, 5, 3, 2}, {20},         {10}};

    // Assert that initial_populations has the right shape.
    if (initial_populations.size() != (int)LctState::InfectionStateBase::Count) {
        mio::log_error("The number of vectors in initial_populations does not match the number of InfectionStates.");
        return 1;
    }
    if ((initial_populations[(int)LctState::InfectionStateBase::Susceptible].size() !=
         LctState::get_num_subcompartments<LctState::InfectionStateBase::Susceptible>()) ||
        (initial_populations[(int)LctState::InfectionStateBase::Exposed].size() !=
         LctState::get_num_subcompartments<LctState::InfectionStateBase::Exposed>()) ||
        (initial_populations[(int)LctState::InfectionStateBase::InfectedNoSymptoms].size() !=
         LctState::get_num_subcompartments<LctState::InfectionStateBase::InfectedNoSymptoms>()) ||
        (initial_populations[(int)LctState::InfectionStateBase::InfectedSymptoms].size() !=
         LctState::get_num_subcompartments<LctState::InfectionStateBase::InfectedSymptoms>()) ||
        (initial_populations[(int)LctState::InfectionStateBase::InfectedSevere].size() !=
         LctState::get_num_subcompartments<LctState::InfectionStateBase::InfectedSevere>()) ||
        (initial_populations[(int)LctState::InfectionStateBase::InfectedCritical].size() !=
         LctState::get_num_subcompartments<LctState::InfectionStateBase::InfectedCritical>()) ||
        (initial_populations[(int)LctState::InfectionStateBase::Recovered].size() !=
         LctState::get_num_subcompartments<LctState::InfectionStateBase::Recovered>()) ||
        (initial_populations[(int)LctState::InfectionStateBase::Dead].size() !=
         LctState::get_num_subcompartments<LctState::InfectionStateBase::Dead>())) {
        mio::log_error("The length of at least one vector in initial_populations does not match the related number of "
                       "subcompartments.");
        return 1;
    }

    // Transfer the initial values in initial_populations to the vector init.
    Eigen::VectorXd init = Eigen::VectorXd::Zero(LctState::Count);
    init[(int)LctState::get_first_index<LctState::InfectionStateBase::Susceptible>()] =
        initial_populations[(int)LctState::InfectionStateBase::Susceptible][0];
    for (unsigned int i = 0; i < LctState::get_num_subcompartments<LctState::InfectionStateBase::Exposed>(); i++) {
        init[(int)LctState::get_first_index<LctState::InfectionStateBase::Exposed>() + i] =
            initial_populations[(int)LctState::InfectionStateBase::Exposed][i];
    }
    for (unsigned int i = 0; i < LctState::get_num_subcompartments<LctState::InfectionStateBase::InfectedNoSymptoms>();
         i++) {
        init[(int)LctState::get_first_index<LctState::InfectionStateBase::InfectedNoSymptoms>() + i] =
            initial_populations[(int)LctState::InfectionStateBase::InfectedNoSymptoms][i];
    }
    for (unsigned int i = 0; i < LctState::get_num_subcompartments<LctState::InfectionStateBase::InfectedSymptoms>();
         i++) {
        init[(int)LctState::get_first_index<LctState::InfectionStateBase::InfectedSymptoms>() + i] =
            initial_populations[(int)LctState::InfectionStateBase::InfectedSymptoms][i];
    }
    for (unsigned int i = 0; i < LctState::get_num_subcompartments<LctState::InfectionStateBase::InfectedSevere>();
         i++) {
        init[(int)LctState::get_first_index<LctState::InfectionStateBase::InfectedSevere>() + i] =
            initial_populations[(int)LctState::InfectionStateBase::InfectedSevere][i];
    }
    for (unsigned int i = 0; i < LctState::get_num_subcompartments<LctState::InfectionStateBase::InfectedCritical>();
         i++) {
        init[(int)LctState::get_first_index<LctState::InfectionStateBase::InfectedCritical>() + i] =
            initial_populations[(int)LctState::InfectionStateBase::InfectedCritical][i];
    }
    init[(int)LctState::get_first_index<LctState::InfectionStateBase::Recovered>()] =
        initial_populations[(int)LctState::InfectionStateBase::Recovered][0];
    init[(int)LctState::get_first_index<LctState::InfectionStateBase::Dead>()] =
        initial_populations[(int)LctState::InfectionStateBase::Dead][0];

    // Initialize model.
    Model model(std::move(init));

    // Set Parameters.
    model.parameters.get<mio::lsecir::TimeExposed>()            = 3.2;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>() = 2;
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
    mio::TimeSeries<ScalarType> result = mio::lsecir::simulate(0, tmax, 0.5, model);
    // Calculate the distribution in the InfectionState%s without subcompartments of the result and print it.
    mio::TimeSeries<ScalarType> population_no_subcompartments = model.calculate_populations(result);
    population_no_subcompartments.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 16, 8);
}