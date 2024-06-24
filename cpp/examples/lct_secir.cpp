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
#include "lct_secir/initializer_flows.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"

#include <vector>

int main()
{
    // Simple example to demonstrate how to run a simulation using an LCT-SECIR model.
    // Parameters, initial values and the number of subcompartments are not meant to represent a realistic scenario.

    using Model    = mio::lsecir::Model<2, 3, 1, 1, 5>;
    using LctState = Model::LctState;

    // Variable defines whether the class Initializer is used to define an initial vector from flows or whether a manually
    // defined initial vector is used to initialize the LCT model.
    bool use_initializer_flows = true;

    ScalarType tmax = 20;

    // Initialize a model.
    Model model(std::move(Eigen::VectorXd::Zero(LctState::Count)));

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

    if (use_initializer_flows) {
        // Example how to use the class Initializer for the definition of an initial vector for the LCT model.
        // Create TimeSeries with num_transitions elements.
        int num_transitions = (int)mio::lsecir::InfectionTransition::Count;
        mio::TimeSeries<ScalarType> flows(num_transitions);

        ScalarType dt = 0.1;

        mio::TimeSeries<ScalarType>::Vector vec_flows(num_transitions);
        vec_flows[(int)mio::lsecir::InfectionTransition::SusceptibleToExposed]                 = 2.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 1.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 4.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::InfectedSymptomsToInfectedSevere]     = 1.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::InfectedSymptomsToRecovered]          = 4.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::InfectedSevereToInfectedCritical]     = 1.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::InfectedSevereToRecovered]            = 1.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::InfectedCriticalToDead]               = 1.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::InfectedCriticalToRecovered]          = 1.0;
        // Add initial time point to time series.
        flows.add_time_point(-110, vec_flows);
        // Add further time points until time 0.
        while (flows.get_last_time() < -1e-10) {
            flows.add_time_point(flows.get_last_time() + dt, vec_flows);
        }

        // Set initialization vector for the LCT model.
        mio::lsecir::Initializer<Model> initializer(std::move(flows), model);
        initializer.set_tol_for_support_max(1e-6);
        auto status = initializer.compute_initialization_vector(1000000, 10, 16000);
        if (status) {
            return 1;
        }
    }
    else {
        // Simple example how to initialize model without flows.
        // Define the initial value vector init with the distribution of the population into subcompartments.
        // This method of defining the vector using a vector of vectors is a bit of overhead, but should remind you how
        // the entries of the initial value vector relate to the defined template parameters of the model or the number of subcompartments.
        // It is also possible to define the initial value vector directly.
        std::vector<std::vector<ScalarType>> initial_populations = {{750}, {30, 20},          {20, 10, 10}, {50},
                                                                    {50},  {10, 10, 5, 3, 2}, {20},         {10}};

        // Assert that initial_populations has the right shape.
        if (initial_populations.size() != (size_t)LctState::InfectionState::Count) {
            mio::log_error(
                "The number of vectors in initial_populations does not match the number of InfectionStates.");
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
            mio::log_error(
                "The length of at least one vector in initial_populations does not match the corresponding number of "
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
        // Set initial values for the model.
        model.set_initial_values(std::move(init));
    }

    // Perform a simulation.
    mio::TimeSeries<ScalarType> result = mio::lsecir::simulate(0, tmax, 0.5, model);
    // Calculate the distribution in the InfectionState%s without subcompartments of the result and print it.
    mio::TimeSeries<ScalarType> population_no_subcompartments = model.calculate_populations(result);
    population_no_subcompartments.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 16, 8);
}