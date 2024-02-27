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

int main()
{
    /** Simple example to demonstrate how to run a simulation using an LCT SECIR model. 
    Parameters, initial values and subcompartments are not meant to represent a realistic scenario. */

    using Model    = mio::lsecir::Model<2, 3, 1, 1, 5>;
    using InfState = Model::InfState;

    ScalarType tmax = 20;

    /* Define the initial value vector init with the distribution of the population into subcompartments.
    This method of defining the vector is a bit of overhead, but should remind you how the entries of the initial 
    value vector relate to the defined template parameters of the model or the number of subcompartments.
    It is also possible to define the initial value vector directly.*/
    // Define initial values for each infection state with the appropriate number of subcompartments.
    ScalarType initial_value_Susceptible                                                           = 750;
    ScalarType initial_value_Exposed[InfState::get_num_subcompartments<InfState::Base::Exposed>()] = {30, 20};
    ScalarType
        initial_value_InfectedNoSymptoms[InfState::get_num_subcompartments<InfState::Base::InfectedNoSymptoms>()] = {
            20, 10, 10};
    ScalarType initial_value_InfectedSymptoms[InfState::get_num_subcompartments<InfState::Base::InfectedSymptoms>()] = {
        50};
    ScalarType initial_value_InfectedSevere[InfState::get_num_subcompartments<InfState::Base::InfectedSevere>()] = {50};
    ScalarType initial_value_InfectedCritical[InfState::get_num_subcompartments<InfState::Base::InfectedCritical>()] = {
        10, 10, 5, 3, 2};
    ScalarType initial_value_Recovered = 20;
    ScalarType initial_value_Dead      = 10;

    // Transfer the initial values to the vector init.
    Eigen::VectorXd init                                           = Eigen::VectorXd::Zero(InfState::Count);
    init[InfState::get_first_index<InfState::Base::Susceptible>()] = initial_value_Susceptible;
    for (unsigned int i = 0; i < InfState::get_num_subcompartments<InfState::Base::Exposed>(); i++) {
        init[InfState::get_first_index<InfState::Base::Exposed>() + i] = initial_value_Exposed[i];
    }
    for (unsigned int i = 0; i < InfState::get_num_subcompartments<InfState::Base::InfectedNoSymptoms>(); i++) {
        init[InfState::get_first_index<InfState::Base::InfectedNoSymptoms>() + i] = initial_value_InfectedNoSymptoms[i];
    }
    for (unsigned int i = 0; i < InfState::get_num_subcompartments<InfState::Base::InfectedSymptoms>(); i++) {
        init[InfState::get_first_index<InfState::Base::InfectedSymptoms>() + i] = initial_value_InfectedSymptoms[i];
    }
    for (unsigned int i = 0; i < InfState::get_num_subcompartments<InfState::Base::InfectedSevere>(); i++) {
        init[InfState::get_first_index<InfState::Base::InfectedSevere>() + i] = initial_value_InfectedSevere[i];
    }
    for (unsigned int i = 0; i < InfState::get_num_subcompartments<InfState::Base::InfectedCritical>(); i++) {
        init[InfState::get_first_index<InfState::Base::InfectedCritical>() + i] = initial_value_InfectedCritical[i];
    }
    init[InfState::get_first_index<InfState::Base::Recovered>()] = initial_value_Recovered;
    init[InfState::get_first_index<InfState::Base::Dead>()]      = initial_value_Dead;

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