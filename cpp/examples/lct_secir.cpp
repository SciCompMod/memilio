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
#include <vector>

int main()
{
    /** Simple example to demonstrate how to run a simulation using an LCT SECIR model. 
    Parameters, initial values and subcompartments are not meant to represent a realistic scenario. */

    mio::lsecir::InfectionState<mio::lsecir::InfectionStateBase, 1, 2, 3, 1, 1, 5, 1, 1> infection_state;

    ScalarType tmax = 20;

    // Define initial distribution of the population in the subcompartments.
    Eigen::VectorXd init(infection_state.get_count());
    init[infection_state.get_firstindex<mio::lsecir::InfectionStateBase::Susceptible>()]            = 750;
    init[infection_state.get_firstindex<mio::lsecir::InfectionStateBase::Exposed>()]                = 30;
    init[infection_state.get_firstindex<mio::lsecir::InfectionStateBase::Exposed>() + 1]            = 20;
    init[infection_state.get_firstindex<mio::lsecir::InfectionStateBase::InfectedNoSymptoms>()]     = 20;
    init[infection_state.get_firstindex<mio::lsecir::InfectionStateBase::InfectedNoSymptoms>() + 1] = 10;
    init[infection_state.get_firstindex<mio::lsecir::InfectionStateBase::InfectedNoSymptoms>() + 2] = 10;
    init[infection_state.get_firstindex<mio::lsecir::InfectionStateBase::InfectedSymptoms>()]       = 50;
    init[infection_state.get_firstindex<mio::lsecir::InfectionStateBase::InfectedSevere>()]         = 50;
    init[infection_state.get_firstindex<mio::lsecir::InfectionStateBase::InfectedCritical>()]       = 10;
    init[infection_state.get_firstindex<mio::lsecir::InfectionStateBase::InfectedCritical>() + 1]   = 10;
    init[infection_state.get_firstindex<mio::lsecir::InfectionStateBase::InfectedCritical>() + 2]   = 5;
    init[infection_state.get_firstindex<mio::lsecir::InfectionStateBase::InfectedCritical>() + 3]   = 3;
    init[infection_state.get_firstindex<mio::lsecir::InfectionStateBase::InfectedCritical>() + 4]   = 2;
    init[infection_state.get_firstindex<mio::lsecir::InfectionStateBase::Recovered>()]              = 20;
    init[infection_state.get_firstindex<mio::lsecir::InfectionStateBase::Dead>()]                   = 10;

    // Initialize model.
    mio::lsecir::Model<2, 3, 1, 1, 5> model(std::move(init));

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