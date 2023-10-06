/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
    /** Simple example to demonstrate how to simulate using an LCT SECIR model. 
    Parameters, initial values and subcompartments are not realistic. */

    // Set vector that specifies the number of subcompartments.
    std::vector<int> num_subcompartments((int)mio::lsecir::InfectionStateBase::Count, 1);
    num_subcompartments[(int)mio::lsecir::InfectionStateBase::Exposed]            = 2;
    num_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] = 3;
    num_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedCritical]   = 5;
    mio::lsecir::InfectionState infectionState(num_subcompartments);

    ScalarType tmax = 20;

    // Define initial distribution of the population in the subcompartments.
    Eigen::VectorXd init(infectionState.get_count());
    init[infectionState.get_firstindex(mio::lsecir::InfectionStateBase::Susceptible)]            = 750;
    init[infectionState.get_firstindex(mio::lsecir::InfectionStateBase::Exposed)]                = 30;
    init[infectionState.get_firstindex(mio::lsecir::InfectionStateBase::Exposed) + 1]            = 20;
    init[infectionState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedNoSymptoms)]     = 20;
    init[infectionState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedNoSymptoms) + 1] = 10;
    init[infectionState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedNoSymptoms) + 2] = 10;
    init[infectionState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedSymptoms)]       = 50;
    init[infectionState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedSevere)]         = 50;
    init[infectionState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical)]       = 10;
    init[infectionState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 1]   = 10;
    init[infectionState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 2]   = 5;
    init[infectionState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 3]   = 3;
    init[infectionState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 4]   = 2;
    init[infectionState.get_firstindex(mio::lsecir::InfectionStateBase::Recovered)]              = 20;
    init[infectionState.get_firstindex(mio::lsecir::InfectionStateBase::Dead)]                   = 10;

    // Initialize model.
    mio::lsecir::Model model(std::move(init), infectionState);

    // Set Parameters.
    model.parameters.get<mio::lsecir::TimeExposed>()            = 2 * 4.2 - 5.2;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>() = 2 * (5.2 - 4.2);
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()   = 5.8;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()     = 9.5;
    // Also possible to change values with setter.
    model.parameters.set<mio::lsecir::TimeInfectedCritical>(7.1);

    model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.05;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::lsecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(5.));

    model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 0.7;
    model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 0.25;
    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.09;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.2;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.25;
    model.parameters.set<mio::lsecir::DeathsPerCritical>(0.3);

    // Perform a simulation.
    mio::TimeSeries<ScalarType> result = mio::lsecir::simulate(0, tmax, 0.5, model);
    // Calculate the distribution in infectionState without subcompartments of the result and print it.
    mio::TimeSeries<ScalarType> populations = model.calculate_populations(result);
    mio::lsecir::print_TimeSeries(populations, model.get_heading_CompartmentsBase());
}