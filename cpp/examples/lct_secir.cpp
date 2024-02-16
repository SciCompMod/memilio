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
#include <vector>

int main()
{
    /** Simple example to demonstrate how to run a simulation using an LCT SECIR model. 
    Parameters, initial values and subcompartments are not meant to represent a realistic scenario. */

    // Variable defines whether the class Initializer is used to define an initial vector from flows or whether a manually
    // defined initial vector is used to initialize the LCT model.
    bool use_initializer_flows = true;

    ScalarType tmax = 20;

    // Set vector that specifies the number of subcompartments.
    std::vector<int> num_subcompartments((int)mio::lsecir::InfectionStateBase::Count, 1);
    num_subcompartments[(int)mio::lsecir::InfectionStateBase::Exposed]            = 2;
    num_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] = 3;
    num_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedCritical]   = 5;
    mio::lsecir::InfectionState infection_state(num_subcompartments);

    // Define parameters used for simulation (and perhaps initialization).
    mio::lsecir::Parameters parameters;
    parameters.get<mio::lsecir::TimeExposed>()            = 3.2;
    parameters.get<mio::lsecir::TimeInfectedNoSymptoms>() = 2;
    parameters.get<mio::lsecir::TimeInfectedSymptoms>()   = 5.8;
    parameters.get<mio::lsecir::TimeInfectedSevere>()     = 9.5;
    // It is also possible to change values with the set function.
    parameters.set<mio::lsecir::TimeInfectedCritical>(7.1);

    parameters.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.05;

    mio::ContactMatrixGroup& contact_matrix = parameters.get<mio::lsecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
    // From SimulationTime 5, the contact pattern is reduced to 30% of the initial value.
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(5.));

    parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 0.7;
    parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 0.25;
    parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.09;
    parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.2;
    parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.25;
    parameters.set<mio::lsecir::DeathsPerCritical>(0.3);

    // Define initial value vector with subcompartments with the method defined in use_initializer_flows.
    Eigen::VectorXd init(infection_state.get_count());

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
        while (flows.get_last_time() < 0) {
            flows.add_time_point(flows.get_last_time() + dt, vec_flows);
        }

        // Get initialization vector for LCT model with subcompartments defined in infection_state.
        mio::lsecir::Initializer initializer(std::move(flows), infection_state, std::move(parameters));
        initializer.set_tol_for_support_max(1e-6);
        init = initializer.compute_initializationvector(1000000, 10, 16000);
    }
    else {
        // Simple example how to initialize model.
        // Define initial distribution of the population in the subcompartments.
        init[infection_state.get_firstindex(mio::lsecir::InfectionStateBase::Susceptible)]            = 750;
        init[infection_state.get_firstindex(mio::lsecir::InfectionStateBase::Exposed)]                = 30;
        init[infection_state.get_firstindex(mio::lsecir::InfectionStateBase::Exposed) + 1]            = 20;
        init[infection_state.get_firstindex(mio::lsecir::InfectionStateBase::InfectedNoSymptoms)]     = 20;
        init[infection_state.get_firstindex(mio::lsecir::InfectionStateBase::InfectedNoSymptoms) + 1] = 10;
        init[infection_state.get_firstindex(mio::lsecir::InfectionStateBase::InfectedNoSymptoms) + 2] = 10;
        init[infection_state.get_firstindex(mio::lsecir::InfectionStateBase::InfectedSymptoms)]       = 50;
        init[infection_state.get_firstindex(mio::lsecir::InfectionStateBase::InfectedSevere)]         = 50;
        init[infection_state.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical)]       = 10;
        init[infection_state.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 1]   = 10;
        init[infection_state.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 2]   = 5;
        init[infection_state.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 3]   = 3;
        init[infection_state.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 4]   = 2;
        init[infection_state.get_firstindex(mio::lsecir::InfectionStateBase::Recovered)]              = 20;
        init[infection_state.get_firstindex(mio::lsecir::InfectionStateBase::Dead)]                   = 10;
    }

    // Initialize model.
    mio::lsecir::Model model(std::move(init), infection_state, std::move(parameters));

    // Perform a simulation.
    mio::TimeSeries<ScalarType> result = mio::lsecir::simulate(0, tmax, 0.5, model);
    // Calculate the distribution in the InfectionState%s without subcompartments of the result and print it.
    mio::TimeSeries<ScalarType> population_no_subcompartments = model.calculate_populations(result);
    population_no_subcompartments.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 16, 8);
}