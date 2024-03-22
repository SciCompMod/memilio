/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Anna Wendler
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
#include "ide_secir/model_ide.h"
#include "ide_secir/parameters.h"
#include "ode_secir/model.h"
#include "infection_state.h"
#include "memilio/config.h"
#include "memilio/utils/logging.h"
#include "memilio/math/eigen.h"

namespace mio
{
namespace isecir
{
void get_flows_from_ode_compartments(mio::osecir::Model& model_ode, mio::TimeSeries<ScalarType> compartments,
                                     mio::TimeSeries<ScalarType>& flows, ScalarType t_window, ScalarType t_max,
                                     ScalarType dt_small, ScalarType dt_big)
{
    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // use scale_timesteps to get from index wrt ODE timestep to index wrt IDE timestep
    // here we assume that we solve the ODE model on a finer scale (or equal scale) than the IDE model
    ScalarType scale_timesteps = dt_big / dt_small;

    // compute index variables with respect to dt_big
    Eigen::Index t_window_index = Eigen::Index(std::ceil(t_window / dt_big));
    Eigen::Index t_max_index    = Eigen::Index(std::ceil(t_max / dt_big));

    Eigen::Index flows_start_index = t_max_index - t_window_index + 1;

    // add time points to TimeSeries containing flows
    for (Eigen::Index i = flows_start_index; i <= t_max_index; i++) {
        flows.add_time_point(i * dt_big, mio::TimeSeries<ScalarType>::Vector::Constant(num_transitions, 0));
        flows.get_last_value()[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)] +=
            compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::Susceptible)] -
            compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::Susceptible)];
    }

    // compute resulting flows as combination of change in compartments and previously computed flows

    // flow from E to C for -global_support_max, ..., 0
    for (int i = flows_start_index; i <= t_max_index; i++) {
        flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] =
            compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::Exposed)] -
            compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::Exposed)] +
            flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)];
    }

    // flow from C to I and from C to R for -global_support_max, ..., 0
    for (int i = flows_start_index; i <= t_max_index; i++) {
        flows[i -
              flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
            (1 - model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]) *
            (compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptoms)] -
             compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptoms)] +
             flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)]);
        flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] =
            (model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]) *
            (compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptoms)] -
             compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptoms)] +
             flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)]);
        ;
    }

    // flow from I to H and from I to R for -global_support_max, ..., 0
    for (int i = flows_start_index; i <= t_max_index; i++) {
        flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)] =
            model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0] *
            (compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] -
             compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] +
             flows[i - flows_start_index]
                  [Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)]);
        flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)] =
            (1 - model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]) *
            (compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] -
             compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] +
             flows[i - flows_start_index]
                  [Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)]);
    }

    // flow from H to U and from H to R for -global_support_max, ..., 0
    for (int i = flows_start_index; i <= t_max_index; i++) {
        flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)] =
            model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0] *
            (compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::InfectedSevere)] -
             compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::InfectedSevere)] +
             flows[i - flows_start_index]
                  [Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)]);
        flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)] =
            (1 - model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]) *
            (compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::InfectedSevere)] -
             compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::InfectedSevere)] +
             flows[i - flows_start_index]
                  [Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)]);
    }

    // flow from U to D and from U to R for -global_support_max, ..., 0
    for (int i = flows_start_index; i <= t_max_index; i++) {
        flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)] =
            model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0] *
            (compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::InfectedCritical)] -
             compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::InfectedCritical)] +
             flows[i - flows_start_index]
                  [Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)]);
        flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)] =
            (1 - model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]) *
            (compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::InfectedCritical)] -
             compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::InfectedCritical)] +
             flows[i - flows_start_index]
                  [Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)]);
    }
}

void compute_initial_flows_for_ide_from_ode(mio::osecir::Model& model_ode, mio::isecir::Model& model_ide,
                                            mio::TimeSeries<ScalarType> secihurd_ode, ScalarType t0_ide,
                                            ScalarType dt_ode, ScalarType dt_ide)
{
    std::cout << "Computing initial flows. \n";

    // get (global) support_max to determine how many flows in the past we have to compute
    ScalarType global_support_max = model_ide.get_global_support_max(dt_ide);

    get_flows_from_ode_compartments(model_ode, secihurd_ode, model_ide.m_transitions, global_support_max, t0_ide,
                                    dt_ode, dt_ide);
}

} // namespace isecir
} // namespace mio