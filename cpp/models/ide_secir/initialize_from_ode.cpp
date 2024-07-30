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
void get_flows_from_ode_compartments(mio::osecir::Model<ScalarType>& model_ode,
                                     mio::TimeSeries<ScalarType> compartments, mio::TimeSeries<ScalarType>& flows,
                                     ScalarType t_max, ScalarType t_window, ScalarType dt_reference,
                                     ScalarType dt_comparison)
{ //TODO: Check that flows is an empty TimeSeries at the beginning.
    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // Use scale_timesteps to get from index wrt ODE timestep to index wrt IDE timestep.
    // Here we assume that we solve the ODE model on a finer scale (or equal scale) than the IDE model.

    if (dt_comparison < 1e-10) {
        dt_comparison = dt_reference;
    }

    ScalarType scale_timesteps = dt_comparison / dt_reference;

    // Compute index variables with respect to dt_comparison.
    Eigen::Index t_window_index = Eigen::Index(std::ceil(t_window / dt_comparison));
    Eigen::Index t_max_index    = Eigen::Index(std::ceil(t_max / dt_comparison));

    Eigen::Index flows_start_index = t_max_index - t_window_index + 1;

    // Add time points to TimeSeries containing flows.
    for (Eigen::Index i = flows_start_index; i <= t_max_index; i++) {
        flows.add_time_point(i * dt_comparison, mio::TimeSeries<ScalarType>::Vector::Constant(num_transitions, 0));
        flows.get_last_value()[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)] +=
            compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::Susceptible)] -
            compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::Susceptible)];
    }

    // compute resulting flows as combination of change in compartments and previously computed flows

    // flow from E to C
    for (int i = flows_start_index; i <= t_max_index; i++) {
        flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] =
            compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::Exposed)] -
            compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::Exposed)] +
            flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)];
    }

    // flow from C to I and from C to R
    for (int i = flows_start_index; i <= t_max_index; i++) {
        flows[i -
              flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
            (1 -
             model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0]) *
            (compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptoms)] -
             compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptoms)] +
             flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)]);
        flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] =
            (model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0]) *
            (compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptoms)] -
             compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptoms)] +
             flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)]);
    }

    // flow from I to H and from I to R
    for (int i = flows_start_index; i <= t_max_index; i++) {
        flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)] =
            model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0] *
            (compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] -
             compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] +
             flows[i - flows_start_index]
                  [Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)]);
        flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)] =
            (1 - model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0]) *
            (compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] -
             compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] +
             flows[i - flows_start_index]
                  [Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)]);
    }

    // flow from H to U and from H to R
    for (int i = flows_start_index; i <= t_max_index; i++) {
        flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)] =
            model_ode.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()[(mio::AgeGroup)0] *
            (compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::InfectedSevere)] -
             compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::InfectedSevere)] +
             flows[i - flows_start_index]
                  [Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)]);
        flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)] =
            (1 - model_ode.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()[(mio::AgeGroup)0]) *
            (compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::InfectedSevere)] -
             compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::InfectedSevere)] +
             flows[i - flows_start_index]
                  [Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)]);
    }

    // flow from U to D and from U to R
    for (int i = flows_start_index; i <= t_max_index; i++) {
        flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)] =
            model_ode.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()[(mio::AgeGroup)0] *
            (compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::InfectedCritical)] -
             compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::InfectedCritical)] +
             flows[i - flows_start_index]
                  [Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)]);
        flows[i - flows_start_index][Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)] =
            (1 - model_ode.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()[(mio::AgeGroup)0]) *
            (compartments[scale_timesteps * (i - 1)][Eigen::Index(mio::osecir::InfectionState::InfectedCritical)] -
             compartments[scale_timesteps * i][Eigen::Index(mio::osecir::InfectionState::InfectedCritical)] +
             flows[i - flows_start_index]
                  [Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)]);
    }
}

void compute_initial_flows_for_ide_from_ode(mio::osecir::Model<ScalarType>& model_ode, mio::isecir::Model& model_ide,
                                            mio::TimeSeries<ScalarType> secihurd_ode, ScalarType t0_ide,
                                            ScalarType dt_ode, ScalarType dt_ide)
{
    std::cout << "Computing initial flows. \n";

    // Use t_window=t0_ide to get flows from t0 onwards.
    get_flows_from_ode_compartments(model_ode, secihurd_ode, model_ide.m_transitions, t0_ide, t0_ide, dt_ode, dt_ide);

    // Correct values in populations.
    if (model_ide.m_populations.get_last_time() != model_ide.m_transitions.get_last_time()) {
        model_ide.m_populations.remove_last_time_point();
        model_ide.m_populations.add_time_point<Eigen::VectorXd>(
            model_ide.m_transitions.get_last_time(),
            TimeSeries<ScalarType>::Vector::Constant((int)InfectionState::Count, 0));
        model_ide.m_populations[0][Eigen::Index(InfectionState::Dead)] =
            secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() -
                         (secihurd_ode.get_last_time() - t0_ide) / dt_ode - 1]
                        [(Eigen::Index)mio::osecir::InfectionState::Dead];
    }
}

} // namespace isecir
} // namespace mio