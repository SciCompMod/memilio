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

void compute_initial_flows_from_ode_compartments(mio::osecir::Model& model_ode, mio::isecir::Model& model_ide,
                                                 mio::TimeSeries<ScalarType> secihurd_ode, ScalarType t0_ide,
                                                 ScalarType dt)
{
    std::cout << "Computing initial flows. \n";
    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // get (global) support_max to determine how many flows in .the past we have to compute
    ScalarType global_support_max         = model_ide.get_global_support_max(dt);
    Eigen::Index global_support_max_index = std::ceil(global_support_max / dt);

    // remove time point
    model_ide.m_transitions.remove_last_time_point();

    ScalarType t0_ide_index = std::ceil(t0_ide / dt);
    unused(secihurd_ode);
    int init_start_index = t0_ide_index - global_support_max_index + 1;
    // flow from S to E for -6*global_support_max, ..., 0 (directly from compartments)
    // add time points to init_transitions here
    for (int i = init_start_index; i <= t0_ide_index; i++) {
        model_ide.m_transitions.add_time_point(i * dt,
                                               mio::TimeSeries<ScalarType>::Vector::Constant(num_transitions, 0));
        model_ide.m_transitions.get_last_value()[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)] =
            secihurd_ode[i - 1][Eigen::Index(mio::osecir::InfectionState::Susceptible)] -
            secihurd_ode[i][Eigen::Index(mio::osecir::InfectionState::Susceptible)];
    }

    // compute resulting flows as combination of change in compartments and previously computed flows
    // Eigen::Index start_shift = t0_ide_index - 6 * global_support_max_index;

    // flow from E to C for -global_support_max, ..., 0
    for (int i = init_start_index; i <= t0_ide_index; i++) {
        model_ide.m_transitions[i - init_start_index]
                               [Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] =
            secihurd_ode[i - 1][Eigen::Index(mio::osecir::InfectionState::Exposed)] -
            secihurd_ode[i][Eigen::Index(mio::osecir::InfectionState::Exposed)] +
            model_ide.m_transitions[i - init_start_index]
                                   [Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)];
    }

    // flow from C to I and from C to R for -global_support_max, ..., 0
    for (int i = init_start_index; i <= t0_ide_index; i++) {
        model_ide.m_transitions[i - init_start_index]
                               [Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
            (1 - model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]) *
            (secihurd_ode[i - 1][Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptoms)] -
             secihurd_ode[i][Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptoms)] +
             model_ide.m_transitions[i - init_start_index]
                                    [Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)]);
        model_ide.m_transitions[i - init_start_index]
                               [Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] =
            (model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]) *
            (secihurd_ode[i - 1][Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptoms)] -
             secihurd_ode[i][Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptoms)] +
             model_ide.m_transitions[i - init_start_index]
                                    [Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)]);
        ;
    }

    // flow from I to H and from I to R for -global_support_max, ..., 0
    for (int i = init_start_index; i <= t0_ide_index; i++) {
        model_ide.m_transitions[i - init_start_index]
                               [Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)] =
            model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0] *
            (secihurd_ode[i - 1][Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] -
             secihurd_ode[i][Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] +
             model_ide.m_transitions[i - init_start_index][Eigen::Index(
                 mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)]);
        model_ide.m_transitions[i - init_start_index]
                               [Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)] =
            (1 - model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]) *
            (secihurd_ode[i - 1][Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] -
             secihurd_ode[i][Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] +
             model_ide.m_transitions[i - init_start_index][Eigen::Index(
                 mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)]);
    }

    // flow from H to U and from H to R for -global_support_max, ..., 0
    for (int i = init_start_index; i <= t0_ide_index; i++) {
        model_ide.m_transitions[i - init_start_index]
                               [Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)] =
            model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0] *
            (secihurd_ode[i - 1][Eigen::Index(mio::osecir::InfectionState::InfectedSevere)] -
             secihurd_ode[i][Eigen::Index(mio::osecir::InfectionState::InfectedSevere)] +
             model_ide.m_transitions[i - init_start_index]
                                    [Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)]);
        model_ide.m_transitions[i - init_start_index]
                               [Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)] =
            (1 - model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]) *
            (secihurd_ode[i - 1][Eigen::Index(mio::osecir::InfectionState::InfectedSevere)] -
             secihurd_ode[i][Eigen::Index(mio::osecir::InfectionState::InfectedSevere)] +
             model_ide.m_transitions[i - init_start_index]
                                    [Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)]);
    }

    // flow from U to D and from U to R for -global_support_max, ..., 0
    for (int i = init_start_index; i <= t0_ide_index; i++) {
        model_ide.m_transitions[i - init_start_index]
                               [Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)] =
            model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0] *
            (secihurd_ode[i - 1][Eigen::Index(mio::osecir::InfectionState::InfectedCritical)] -
             secihurd_ode[i][Eigen::Index(mio::osecir::InfectionState::InfectedCritical)] +
             model_ide.m_transitions[i - init_start_index]
                                    [Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)]);
        model_ide.m_transitions[i - init_start_index]
                               [Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)] =
            (1 - model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]) *
            (secihurd_ode[i - 1][Eigen::Index(mio::osecir::InfectionState::InfectedCritical)] -
             secihurd_ode[i][Eigen::Index(mio::osecir::InfectionState::InfectedCritical)] +
             model_ide.m_transitions[i - init_start_index]
                                    [Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)]);
    }
}

} // namespace isecir
} // namespace mio