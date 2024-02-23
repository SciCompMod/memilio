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

#include "memilio/config.h"
#include "memilio/utils/compiler_diagnostics.h"

//see below for line that causes this warning
GCC_CLANG_DIAGNOSTIC(push)
GCC_CLANG_DIAGNOSTIC(ignored "-Wmaybe-uninitialized")

#include "ide_secir/initialflows.h"
#include "memilio/utils/logging.h"
#include <iterator>

#ifdef MEMILIO_HAS_JSONCPP

#include "memilio/io/epi_data.h"
#include "memilio/utils/memory.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/utils/stl_util.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/date.h"

#include <json/json.h>

#include <boost/filesystem.hpp>

#include <numeric>
#include <vector>
#include <iostream>
#include <string>
#include <random>
#include <fstream>

namespace mio
{
namespace isecir
{
ScalarType compute_global_support_max(Model model, ScalarType dt)
{
    // for simplicity use this function for now to avoid having to deal with different support_max values for the different transitions
    // depending on the variance of the support_max values it woukd make more sense to treat each transition individually
    ScalarType global_support_max = std::max(
        {model.parameters.get<TransitionDistributions>()[(int)InfectionTransition::ExposedToInfectedNoSymptoms]
             .get_support_max(dt),
         model.parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
             .get_support_max(dt),
         model.parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedNoSymptomsToRecovered]
             .get_support_max(dt),
         model.parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSymptomsToInfectedSevere]
             .get_support_max(dt),
         model.parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSymptomsToRecovered]
             .get_support_max(dt),
         model.parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSevereToInfectedCritical]
             .get_support_max(dt),
         model.parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSevereToRecovered]
             .get_support_max(dt),
         model.parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedCriticalToDead]
             .get_support_max(dt),
         model.parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedCriticalToRecovered]
             .get_support_max(dt)});

    return global_support_max;
}

//TODO: do this in StateAgeFunction or is this only needed here? Or should we take these values from the literature/ consider them as known
// in the simulation?
ScalarType compute_mean(int idx_CurrentFlow)
{
    // ScalarType mean{};
    unused(idx_CurrentFlow);
    // use dummy value = 1 for now
    return 1.;
}

// Here we want to implement how to compute previous flow directly from RKI data by using the mean duration spent in
// some InfectionState
void compute_flows_with_mean(Model& model, Eigen::Index idx_CurrentFlow, Eigen::Index idx_OutgoingFlow, ScalarType dt,
                             Eigen::Index time_series_index)
{
    // compute mean
    ScalarType mean = compute_mean(idx_CurrentFlow);

    // compute how many time points we need to shift values to shift values using the mean
    Eigen::Index mean_index = (Eigen::Index)std::round(mean / dt);

    model.m_transitions.get_value(time_series_index)[Eigen::Index(idx_CurrentFlow)] =
        model.m_transitions.get_value(time_series_index + mean_index)[Eigen::Index(idx_OutgoingFlow)];
}

// TODO: this formula can lead to negative flows!!! Not used at the moment
// TODO: Think about good name for CurrentFlow and OutgoingFlow; also think about if we should use them as int
// or Eigen::Index in argument
// Here we want to implement how to compute the previous flow
void compute_previous_flows(Model& model, Eigen::Index idx_CurrentFlow, Eigen::Index idx_OutgoingFlow,
                            Eigen::Index time_series_index, ScalarType dt)
{
    // TODO: make this formula more similar to compute_flows regarding indices or is this more confusing?
    // TODO: case distinction for gamma'(dt) = 0 and mu = 0

    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(
        model.parameters.get<TransitionDistributions>()[idx_CurrentFlow].get_support_max(dt) / dt);

    ScalarType sum = 0.;
    for (Eigen::Index i = time_series_index + 1 - calc_time_index; i < time_series_index; i++) {

        ScalarType state_age = (time_series_index + 1 - i) * dt;

        // backward difference scheme to approximate first derivative of TransitionDistribution
        sum += (model.parameters.get<TransitionDistributions>()[idx_OutgoingFlow].eval(state_age) -
                model.parameters.get<TransitionDistributions>()[idx_OutgoingFlow].eval(state_age - dt)) /
               dt * model.m_transitions.get_value(i)[idx_CurrentFlow];
    }

    // std::cout << "Sum: " << sum << "\n";

    model.m_transitions.get_value(time_series_index)[Eigen::Index(idx_CurrentFlow)] =
        (-1 / (model.parameters.get<TransitionDistributions>()[idx_OutgoingFlow].eval(dt) -
               model.parameters.get<TransitionDistributions>()[idx_OutgoingFlow].eval(0))) *
        (model.m_transitions.get_value(time_series_index)[Eigen::Index(idx_OutgoingFlow)] /
             (dt * model.parameters.get<TransitionProbabilities>()[idx_OutgoingFlow]) +
         sum);
}

// we assume that we start the simulation at time 0 and want to compute the necessary flows
// in the past for the initialization of the model
void set_initial_flows(Model& model, ScalarType dt, ScalarType rki_cases_dummy, ScalarType rki_deaths_dummy)
{
    unused(rki_deaths_dummy);
    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // get (global) support_max to determine how many flows in the past we have to compute
    ScalarType global_support_max         = compute_global_support_max(model, dt);
    Eigen::Index global_support_max_index = std::ceil(global_support_max / dt);

    // remove time point that was needed for initialization of model
    // TODO: find out of there is a better way to handle that
    // maybe make flows an optional argument in construction of Model and then check in check constraints that there were flows set
    model.m_transitions.remove_last_time_point();

    // we already know the flows from C to I for a sufficient number of time points in the past

    // now we first want to compute the flows that are coming later, i.e. flows starting at I or at a later compartment

    // define start shift as the number of time pints in the past where we start the computation
    // this is the first timepoint in the m_transitions TimeSeries and we want to store the following flows at the right time point
    Eigen::Index start_shift = -4 * global_support_max_index;
    unused(start_shift);

    // TODO: we need values for t>0 so that we can shift values accordingly, remove them before starting our model?
    // write rki data into time series for transition from InfectedNoSymptoms to InfectedSymptoms
    for (int i = start_shift + 1; i <= 2 * global_support_max_index; i++) {
        // add time point
        model.m_transitions.add_time_point(i * dt, mio::TimeSeries<ScalarType>::Vector::Constant(num_transitions, 0.));
        // C to I
        model.m_transitions
            .get_last_value()[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
            rki_cases_dummy;
    }

    // Compute remaining flows
    // E to C
    // TODO: think about when to start for loops, this should also depend on the mean and not only global_max_support

    // start with computing necessary values using mean
    for (Eigen::Index i = start_shift + 1; i <= global_support_max_index; i++) {

        compute_flows_with_mean(model, Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms),
                                Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms),
                                dt, i - start_shift - 1);
    }

    // // compute remaining flows from E to C using formula based on flow, see Overleaf
    // for (Eigen::Index i = start_shift + 2 * global_support_max_index + 1; i <= 0; i++) {

    //     compute_previous_flows(model, Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms),
    //                            Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms),
    //                            i - start_shift - 1, dt);
    // }

    // S to E
    // start with computing necessary values using mean
    for (Eigen::Index i = start_shift + 1; i <= 0; i++) {
        std::cout << "Current time index: " << i << "\n";
        std::cout << "numtimepoints: " << model.m_transitions.get_num_time_points() << "\n";
        std::cout << "Index in TimeSeries: " << i - start_shift << std::endl;

        compute_flows_with_mean(model, Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed),
                                Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms), dt,
                                i - start_shift - 1);
    }

    // // compute remaining flows from S to E using formula based on flow, see Overleaf
    // // for (Eigen::Index i = start_shift + global_support_max_index + 1; i <= 0; i++) {
    // for (Eigen::Index i = start_shift + 1; i <= 0; i++) {

    //     compute_previous_flows(model, Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed),
    //                            Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms),
    //                            i - start_shift - 1, dt);
    // }

    // I to H for -3*global_support_max, ..., 0
    for (int i = -3 * global_support_max_index + 1; i <= 0; i++) {
        // I to H
        model.compute_flow(int(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere),
                           Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt, true,
                           i - start_shift);
    }

    // H to U for -2*global_support_max, ..., 0
    for (int i = -2 * global_support_max_index + 1; i <= 0; i++) {
        // H to U
        model.compute_flow((int)InfectionTransition::InfectedSevereToInfectedCritical,
                           Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere), dt, true,
                           i - start_shift);
    }

    // C, I, H, U to R and U to D for -1*global_support_max, ..., 0
    for (int i = -global_support_max_index + 1; i <= 0; i++) {
        // C to R
        model.compute_flow((int)InfectionTransition::InfectedNoSymptomsToRecovered,
                           Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms), dt, true, i - start_shift);
        // I to R
        model.compute_flow((int)InfectionTransition::InfectedSymptomsToRecovered,
                           Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt, true,
                           i - start_shift);
        // H to R
        model.compute_flow((int)InfectionTransition::InfectedSevereToRecovered,
                           Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere), dt, true,
                           i - start_shift);
        // U to R
        model.compute_flow((int)InfectionTransition::InfectedCriticalToRecovered,
                           Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical), dt, true,
                           i - start_shift);
        // U to D
        model.compute_flow((int)InfectionTransition::InfectedCriticalToDead,
                           Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical), dt, true,
                           i - start_shift);
    }
}
} // namespace isecir
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP