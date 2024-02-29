/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Anna Wendler, Lena Ploetzke
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

#include "ide_secir/parameters_io.h"
#include "memilio/config.h"
#include "memilio/utils/compiler_diagnostics.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "memilio/io/epi_data.h"

#include "ide_secir/model.h"
#include "ide_secir/infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/io/io.h"
#include "memilio/utils/date.h"

#include <string>
#include <iostream>
#include <cmath>

namespace mio
{
namespace isecir
{
/* 
    TODO: -implement a function like "check_quality" or constraints, to check eg if all the flows are nonnegative (and maybe set them to 0 otherwise)
    - Discuss how to solve the problem with the mean or previous flow.
    */

//TODO: do this in StateAgeFunction or is this only needed here? Or should we take these values from the literature/ consider them as known
// in the simulation?
ScalarType compute_mean(Eigen::Index idx_CurrentFlow)
{
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
    Eigen::Index mean_index = (Eigen::Index)(std::round(mean / dt));

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
IOResult<void> set_initial_flows(Model& model, ScalarType dt, std::string const& path, Date date,
                                 ScalarType scale_confirmed_cases)
{
    //--- Preparations ---
    // Try to get rki data from path.
    BOOST_OUTCOME_TRY(rki_data, mio::read_confirmed_cases_noage(path));
    auto max_date_entry = std::max_element(rki_data.begin(), rki_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == rki_data.end()) {
        log_error("RKI data file is empty.");
        return failure(StatusCode::InvalidFileFormat, path + ", file is empty.");
    }
    auto max_date = max_date_entry->date;
    if (max_date < date) {
        log_error("Specified date does not exist in RKI data.");
        return failure(StatusCode::OutOfRange, path + ", specified date does not exist in RKI data.");
    }

    // Get (global) support_max to determine how many flows in the past we have to compute
    ScalarType global_support_max         = model.get_global_support_max(dt);
    Eigen::Index global_support_max_index = Eigen::Index(std::ceil(global_support_max / dt));

    // m_transitions should be empty at the beginning.
    if (model.m_transitions.get_num_time_points() > 0) {
        model.m_transitions = TimeSeries<ScalarType>((int)InfectionTransition::Count);
    }

    // The first time we need is -4 * global_support_max
    Eigen::Index start_shift = -4 * global_support_max_index;
    // Create TimeSeries with zeros.
    for (Eigen::Index i = start_shift; i <= 2 * global_support_max_index; i++) {
        // Add time point.
        model.m_transitions.add_time_point(
            i * dt, TimeSeries<ScalarType>::Vector::Constant((int)InfectionTransition::Count, 0.));
    }

    //--- Calculate the flow InfectedNoSymptomsToInfectedSymptoms using the RKI data and store in the m_transitions object.---
    ScalarType min_offset_needed = std::ceil(
        model.m_transitions.get_time(0) -
        1); //Need -1 if first time point is integer and just the floor value if not, therefore use ceil and -1
    ScalarType max_offset_needed = std::ceil(model.m_transitions.get_last_time());

    bool min_offset_needed_avail = false;
    bool max_offset_needed_avail = false;
    // Go through the entries of rki_data and check if date is needed for calculation. Confirmed cases are scaled.
    Eigen::Index dummy_idx = 0;
    for (auto&& entry : rki_data) {
        auto offset = get_offset_in_days(entry.date, date);
        if (!(offset < min_offset_needed) && !(offset > max_offset_needed)) {
            if (offset == min_offset_needed) {
                min_offset_needed_avail = true;
            }
            dummy_idx = Eigen::Index(std::ceil((offset - model.m_transitions.get_time(0)) / dt));
            if (dummy_idx >= 0) {
                model
                    .m_transitions[dummy_idx][Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
                    scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == max_offset_needed) {
                max_offset_needed_avail = true;
            }
            if (offset == 0) {
                model.m_populations[0][(int)InfectionState::Dead] = entry.num_deaths;
                model.m_total_confirmed_cases                     = scale_confirmed_cases * entry.num_confirmed;
            }
        }
    }

    if (!max_offset_needed_avail || !min_offset_needed_avail) {
        log_error("Necessary range of dates needed to compute initial values does not exist in RKI data.");
        return failure(StatusCode::OutOfRange, path + ", necessary range of dates does not exist in RKI data.");
    }

    // ----------------
    // The first time we need is -4 * global_support_max

    // we already know the flows from C to I for a sufficient number of time points in the past

    // now we first want to compute the flows that are coming later, i.e. flows starting at I or at a later compartment

    // define start shift as the number of time pints in the past where we start the computation
    // this is the first timepoint in the m_transitions TimeSeries and we want to store the following flows at the right time point

    // TODO: we need values for t>0 so that we can shift values accordingly, remove them before starting our model?
    // TODO: write rki data into time series for transition from InfectedNoSymptoms to InfectedSymptoms
    /*for (Eigen::Index i = start_shift + 1; i <= 2 * global_support_max_index; i++) {
        // add time point
        model.m_transitions.add_time_point(
            i * dt, TimeSeries<ScalarType>::Vector::Constant((int)InfectionTransition::Count, 0.));
        // C to I
        model.m_transitions.get_last_value()[Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
            rki_cases_dummy;
    }

    // Compute remaining flows
    // E to C
    // TODO: think about when to start for loops, this should also depend on the mean and not only global_max_support

    // start with computing necessary values using mean
    for (Eigen::Index i = start_shift + 1; i <= global_support_max_index; i++) {

        compute_flows_with_mean(model, Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms),
                                Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt,
                                i - start_shift - 1);
    }

    // // compute remaining flows from E to C using formula based on flow, see Overleaf
    // for (Eigen::Index i = start_shift + 2 * global_support_max_index + 1; i <= 0; i++) {

    //     compute_previous_flows(model, Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms),
    //                            Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms),
    //                            i - start_shift - 1, dt);
    // }

    // S to E
    // start with computing necessary values using mean
    for (Eigen::Index i = start_shift + 1; i <= 0; i++) {*/
    /*std::cout << "Current time index: " << i << "\n";
        std::cout << "numtimepoints: " << model.m_transitions.get_num_time_points() << "\n";
        std::cout << "Index in TimeSeries: " << i - start_shift << std::endl;*/

    /* compute_flows_with_mean(model, Eigen::Index(InfectionTransition::SusceptibleToExposed),
                                Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms), dt,
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
    for (Eigen::Index i = -3 * global_support_max_index + 1; i <= 0; i++) {
        // I to H
        model.compute_flow(int(InfectionTransition::InfectedSymptomsToInfectedSevere),
                           Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt, true,
                           i - start_shift);
    }

    // H to U for -2*global_support_max, ..., 0
    for (Eigen::Index i = -2 * global_support_max_index + 1; i <= 0; i++) {
        // H to U
        model.compute_flow((int)InfectionTransition::InfectedSevereToInfectedCritical,
                           Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere), dt, true,
                           i - start_shift);
    }

    // C, I, H, U to R and U to D for -1*global_support_max, ..., 0
    for (Eigen::Index i = -global_support_max_index + 1; i <= 0; i++) {
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
*/
    // At the end of the computation delete all timepoints in the future in m_transitions.
    // TODO: maybe change this and just keep from -max_supp until 0. Too many negative flows are unneccessary but not a problem, too many positive points are a problem.
    while (model.m_transitions.get_last_time() > 0) {
        model.m_transitions.remove_last_time_point();
    }
    return mio::success();
}

} // namespace isecir
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP