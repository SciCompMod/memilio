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
    TODO talk: -implement a function like "check_quality" or constraints, to check eg if all the flows are nonnegative (and maybe set them to 0 otherwise)?
    - Discuss how to solve the problem with the mean or previous flow.
    - compute_flow in model.h: one index is eigen, one int...
    - Seperate function for reading data into flow?
    - GROSSES Problem: ist nicht unwahrscheinlich, dass gamma(t_1) null ist..
    */

// TODO: do this in StateAgeFunction or is this only needed here? Or should we take these values from the literature/ consider them as known
// in the simulation?
ScalarType compute_mean(Eigen::Index idx_CurrentFlow)
{
    unused(idx_CurrentFlow);
    // use dummy value = 1 for now
    return 1.;
}

// Here we want to implement how to compute previous flow directly from RKI data by using the mean duration spent in
// some InfectionState
void compute_flows_with_mean(Model& model, Eigen::Index idx_InfectionTransitions, Eigen::Index idx_OutgoingFlow,
                             ScalarType dt, Eigen::Index current_time_index)
{
    // compute mean
    ScalarType mean = compute_mean(idx_InfectionTransitions);

    // compute how many time points we need to shift values to shift values using the mean
    Eigen::Index mean_index = (Eigen::Index)(std::round(mean / dt));

    model.m_transitions.get_value(current_time_index)[Eigen::Index(idx_InfectionTransitions)] =
        (1 / model.parameters.get<TransitionProbabilities>()[idx_OutgoingFlow]) *
        model.m_transitions.get_value(current_time_index + mean_index)[Eigen::Index(idx_OutgoingFlow)];
}

// TODO: this formula can lead to negative flows!!! Not used at the moment
// TODO: Think about good name for CurrentFlow and OutgoingFlow; also think about if we should use them as int
// or Eigen::Index in argument
// Here we want to implement how to compute the previous flow
void compute_previous_flows(Model& model, Eigen::Index idx_CurrentFlow, Eigen::Index idx_OutgoingFlow,
                            Eigen::Index current_time_index, ScalarType dt)
{
    // TODO: make this formula more similar to compute_flows regarding indices or is this more confusing?
    // TODO: case distinction for gamma'(dt) = 0 and mu = 0

    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(
        model.parameters.get<TransitionDistributions>()[idx_CurrentFlow].get_support_max(dt) / dt);

    ScalarType sum = 0.;
    for (Eigen::Index i = current_time_index + 1 - calc_time_index; i < current_time_index; i++) {

        ScalarType state_age = (current_time_index + 1 - i) * dt;

        // backward difference scheme to approximate first derivative of TransitionDistribution
        sum += (model.parameters.get<TransitionDistributions>()[idx_OutgoingFlow].eval(state_age) -
                model.parameters.get<TransitionDistributions>()[idx_OutgoingFlow].eval(state_age - dt)) /
               dt * model.m_transitions.get_value(i)[idx_CurrentFlow];
    }

    std::cout << "Sum: " << sum << "\n";

    model.m_transitions.get_value(current_time_index)[Eigen::Index(idx_CurrentFlow)] =
        (-dt / (model.parameters.get<TransitionDistributions>()[idx_OutgoingFlow].eval(dt) -
                model.parameters.get<TransitionDistributions>()[idx_OutgoingFlow].eval(0))) *
        (model.m_transitions.get_value(current_time_index)[Eigen::Index(idx_OutgoingFlow)] /
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

    // The first time we need is -4 * global_support_max.
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
    // Define variables to store the first and the last index of the TimeSeries where the entry of rki_data is potentially needed.
    Eigen::Index idx_needed_first = 0;
    Eigen::Index idx_needed_last  = 0;
    ScalarType time_idx           = 0;
    for (auto&& entry : rki_data) {
        int offset = get_offset_in_days(entry.date, date);
        if (!(offset < min_offset_needed) && !(offset > max_offset_needed)) {
            if (offset == min_offset_needed) {
                min_offset_needed_avail = true;
            }
            // Smallest index for which the entry is needed is when offset=std::ceil(time_idx).
            idx_needed_first =
                Eigen::Index(std::max(std::floor((offset - model.m_transitions.get_time(0) - 1) / dt), 0.));
            // Biggest index for which the entry is needed is when offset=std::floor(time_idx-dt).
            idx_needed_last = Eigen::Index(std::min(std::ceil((offset - model.m_transitions.get_time(0) + 1) / dt),
                                                    double(model.m_transitions.get_num_time_points() - 1)));

            for (Eigen::Index i = idx_needed_first; i <= idx_needed_last; i++) {
                time_idx = model.m_transitions.get_time(i);
                if (offset == int(std::floor(time_idx))) {
                    model.m_transitions[i][Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] +=
                        (1 - (time_idx - std::floor(time_idx))) * scale_confirmed_cases * entry.num_confirmed;
                }
                if (offset == int(std::ceil(time_idx))) {
                    model.m_transitions[i][Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] +=
                        (time_idx - std::floor(time_idx)) * scale_confirmed_cases * entry.num_confirmed;
                }
                if (offset == int(std::floor(time_idx - dt))) {
                    model.m_transitions[i][Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] -=
                        (1 - (time_idx - dt - std::floor(time_idx - dt))) * scale_confirmed_cases * entry.num_confirmed;
                }
                if (offset == int(std::ceil(time_idx - dt))) {
                    model.m_transitions[i][Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] -=
                        (time_idx - dt - std::floor(time_idx - dt)) * scale_confirmed_cases * entry.num_confirmed;
                }
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

    //--- Calculate the flows "after" InfectedNoSymptomsToInfectedSymptoms. ---
    // I to H for -3 * global_support_max, ..., 0
    for (Eigen::Index i = -3 * global_support_max_index; i <= 0; i++) {
        model.compute_flow(int(InfectionTransition::InfectedSymptomsToInfectedSevere),
                           Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt,
                           i - start_shift);
    }
    // H to U for -2 * global_support_max, ..., 0
    for (Eigen::Index i = -2 * global_support_max_index; i <= 0; i++) {
        model.compute_flow((int)InfectionTransition::InfectedSevereToInfectedCritical,
                           Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere), dt, i - start_shift);
    }
    // I, H, U to R and U to D for -global_support_max, ..., 0
    for (Eigen::Index i = -global_support_max_index; i <= 0; i++) {
        // I to R
        model.compute_flow((int)InfectionTransition::InfectedSymptomsToRecovered,
                           Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt,
                           i - start_shift);
        // H to R
        model.compute_flow((int)InfectionTransition::InfectedSevereToRecovered,
                           Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere), dt, i - start_shift);
        // U to R
        model.compute_flow((int)InfectionTransition::InfectedCriticalToRecovered,
                           Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical), dt, i - start_shift);
        // U to D
        model.compute_flow((int)InfectionTransition::InfectedCriticalToDead,
                           Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical), dt, i - start_shift);
    }

    //--- Calculate the remaining flows. ---
    // E to C for -global_support_max, ..., 0
    for (Eigen::Index i = -global_support_max_index; i <= 0; i++) {
        compute_flows_with_mean(model, Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms),
                                Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt,
                                i - start_shift - 1
        // Eigen::Index idx_CurrentFlow, Eigen::Index idx_OutgoingFlow, ScalarType dt,
        // Eigen::Index current_time_index)
        //model.compute_flow((int)InfectionTransition::InfectedNoSymptomsToRecovered,
        //             Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms), dt, i - start_shift);
    }

    //------------
    /*// TODO: think about when to start for loops, this should also depend on the mean and not only global_max_support
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

    // S to E for -global_support_max, ..., 0
    // start with computing necessary values using mean
    for (Eigen::Index i = start_shift + 1; i <= 0; i++) {

        compute_flows_with_mean(model, Eigen::Index(InfectionTransition::SusceptibleToExposed),
                                Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms), dt,
                                i - start_shift - 1);
    }

    // // compute remaining flows from S to E using formula based on flow, see Overleaf
    // // for (Eigen::Index i = start_shift + global_support_max_index + 1; i <= 0; i++) {
    // for (Eigen::Index i = start_shift + 1; i <= 0; i++) {

    //     compute_previous_flows(model, Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed),
    //                            Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms),
    //                            i - start_shift - 1, dt);
    // }*/

    // C to R for -global_support_max, ..., 0
    // If we previously calculated the transition from E to C, we can calculate this transition using the standard formula.
    for (Eigen::Index i = -global_support_max_index; i <= 0; i++) {
        model.compute_flow((int)InfectionTransition::InfectedNoSymptomsToRecovered,
                           Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms), dt, i - start_shift);
    }

    // At the end of the calculation, delete all time points that are not required for the simulation.
    auto transition_copy(model.m_transitions);
    model.m_transitions = TimeSeries<ScalarType>((int)InfectionTransition::Count);
    for (int i = -global_support_max_index; i <= 0; i++) {
        model.m_transitions.add_time_point(i * dt, transition_copy.get_value(i - start_shift));
    }

    return mio::success();
}

} // namespace isecir
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP