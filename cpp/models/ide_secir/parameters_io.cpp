/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Lena Ploetzke, Anna Wendler
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
#include "memilio/epidemiology/age_group.h"
#include <cstddef>
#include <iostream>

#ifdef MEMILIO_HAS_JSONCPP

#include "ide_secir/model.h"
#include "ide_secir/infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"
#include "memilio/utils/date.h"

#include <string>
#include <cmath>

namespace mio
{
namespace isecir
{

IOResult<void> set_initial_flows(Model& model, ScalarType dt, std::string const& path, Date date,
                                 ScalarType scale_confirmed_cases)
{
    //--- Preparations ---
    // Try to get RKI data from path.
    BOOST_OUTCOME_TRY(auto&& rki_data, mio::read_confirmed_cases_data(path));
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
    auto min_date_entry = std::min_element(rki_data.begin(), rki_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    auto min_date       = min_date_entry->date;

    // Get (global) support_max to determine how many flows in the past we have to compute.
    ScalarType global_support_max         = model.get_global_support_max(dt);
    Eigen::Index global_support_max_index = Eigen::Index(std::ceil(global_support_max / dt));

    // Get the number of AgeGroups.
    const size_t num_age_groups = ConfirmedCasesDataEntry::age_group_names.size();

    // transitions should be empty at the beginning.

    if (model.transitions.get_num_time_points() > 0) {
        model.transitions = TimeSeries<ScalarType>(Eigen::Index(InfectionTransition::Count) * num_age_groups);
    }
    if (model.populations.get_time(0) != 0) {
        model.populations.remove_last_time_point();
        model.populations.add_time_point<Eigen::VectorXd>(
            0, TimeSeries<ScalarType>::Vector::Constant((int)InfectionState::Count * num_age_groups, 0));
    }

    // The first time we need is -4 * global_support_max.
    Eigen::Index start_shift = 4 * global_support_max_index;
    // The last time needed is dependent on the mean stay time in the Exposed compartment and
    // the mean stay time of asymptomatic individuals in InfectedNoSymptoms.
    // The mean stay time in a compartment may be dependent on the AgeGroup.
    CustomIndexArray<ScalarType, AgeGroup> mean_ExposedToInfectedNoSymptoms =
        CustomIndexArray<ScalarType, AgeGroup>(AgeGroup(num_age_groups), 0.);
    CustomIndexArray<ScalarType, AgeGroup> mean_InfectedNoSymptomsToInfectedSymptoms =
        CustomIndexArray<ScalarType, AgeGroup>(AgeGroup(num_age_groups), 0.);
    CustomIndexArray<ScalarType, AgeGroup> mean_InfectedSymptomsToInfectedSevere =
        CustomIndexArray<ScalarType, AgeGroup>(AgeGroup(num_age_groups), 0.);
    CustomIndexArray<ScalarType, AgeGroup> mean_InfectedSevereToInfectedCritical =
        CustomIndexArray<ScalarType, AgeGroup>(AgeGroup(num_age_groups), 0.);
    CustomIndexArray<ScalarType, AgeGroup> mean_InfectedCriticalToDead =
        CustomIndexArray<ScalarType, AgeGroup>(AgeGroup(num_age_groups), 0.);
    Eigen::Index last_time_index_needed = 0;

    for (AgeGroup group = AgeGroup(0); group < AgeGroup(num_age_groups); group++) {
        // Set the Dead compartment to 0 so that RKI data can be added correctly.
        int Di                   = model.get_state_flat_index(Eigen::Index(InfectionState::Dead), group);
        model.populations[0][Di] = 0;

        mean_ExposedToInfectedNoSymptoms[group] =
            model.parameters
                .get<TransitionDistributions>()[group][Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms)]
                .get_mean(dt);
        mean_InfectedNoSymptomsToInfectedSymptoms[group] =
            model.parameters
                .get<TransitionDistributions>()[group]
                                               [Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)]
                .get_mean(dt);
        mean_InfectedSymptomsToInfectedSevere[group] =
            model.parameters
                .get<TransitionDistributions>()[group]
                                               [Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere)]
                .get_mean(dt);
        mean_InfectedSevereToInfectedCritical[group] =
            model.parameters
                .get<TransitionDistributions>()[group]
                                               [Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical)]
                .get_mean(dt);
        mean_InfectedCriticalToDead[group] =
            model.parameters
                .get<TransitionDistributions>()[group][Eigen::Index(InfectionTransition::InfectedCriticalToDead)]
                .get_mean(dt);
        if (last_time_index_needed <
            Eigen::Index(std::ceil(
                (mean_ExposedToInfectedNoSymptoms[group] + mean_InfectedNoSymptomsToInfectedSymptoms[group]) / dt))) {
            last_time_index_needed = Eigen::Index(std::ceil(
                (mean_ExposedToInfectedNoSymptoms[group] + mean_InfectedNoSymptomsToInfectedSymptoms[group]) / dt));
        }
    }

    // Create TimeSeries with zeros. The index of time zero is start_shift.
    for (Eigen::Index i = -start_shift; i <= last_time_index_needed; i++) {
        // Add time point.
        model.transitions.add_time_point(
            i * dt, TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count * num_age_groups, 0.));
    }

    // Setting the entries in m_total_confirmed_cases to zero before overwriting it with the RKI data.
    model.total_confirmed_cases = CustomIndexArray<ScalarType, AgeGroup>(AgeGroup(num_age_groups), 0.);
    //--- Calculate the flow InfectedNoSymptomsToInfectedSymptoms using the RKI data and store in the transitions object.---
    ScalarType min_offset_needed = std::ceil(
        model.transitions.get_time(0) -
        1); // Need -1 if first time point is integer and just the floor value if not, therefore use ceil and -1
    ScalarType max_offset_needed = std::ceil(model.transitions.get_last_time());
    bool min_offset_needed_avail = false;
    bool max_offset_needed_avail = false;

    // Go through the entries of rki_data and check if date is needed for calculation. Confirmed cases are scaled.
    // Define dummy variables to store the first and the last index of the TimeSeries where the considered entry of
    // rki_data is potentially needed.
    Eigen::Index idx_needed_first = 0;
    Eigen::Index idx_needed_last  = 0;
    ScalarType time_idx           = 0;
    for (auto&& entry : rki_data) {
        int offset     = get_offset_in_days(entry.date, date);
        AgeGroup group = entry.age_group;
        if ((offset >= min_offset_needed) && (offset <= max_offset_needed)) {
            if (offset == min_offset_needed) {
                min_offset_needed_avail = true;
            }
            if (offset == max_offset_needed) {
                max_offset_needed_avail = true;
            }
            // Smallest index for which the entry is needed.
            idx_needed_first =
                Eigen::Index(std::max(std::floor((offset - model.transitions.get_time(0) - 1) / dt), 0.));
            // Biggest index for which the entry is needed.
            idx_needed_last = Eigen::Index(std::min(std::ceil((offset - model.transitions.get_time(0) + 1) / dt),
                                                    double(model.transitions.get_num_time_points() - 1)));

            int INStISyi = model.get_transition_flat_index(
                Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), group);

            for (Eigen::Index i = idx_needed_first; i <= idx_needed_last; i++) {

                time_idx = model.transitions.get_time(i);
                if (offset == int(std::floor(time_idx))) {
                    model.transitions[i][INStISyi] +=
                        (1 - (time_idx - std::floor(time_idx))) * scale_confirmed_cases * entry.num_confirmed;
                }
                if (offset == int(std::ceil(time_idx))) {
                    model.transitions[i][INStISyi] +=
                        (time_idx - std::floor(time_idx)) * scale_confirmed_cases * entry.num_confirmed;
                }
                if (offset == int(std::floor(time_idx - dt))) {
                    model.transitions[i][INStISyi] -=
                        (1 - (time_idx - dt - std::floor(time_idx - dt))) * scale_confirmed_cases * entry.num_confirmed;
                }
                if (offset == int(std::ceil(time_idx - dt))) {
                    model.transitions[i][INStISyi] -=
                        (time_idx - dt - std::floor(time_idx - dt)) * scale_confirmed_cases * entry.num_confirmed;
                }
            }

            // Compute Dead by shifting RKI data according to mean stay times.
            // This is done because the RKI reports death with the date of positive test instead of the date of deaths.
            int Di = model.get_state_flat_index(Eigen::Index(InfectionState::Dead), group);
            if (offset ==
                std::floor(-mean_InfectedSymptomsToInfectedSevere[group] -
                           mean_InfectedSevereToInfectedCritical[group] - mean_InfectedCriticalToDead[group])) {
                model.populations[0][Di] +=
                    (1 -
                     (-mean_InfectedSymptomsToInfectedSevere[group] - mean_InfectedSevereToInfectedCritical[group] -
                      mean_InfectedCriticalToDead[group] -
                      std::floor(-mean_InfectedSymptomsToInfectedSevere[group] -
                                 mean_InfectedSevereToInfectedCritical[group] - mean_InfectedCriticalToDead[group]))) *
                    entry.num_deaths;
            }
            if (offset ==
                std::ceil(-mean_InfectedSymptomsToInfectedSevere[group] - mean_InfectedSevereToInfectedCritical[group] -
                          mean_InfectedCriticalToDead[group])) {
                model.populations[0][Di] +=
                    (-mean_InfectedSymptomsToInfectedSevere[group] - mean_InfectedSevereToInfectedCritical[group] -
                     mean_InfectedCriticalToDead[group] -
                     std::floor(-mean_InfectedSymptomsToInfectedSevere[group] -
                                mean_InfectedSevereToInfectedCritical[group] - mean_InfectedCriticalToDead[group])) *
                    entry.num_deaths;
            }

            if (offset == 0) {
                model.total_confirmed_cases[group] = scale_confirmed_cases * entry.num_confirmed;
            }
        }
    }

    if (!max_offset_needed_avail) {
        log_error("Necessary range of dates needed to compute initial values does not exist in RKI data.");
        return failure(StatusCode::OutOfRange, path + ", necessary range of dates does not exist in RKI data.");
    }
    if (!min_offset_needed_avail) {
        std::string min_date_string =
            std::to_string(min_date.day) + "." + std::to_string(min_date.month) + "." + std::to_string(min_date.year);
        // Get first date that is needed.
        mio::Date min_offset_date          = offset_date_by_days(date, int(min_offset_needed));
        std::string min_offset_date_string = std::to_string(min_offset_date.day) + "." +
                                             std::to_string(min_offset_date.month) + "." +
                                             std::to_string(min_offset_date.year);
        log_warning("RKI data is needed from " + min_offset_date_string +
                    " to compute initial values. RKI data is only available from " + min_date_string +
                    ". Missing dates were set to 0.");
    }

    //--- Calculate the flows "after" InfectedNoSymptomsToInfectedSymptoms. ---
    // Set m_support_max_vector and m_derivative_vector in the model which is needed for the following computations.
    model.set_transitiondistributions_support_max(dt);
    model.set_transitiondistributions_derivative(dt);

    for (AgeGroup group = AgeGroup(0); group < AgeGroup(num_age_groups); group++) {
        //--- Calculate the flows "after" InfectedNoSymptomsToInfectedSymptoms. ---
        // Compute flow InfectedSymptomsToInfectedSevere for -3 * global_support_max, ..., 0.
        for (Eigen::Index i = -3 * global_support_max_index; i <= 0; i++) {
            model.compute_flow(Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere),
                               Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt,
                               i + start_shift, group);
        }
        // Compute flow InfectedSevereToInfectedCritical for -2 * global_support_max, ..., 0.
        for (Eigen::Index i = -2 * global_support_max_index; i <= 0; i++) {
            model.compute_flow(Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical),
                               Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere), dt, i + start_shift,
                               group);
        }
        // Compute flows from InfectedSymptoms, InfectedSevere and InfectedCritical to Recovered and
        // flow InfectedCriticalToDead for -global_support_max, ..., 0.
        for (Eigen::Index i = -global_support_max_index; i <= 0; i++) {
            // Compute flow InfectedSymptomsToRecovered.
            model.compute_flow(Eigen::Index(InfectionTransition::InfectedSymptomsToRecovered),
                               Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt,
                               i + start_shift, group);
            // Compute flow InfectedSevereToRecovered.
            model.compute_flow(Eigen::Index(InfectionTransition::InfectedSevereToRecovered),
                               Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere), dt, i + start_shift,
                               group);
            // Compute flow InfectedCriticalToRecovered.
            model.compute_flow(Eigen::Index(InfectionTransition::InfectedCriticalToRecovered),
                               Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical), dt, i + start_shift,
                               group);
            // Compute flow InfectedCriticalToDead.
            model.compute_flow(Eigen::Index(InfectionTransition::InfectedCriticalToDead),
                               Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical), dt, i + start_shift,
                               group);
        }

        //--- Calculate the remaining flows. ---
        // Compute flow ExposedToInfectedNoSymptoms for -2 * global_support_max, ..., 0.
        // Use mean value of the TransitionDistribution InfectedNoSymptomsToInfectedSymptoms for the calculation.

        Eigen::Index index_shift_mean = Eigen::Index(std::round(mean_InfectedNoSymptomsToInfectedSymptoms[group] / dt));
        int EtINSi =
            model.get_transition_flat_index(Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms), group);
        int INStISyi = model.get_transition_flat_index(
            Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), group);
        for (Eigen::Index i = -2 * global_support_max_index; i <= 0; i++) {
            model.transitions[i + start_shift][EtINSi] =
                (1 / model.parameters.get<TransitionProbabilities>()[group][Eigen::Index(
                         InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)]) *
                model.transitions[i + start_shift + index_shift_mean][INStISyi];
        }

        // Compute flow SusceptibleToExposed for -global_support_max, ..., 0.
        // Use mean values of the TransitionDistribution ExposedToInfectedNoSymptoms and of the
        // TransitionDistribution InfectedNoSymptomsToInfectedSymptoms for the calculation.
        index_shift_mean = Eigen::Index(std::round(
            (mean_ExposedToInfectedNoSymptoms[group] + mean_InfectedNoSymptomsToInfectedSymptoms[group]) / dt));
        int StEi = model.get_transition_flat_index(Eigen::Index(InfectionTransition::SusceptibleToExposed), group);

        for (Eigen::Index i = -global_support_max_index; i <= 0; i++) {
            model.transitions[i + start_shift][StEi] =
                (1 / model.parameters.get<TransitionProbabilities>()[group][Eigen::Index(
                         InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)]) *
                model.transitions[i + start_shift + index_shift_mean][INStISyi];
        }

        // InfectedNoSymptomsToRecovered for -global_support_max, ..., 0.
        // If we previously calculated the transition ExposedToInfectedNoSymptoms, we can calculate this transition
        // using the standard formula.
        for (Eigen::Index i = -global_support_max_index; i <= 0; i++) {
            model.compute_flow(Eigen::Index(InfectionTransition::InfectedNoSymptomsToRecovered),
                               Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms), dt, i + start_shift,
                               group);
        }
    }

    // At the end of the calculation, delete all time points that are not required for the simulation.
    auto transitions_copy(model.transitions);
    model.transitions = TimeSeries<ScalarType>(Eigen::Index(InfectionTransition::Count) * num_age_groups);
    for (Eigen::Index i = -global_support_max_index; i <= 0; i++) {
        model.transitions.add_time_point(i * dt, transitions_copy.get_value(i + start_shift));
    }

    return mio::success();
}

} // namespace isecir
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP
