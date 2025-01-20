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
#ifndef IDE_INITIALFLOWS_H
#define IDE_INITIALFLOWS_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "ide_secir/model.h"
#include "ide_secir/infection_state.h"

#include "memilio/epidemiology/age_group.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"
#include "memilio/utils/date.h"
#include "memilio/utils/logging.h"

#include <string>

namespace mio
{
namespace isecir
{

/**
* @brief Computes a TimeSeries of flows to provide initial data for an IDE-SECIR model with data from RKI.
*   
* The flows InfectedNoSymptomsToInfectedSymptoms are calculated using the confirmed cases in the RKI data.
* If necessary, the RKI data are linearly interpolated within one day.
* The RKI data should contain data for each needed day without division of age groups, the completeness of the dates is not verified.
* Data can be downloaded e.g. with the file pycode/memilio-epidata/memilio/epidata/getCaseData.py, 
* which creates a file named cases_all_germany.json or a similar name. One should set impute_dates=True so that missing dates are imputed.
*
* The flows InfectedSymptomsToInfectedSevere, InfectedSymptomsToRecovered, InfectedSevereToInfectedCritical,
* InfectedSevereToRecovered, InfectedCriticalToDead and InfectedCriticalToRecovered can then be calculated using
* the InfectedNoSymptomsToInfectedSymptoms flow with the standard formula from the IDE model.
* The ExposedToInfectedNoSymptoms and InfectedNoSymptomsToInfectedSymptoms flows are calculated 
* using the means of the respective TransitionDistribution. 
* The flow InfectedNoSymptomsToInfectedSymptoms is calculated with the standard formula from the IDE model
* using the results for ExposedToInfectedNoSymptoms.
*
* The number of deaths used in the model is computed by shifting the reported RKI data according to the mean values of the transitions 
* InfectedSymptomsToInfectedSevere, InfectedSevereToInfectedCritical and InfectedCriticalToDead.
* We also set the number of total confirmed cases in the model. 
* Therefore the initialization method using the total confirmed cases is used in the model. See also the documentation of the model class.
* 
* The start date of the model simulation is set to t0=0.
*
* @param[in, out] model The model for which the initial flows should be computed.
* @param[in] dt Time step size.
* @param[in] rki_data Vector containing RKI data.
* @param[in] date The start date of the simulation and the last time point of the TimeSeries used for initialization.
* @param[in] scale_confirmed_cases Vector with factor(s for each age group) by which to scale the confirmed cases of 
*   rki_data to consider unreported cases.
* @tparam EntryType is expected to be ConfirmedCasesNoAgeEntry for data that is not age resolved and 
*   ConfirmedCasesDataEntry for age resolved data. See also epi_data.h.
* @returns Any io errors that happen during reading of the files.
*/

template <typename EntryType>
IOResult<void> set_initial_flows(Model& model, const ScalarType dt, const std::vector<EntryType> rki_data,
                                 const Date date, const CustomIndexArray<ScalarType, AgeGroup> scale_confirmed_cases)
{
    // Check if scale_confirmed_cases has the right size (= number of age groups).
    assert(model.get_num_agegroups() == (size_t)scale_confirmed_cases.size());
    // Check if the correct EntryType was used.
    if constexpr (std::is_same_v<EntryType, ConfirmedCasesDataEntry>) {
        assert(model.get_num_agegroups() == (size_t)EntryType::age_group_names.size());
    }
    else {
        assert(model.get_num_agegroups() == 1);
    }

    //--- Preparations ---
    auto max_date_entry = std::max_element(rki_data.begin(), rki_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == rki_data.end()) {
        log_error("RKI data file is empty.");
        return failure(StatusCode::InvalidFileFormat, "RKI data file is empty.");
    }
    auto max_date = max_date_entry->date;
    if (max_date < date) {
        log_error("Specified date does not exist in RKI data.");
        return failure(StatusCode::OutOfRange, "Specified date does not exist in RKI data.");
    }

    // Get (global) support_max to determine how many flows in the past we have to compute.
    ScalarType global_support_max         = model.get_global_support_max(dt);
    Eigen::Index global_support_max_index = Eigen::Index(std::ceil(global_support_max / dt));

    // Get the number of AgeGroups.
    const size_t num_age_groups = model.get_num_agegroups();

    // m_transitions should be empty at the beginning.
    if (model.m_transitions.get_num_time_points() > 0) {
        model.m_transitions = TimeSeries<ScalarType>(Eigen::Index(InfectionTransition::Count) * num_age_groups);
    }
    if (model.m_populations.get_time(0) != 0) {
        model.m_populations.remove_last_time_point();
        model.m_populations.add_time_point<Eigen::VectorXd>(
            0, TimeSeries<ScalarType>::Vector::Constant((int)InfectionState::Count * num_age_groups, 0));
    }

    // The first time we need is -4 * global_support_max because we need values for
    // InfectedNoSymptomsToInfectedSymptoms on this time window to compute all consecutive transitions on the time
    // window from -global_support_max to 0.
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
        int Di                     = model.get_state_flat_index(Eigen::Index(InfectionState::Dead), group);
        model.m_populations[0][Di] = 0;

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
        model.m_transitions.add_time_point(
            i * dt, TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count * num_age_groups, 0.));
    }

    // Setting the entries in m_total_confirmed_cases to zero before overwriting it with the RKI data.
    model.m_total_confirmed_cases = CustomIndexArray<ScalarType, AgeGroup>(AgeGroup(num_age_groups), 0.);
    //--- Calculate the flow InfectedNoSymptomsToInfectedSymptoms using the RKI data and store in the m_transitions object.---
    ScalarType min_offset_needed = std::ceil(
        model.m_transitions.get_time(0) -
        1); // Need -1 if first time point is integer and just the floor value if not, therefore use ceil and -1
    ScalarType max_offset_needed = std::ceil(model.m_transitions.get_last_time());
    bool min_offset_needed_avail = false;
    bool max_offset_needed_avail = false;

    // Go through the entries of rki_data and check if date is needed for calculation. Confirmed cases are scaled.
    // Define dummy variables to store the first and the last index of the TimeSeries where the considered entry of
    // rki_data is potentially needed.
    Eigen::Index idx_needed_first = 0;
    Eigen::Index idx_needed_last  = 0;
    ScalarType time_idx           = 0;

    for (auto&& entry : rki_data) {
        int offset = get_offset_in_days(entry.date, date);

        // Get the index regarding the age group.
        // If we don't have age resolution and use EntryType=ConfirmedCasesNoAge, the index is set to 1.
        // If we consider multiple age groups and use EntryType=ConfirmedCasesDataEntry, it is determined accordingly.
        AgeGroup group = AgeGroup(0);
        if constexpr (std::is_same_v<EntryType, ConfirmedCasesDataEntry>) {
            group = entry.age_group;
        }

        if ((offset >= min_offset_needed) && (offset <= max_offset_needed)) {
            if (offset == min_offset_needed) {
                min_offset_needed_avail = true;
            }
            if (offset == max_offset_needed) {
                max_offset_needed_avail = true;
            }
            // Smallest index for which the entry is needed.
            idx_needed_first =
                Eigen::Index(std::max(std::floor((offset - model.m_transitions.get_time(0) - 1) / dt), 0.));
            // Biggest index for which the entry is needed.
            idx_needed_last = Eigen::Index(std::min(std::ceil((offset - model.m_transitions.get_time(0) + 1) / dt),
                                                    double(model.m_transitions.get_num_time_points() - 1)));

            int INStISyi = model.get_transition_flat_index(
                Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), group);

            for (Eigen::Index i = idx_needed_first; i <= idx_needed_last; i++) {

                time_idx = model.m_transitions.get_time(i);
                if (offset == int(std::floor(time_idx))) {
                    model.m_transitions[i][INStISyi] +=
                        (1 - (time_idx - std::floor(time_idx))) * scale_confirmed_cases[group] * entry.num_confirmed;
                }
                if (offset == int(std::ceil(time_idx))) {
                    model.m_transitions[i][INStISyi] +=
                        (time_idx - std::floor(time_idx)) * scale_confirmed_cases[group] * entry.num_confirmed;
                }
                if (offset == int(std::floor(time_idx - dt))) {
                    model.m_transitions[i][INStISyi] -= (1 - (time_idx - dt - std::floor(time_idx - dt))) *
                                                        scale_confirmed_cases[group] * entry.num_confirmed;
                }
                if (offset == int(std::ceil(time_idx - dt))) {
                    model.m_transitions[i][INStISyi] -= (time_idx - dt - std::floor(time_idx - dt)) *
                                                        scale_confirmed_cases[group] * entry.num_confirmed;
                }
            }

            // Compute Dead by shifting RKI data according to mean stay times.
            // This is done because the RKI reports death with the date of positive test instead of the date of deaths.
            int Di = model.get_state_flat_index(Eigen::Index(InfectionState::Dead), group);
            if (offset ==
                std::floor(-mean_InfectedSymptomsToInfectedSevere[group] -
                           mean_InfectedSevereToInfectedCritical[group] - mean_InfectedCriticalToDead[group])) {
                model.m_populations[0][Di] +=
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
                model.m_populations[0][Di] +=
                    (-mean_InfectedSymptomsToInfectedSevere[group] - mean_InfectedSevereToInfectedCritical[group] -
                     mean_InfectedCriticalToDead[group] -
                     std::floor(-mean_InfectedSymptomsToInfectedSevere[group] -
                                mean_InfectedSevereToInfectedCritical[group] - mean_InfectedCriticalToDead[group])) *
                    entry.num_deaths;
            }

            if (offset == 0) {
                model.m_total_confirmed_cases[group] = scale_confirmed_cases[group] * entry.num_confirmed;
            }
        }
    }

    if (!max_offset_needed_avail) {
        log_error("Necessary range of dates needed to compute initial values does not exist in RKI data.");
        return failure(StatusCode::OutOfRange, "Necessary range of dates does not exist in RKI data.");
    }
    if (!min_offset_needed_avail) {
        auto min_date_entry = std::min_element(rki_data.begin(), rki_data.end(), [](auto&& a, auto&& b) {
            return a.date < b.date;
        });
        auto min_date       = min_date_entry->date;

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

    //--- Calculate the flows "after" InfectedNoSymptomsToInfectedSymptoms (that were set above using rki_data). ---
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
            model.m_transitions[i + start_shift][EtINSi] =
                (1 / model.parameters.get<TransitionProbabilities>()[group][Eigen::Index(
                         InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)]) *
                model.m_transitions[i + start_shift + index_shift_mean][INStISyi];
        }

        // Compute flow SusceptibleToExposed for -global_support_max, ..., 0.
        // Use mean values of the TransitionDistribution ExposedToInfectedNoSymptoms and of the
        // TransitionDistribution InfectedNoSymptomsToInfectedSymptoms for the calculation.
        index_shift_mean = Eigen::Index(std::round(
            (mean_ExposedToInfectedNoSymptoms[group] + mean_InfectedNoSymptomsToInfectedSymptoms[group]) / dt));
        int StEi = model.get_transition_flat_index(Eigen::Index(InfectionTransition::SusceptibleToExposed), group);

        for (Eigen::Index i = -global_support_max_index; i <= 0; i++) {
            model.m_transitions[i + start_shift][StEi] =
                (1 / model.parameters.get<TransitionProbabilities>()[group][Eigen::Index(
                         InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)]) *
                model.m_transitions[i + start_shift + index_shift_mean][INStISyi];
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
    auto transition_copy(model.m_transitions);
    model.m_transitions = TimeSeries<ScalarType>(Eigen::Index(InfectionTransition::Count) * num_age_groups);
    for (Eigen::Index i = -global_support_max_index; i <= 0; i++) {
        model.m_transitions.add_time_point(i * dt, transition_copy.get_value(i + start_shift));
    }

    return mio::success();
}

} // namespace isecir
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // IDE_INITIALFLOWS_H
