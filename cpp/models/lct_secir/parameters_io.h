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

#ifndef LCTSECIR_PARAMETERS_IO_H
#define LCTSECIR_PARAMETERS_IO_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "lct_secir/model.h"
#include "lct_secir/infection_state.h"
#include "lct_secir/parameters.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"
#include "memilio/utils/date.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/metaprogramming.h"
#include "memilio/math/eigen.h"
#include "memilio/epidemiology/lct_populations.h"

#include <string>
namespace mio
{
namespace lsecir
{

namespace details
{ // Use namespace to hide functions that are not intended to be used outside this file.

template <class EntryType>
using has_age_member_t = decltype(EntryType::age_group_names);
/**
 * @brief Check whether a check_constraints function exists.
 * @tparam The type to check for the existence of the member function.
 */
template <class EntryType>
using is_age_resolved_entry_type = is_expression_valid<has_age_member_t, EntryType>;

/**
* @brief Processes one entry of an RKI data set for the definition of an initial value vector for an LCT model.
*   
* Takes one entry of an RKI data vector and changes the value of an initial value vector accordingly.
* @param[in, out] init The initial value vector (updated according to the data entry).
* @param[in] entry The entry of the RKI data set.
* @param[in] offset The offset between the date of the entry and the date for which the 
*   initial value vector is calculated.
* @param[in] staytimes A vector of the average time spent in each compartment defined in InfectionState,
*    as given by the model (for the group under consideration) for which the initial value vector is calculated.
* @param[in] inv_prob_SymptomsPerNoSymptoms The inverse of the probability InfectedSymptomsPerInfectedNoSymptoms as
*   given by the model (for the group under consideration) for which the initial value vector is calculated.
* @param[in] severePerInfectedSymptoms Probability as given by the model (for the group under consideration)
*    for which the initial value vector is calculated.
* @param[in] criticalPerSevere Probability as given by the model (for the group under consideration)
*    for which the initial value vector is calculated.
* @param[in] scale_confirmed_cases Factor by which to scale the confirmed cases of RKI data to consider unreported cases.
* @tparam LctState The LctState of the group under consideration.
* @tparam EntryType The type of the data entry of the RKI data.
*/
template <class Populations, class EntryType, size_t Group>
void process_entry(Populations& populations, const EntryType& entry, int offset,
                   const std::vector<ScalarType> staytimes, ScalarType inv_prob_SymptomsPerNoSymptoms,
                   ScalarType severePerInfectedSymptoms, ScalarType criticalPerSevere, ScalarType scale_confirmed_cases)
{
    using LctStateGroup      = type_at_index_t<Group, typename Populations::LctStatesGroups>;
    size_t first_index_group = populations.template get_first_index_of_group<Group>();

    // Add confirmed cases at date to compartment Recovered.
    if (offset == 0) {
        populations[first_index_group + LctStateGroup::template get_first_index<InfectionState::Recovered>()] +=
            scale_confirmed_cases * entry.num_confirmed;
    }

    // Compute initial values for compartment InfectedNoSymptoms.
    if (offset >= 0 && offset <= std::ceil(staytimes[(size_t)InfectionState::InfectedNoSymptoms])) {
        // Mean stay time in each subcompartment.
        ScalarType time_INS_per_subcomp =
            staytimes[(size_t)InfectionState::InfectedNoSymptoms] /
            (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>();
        // Index of the last subcompartment of InfectedNoSymptoms.
        size_t last_index_INS =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedNoSymptoms>() +
            LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>() - 1;
        for (int i = 0; i < (int)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>();
             i++) {
            if (offset == std::floor(i * time_INS_per_subcomp)) {
                populations[last_index_INS - i] -=
                    (1 - (i * time_INS_per_subcomp - std::floor(i * time_INS_per_subcomp))) *
                    inv_prob_SymptomsPerNoSymptoms * scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == std::ceil(i * time_INS_per_subcomp)) {
                populations[last_index_INS - i] -= (i * time_INS_per_subcomp - std::floor(i * time_INS_per_subcomp)) *
                                                   inv_prob_SymptomsPerNoSymptoms * scale_confirmed_cases *
                                                   entry.num_confirmed;
            }
            if (offset == std::floor((i + 1) * time_INS_per_subcomp)) {
                populations[last_index_INS - i] +=
                    (1 - ((i + 1) * time_INS_per_subcomp - std::floor((i + 1) * time_INS_per_subcomp))) *
                    inv_prob_SymptomsPerNoSymptoms * scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == std::ceil((i + 1) * time_INS_per_subcomp)) {
                populations[last_index_INS - i] +=
                    ((i + 1) * time_INS_per_subcomp - std::floor((i + 1) * time_INS_per_subcomp)) *
                    inv_prob_SymptomsPerNoSymptoms * scale_confirmed_cases * entry.num_confirmed;
            }
        }
    }

    // Compute initial values for compartment Exposed.
    if (offset >= std::floor(staytimes[(size_t)InfectionState::InfectedNoSymptoms]) &&
        offset <= std::ceil(staytimes[(size_t)InfectionState::Exposed] +
                            staytimes[(size_t)InfectionState::InfectedNoSymptoms])) {
        // Mean stay time in each subcompartment.
        ScalarType time_Exposed_per_subcomp =
            staytimes[(size_t)InfectionState::Exposed] /
            (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::Exposed>();
        // Index of the last subcompartment of Exposed.
        size_t last_index_Exposed = first_index_group +
                                    LctStateGroup::template get_first_index<InfectionState::Exposed>() +
                                    LctStateGroup::template get_num_subcompartments<InfectionState::Exposed>() - 1;
        for (int i = 0; i < (int)LctStateGroup::template get_num_subcompartments<InfectionState::Exposed>(); i++) {
            if (offset ==
                std::floor(staytimes[(size_t)InfectionState::InfectedNoSymptoms] + i * time_Exposed_per_subcomp)) {
                populations[last_index_Exposed - i] -=
                    (1 - (staytimes[(size_t)InfectionState::InfectedNoSymptoms] + i * time_Exposed_per_subcomp -
                          std::floor(staytimes[(size_t)InfectionState::InfectedNoSymptoms] +
                                     i * time_Exposed_per_subcomp))) *
                    inv_prob_SymptomsPerNoSymptoms * scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset ==
                std::ceil(staytimes[(size_t)InfectionState::InfectedNoSymptoms] + i * time_Exposed_per_subcomp)) {
                populations[last_index_Exposed - i] -=
                    (staytimes[(size_t)InfectionState::InfectedNoSymptoms] + i * time_Exposed_per_subcomp -
                     std::floor(staytimes[(size_t)InfectionState::InfectedNoSymptoms] + i * time_Exposed_per_subcomp)) *
                    inv_prob_SymptomsPerNoSymptoms * scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == std::floor(staytimes[(size_t)InfectionState::InfectedNoSymptoms] +
                                     (i + 1) * time_Exposed_per_subcomp)) {
                populations[last_index_Exposed - i] +=
                    (1 - (staytimes[(size_t)InfectionState::InfectedNoSymptoms] + (i + 1) * time_Exposed_per_subcomp -
                          std::floor(staytimes[(size_t)InfectionState::InfectedNoSymptoms] +
                                     (i + 1) * time_Exposed_per_subcomp))) *
                    inv_prob_SymptomsPerNoSymptoms * scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset ==
                std::ceil(staytimes[(size_t)InfectionState::InfectedNoSymptoms] + (i + 1) * time_Exposed_per_subcomp)) {
                populations[last_index_Exposed - i] +=
                    (staytimes[(size_t)InfectionState::InfectedNoSymptoms] + (i + 1) * time_Exposed_per_subcomp -
                     std::floor(staytimes[(size_t)InfectionState::InfectedNoSymptoms] +
                                (i + 1) * time_Exposed_per_subcomp)) *
                    inv_prob_SymptomsPerNoSymptoms * scale_confirmed_cases * entry.num_confirmed;
            }
        }
    }

    // Compute initial values for compartment InfectedSymptoms.
    if (offset >= std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms]) && offset <= 0) {
        // Mean stay time in each subcompartment.
        ScalarType time_IS_per_subcomp =
            staytimes[(size_t)InfectionState::InfectedSymptoms] /
            (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms>();
        // Index of the first subcompartment of InfectedSymptoms.
        size_t first_index_IS =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedSymptoms>();
        for (int i = 0; i < (int)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms>();
             i++) {
            if (offset == std::floor(-time_IS_per_subcomp * (i + 1))) {
                populations[first_index_IS + i] -=
                    (1 - (-(i + 1.) * time_IS_per_subcomp - std::floor(-(i + 1) * time_IS_per_subcomp))) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == std::ceil(-time_IS_per_subcomp * (i + 1))) {
                populations[first_index_IS + i] -=
                    (-(i + 1) * time_IS_per_subcomp - std::floor(-(i + 1) * time_IS_per_subcomp)) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == std::floor(-time_IS_per_subcomp * i)) {
                populations[first_index_IS + i] +=
                    (1 - (-i * time_IS_per_subcomp - std::floor(-i * time_IS_per_subcomp))) * scale_confirmed_cases *
                    entry.num_confirmed;
            }
            if (offset == std::ceil(-time_IS_per_subcomp * i)) {
                populations[first_index_IS + i] += (-i * time_IS_per_subcomp - std::floor(-i * time_IS_per_subcomp)) *
                                                   scale_confirmed_cases * entry.num_confirmed;
            }
        }
    }

    // Compute initial values for compartment InfectedSevere.
    if (offset >= std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                             staytimes[(size_t)InfectionState::InfectedSevere]) &&
        offset <= std::ceil(-staytimes[(size_t)InfectionState::InfectedSymptoms])) {
        // Mean stay time in each subcompartment.
        ScalarType time_ISev_per_subcomp =
            staytimes[(size_t)InfectionState::InfectedSevere] /
            (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere>();
        // Index of the first subcompartment of InfectedSevere.
        size_t first_index_ISev =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedSevere>();
        for (int i = 0; i < (int)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere>();
             i++) {
            if (offset ==
                std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] - time_ISev_per_subcomp * (i + 1))) {
                populations[first_index_ISev + i] -=
                    severePerInfectedSymptoms *
                    (1 - (-staytimes[(size_t)InfectionState::InfectedSymptoms] - (i + 1) * time_ISev_per_subcomp -
                          std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                     (i + 1) * time_ISev_per_subcomp))) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset ==
                std::ceil(-staytimes[(size_t)InfectionState::InfectedSymptoms] - time_ISev_per_subcomp * (i + 1))) {
                populations[first_index_ISev + i] -=
                    severePerInfectedSymptoms *
                    (-staytimes[(size_t)InfectionState::InfectedSymptoms] - (i + 1) * time_ISev_per_subcomp -
                     std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                (i + 1) * time_ISev_per_subcomp)) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset ==
                std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] - time_ISev_per_subcomp * i)) {
                populations[first_index_ISev + i] +=
                    severePerInfectedSymptoms *
                    (1 -
                     (-staytimes[(size_t)InfectionState::InfectedSymptoms] - i * time_ISev_per_subcomp -
                      std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] - i * time_ISev_per_subcomp))) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == std::ceil(-staytimes[(size_t)InfectionState::InfectedSymptoms] - time_ISev_per_subcomp * i)) {
                populations[first_index_ISev + i] +=
                    severePerInfectedSymptoms *
                    (-staytimes[(size_t)InfectionState::InfectedSymptoms] - i * time_ISev_per_subcomp -
                     std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] - i * time_ISev_per_subcomp)) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
        }
    }

    // Compute initial values for compartment InfectedCritical.
    if (offset >= std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                             staytimes[(size_t)InfectionState::InfectedSevere] -
                             staytimes[(size_t)InfectionState::InfectedCritical]) &&
        offset <= std::ceil(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                            staytimes[(size_t)InfectionState::InfectedSevere])) {
        // Mean stay time in each subcompartment.
        ScalarType time_ICri_per_subcomp =
            staytimes[(size_t)InfectionState::InfectedCritical] /
            (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical>();
        // Index of the first subcompartment of InfectedCritical.
        size_t first_index_ICri =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedCritical>();
        for (int i = 0; i < (int)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical>();
             i++) {
            if (offset ==
                std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                           staytimes[(size_t)InfectionState::InfectedSevere] - time_ICri_per_subcomp * (i + 1))) {
                populations[first_index_ICri + i] -=
                    severePerInfectedSymptoms * criticalPerSevere *
                    (1 - (-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                          staytimes[(size_t)InfectionState::InfectedSevere] - (i + 1) * time_ICri_per_subcomp -
                          std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                     staytimes[(size_t)InfectionState::InfectedSevere] -
                                     (i + 1) * time_ICri_per_subcomp))) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset ==
                std::ceil(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                          staytimes[(size_t)InfectionState::InfectedSevere] - time_ICri_per_subcomp * (i + 1))) {
                populations[first_index_ICri + i] -=
                    severePerInfectedSymptoms * criticalPerSevere *
                    (-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                     staytimes[(size_t)InfectionState::InfectedSevere] - (i + 1) * time_ICri_per_subcomp -
                     std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                staytimes[(size_t)InfectionState::InfectedSevere] - (i + 1) * time_ICri_per_subcomp)) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                     staytimes[(size_t)InfectionState::InfectedSevere] - time_ICri_per_subcomp * i)) {
                populations[first_index_ICri + i] +=
                    severePerInfectedSymptoms * criticalPerSevere *
                    (1 - (-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                          staytimes[(size_t)InfectionState::InfectedSevere] - i * time_ICri_per_subcomp -
                          std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                     staytimes[(size_t)InfectionState::InfectedSevere] - i * time_ICri_per_subcomp))) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == std::ceil(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                    staytimes[(size_t)InfectionState::InfectedSevere] - time_ICri_per_subcomp * i)) {
                populations[first_index_ICri + i] +=
                    severePerInfectedSymptoms * criticalPerSevere *
                    (-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                     staytimes[(size_t)InfectionState::InfectedSevere] - i * time_ICri_per_subcomp -
                     std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                staytimes[(size_t)InfectionState::InfectedSevere] - i * time_ICri_per_subcomp)) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
        }
    }

    // Compute Dead.
    if (offset == std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                             staytimes[(size_t)InfectionState::InfectedSevere] -
                             staytimes[(size_t)InfectionState::InfectedCritical])) {
        populations[first_index_group + LctStateGroup::template get_first_index<InfectionState::Dead>()] +=
            (1 -
             (-staytimes[(size_t)InfectionState::InfectedSymptoms] - staytimes[(size_t)InfectionState::InfectedSevere] -
              staytimes[(size_t)InfectionState::InfectedCritical] -
              std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                         staytimes[(size_t)InfectionState::InfectedSevere] -
                         staytimes[(size_t)InfectionState::InfectedCritical]))) *
            entry.num_deaths;
    }
    if (offset == std::ceil(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                            staytimes[(size_t)InfectionState::InfectedSevere] -
                            staytimes[(size_t)InfectionState::InfectedCritical])) {
        populations[first_index_group + LctStateGroup::template get_first_index<InfectionState::Dead>()] +=
            (-staytimes[(size_t)InfectionState::InfectedSymptoms] - staytimes[(size_t)InfectionState::InfectedSevere] -
             staytimes[(size_t)InfectionState::InfectedCritical] -
             std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                        staytimes[(size_t)InfectionState::InfectedSevere] -
                        staytimes[(size_t)InfectionState::InfectedCritical])) *
            entry.num_deaths;
    }
}

/**
* @brief Computes an initialization vector for an age-resolved LCT model with case data from RKI with age groups.
*   
* Please use the set_initial_data_from_confirmed_cases() function, which calls this function automatically!
* Calculates an initial value vector for an LCT model and updates the initial value vector in the model.
* The number of groups in the LCT model has to be the same as the number of RKI age groups. 
* For the computation, expected stay times in the subcompartments are used. To calculate the initial values, 
* we assume for simplicity that individuals stay in the subcompartment for exactly the expected time.
* The RKI data are linearly interpolated within one day to match the expected stay time in a subcompartment.
* The RKI data should contain data for each required day with division in age groups, the completeness 
* of the dates is not verified.
* Data can be downloaded e.g. with the file pycode/memilio-epidata/memilio/epidata/getCaseData.py, 
* which creates a file named cases_all_age.json or a similar name. 
* One should set impute_dates=True so that missing dates are imputed.
*
* @param[in, out] model The model for which the inital data should be computed and set.
* @param[in] rki_data Vector with the RKI data.
* @param[in] date Date for which the initial values should be computed. date is the start date of the simulation.
* @param[in] total_population Total size of the population for every age group. 
* @param[in] scale_confirmed_cases Factor for each age group by which to scale the confirmed cases of the rki data to consider 
*       unreported cases.
* @tparam Model is expected to be an LCT-SECIR model defined in models/lct_secir/model.h.
* @tparam Group The age group for which the initial values should be calculated. The function is called recursively 
*   such that the initial values are calculated for every age group if Group is zero at the beginning.
* @returns Any io errors that happen during reading of the files and the calculation.
*/
template <class Populations, class EntryType, size_t Group = 0>
IOResult<void> set_initial_data_from_confirmed_cases(Populations& populations, const std::vector<EntryType>& rki_data,
                                                     const Parameters& parameters, Date date,
                                                     const std::vector<ScalarType>& total_population,
                                                     const std::vector<ScalarType>& scale_confirmed_cases)
{
    static_assert((Group < Populations::num_groups) && (Group >= 0), "The template parameter Group should be valid.");
    using LctStateGroup      = type_at_index_t<Group, typename Populations::LctStatesGroups>;
    size_t first_index_group = populations.template get_first_index_of_group<Group>();

    // Define variables for parameters.
    std::vector<ScalarType> staytimes((size_t)InfectionState::Count, -1.);
    staytimes[(size_t)InfectionState::Exposed]            = parameters.template get<TimeExposed>()[Group];
    staytimes[(size_t)InfectionState::InfectedNoSymptoms] = parameters.template get<TimeInfectedNoSymptoms>()[Group];
    staytimes[(size_t)InfectionState::InfectedSymptoms]   = parameters.template get<TimeInfectedSymptoms>()[Group];
    staytimes[(size_t)InfectionState::InfectedSevere]     = parameters.template get<TimeInfectedSevere>()[Group];
    staytimes[(size_t)InfectionState::InfectedCritical]   = parameters.template get<TimeInfectedCritical>()[Group];
    ScalarType inv_prob_SymptomsPerNoSymptoms =
        1 / (1 - parameters.template get<RecoveredPerInfectedNoSymptoms>()[Group]);
    ScalarType prob_SeverePerInfectedSymptoms = parameters.template get<SeverePerInfectedSymptoms>()[Group];
    ScalarType prob_CriticalPerSevere         = parameters.template get<CriticalPerSevere>()[Group];

    ScalarType min_offset_needed = std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                              staytimes[(size_t)InfectionState::InfectedSevere] -
                                              staytimes[(size_t)InfectionState::InfectedCritical]);
    ScalarType max_offset_needed =
        std::ceil(staytimes[(size_t)InfectionState::Exposed] + staytimes[(size_t)InfectionState::InfectedNoSymptoms]);

    bool min_offset_needed_avail = false;
    bool max_offset_needed_avail = false;

    // Go through the entries of rki_data and check if the entry has the age_group Group and if the date is
    // needed for calculation. Confirmed cases are scaled by scale_confirmed_cases.
    for (auto&& entry : rki_data) {
        if constexpr (is_age_resolved_entry_type<EntryType>::value) {
            if ((size_t)entry.age_group == Group) {
                int offset = get_offset_in_days(entry.date, date);
                if ((offset >= min_offset_needed) && (offset <= max_offset_needed)) {
                    if (offset == max_offset_needed) {
                        max_offset_needed_avail = true;
                    }
                    if (offset == min_offset_needed) {
                        min_offset_needed_avail = true;
                    }
                    process_entry<Populations, EntryType, Group>(
                        populations, entry, offset, staytimes, inv_prob_SymptomsPerNoSymptoms,
                        prob_SeverePerInfectedSymptoms, prob_CriticalPerSevere, scale_confirmed_cases[Group]);
                }
            }
        }
        else {
            int offset = get_offset_in_days(entry.date, date);
            if ((offset >= min_offset_needed) && (offset <= max_offset_needed)) {
                if (offset == max_offset_needed) {
                    max_offset_needed_avail = true;
                }
                if (offset == min_offset_needed) {
                    min_offset_needed_avail = true;
                }
                process_entry<Populations, EntryType, Group>(
                    populations, entry, offset, staytimes, inv_prob_SymptomsPerNoSymptoms,
                    prob_SeverePerInfectedSymptoms, prob_CriticalPerSevere, scale_confirmed_cases[Group]);
            }
        }
    }

    // Compute Recovered.
    populations[first_index_group + LctStateGroup::template get_first_index<InfectionState::Recovered>()] -=
        populations.get_compartments()
            .segment(first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedSymptoms>(),
                     LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms>() +
                         LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere>() +
                         LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical>())
            .sum();
    populations[first_index_group + LctStateGroup::template get_first_index<InfectionState::Recovered>()] -=
        populations[first_index_group + LctStateGroup::template get_first_index<InfectionState::Dead>()];

    // Compute Susceptibles.
    populations[first_index_group + LctStateGroup::template get_first_index<InfectionState::Susceptible>()] =
        total_population[Group] -
        populations.get_compartments().segment(first_index_group + 1, LctStateGroup::Count - 1).sum();

    // Check whether all necessary dates are available.
    if (!max_offset_needed_avail || !min_offset_needed_avail) {
        log_error(
            "Necessary range of dates needed to compute initial values does not exist in RKI data for group {:d}.",
            Group);
        return failure(StatusCode::OutOfRange,
                       "Necessary range of dates needed to compute initial values does not exist in RKI data.");
    }

    if constexpr (Group + 1 < Populations::num_groups) {
        return set_initial_data_from_confirmed_cases<Populations, EntryType, Group + 1>(
            populations, rki_data, parameters, date, total_population, scale_confirmed_cases);
    }
    else {
        return mio::success();
    }
}

} // namespace details

/**
* @brief Computes an initialization vector for an LCT model with case data from RKI.
*   
* Use just one group in the model definition to not divide between age groups.
* Otherwise, the number of groups has to match the number of RKI age groups.
* The function calculates an initial value vector for an LCT model and updates the initial value vector in the model.
* For the computation expected stay times in the subcompartments are used. To calculate the initial values, 
* we assume for simplicity that individuals stay in the subcompartment for exactly the expected time.
* The RKI data are linearly interpolated within one day to match the expected stay time in a subcompartment.
* The RKI data should contain data for each needed day with or without division of age groups, 
*   the completeness of the dates is not verified.
* Data can be downloaded e.g. with the file pycode/memilio-epidata/memilio/epidata/getCaseData.py, which creates files
* named e.g. cases_all_germany.json for no groups or cases_all_age.json with division in age groups or similar names.
* One should set impute_dates=True so that missing dates are imputed.
* The data and the number of entries in the total_population and scale_confirmed_cases vectors have to match the 
*   number of groups used in Model.
*
* @param[in, out] model The model for which the inital data should be computed and set.
* @param[in] path Path to the RKI data file.
* @param[in] date Date for which the initial values should be computed. date is the start date of the simulation.
* @param[in] total_population Total size of the population of Germany or of every age group. 
* @param[in] scale_confirmed_cases Factor(s for each age group) by which to scale the confirmed cases of the rki data 
*   to consider unreported cases.
* @tparam Model is expected to be an LCT-SECIR model defined in models/lct_secir/model.h.
* @returns Any io errors that happen during reading of the files.
*/
template <class Populations, class EntryType>
IOResult<void> set_initial_data_from_confirmed_cases(const std::vector<EntryType>& rki_data, Populations& populations,
                                                     const Parameters& parameters, const Date date,
                                                     const std::vector<ScalarType>& total_population,
                                                     const std::vector<ScalarType>& scale_confirmed_cases)
{
    assert(total_population.size() == Populations::num_groups);
    assert(scale_confirmed_cases.size() == Populations::num_groups);
    if constexpr (Populations::num_groups > 1) {
        static_assert(details::is_age_resolved_entry_type<EntryType>::value);
        assert(ConfirmedCasesDataEntry::age_group_names.size() == Populations::num_groups);
    }
    auto max_date_entry = std::max_element(rki_data.begin(), rki_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == rki_data.end()) {
        log_error("RKI data file is empty.");
        return failure(StatusCode::InvalidFileFormat, "RKI data is empty.");
    }
    auto max_date = max_date_entry->date;
    if (max_date < date) {
        log_error("Specified date does not exist in RKI data.");
        return failure(StatusCode::OutOfRange, "Specified date does not exist in RKI data.");
    }
    return details::set_initial_data_from_confirmed_cases<Populations, EntryType>(
        populations, rki_data, parameters, date, total_population, scale_confirmed_cases);
}

} // namespace lsecir
} // namespace mio
#endif // MEMILIO_HAS_JSONCPP

#endif // LCTSECIR_PARAMETERS_IO_H
