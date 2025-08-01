/*
* Copyright (C) 2020-2025 MEmilio
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
#include "memilio/math/floating_point.h"
#include "memilio/epidemiology/lct_populations.h"

#include <string>
#include <vector>
#include <type_traits>
namespace mio
{
namespace lsecir
{

namespace details
{ // Use namespace to hide functions that are not intended to be used outside this file.

/**
* @brief Processes one entry of an RKI data set for the definition of an initial value vector for an LCT population.
*
* Takes one entry of an RKI data vector and changes the value in populations accordingly.
* This function provides sub-functionality of the set_initial_values_from_confirmed_cases() function.
*
* @param[out] populations The populations for which the inital data should be computed and set.
* @param[in] entry The entry of the RKI data set.
* @param[in] offset The offset between the date of the entry and the date for which the
*   initial value vector is calculated.
* @param[in] staytimes A vector of the average time spent in each compartment defined in InfectionState
*    (for the group under consideration) for which the initial value vector is calculated.
* @param[in] inv_prob_SymptomsPerNoSymptoms The inverse of the probability InfectedSymptomsPerInfectedNoSymptoms
*    (for the group under consideration) for which the initial value vector is calculated.
* @param[in] severePerInfectedSymptoms Probability (for the group under consideration)
*    for which the initial value vector is calculated.
* @param[in] criticalPerSevere Probability (for the group under consideration)
*    for which the initial value vector is calculated.
* @param[in] scale_confirmed_cases Factor by which to scale the confirmed cases of RKI data to consider unreported cases.
* @tparam Populations is expected to be an LctPopulations defined in epidemiology/lct_populations.
*   This defines the number of age groups and the numbers of subcompartments.
* @tparam EntryType The type of the data entry of the RKI data.
* @tparam Group The age group of the entry the should be processed.
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

    // Compute Dead by shifting RKI data according to mean stay times.
    // This is because the RKI reports deaths with the date of the positive test, not the date of death.
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
* @brief Computes an initialization vector for an LCT population with case data from RKI recursively for each age group
*    (or for one age group in the case without age resolution).
*
* Please use the set_initial_values_from_reported_data() function, which calls this function automatically!
* This function calculates a segment referring to the defined age group of the initial value vector with
*   subcompartments using the rki_data and the parameters.
* The values for the whole initial value vector stored in populations are calculated recursively.
*
* @param[in] rki_data Vector with the RKI data.
* @param[out] populations The populations for which the inital data should be computed and set.
* @param[in] parameters The parameters that should be used to calculate the initial values.
*   Probabilities and mean stay times are used.
* @param[in] date Date for which the initial values should be computed. date is the start date of the simulation.
* @param[in] total_population Total size of the population of Germany or of every age group.
* @param[in] scale_confirmed_cases Factor(s for each age group) by which to scale the confirmed cases of the rki data
*   to consider unreported cases.
* @tparam Populations is expected to be an LctPopulations defined in epidemiology/lct_populations.
*   This defines the number of age groups and the numbers of subcompartments.
* @tparam EntryType is expected to be ConfirmedCasesNoAgeEntry for data that is not age resolved and
*   ConfirmedCasesDataEntry for age resolved data. See also epi_data.h.
* @tparam Group The age group for which the initial values should be calculated. The function is called recursively
*   such that the initial values are calculated for every age group if Group is zero at the beginning.
* @returns Any io errors that happen during data processing.
*/
template <class Populations, class EntryType, size_t Group = 0>
IOResult<void> set_initial_values_from_confirmed_cases(Populations& populations, const std::vector<EntryType>& rki_data,
                                                       const Parameters<ScalarType>& parameters, const Date date,
                                                       const std::vector<ScalarType>& total_population,
                                                       const std::vector<ScalarType>& scale_confirmed_cases)
{
    static_assert((Group < Populations::num_groups) && (Group >= 0), "The template parameter Group should be valid.");
    using LctStateGroup      = type_at_index_t<Group, typename Populations::LctStatesGroups>;
    size_t first_index_group = populations.template get_first_index_of_group<Group>();

    // Define variables for parameters.
    std::vector<ScalarType> staytimes((size_t)InfectionState::Count, -1.);
    staytimes[(size_t)InfectionState::Exposed] = parameters.template get<TimeExposed<ScalarType>>()[Group];
    staytimes[(size_t)InfectionState::InfectedNoSymptoms] =
        parameters.template get<TimeInfectedNoSymptoms<ScalarType>>()[Group];
    staytimes[(size_t)InfectionState::InfectedSymptoms] =
        parameters.template get<TimeInfectedSymptoms<ScalarType>>()[Group];
    staytimes[(size_t)InfectionState::InfectedSevere] =
        parameters.template get<TimeInfectedSevere<ScalarType>>()[Group];
    staytimes[(size_t)InfectionState::InfectedCritical] =
        parameters.template get<TimeInfectedCritical<ScalarType>>()[Group];
    ScalarType inv_prob_SymptomsPerNoSymptoms =
        1 / (1 - parameters.template get<RecoveredPerInfectedNoSymptoms<ScalarType>>()[Group]);
    ScalarType prob_SeverePerInfectedSymptoms = parameters.template get<SeverePerInfectedSymptoms<ScalarType>>()[Group];
    ScalarType prob_CriticalPerSevere         = parameters.template get<CriticalPerSevere<ScalarType>>()[Group];

    ScalarType min_offset_needed = std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                              staytimes[(size_t)InfectionState::InfectedSevere] -
                                              staytimes[(size_t)InfectionState::InfectedCritical]);
    ScalarType max_offset_needed =
        std::ceil(staytimes[(size_t)InfectionState::Exposed] + staytimes[(size_t)InfectionState::InfectedNoSymptoms]);

    bool min_offset_needed_avail = false;
    bool max_offset_needed_avail = false;

    // Go through the entries of rki_data and check if the entry is age resolved and is referring to
    // the age_group Group in the case with age resolution. If the date is
    // needed for calculation, another function to handle the entry is called. Confirmed cases are scaled by scale_confirmed_cases.
    for (auto&& entry : rki_data) {
        if constexpr (std::is_same_v<EntryType, ConfirmedCasesDataEntry>) {
            if (!((size_t)entry.age_group == Group)) {
                continue;
            }
        }
        int offset = get_offset_in_days(entry.date, date);
        if ((offset >= min_offset_needed) && (offset <= max_offset_needed)) {
            if (offset == max_offset_needed) {
                max_offset_needed_avail = true;
            }
            if (offset == min_offset_needed) {
                min_offset_needed_avail = true;
            }
            process_entry<Populations, EntryType, Group>(populations, entry, offset, staytimes,
                                                         inv_prob_SymptomsPerNoSymptoms, prob_SeverePerInfectedSymptoms,
                                                         prob_CriticalPerSevere, scale_confirmed_cases[Group]);
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

    // Check if all values for populations are valid.
    for (size_t i = first_index_group; i < LctStateGroup::Count; i++) {
        if (floating_point_less<ScalarType>((ScalarType)populations[i], 0., Limits<ScalarType>::zero_tolerance())) {
            log_error("Something went wrong in the initialization of group {:d}. At least one entry is negative.",
                      Group);
            return failure(StatusCode::InvalidValue,
                           "Something went wrong in the initialization as at least one entry is negative.");
        }
    }

    if constexpr (Group + 1 < Populations::num_groups) {
        return set_initial_values_from_confirmed_cases<Populations, EntryType, Group + 1>(
            populations, rki_data, parameters, date, total_population, scale_confirmed_cases);
    }
    else {
        return success();
    }
}

/**
* @brief Computes the total number of patients in Intensive Care Units (in all groups and subcompartments)
*       in the provided Population.
*
* This function calculates the total number of individuals within the compartment InfectedCritical
* for the Populations-data provided, irrespective of their subcompartment or group.
* This total number can be used to scale the entries so that the total number in InfectedCritical is equal to
* the number of ICU patients reported in the DIVI data.
* Please use the set_initial_values_from_reported_data() function, which calls this function automatically!
*
* @param[in] populations The populations for which the total number in InfectedCritical should be computed.
* @tparam Populations is expected to be an LctPopulations defined in epidemiology/lct_populations.
*   This defines the number of age groups and the numbers of subcompartments.
* @tparam Group The age group for which the total number should be calculated. The function is called recursively
*   such that the total number in InfectedCritical within all groups is calculated if Group is zero at the beginning.
* @returns The total number of patients in Intensive Care Units (in all groups and subcompartments).
*/
template <class Populations, size_t Group = 0>
ScalarType get_total_InfectedCritical_from_populations(const Populations& populations)
{
    using LctStateGroup      = type_at_index_t<Group, typename Populations::LctStatesGroups>;
    size_t first_index_group = populations.template get_first_index_of_group<Group>();

    ScalarType infectedCritical_Group =
        populations.get_compartments()
            .segment(first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedCritical>(),
                     LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical>())
            .sum();

    if constexpr (Group + 1 < Populations::num_groups) {
        return infectedCritical_Group +
               get_total_InfectedCritical_from_populations<Populations, Group + 1>(populations);
    }
    else {
        return infectedCritical_Group;
    }
}

/**
* @brief Extract the reported number of patients in ICU for a specific date from DIVI data.
*
* @param[in] divi_data Vector with reported DIVI data.
* @param[in] date Date for which the reported number of patients in ICU should be extracted.
* @returns The reported number of patients in ICU or any io errors that happen during data processing.
*/
IOResult<ScalarType> get_icu_from_divi_data(const std::vector<DiviEntry>& divi_data, const Date date)
{
    for (auto&& entry : divi_data) {
        int offset = get_offset_in_days(entry.date, date);
        if (offset == 0) {
            return success(entry.num_icu);
        }
    }
    log_error("Specified date does not exist in DIVI data.");
    return failure(StatusCode::OutOfRange, "Specified date does not exist in DIVI data.");
}

/**
* @brief Rescales the entries for InfectedCritical in populations such that the total number
*   equals the reported number.
*
* This function rescales the entries for InfectedCritical in the given population for every group and subcompartment
* such that the total number in all InfectedCritical compartments equals the reported number infectedCritical_reported.
*
* If the total number of individuals in InfectedCritical in populations is zero and the reported number is not,
* the reported number is distributed uniformly across the groups.
* Within the groups, the number is distributed uniformly to the subcompartments.
* Note that especially the uniform distribution across groups is not necessarily realistic,
* because the need for intensive care can differ by group.
*
* @param[in,out] populations The populations for which the entries of the InfectedCritical compartments are rescaled.
* @param[in] infectedCritical_reported The reported number for patients in ICU. The total number of individuals in the
*   InfectedCritical compartment in populations will be equal to this value afterward.
*   You can calculate this value with the get_icu_from_divi_data() function.
* @param[in] infectedCritical_populations The current total number of individuals in the InfectedCritical compartment
*   in populations. You can calculate this value with the get_total_InfectedCritical_from_populations() function.
* @tparam Populations is expected to be an LctPopulations defined in epidemiology/lct_populations.
*   This defines the number of age groups and the numbers of subcompartments.
* @tparam Group The age group for which the entries of InfectedCritical should be scaled.
*   The function is called recursively for the groups. The total number in the InfectedCritical compartments is only
*   equal to infectedCritical_reported after the function call if Group is set to zero in the beginning.
* @returns Any io errors that happen during the scaling.
*/
template <class Populations, size_t Group = 0>
IOResult<void> rescale_to_divi_data(Populations& populations, const ScalarType infectedCritical_reported,
                                    const ScalarType infectedCritical_populations)
{
    if (floating_point_less<ScalarType>(infectedCritical_reported, 0., Limits<ScalarType>::zero_tolerance())) {
        log_error("The provided reported number of InfectedCritical is negative. Please check the data.");
        return failure(StatusCode::InvalidValue,
                       "The provided reported number of InfectedCritical is negative. Please check the data.");
    }

    using LctStateGroup      = type_at_index_t<Group, typename Populations::LctStatesGroups>;
    size_t first_index_group = populations.template get_first_index_of_group<Group>();

    if (floating_point_equal<ScalarType>(infectedCritical_populations, 0., Limits<ScalarType>::zero_tolerance())) {
        if (!(floating_point_equal<ScalarType>(infectedCritical_reported, 0., Limits<ScalarType>::zero_tolerance()))) {
            log_info("The calculated number of patients in intensive care is zero, although the reported number is "
                     "not. The reported number is uniformly distributed across groups and subcompartments. Note "
                     "that this is not necessarily realistic.");
            size_t num_InfectedCritical =
                LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical>();
            ScalarType num_age_groups = (ScalarType)Populations::num_groups;
            // Distribute reported number uniformly to age groups and subcompartments.
            for (size_t subcompartment = 0;
                 subcompartment < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical>();
                 subcompartment++) {
                populations[first_index_group +
                            LctStateGroup::template get_first_index<InfectionState::InfectedCritical>() +
                            subcompartment] =
                    infectedCritical_reported / ((ScalarType)num_InfectedCritical * num_age_groups);
            }
            // Adjust Recovered compartment.
            populations[first_index_group + LctStateGroup::template get_first_index<InfectionState::Recovered>()] -=
                infectedCritical_reported / num_age_groups;
            // Number of Susceptibles is not affected because Recovered is adjusted accordingly.
        }
    }
    else {
        // Adjust number of Recovered by adding the old number in InfectedCritical
        // and subtracting the new number (= scaling_factor * old number).
        ScalarType scaling_factor_infectedCritical = infectedCritical_reported / infectedCritical_populations;
        populations[first_index_group + LctStateGroup::template get_first_index<InfectionState::Recovered>()] +=
            (1 - scaling_factor_infectedCritical) *
            populations.get_compartments()
                .segment(first_index_group +
                             LctStateGroup::template get_first_index<InfectionState::InfectedCritical>(),
                         LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical>())
                .sum();
        // Adjust InfectedCritical.
        for (size_t subcompartment = 0;
             subcompartment < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical>();
             subcompartment++) {
            populations[first_index_group +
                        LctStateGroup::template get_first_index<InfectionState::InfectedCritical>() + subcompartment] *=
                scaling_factor_infectedCritical;
        }
        // Number of Susceptibles is not affected because Recovered is adjusted accordingly.
    }
    if (floating_point_less<ScalarType>(
            (ScalarType)
                populations[first_index_group + LctStateGroup::template get_first_index<InfectionState::Recovered>()],
            0., Limits<ScalarType>::zero_tolerance())) {
        log_error(
            "Scaling with reported DIVI data led to a negative entry in the Recovered compartment for group {:d}.",
            Group);
        return failure(StatusCode::InvalidValue,
                       "Scaling with reported DIVI data led to a negative entry in a Recovered compartment.");
    }
    if constexpr (Group + 1 < Populations::num_groups) {
        return rescale_to_divi_data<Populations, Group + 1>(populations, infectedCritical_reported,
                                                            infectedCritical_populations);
    }
    else {
        return success();
    }
}
} // namespace details

/**
* @brief Computes an initialization vector for an LCT population with case data from RKI (and possibly DIVI data).
*
* Use just one group in the definition of the populations to not divide between age groups.
* Otherwise, the number of groups has to match the number of RKI age groups.
* The function calculates an initial value vector referring to an LCT population and updates the initial value vector
* in the populations class.
* For the computation expected stay times in the subcompartments defined in the parameters variable are used.
* To calculate the initial values, we assume for simplicity that individuals stay in the subcompartment
* for exactly the expected time.
* The RKI data are linearly interpolated within one day to match the expected stay time in a subcompartment.
* The RKI data should contain data for each needed day with or without division of age groups,
*   the completeness of the dates is not verified.
* Data can be downloaded e.g. with the file pycode/memilio-epidata/memilio/epidata/getCaseData.py, which creates files
* named e.g. cases_all_germany.json for no groups or cases_all_age.json with division in age groups or similar names.
* One should set impute_dates=True so that missing dates are imputed.
* To read the data into a vector, use the functionality from epi_data.h.
* The data and the number of entries in the total_population and scale_confirmed_cases vectors have to match the
*   number of groups used in Populations.
*
* Additionally, one can scale the result from the calculation with the RKI data to match the reported number of
* patients in ICUs. The patient numbers are provided by DIVI and can be downloaded e.g. using
* pycode/memilio-epidata/memilio/epidata/getDIVIData.py (One should also set impute_dates=True so that missing dates
* are imputed.). Again, to read the data into a vector, use the functionality from epi_data.h.
*
* @param[in] rki_data Vector with the RKI data.
* @param[out] populations The populations for which the inital data should be computed and set.
* @param[in] parameters The parameters that should be used to calculate the initial values.
*   Probabilities and mean stay times are used.
* @param[in] date Date for which the initial values should be computed. date is the start date of the simulation.
* @param[in] total_population Total size of the population of Germany or of every age group.
* @param[in] scale_confirmed_cases Factor(s for each age group) by which to scale the confirmed cases of the rki data
*   to consider unreported cases.
* @param[in] divi_data Vector with DIVI data used to scale the number of individuals in the InfectedCritical
*   compartments in populations so that the total number match the reported number.
*   For the default value (an empty vector), the calculated populations using the RKI data is not scaled.
* @tparam Populations is expected to be an LctPopulations defined in epidemiology/lct_populations.
*   This defines the number of age groups and the numbers of subcompartments.
* @tparam EntryType is expected to be ConfirmedCasesNoAgeEntry for data that is not age resolved and
*   ConfirmedCasesDataEntry for age resolved data. See also epi_data.h.
* @returns Any io errors that happen during data processing.
*/
template <class Populations, class EntryType>
IOResult<void> set_initial_values_from_reported_data(const std::vector<EntryType>& rki_data, Populations& populations,
                                                     const Parameters<ScalarType>& parameters, const Date date,
                                                     const std::vector<ScalarType>& total_population,
                                                     const std::vector<ScalarType>& scale_confirmed_cases,
                                                     const std::vector<DiviEntry>& divi_data = std::vector<DiviEntry>())
{ // Check if the inputs are matching.
    assert(total_population.size() == Populations::num_groups);
    assert(scale_confirmed_cases.size() == Populations::num_groups);
    if constexpr (Populations::num_groups > 1) {
        static_assert(std::is_same_v<EntryType, ConfirmedCasesDataEntry>);
        assert(ConfirmedCasesDataEntry::age_group_names.size() == Populations::num_groups);
    }
    else {
        static_assert(std::is_same_v<EntryType, ConfirmedCasesNoAgeEntry>);
    }
    // Check if RKI data vector is valid.
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
    // Initially set populations to zero.
    for (size_t i = 0; i < populations.get_num_compartments(); i++) {
        populations[i] = 0;
    }
    // Set populations using the RKI data.
    IOResult<void> ioresult_confirmedcases = details::set_initial_values_from_confirmed_cases<Populations, EntryType>(
        populations, rki_data, parameters, date, total_population, scale_confirmed_cases);
    if (!(ioresult_confirmedcases)) {
        return ioresult_confirmedcases;
    }

    // Check if DIVI data is provided and scale the result in populations accordingly.
    if (!divi_data.empty()) {
        ScalarType infectedCritical_populations =
            details::get_total_InfectedCritical_from_populations<Populations>(populations);
        auto infectedCritical_reported = details::get_icu_from_divi_data(divi_data, date);
        if (!(infectedCritical_reported)) {
            return infectedCritical_reported.error();
        }
        return details::rescale_to_divi_data<Populations>(populations, infectedCritical_reported.value(),
                                                          infectedCritical_populations);
    }

    return success();
}

} // namespace lsecir
} // namespace mio
#endif // MEMILIO_HAS_JSONCPP

#endif // LCTSECIR_PARAMETERS_IO_H
