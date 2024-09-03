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
#include "ode_secir/infection_state.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "lct_secir/model.h"
#include "lct_secir/infection_state.h"
#include "lct_secir/parameters.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"
#include "memilio/utils/date.h"
#include "memilio/utils/logging.h"
#include "memilio/math/eigen.h"

#include <string>
namespace mio
{
namespace lsecir
{

namespace details
{ // Use namespace to hide functions that are not intended to be used outside this file.

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
* @param[in] scale_confirmed_cases Factor by which to scale the confirmed cases of rki data to consider unreported cases.
* @tparam LctState The LctState of the group under consideration.
* @tparam EntryType The type of the data entry of the RKI data.
*/
template <class LctState, class EntryType>
void process_entry(Eigen::VectorXd& init, const EntryType& entry, int offset, const std::vector<ScalarType> staytimes,
                   ScalarType inv_prob_SymptomsPerNoSymptoms, ScalarType severePerInfectedSymptoms,
                   ScalarType criticalPerSevere, ScalarType scale_confirmed_cases)
{

    // Add confirmed cases at date to compartment Recovered.
    if (offset == 0) {
        init[LctState::template get_first_index<InfectionState::Recovered>()] +=
            scale_confirmed_cases * entry.num_confirmed;
    }

    // Compute initial values for compartment InfectedNoSymptoms.
    if (offset >= 0 && offset <= std::ceil(staytimes[(size_t)InfectionState::InfectedNoSymptoms])) {
        // Mean stay time in each subcompartment.
        ScalarType time_INS_per_subcomp =
            staytimes[(size_t)InfectionState::InfectedNoSymptoms] /
            (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>();
        // Index of the last subcompartment of InfectedNoSymptoms.
        size_t last_index_INS = LctState::template get_first_index<InfectionState::InfectedNoSymptoms>() +
                                LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>() - 1;
        for (int i = 0; i < (int)LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>();
             i++) {
            if (offset == std::floor(i * time_INS_per_subcomp)) {
                init[last_index_INS - i] -= (1 - (i * time_INS_per_subcomp - std::floor(i * time_INS_per_subcomp))) *
                                            inv_prob_SymptomsPerNoSymptoms * scale_confirmed_cases *
                                            entry.num_confirmed;
            }
            if (offset == std::ceil(i * time_INS_per_subcomp)) {
                init[last_index_INS - i] -= (i * time_INS_per_subcomp - std::floor(i * time_INS_per_subcomp)) *
                                            inv_prob_SymptomsPerNoSymptoms * scale_confirmed_cases *
                                            entry.num_confirmed;
            }
            if (offset == std::floor((i + 1) * time_INS_per_subcomp)) {
                init[last_index_INS - i] +=
                    (1 - ((i + 1) * time_INS_per_subcomp - std::floor((i + 1) * time_INS_per_subcomp))) *
                    inv_prob_SymptomsPerNoSymptoms * scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == std::ceil((i + 1) * time_INS_per_subcomp)) {
                init[last_index_INS - i] +=
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
            (ScalarType)LctState::template get_num_subcompartments<InfectionState::Exposed>();
        // Index of the last subcompartment of Exposed.
        size_t last_index_Exposed = LctState::template get_first_index<InfectionState::Exposed>() +
                                    LctState::template get_num_subcompartments<InfectionState::Exposed>() - 1;
        for (int i = 0; i < (int)LctState::template get_num_subcompartments<InfectionState::Exposed>(); i++) {
            if (offset ==
                std::floor(staytimes[(size_t)InfectionState::InfectedNoSymptoms] + i * time_Exposed_per_subcomp)) {
                init[last_index_Exposed - i] -=
                    (1 - (staytimes[(size_t)InfectionState::InfectedNoSymptoms] + i * time_Exposed_per_subcomp -
                          std::floor(staytimes[(size_t)InfectionState::InfectedNoSymptoms] +
                                     i * time_Exposed_per_subcomp))) *
                    inv_prob_SymptomsPerNoSymptoms * scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset ==
                std::ceil(staytimes[(size_t)InfectionState::InfectedNoSymptoms] + i * time_Exposed_per_subcomp)) {
                init[last_index_Exposed - i] -=
                    (staytimes[(size_t)InfectionState::InfectedNoSymptoms] + i * time_Exposed_per_subcomp -
                     std::floor(staytimes[(size_t)InfectionState::InfectedNoSymptoms] + i * time_Exposed_per_subcomp)) *
                    inv_prob_SymptomsPerNoSymptoms * scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == std::floor(staytimes[(size_t)InfectionState::InfectedNoSymptoms] +
                                     (i + 1) * time_Exposed_per_subcomp)) {
                init[last_index_Exposed - i] +=
                    (1 - (staytimes[(size_t)InfectionState::InfectedNoSymptoms] + (i + 1) * time_Exposed_per_subcomp -
                          std::floor(staytimes[(size_t)InfectionState::InfectedNoSymptoms] +
                                     (i + 1) * time_Exposed_per_subcomp))) *
                    inv_prob_SymptomsPerNoSymptoms * scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset ==
                std::ceil(staytimes[(size_t)InfectionState::InfectedNoSymptoms] + (i + 1) * time_Exposed_per_subcomp)) {
                init[last_index_Exposed - i] +=
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
            (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>();
        // Index of the first subcompartment of InfectedSymptoms.
        size_t first_index_IS = LctState::template get_first_index<InfectionState::InfectedSymptoms>();
        for (int i = 0; i < (int)LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>(); i++) {
            if (offset == std::floor(-time_IS_per_subcomp * (i + 1))) {
                init[first_index_IS + i] -=
                    (1 - (-(i + 1.) * time_IS_per_subcomp - std::floor(-(i + 1) * time_IS_per_subcomp))) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == std::ceil(-time_IS_per_subcomp * (i + 1))) {
                init[first_index_IS + i] -=
                    (-(i + 1) * time_IS_per_subcomp - std::floor(-(i + 1) * time_IS_per_subcomp)) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == std::floor(-time_IS_per_subcomp * i)) {
                init[first_index_IS + i] += (1 - (-i * time_IS_per_subcomp - std::floor(-i * time_IS_per_subcomp))) *
                                            scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == std::ceil(-time_IS_per_subcomp * i)) {
                init[first_index_IS + i] += (-i * time_IS_per_subcomp - std::floor(-i * time_IS_per_subcomp)) *
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
            (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedSevere>();
        // Index of the first subcompartment of InfectedSevere.
        size_t first_index_ISev = LctState::template get_first_index<InfectionState::InfectedSevere>();
        for (int i = 0; i < (int)LctState::template get_num_subcompartments<InfectionState::InfectedSevere>(); i++) {
            if (offset ==
                std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] - time_ISev_per_subcomp * (i + 1))) {
                init[first_index_ISev + i] -=
                    severePerInfectedSymptoms *
                    (1 - (-staytimes[(size_t)InfectionState::InfectedSymptoms] - (i + 1) * time_ISev_per_subcomp -
                          std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                     (i + 1) * time_ISev_per_subcomp))) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset ==
                std::ceil(-staytimes[(size_t)InfectionState::InfectedSymptoms] - time_ISev_per_subcomp * (i + 1))) {
                init[first_index_ISev + i] -=
                    severePerInfectedSymptoms *
                    (-staytimes[(size_t)InfectionState::InfectedSymptoms] - (i + 1) * time_ISev_per_subcomp -
                     std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                (i + 1) * time_ISev_per_subcomp)) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset ==
                std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] - time_ISev_per_subcomp * i)) {
                init[first_index_ISev + i] +=
                    severePerInfectedSymptoms *
                    (1 -
                     (-staytimes[(size_t)InfectionState::InfectedSymptoms] - i * time_ISev_per_subcomp -
                      std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] - i * time_ISev_per_subcomp))) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == std::ceil(-staytimes[(size_t)InfectionState::InfectedSymptoms] - time_ISev_per_subcomp * i)) {
                init[first_index_ISev + i] +=
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
            (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedCritical>();
        // Index of the first subcompartment of InfectedCritical.
        size_t first_index_ICri = LctState::template get_first_index<InfectionState::InfectedCritical>();
        for (int i = 0; i < (int)LctState::template get_num_subcompartments<InfectionState::InfectedCritical>(); i++) {
            if (offset ==
                std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                           staytimes[(size_t)InfectionState::InfectedSevere] - time_ICri_per_subcomp * (i + 1))) {
                init[first_index_ICri + i] -=
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
                init[first_index_ICri + i] -=
                    severePerInfectedSymptoms * criticalPerSevere *
                    (-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                     staytimes[(size_t)InfectionState::InfectedSevere] - (i + 1) * time_ICri_per_subcomp -
                     std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                staytimes[(size_t)InfectionState::InfectedSevere] - (i + 1) * time_ICri_per_subcomp)) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                     staytimes[(size_t)InfectionState::InfectedSevere] - time_ICri_per_subcomp * i)) {
                init[first_index_ICri + i] +=
                    severePerInfectedSymptoms * criticalPerSevere *
                    (1 - (-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                          staytimes[(size_t)InfectionState::InfectedSevere] - i * time_ICri_per_subcomp -
                          std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                     staytimes[(size_t)InfectionState::InfectedSevere] - i * time_ICri_per_subcomp))) *
                    scale_confirmed_cases * entry.num_confirmed;
            }
            if (offset == std::ceil(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                    staytimes[(size_t)InfectionState::InfectedSevere] - time_ICri_per_subcomp * i)) {
                init[first_index_ICri + i] +=
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
        init[LctState::template get_first_index<InfectionState::Dead>()] +=
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
        init[LctState::template get_first_index<InfectionState::Dead>()] +=
            (-staytimes[(size_t)InfectionState::InfectedSymptoms] - staytimes[(size_t)InfectionState::InfectedSevere] -
             staytimes[(size_t)InfectionState::InfectedCritical] -
             std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                        staytimes[(size_t)InfectionState::InfectedSevere] -
                        staytimes[(size_t)InfectionState::InfectedCritical])) *
            entry.num_deaths;
    }
}

/**
* @brief Computes an initialization vector for an LCT model with case data from RKI without age groups.
*   
* Please use the set_initial_data_from_confirmed_cases() function, which calls this function automatically!
* Calculates an initial value vector for an LCT model with one group and updates the initial value vector in the model.
* For the computation expected stay times in the subcompartments are used. To calculate the initial values, 
* we assume for simplicity that individuals stay in the subcompartment for exactly the expected time.
* The RKI data are linearly interpolated within one day to match the expected stay time in a subcompartment.
* The RKI data should contain data for each needed day without division of age groups, the completeness of the dates is
* not verified.
* Data can be downloaded e.g. with the file pycode/memilio-epidata/memilio/epidata/getCaseData.py, 
* which creates a file named cases_all_germany.json or a similar name. One should set impute_dates=True such 
* that missing dates are imputed.
*
* @param[in, out] model The model for which the inital data should be computed and set.
* @param[in] path Path to the RKI file.
* @param[in] date Date for which the initial values should be computed. date is the start date of the simulation.
* @param[in] total_population Total size of the population of Germany. 
*       The sum of all values in the vector of subcompartments will be equal to that value.
* @param[in] scale_confirmed_cases Factor by which to scale the confirmed cases of rki data to consider unreported cases.
* @tparam Model is expected to be an LCT-SECIR model defined in models/lct_secir/model.h with one group.
* @returns Any io errors that happen during reading of the files and the calculation.
*/
template <class Model>
IOResult<void> set_initial_data_from_confirmed_cases_notageres(Model& model, const std::string& path, Date date,
                                                               ScalarType total_population,
                                                               ScalarType scale_confirmed_cases)
{
    using LctState = typename std::tuple_element_t<0, typename Model::tupleLctStates>;

    // Try to get rki data from path.
    BOOST_OUTCOME_TRY(auto&& rki_data, mio::read_confirmed_cases_noage(path));
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
    // Compute initial values for all subcompartments.
    Eigen::VectorXd init = Eigen::VectorXd::Zero(LctState::Count);

    // Define variables for parameters.
    std::vector<ScalarType> staytimes((size_t)InfectionState::Count, -1.);
    staytimes[(size_t)InfectionState::Exposed]            = model.parameters.template get<TimeExposed>()[0];
    staytimes[(size_t)InfectionState::InfectedNoSymptoms] = model.parameters.template get<TimeInfectedNoSymptoms>()[0];
    staytimes[(size_t)InfectionState::InfectedSymptoms]   = model.parameters.template get<TimeInfectedSymptoms>()[0];
    staytimes[(size_t)InfectionState::InfectedSevere]     = model.parameters.template get<TimeInfectedSevere>()[0];
    staytimes[(size_t)InfectionState::InfectedCritical]   = model.parameters.template get<TimeInfectedCritical>()[0];
    ScalarType inv_prob_SymptomsPerNoSymptoms =
        1 / (1 - model.parameters.template get<RecoveredPerInfectedNoSymptoms>()[0]);
    ScalarType prob_SeverePerInfectedSymptoms = model.parameters.template get<SeverePerInfectedSymptoms>()[0];
    ScalarType prob_CriticalPerSevere         = model.parameters.template get<CriticalPerSevere>()[0];

    ScalarType min_offset_needed = std::floor(-staytimes[(size_t)InfectionState::InfectedSymptoms] -
                                              staytimes[(size_t)InfectionState::InfectedSevere] -
                                              staytimes[(size_t)InfectionState::InfectedCritical]);
    ScalarType max_offset_needed =
        std::ceil(staytimes[(size_t)InfectionState::Exposed] + staytimes[(size_t)InfectionState::InfectedNoSymptoms]);

    bool min_offset_needed_avail = false;
    bool max_offset_needed_avail = false;

    // Go through the entries of rki_data and check if date is needed for calculation. Confirmed cases are scaled by scale_confirmed_cases.
    for (auto&& entry : rki_data) {
        int offset = get_offset_in_days(entry.date, date);
        if ((offset >= min_offset_needed) && (offset <= max_offset_needed)) {
            if (offset == max_offset_needed) {
                max_offset_needed_avail = true;
            }
            if (offset == min_offset_needed) {
                min_offset_needed_avail = true;
            }
            process_entry<LctState, ConfirmedCasesNoAgeEntry>(
                init, entry, offset, staytimes, inv_prob_SymptomsPerNoSymptoms, prob_SeverePerInfectedSymptoms,
                prob_CriticalPerSevere, scale_confirmed_cases);
        }
    }

    // Compute Recovered.
    init[LctState::template get_first_index<InfectionState::Recovered>()] -=
        init.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>() +
                         LctState::template get_num_subcompartments<InfectionState::InfectedSevere>() +
                         LctState::template get_num_subcompartments<InfectionState::InfectedCritical>())
            .sum();
    init[LctState::template get_first_index<InfectionState::Recovered>()] -=
        init[LctState::template get_first_index<InfectionState::Dead>()];

    // Compute Susceptibles.
    init[LctState::template get_first_index<InfectionState::Susceptible>()] =
        total_population - init.segment(1, LctState::Count - 1).sum();

    // Check whether all necessary dates are available.
    if (!max_offset_needed_avail || !min_offset_needed_avail) {
        log_error("Necessary range of dates needed to compute initial values does not exist in RKI data.");
        return failure(StatusCode::OutOfRange, path + ", necessary range of dates does not exist in RKI data.");
    }

    // Set computed initial value vector.
    for (size_t i = 0; i < LctState::Count; i++) {
        model.populations[i] = init[i];
    }

    if (model.check_constraints()) {
        return failure(StatusCode::InvalidValue,
                       "One or more constraint(s) of the model are not fulfilled. Please check the model setup. "
                       "Otherwise, the data are not fitting the input variables.");
    }
    else {
        return mio::success();
    }
}

/**
* @brief Computes an initialization vector for an age-resolved LCT model with case data from RKI with age groups.
*   
* Please use the set_initial_data_from_confirmed_cases() function, which calls this function automatically!
* Calculates an initial value vector for an LCT model and updates the initial value vector in the model.
* The number of groups in the LCT model has to be the same as the number of RKI age groups. 
* For the computation expected stay times in the subcompartments are used. To calculate the initial values, 
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
* @tparam Group The age group for which the initial values should be calculted. The function is called recursively 
*   such that the initial values are calculated for every age group if Group is zero at the beginning.
* @returns Any io errors that happen during reading of the files and the calculation.
*/
template <class Model, size_t Group = 0, std::enable_if_t<(Group < Model::m_groups) && (Group >= 0), bool> = true>
IOResult<void> set_initial_data_from_confirmed_cases_ageres(Model& model,
                                                            const std::vector<ConfirmedCasesDataEntry>& rki_data,
                                                            Date date, const std::vector<ScalarType>& total_population,
                                                            const std::vector<ScalarType>& scale_confirmed_cases)
{
    using LctState = typename std::tuple_element_t<Group, typename Model::tupleLctStates>;

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

    // Compute initial values for all subcompartments of the group.
    Eigen::VectorXd init = Eigen::VectorXd::Zero(LctState::Count);

    // Define variables for parameters.
    std::vector<ScalarType> staytimes((size_t)InfectionState::Count, -1.);
    staytimes[(size_t)InfectionState::Exposed] = model.parameters.template get<TimeExposed>()[Group];
    staytimes[(size_t)InfectionState::InfectedNoSymptoms] =
        model.parameters.template get<TimeInfectedNoSymptoms>()[Group];
    staytimes[(size_t)InfectionState::InfectedSymptoms] = model.parameters.template get<TimeInfectedSymptoms>()[Group];
    staytimes[(size_t)InfectionState::InfectedSevere]   = model.parameters.template get<TimeInfectedSevere>()[Group];
    staytimes[(size_t)InfectionState::InfectedCritical] = model.parameters.template get<TimeInfectedCritical>()[Group];
    ScalarType inv_prob_SymptomsPerNoSymptoms =
        1 / (1 - model.parameters.template get<RecoveredPerInfectedNoSymptoms>()[Group]);
    ScalarType prob_SeverePerInfectedSymptoms = model.parameters.template get<SeverePerInfectedSymptoms>()[Group];
    ScalarType prob_CriticalPerSevere         = model.parameters.template get<CriticalPerSevere>()[Group];

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
        if ((size_t)entry.age_group == Group) {
            int offset = get_offset_in_days(entry.date, date);
            if ((offset >= min_offset_needed) && (offset <= max_offset_needed)) {
                if (offset == max_offset_needed) {
                    max_offset_needed_avail = true;
                }
                if (offset == min_offset_needed) {
                    min_offset_needed_avail = true;
                }
                process_entry<LctState, ConfirmedCasesDataEntry>(
                    init, entry, offset, staytimes, inv_prob_SymptomsPerNoSymptoms, prob_SeverePerInfectedSymptoms,
                    prob_CriticalPerSevere, scale_confirmed_cases[Group]);
            }
        }
    }

    // Compute Recovered.
    init[LctState::template get_first_index<InfectionState::Recovered>()] -=
        init.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>() +
                         LctState::template get_num_subcompartments<InfectionState::InfectedSevere>() +
                         LctState::template get_num_subcompartments<InfectionState::InfectedCritical>())
            .sum();
    init[LctState::template get_first_index<InfectionState::Recovered>()] -=
        init[LctState::template get_first_index<InfectionState::Dead>()];

    // Compute Susceptibles.
    init[LctState::template get_first_index<InfectionState::Susceptible>()] =
        total_population[Group] - init.segment(1, LctState::Count - 1).sum();

    // Check whether all necessary dates are available.
    if (!max_offset_needed_avail || !min_offset_needed_avail) {
        log_error("Necessary range of dates needed to compute initial values does not exist in RKI data.");
        return failure(StatusCode::OutOfRange,
                       "Necessary range of dates needed to compute initial values does not exist in RKI data.");
    }

    // Set computed initial value vector.
    size_t first_index = model.populations.template get_first_index_group<Group>();
    for (size_t i = 0; i < LctState::Count; i++) {
        model.populations[first_index + i] = init[i];
    }
    if constexpr (Group + 1 < Model::m_groups) {
        return set_initial_data_from_confirmed_cases_ageres<Model, Group + 1>(model, rki_data, date, total_population,
                                                                              scale_confirmed_cases);
    }
    else {
        if (model.check_constraints()) {
            return failure(StatusCode::InvalidValue,
                           "One or more constraint(s) of the model are not fulfilled. Please check the model setup. "
                           "Otherwise, the data are not fitting the input variables.");
        }
        else {
            return mio::success();
        }
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
template <class Model>
IOResult<void> set_initial_data_from_confirmed_cases(Model& model, const std::string& path, const Date date,
                                                     const std::vector<ScalarType>& total_population,
                                                     const std::vector<ScalarType>& scale_confirmed_cases)
{

    assert(total_population.size() == Model::m_groups);
    assert(scale_confirmed_cases.size() == Model::m_groups);
    if constexpr (Model::m_groups > 1) {
        assert(ConfirmedCasesDataEntry::age_group_names.size() == Model::m_groups);
        BOOST_OUTCOME_TRY(auto&& rki_data, mio::read_confirmed_cases_data(path));
        return details::set_initial_data_from_confirmed_cases_ageres<Model>(model, rki_data, date, total_population,
                                                                            scale_confirmed_cases);
    }
    else {
        return details::set_initial_data_from_confirmed_cases_notageres<Model>(model, path, date, total_population[0],
                                                                               scale_confirmed_cases[0]);
    }
}

} // namespace lsecir
} // namespace mio
#endif // MEMILIO_HAS_JSONCPP

#endif // LCTSECIR_PARAMETERS_IO_H
