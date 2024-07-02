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
#include "memilio/math/eigen.h"

#include <string>
#include <iostream>
namespace mio
{
namespace lsecir
{

/**
* @brief Computes an initialization vector for a LCT model with case data from RKI.
*   
* Calculates an initial value vector for an LCT model and updates the initial value vector in the model.
* For the computation expected stay times in the subcompartments are used. To calculate the initial values, 
* we assume for simplicity that individuals stay in the subcompartment for exactly the expected time.
* The RKI data are linearly interpolated within one day to match the expected stay time in a subcompartment.
* The RKI data should contain data for each needed day without division of age groups, the completeness of the dates is not verified.
* Data can be downloaded e.g. with the file pycode/memilio-epidata/memilio/epidata/getCaseData.py, 
* which creates a file named cases_all_germany.json or a similar name. One should set impute_dates=True so that missing dates are imputed.
*
* @param[in, out] model The model for which the inital data should be computed and set.
* @param[in] path Path to the RKI file.
* @param[in] date Date for which the initial values should be computed. date is the start date of the simulation.
* @param[in] total_population Total size of the population of Germany. 
*       The sum of all values in the vector of subcompartments will be equal to that value.
* @param[in] scale_confirmed_cases Factor by which to scale the confirmed cases of rki data to consider unreported cases.
* @tparam Model is expected to be an LCT-SECIR model defined in models/lct_secir/model.h.
* @returns Any io errors that happen during reading of the files.
*/
template <class Model>
IOResult<void> set_initial_data_from_confirmed_cases(Model& model, const std::string& path, Date date,
                                                     ScalarType total_population, ScalarType scale_confirmed_cases = 1.)
{
    using LctState = typename Model::LctState;

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
    // Define variables for parameters that are often needed.

    ScalarType timeExposed            = model.parameters.template get<TimeExposed>();
    ScalarType timeInfectedNoSymptoms = model.parameters.template get<TimeInfectedNoSymptoms>();
    ScalarType timeInfectedSymptoms   = model.parameters.template get<TimeInfectedSymptoms>();
    ScalarType timeInfectedSevere     = model.parameters.template get<TimeInfectedSevere>();
    ScalarType timeInfectedCritical   = model.parameters.template get<TimeInfectedCritical>();
    ScalarType scale_asymptomatic     = 1 / (1 - model.parameters.template get<RecoveredPerInfectedNoSymptoms>());

    ScalarType min_offset_needed = std::floor(-timeInfectedSymptoms - timeInfectedSevere - timeInfectedCritical);
    ScalarType max_offset_needed = std::ceil(timeExposed + timeInfectedNoSymptoms);

    bool min_offset_needed_avail = false;
    bool max_offset_needed_avail = false;

    // Go through the entries of rki_data and check if date is needed for calculation. Confirmed cases are scaled by scale_confirmed_cases.
    for (auto&& entry : rki_data) {
        int offset = get_offset_in_days(entry.date, date);
        if ((offset >= min_offset_needed) && (offset <= max_offset_needed)) {
            // Add confirmed cases at date to compartment Recovered.
            if (offset == 0) {
                init[LctState::template get_first_index<InfectionState::Recovered>()] +=
                    scale_confirmed_cases * entry.num_confirmed;
            }

            // Compute initial values for compartment InfectedNoSymptoms.
            if (offset >= 0 && offset <= std::ceil(timeInfectedNoSymptoms)) {
                // Mean stay time in each subcompartment.
                ScalarType timeInfectedNoSymptoms_i =
                    timeInfectedNoSymptoms /
                    (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>();
                // Index of the last subcompartment of InfectedNoSymptoms.
                size_t idxInfectedNoSymptoms_last =
                    LctState::template get_first_index<InfectionState::InfectedNoSymptoms>() +
                    LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>() - 1;
                for (int i = 0;
                     i < (int)LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>(); i++) {
                    if (offset == std::floor(i * timeInfectedNoSymptoms_i)) {
                        init[idxInfectedNoSymptoms_last - i] -=
                            (1 - (i * timeInfectedNoSymptoms_i - std::floor(i * timeInfectedNoSymptoms_i))) *
                            scale_asymptomatic * scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(i * timeInfectedNoSymptoms_i)) {
                        init[idxInfectedNoSymptoms_last - i] -=
                            (i * timeInfectedNoSymptoms_i - std::floor(i * timeInfectedNoSymptoms_i)) *
                            scale_asymptomatic * scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::floor((i + 1) * timeInfectedNoSymptoms_i)) {
                        init[idxInfectedNoSymptoms_last - i] +=
                            (1 -
                             ((i + 1) * timeInfectedNoSymptoms_i - std::floor((i + 1) * timeInfectedNoSymptoms_i))) *
                            scale_asymptomatic * scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil((i + 1) * timeInfectedNoSymptoms_i)) {
                        init[idxInfectedNoSymptoms_last - i] +=
                            ((i + 1) * timeInfectedNoSymptoms_i - std::floor((i + 1) * timeInfectedNoSymptoms_i)) *
                            scale_asymptomatic * scale_confirmed_cases * entry.num_confirmed;
                    }
                }
            }

            // Compute initial values for compartment Exposed.
            if (offset >= std::floor(timeInfectedNoSymptoms) && offset <= max_offset_needed) {
                // Mean stay time in each subcompartment.
                ScalarType timeExposed_i =
                    timeExposed / (ScalarType)LctState::template get_num_subcompartments<InfectionState::Exposed>();
                // Index of the last subcompartment of Exposed.
                size_t idxExposed_last = LctState::template get_first_index<InfectionState::Exposed>() +
                                         LctState::template get_num_subcompartments<InfectionState::Exposed>() - 1;
                for (int i = 0; i < (int)LctState::template get_num_subcompartments<InfectionState::Exposed>(); i++) {
                    if (offset == std::floor(timeInfectedNoSymptoms + i * timeExposed_i)) {
                        init[idxExposed_last - i] -= (1 - (timeInfectedNoSymptoms + i * timeExposed_i -
                                                           std::floor(timeInfectedNoSymptoms + i * timeExposed_i))) *
                                                     scale_asymptomatic * scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(timeInfectedNoSymptoms + i * timeExposed_i)) {
                        init[idxExposed_last - i] -= (timeInfectedNoSymptoms + i * timeExposed_i -
                                                      std::floor(timeInfectedNoSymptoms + i * timeExposed_i)) *
                                                     scale_asymptomatic * scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::floor(timeInfectedNoSymptoms + (i + 1) * timeExposed_i)) {
                        init[idxExposed_last - i] +=
                            (1 - (timeInfectedNoSymptoms + (i + 1) * timeExposed_i -
                                  std::floor(timeInfectedNoSymptoms + (i + 1) * timeExposed_i))) *
                            scale_asymptomatic * scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(timeInfectedNoSymptoms + (i + 1) * timeExposed_i)) {
                        init[idxExposed_last - i] += (timeInfectedNoSymptoms + (i + 1) * timeExposed_i -
                                                      std::floor(timeInfectedNoSymptoms + (i + 1) * timeExposed_i)) *
                                                     scale_asymptomatic * scale_confirmed_cases * entry.num_confirmed;
                    }
                }
            }

            // Compute initial values for compartment InfectedSymptoms.
            if (offset >= std::floor(-timeInfectedSymptoms) && offset <= 0) {
                // Mean stay time in each subcompartment.
                ScalarType timeInfectedSymptoms_i =
                    timeInfectedSymptoms /
                    (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>();
                // Index of the first subcompartment of InfectedSymptoms.
                size_t idxInfectedSymptoms_first =
                    LctState::template get_first_index<InfectionState::InfectedSymptoms>();
                for (int i = 0; i < (int)LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>();
                     i++) {
                    if (offset == std::floor(-timeInfectedSymptoms_i * (i + 1))) {
                        init[idxInfectedSymptoms_first + i] -=
                            (1 - (-(i + 1.) * timeInfectedSymptoms_i - std::floor(-(i + 1) * timeInfectedSymptoms_i))) *
                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(-timeInfectedSymptoms_i * (i + 1))) {
                        init[idxInfectedSymptoms_first + i] -=
                            (-(i + 1) * timeInfectedSymptoms_i - std::floor(-(i + 1) * timeInfectedSymptoms_i)) *
                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::floor(-timeInfectedSymptoms_i * i)) {
                        init[idxInfectedSymptoms_first + i] +=
                            (1 - (-i * timeInfectedSymptoms_i - std::floor(-i * timeInfectedSymptoms_i))) *
                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(-timeInfectedSymptoms_i * i)) {
                        init[idxInfectedSymptoms_first + i] +=
                            (-i * timeInfectedSymptoms_i - std::floor(-i * timeInfectedSymptoms_i)) *
                            scale_confirmed_cases * entry.num_confirmed;
                    }
                }
            }

            // Compute initial values for compartment InfectedSevere.
            if (offset >= std::floor(-timeInfectedSymptoms - timeInfectedSevere) &&
                offset <= std::ceil(-timeInfectedSymptoms)) {
                // Mean stay time in each subcompartment.
                ScalarType timeInfectedSevere_i =
                    timeInfectedSevere /
                    (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedSevere>();
                // Transmission probability.
                ScalarType prob_SeverePerInfectedSymptoms = model.parameters.template get<SeverePerInfectedSymptoms>();
                // Index of the first subcompartment of InfectedSevere.
                size_t idxInfectedSevere_first = LctState::template get_first_index<InfectionState::InfectedSevere>();
                for (int i = 0; i < (int)LctState::template get_num_subcompartments<InfectionState::InfectedSevere>();
                     i++) {
                    if (offset == std::floor(-timeInfectedSymptoms - timeInfectedSevere_i * (i + 1))) {
                        init[idxInfectedSevere_first + i] -=
                            prob_SeverePerInfectedSymptoms *
                            (1 - (-timeInfectedSymptoms - (i + 1) * timeInfectedSevere_i -
                                  std::floor(-timeInfectedSymptoms - (i + 1) * timeInfectedSevere_i))) *
                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(-timeInfectedSymptoms - timeInfectedSevere_i * (i + 1))) {
                        init[idxInfectedSevere_first + i] -=
                            prob_SeverePerInfectedSymptoms *
                            (-timeInfectedSymptoms - (i + 1) * timeInfectedSevere_i -
                             std::floor(-timeInfectedSymptoms - (i + 1) * timeInfectedSevere_i)) *
                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::floor(-timeInfectedSymptoms - timeInfectedSevere_i * i)) {
                        init[idxInfectedSevere_first + i] +=
                            prob_SeverePerInfectedSymptoms *
                            (1 - (-timeInfectedSymptoms - i * timeInfectedSevere_i -
                                  std::floor(-timeInfectedSymptoms - i * timeInfectedSevere_i))) *
                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(-timeInfectedSymptoms - timeInfectedSevere_i * i)) {
                        init[idxInfectedSevere_first + i] +=
                            prob_SeverePerInfectedSymptoms *
                            (-timeInfectedSymptoms - i * timeInfectedSevere_i -
                             std::floor(-timeInfectedSymptoms - i * timeInfectedSevere_i)) *
                            scale_confirmed_cases * entry.num_confirmed;
                    }
                }
            }

            // Compute initial values for compartment InfectedCritical.
            if (offset >= min_offset_needed && offset <= std::ceil(-timeInfectedSymptoms - timeInfectedSevere)) {
                // Mean stay time in each subcompartment.
                ScalarType timeInfectedCritical_i =
                    timeInfectedCritical /
                    (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedCritical>();
                // Transmission probabilities.
                ScalarType prob_SeverePerInfectedSymptoms = model.parameters.template get<SeverePerInfectedSymptoms>();
                ScalarType prob_CriticalPerSevere         = model.parameters.template get<CriticalPerSevere>();
                // Index of the first subcompartment of InfectedCritical.
                size_t idxInfectedCritical_first =
                    LctState::template get_first_index<InfectionState::InfectedCritical>();
                for (int i = 0; i < (int)LctState::template get_num_subcompartments<InfectionState::InfectedCritical>();
                     i++) {
                    if (offset ==
                        std::floor(-timeInfectedSymptoms - timeInfectedSevere - timeInfectedCritical_i * (i + 1))) {
                        init[idxInfectedCritical_first + i] -=
                            prob_SeverePerInfectedSymptoms * prob_CriticalPerSevere *
                            (1 - (-timeInfectedSymptoms - timeInfectedSevere - (i + 1) * timeInfectedCritical_i -
                                  std::floor(-timeInfectedSymptoms - timeInfectedSevere -
                                             (i + 1) * timeInfectedCritical_i))) *
                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset ==
                        std::ceil(-timeInfectedSymptoms - timeInfectedSevere - timeInfectedCritical_i * (i + 1))) {
                        init[idxInfectedCritical_first + i] -=
                            prob_SeverePerInfectedSymptoms * prob_CriticalPerSevere *
                            (-timeInfectedSymptoms - timeInfectedSevere - (i + 1) * timeInfectedCritical_i -
                             std::floor(-timeInfectedSymptoms - timeInfectedSevere -
                                        (i + 1) * timeInfectedCritical_i)) *
                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::floor(-timeInfectedSymptoms - timeInfectedSevere - timeInfectedCritical_i * i)) {
                        init[idxInfectedCritical_first + i] +=
                            prob_SeverePerInfectedSymptoms * prob_CriticalPerSevere *
                            (1 -
                             (-timeInfectedSymptoms - timeInfectedSevere - i * timeInfectedCritical_i -
                              std::floor(-timeInfectedSymptoms - timeInfectedSevere - i * timeInfectedCritical_i))) *
                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(-timeInfectedSymptoms - timeInfectedSevere - timeInfectedCritical_i * i)) {
                        init[idxInfectedCritical_first + i] +=
                            prob_SeverePerInfectedSymptoms * prob_CriticalPerSevere *
                            (-timeInfectedSymptoms - timeInfectedSevere - i * timeInfectedCritical_i -
                             std::floor(-timeInfectedSymptoms - timeInfectedSevere - i * timeInfectedCritical_i)) *
                            scale_confirmed_cases * entry.num_confirmed;
                    }
                }
            }

            // Compute Dead.
            if (offset == min_offset_needed) {
                min_offset_needed_avail = true;
                init[LctState::template get_first_index<InfectionState::Dead>()] +=
                    (1 - (-timeInfectedSymptoms - timeInfectedSevere - timeInfectedCritical -
                          std::floor(-timeInfectedSymptoms - timeInfectedSevere - timeInfectedCritical))) *
                    entry.num_deaths;
            }
            if (offset == std::ceil(-timeInfectedSymptoms - timeInfectedSevere - timeInfectedCritical)) {
                init[LctState::template get_first_index<InfectionState::Dead>()] +=
                    (-timeInfectedSymptoms - timeInfectedSevere - timeInfectedCritical -
                     std::floor(-timeInfectedSymptoms - timeInfectedSevere - timeInfectedCritical)) *
                    entry.num_deaths;
            }
            if (offset == max_offset_needed) {
                max_offset_needed_avail = true;
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
        total_population - init.segment(1, LctState::Count - 1).sum();

    // Check whether all necessary dates are available.
    if (!max_offset_needed_avail || !min_offset_needed_avail) {
        log_error("Necessary range of dates needed to compute initial values does not exist in RKI data.");
        return failure(StatusCode::OutOfRange, path + ", necessary range of dates does not exist in RKI data.");
    }

    // Set computed initial value vector.
    for (size_t i = 0; i < LctState::Count; i++) {
        model.populations[mio::Index<LctState>(i)] = init[i];
    }

    return mio::success();
}

} // namespace lsecir
} // namespace mio
#endif // MEMILIO_HAS_JSONCPP

#endif // LCTSECIR_PARAMETERS_IO_H
