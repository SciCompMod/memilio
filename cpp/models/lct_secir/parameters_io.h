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
* Computes an initial value vector for an LCT model with the defined infectionState and the parameters. 
* For the computation expected stay times in the subcompartments are used. To calculate the initial values, 
* we assume for simplicity that individuals stay in the subcompartment for exactly the expected time.
* The RKI data are linearly interpolated within one day.
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
* @returns Any io errors that happen during reading of the files.
*/
template <class Model>
IOResult<void> set_initial_data_from_confirmed_cases(Model& model, const std::string& path, Date date,
                                                     ScalarType total_population, ScalarType scale_confirmed_cases = 1.)
{
    using LctState = typename Model::LctState;

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
    // Compute initial values for all subcompartments.
    Eigen::VectorXd init = Eigen::VectorXd::Zero(LctState::Count);
    // Define variables for parameters that are often needed.
    Parameters parameters = model.parameters;

    ScalarType timeExposed            = parameters.get<TimeExposed>();
    ScalarType timeInfectedNoSymptoms = parameters.get<TimeInfectedNoSymptoms>();
    ScalarType timeInfectedSymptoms   = parameters.get<TimeInfectedSymptoms>();
    ScalarType timeInfectedSevere     = parameters.get<TimeInfectedSevere>();
    ScalarType timeInfectedCritical   = parameters.get<TimeInfectedCritical>();
    ScalarType scale_asymptomatic     = 1 / (1 - parameters.get<RecoveredPerInfectedNoSymptoms>());

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
                ScalarType TCi = timeInfectedNoSymptoms /
                                 LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>();
                int idx_Cn = LctState::template get_first_index<InfectionState::InfectedNoSymptoms>() +
                             LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>() - 1;
                for (int i = 0;
                     i < (int)LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>(); i++) {
                    if (offset == std::floor(i * TCi)) {
                        init[idx_Cn - i] -= (1 - (i * TCi - std::floor(i * TCi))) * scale_asymptomatic *
                                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(i * TCi)) {
                        init[idx_Cn - i] -= (i * TCi - std::floor(i * TCi)) * scale_asymptomatic *
                                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::floor((i + 1) * TCi)) {
                        init[idx_Cn - i] += (1 - ((i + 1) * TCi - std::floor((i + 1) * TCi))) * scale_asymptomatic *
                                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil((i + 1) * TCi)) {
                        init[idx_Cn - i] += ((i + 1) * TCi - std::floor((i + 1) * TCi)) * scale_asymptomatic *
                                            scale_confirmed_cases * entry.num_confirmed;
                    }
                }
            }

            // Compute initial values for compartment Exposed.
            if (offset >= std::floor(timeInfectedNoSymptoms) && offset <= max_offset_needed) {
                ScalarType TEi = timeExposed / LctState::template get_num_subcompartments<InfectionState::Exposed>();
                int idx_En     = LctState::template get_first_index<InfectionState::Exposed>() +
                             LctState::template get_num_subcompartments<InfectionState::Exposed>() - 1;
                for (int i = 0; i < (int)LctState::template get_num_subcompartments<InfectionState::Exposed>(); i++) {
                    if (offset == std::floor(timeInfectedNoSymptoms + i * TEi)) {
                        init[idx_En - i] -=
                            (1 - (timeInfectedNoSymptoms + i * TEi - std::floor(timeInfectedNoSymptoms + i * TEi))) *
                            scale_asymptomatic * scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(timeInfectedNoSymptoms + i * TEi)) {
                        init[idx_En - i] -=
                            (timeInfectedNoSymptoms + i * TEi - std::floor(timeInfectedNoSymptoms + i * TEi)) *
                            scale_asymptomatic * scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::floor(timeInfectedNoSymptoms + (i + 1) * TEi)) {
                        init[idx_En - i] += (1 - (timeInfectedNoSymptoms + (i + 1) * TEi -
                                                  std::floor(timeInfectedNoSymptoms + (i + 1) * TEi))) *
                                            scale_asymptomatic * scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(timeInfectedNoSymptoms + (i + 1) * TEi)) {
                        init[idx_En - i] += (timeInfectedNoSymptoms + (i + 1) * TEi -
                                             std::floor(timeInfectedNoSymptoms + (i + 1) * TEi)) *
                                            scale_asymptomatic * scale_confirmed_cases * entry.num_confirmed;
                    }
                }
            }

            // Compute initial values for compartment InfectedSymptoms.
            if (offset >= std::floor(-timeInfectedSymptoms) && offset <= 0) {
                ScalarType TIi =
                    timeInfectedSymptoms /
                    (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>();
                int idx_I1 = LctState::template get_first_index<InfectionState::InfectedSymptoms>();
                for (int i = 0; i < (int)LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>();
                     i++) {
                    if (offset == std::floor(-TIi * (i + 1))) {
                        init[idx_I1 + i] -= (1 - (-(i + 1.) * TIi - std::floor(-(i + 1) * TIi))) *
                                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(-TIi * (i + 1))) {
                        init[idx_I1 + i] -=
                            (-(i + 1) * TIi - std::floor(-(i + 1) * TIi)) * scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::floor(-TIi * i)) {
                        init[idx_I1 + i] += (1 - (-1 * i * TIi - std::floor(-1 * i * TIi))) * scale_confirmed_cases *
                                            entry.num_confirmed;
                    }
                    if (offset == std::ceil(-TIi * i)) {
                        init[idx_I1 + i] +=
                            (-i * TIi - std::floor(-i * TIi)) * scale_confirmed_cases * entry.num_confirmed;
                    }
                }
            }

            // Compute initial values for compartment InfectedSevere.
            if (offset >= std::floor(-timeInfectedSymptoms - timeInfectedSevere) &&
                offset <= std::ceil(-timeInfectedSymptoms)) {
                ScalarType THi =
                    timeInfectedSevere /
                    (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedSevere>();
                ScalarType mu_IH = parameters.get<SeverePerInfectedSymptoms>();
                int idx_H1       = LctState::template get_first_index<InfectionState::InfectedSevere>();
                for (int i = 0; i < (int)LctState::template get_num_subcompartments<InfectionState::InfectedSevere>();
                     i++) {
                    if (offset == std::floor(-timeInfectedSymptoms - THi * (i + 1))) {
                        init[idx_H1 + i] -= mu_IH *
                                            (1 - (-timeInfectedSymptoms - (i + 1) * THi -
                                                  std::floor(-timeInfectedSymptoms - (i + 1) * THi))) *
                                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(-timeInfectedSymptoms - THi * (i + 1))) {
                        init[idx_H1 + i] -= mu_IH *
                                            (-timeInfectedSymptoms - (i + 1) * THi -
                                             std::floor(-timeInfectedSymptoms - (i + 1) * THi)) *
                                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::floor(-timeInfectedSymptoms - THi * i)) {
                        init[idx_H1 + i] +=
                            mu_IH *
                            (1 - (-timeInfectedSymptoms - i * THi - std::floor(-timeInfectedSymptoms - i * THi))) *
                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(-timeInfectedSymptoms - THi * i)) {
                        init[idx_H1 + i] +=
                            mu_IH * (-timeInfectedSymptoms - i * THi - std::floor(-timeInfectedSymptoms - i * THi)) *
                            scale_confirmed_cases * entry.num_confirmed;
                    }
                }
            }

            // Compute initial values for compartment InfectedCritical.
            if (offset >= min_offset_needed && offset <= std::ceil(-timeInfectedSymptoms - timeInfectedSevere)) {
                ScalarType TUi = timeInfectedCritical /
                                 LctState::template get_num_subcompartments<InfectionState::InfectedCritical>();
                ScalarType mu_IH = parameters.get<SeverePerInfectedSymptoms>();
                ScalarType mu_HU = parameters.get<CriticalPerSevere>();
                int idx_U1       = LctState::template get_first_index<InfectionState::InfectedCritical>();
                for (int i = 0; i < (int)LctState::template get_num_subcompartments<InfectionState::InfectedCritical>();
                     i++) {
                    if (offset == std::floor(-timeInfectedSymptoms - timeInfectedSevere - TUi * (i + 1))) {
                        init[idx_U1 + i] -=
                            mu_IH * mu_HU *
                            (1 - (-timeInfectedSymptoms - timeInfectedSevere - (i + 1) * TUi -
                                  std::floor(-timeInfectedSymptoms - timeInfectedSevere - (i + 1) * TUi))) *
                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(-timeInfectedSymptoms - timeInfectedSevere - TUi * (i + 1))) {
                        init[idx_U1 + i] -= mu_IH * mu_HU *
                                            (-timeInfectedSymptoms - timeInfectedSevere - (i + 1) * TUi -
                                             std::floor(-timeInfectedSymptoms - timeInfectedSevere - (i + 1) * TUi)) *
                                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::floor(-timeInfectedSymptoms - timeInfectedSevere - TUi * i)) {
                        init[idx_U1 + i] += mu_IH * mu_HU *
                                            (1 - (-timeInfectedSymptoms - timeInfectedSevere - i * TUi -
                                                  std::floor(-timeInfectedSymptoms - timeInfectedSevere - i * TUi))) *
                                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(-timeInfectedSymptoms - timeInfectedSevere - TUi * i)) {
                        init[idx_U1 + i] += mu_IH * mu_HU *
                                            (-timeInfectedSymptoms - timeInfectedSevere - i * TUi -
                                             std::floor(-timeInfectedSymptoms - timeInfectedSevere - i * TUi)) *
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

    // Compute Susceptibles
    init[LctState::template get_first_index<InfectionState::Susceptible>()] =
        total_population - init.segment(1, LctState::Count - 1).sum();
    if (!max_offset_needed_avail || !min_offset_needed_avail) {
        log_error("Necessary range of dates needed to compute initial values does not exist in RKI data.");
        return failure(StatusCode::OutOfRange, path + ", necessary range of dates does not exist in RKI data.");
    }

    model.set_initial_values(std::move(init));
    return mio::success();
}

} // namespace lsecir
} // namespace mio
#endif // MEMILIO_HAS_JSONCPP

#endif // LCTSECIR_PARAMETERS_IO_H
