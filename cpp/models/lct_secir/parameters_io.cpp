/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

#include "lct_secir/parameters_io.h"
#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

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

IOResult<Eigen::VectorXd> get_initial_data_from_file(std::string const& path, Date date, InfectionState infectionState,
                                                     Parameters&& parameters, ScalarType total_population,
                                                     ScalarType scale_confirmed_cases)
{
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
        log_error("Specified date does not exist in RKI data");
        return failure(StatusCode::OutOfRange, path + ", specified date does not exist in RKI data.");
    }
    // Compute initial values for all subcompartments.
    Eigen::VectorXd init = Eigen::VectorXd::Zero(infectionState.get_count());
    // TODO: was soll die 6 da?
    //auto days_surplus = std::min(get_offset_in_days(max_date, date) - 6, 0);
    // Define variables for parameters that are often needed.
    ScalarType timeExposed            = parameters.get<TimeExposed>();
    ScalarType timeInfectedNoSymptoms = parameters.get<TimeInfectedNoSymptoms>();
    ScalarType timeInfectedSymptoms   = parameters.get<TimeInfectedSymptoms>();
    ScalarType timeInfectedSevere     = parameters.get<TimeInfectedSevere>();
    ScalarType timeInfectedCritical   = parameters.get<TimeInfectedCritical>();
    ScalarType scale_mu_CR            = 1 / (1 - parameters.get<RecoveredPerInfectedNoSymptoms>());

    ScalarType min_offset_needed = std::floor(-timeInfectedSymptoms - timeInfectedSevere - timeInfectedCritical);
    ScalarType max_offset_needed = std::ceil(timeExposed + timeInfectedNoSymptoms);

    // Go through the entries of rki_data and check if date is needed for calculation. Confirmed cases are scaled.
    for (auto&& entry : rki_data) {
        auto offset = get_offset_in_days(entry.date, date);
        if (!(offset < min_offset_needed) || !(offset > max_offset_needed)) {
            // Add confirmed cases at date to compartment Recovered.
            if (offset == 0) {
                init[infectionState.get_firstindex(InfectionStateBase::Recovered)] +=
                    scale_confirmed_cases * entry.num_confirmed;
            }

            // Compute initial values for compartment InfectedNoSymptoms.
            if (offset >= 0 && offset <= std::ceil(timeInfectedNoSymptoms)) {
                ScalarType TCi =
                    timeInfectedNoSymptoms / infectionState.get_number(InfectionStateBase::InfectedNoSymptoms);
                int idx_Cn = infectionState.get_firstindex(InfectionStateBase::InfectedNoSymptoms) +
                             infectionState.get_number(InfectionStateBase::InfectedNoSymptoms) - 1;
                for (int i = 0; i < infectionState.get_number(InfectionStateBase::InfectedNoSymptoms); i++) {
                    if (offset == std::floor(i * TCi)) {
                        init[idx_Cn - i] -= (1 - (i * TCi - std::floor(i * TCi))) * scale_mu_CR *
                                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(i * TCi)) {
                        init[idx_Cn - i] -=
                            (i * TCi - std::floor(i * TCi)) * scale_mu_CR * scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::floor((i + 1) * TCi)) {
                        init[idx_Cn - i] += (1 - ((i + 1) * TCi - std::floor((i + 1) * TCi))) * scale_mu_CR *
                                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil((i + 1) * TCi)) {
                        init[idx_Cn - i] += ((i + 1) * TCi - std::floor((i + 1) * TCi)) * scale_mu_CR *
                                            scale_confirmed_cases * entry.num_confirmed;
                    }
                }
            }

            // Compute initial values for compartment Exposed.
            if (offset >= std::floor(timeInfectedNoSymptoms) && offset <= max_offset_needed) {
                ScalarType TEi = timeExposed / infectionState.get_number(InfectionStateBase::Exposed);
                int idx_En     = infectionState.get_firstindex(InfectionStateBase::Exposed) +
                             infectionState.get_number(InfectionStateBase::Exposed) - 1;
                for (int i = 0; i < infectionState.get_number(InfectionStateBase::Exposed); i++) {
                    if (offset == std::floor(timeInfectedNoSymptoms + i * TEi)) {
                        init[idx_En - i] -=
                            (1 - (timeInfectedNoSymptoms + i * TEi - std::floor(timeInfectedNoSymptoms + i * TEi))) *
                            scale_mu_CR * scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(timeInfectedNoSymptoms + i * TEi)) {
                        init[idx_En - i] -=
                            (timeInfectedNoSymptoms + i * TEi - std::floor(timeInfectedNoSymptoms + i * TEi)) *
                            scale_mu_CR * scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::floor(timeInfectedNoSymptoms + (i + 1) * TEi)) {
                        init[idx_En - i] += (1 - (timeInfectedNoSymptoms + (i + 1) * TEi -
                                                  std::floor(timeInfectedNoSymptoms + (i + 1) * TEi))) *
                                            scale_mu_CR * scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(timeInfectedNoSymptoms + (i + 1) * TEi)) {
                        init[idx_En - i] += (timeInfectedNoSymptoms + (i + 1) * TEi -
                                             std::floor(timeInfectedNoSymptoms + (i + 1) * TEi)) *
                                            scale_mu_CR * scale_confirmed_cases * entry.num_confirmed;
                    }
                }
            }

            // Compute initial values for compartment InfectedSymptoms.
            if (offset >= std::floor(-timeInfectedSymptoms) && offset <= 0) {
                ScalarType TIi = timeInfectedSymptoms / infectionState.get_number(InfectionStateBase::InfectedSymptoms);
                int idx_I1     = infectionState.get_firstindex(InfectionStateBase::InfectedSymptoms);
                for (int i = 0; i < infectionState.get_number(InfectionStateBase::InfectedSymptoms); i++) {
                    if (offset == std::floor(-TIi * (i + 1))) {
                        init[idx_I1 + i] -= (1 - (-(i + 1) * TIi - std::floor(-(i + 1) * TIi))) *
                                            scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::ceil(-TIi * (i + 1))) {
                        init[idx_I1 + i] -=
                            (-(i + 1) * TIi - std::floor(-(i + 1) * TIi)) * scale_confirmed_cases * entry.num_confirmed;
                    }
                    if (offset == std::floor(-TIi * i)) {
                        init[idx_I1 + i] +=
                            (1 - (-i * TIi - std::floor(-i * TIi))) * scale_confirmed_cases * entry.num_confirmed;
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
                ScalarType THi   = timeInfectedSevere / infectionState.get_number(InfectionStateBase::InfectedSevere);
                ScalarType mu_IH = parameters.get<SeverePerInfectedSymptoms>();
                int idx_H1       = infectionState.get_firstindex(InfectionStateBase::InfectedSevere);
                for (int i = 0; i < infectionState.get_number(InfectionStateBase::InfectedSevere); i++) {
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
                ScalarType TUi = timeInfectedCritical / infectionState.get_number(InfectionStateBase::InfectedCritical);
                ScalarType mu_IH = parameters.get<SeverePerInfectedSymptoms>();
                ScalarType mu_HU = parameters.get<CriticalPerSevere>();
                int idx_U1       = infectionState.get_firstindex(InfectionStateBase::InfectedCritical);
                for (int i = 0; i < infectionState.get_number(InfectionStateBase::InfectedCritical); i++) {
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
                init[infectionState.get_firstindex(InfectionStateBase::Dead)] +=
                    (1 - (-timeInfectedSymptoms - timeInfectedSevere - timeInfectedCritical -
                          std::floor(-timeInfectedSymptoms - timeInfectedSevere - timeInfectedCritical))) *
                    entry.num_deaths;
            }
            if (offset == std::ceil(-timeInfectedSymptoms - timeInfectedSevere - timeInfectedCritical)) {
                init[infectionState.get_firstindex(InfectionStateBase::Dead)] +=
                    (-timeInfectedSymptoms - timeInfectedSevere - timeInfectedCritical -
                     std::floor(-timeInfectedSymptoms - timeInfectedSevere - timeInfectedCritical)) *
                    entry.num_deaths;
            }
        }
    }

    // Compute Recovered.
    init[infectionState.get_firstindex(InfectionStateBase::Recovered)] -=
        init.segment(infectionState.get_firstindex(InfectionStateBase::InfectedSymptoms),
                     infectionState.get_number(InfectionStateBase::InfectedSymptoms) +
                         infectionState.get_number(InfectionStateBase::InfectedSevere) +
                         infectionState.get_number(InfectionStateBase::InfectedCritical))
            .sum();
    init[infectionState.get_firstindex(InfectionStateBase::Recovered)] -=
        init[infectionState.get_firstindex(InfectionStateBase::Dead)];

    // Compute Susceptibles
    init[infectionState.get_firstindex(InfectionStateBase::Susceptible)] =
        total_population - init.segment(1, infectionState.get_count() - 1).sum();

    return init;
}

} // namespace lsecir
} // namespace mio
#endif // MEMILIO_HAS_JSONCPP