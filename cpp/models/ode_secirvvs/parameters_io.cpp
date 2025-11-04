/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Wadim Koslow, Daniel Abele, Martin J. Kühn
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
#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "ode_secirvvs/parameters_io.h"

namespace mio
{
namespace osecirvvs
{
namespace details
{

IOResult<void> compute_confirmed_cases_data(
    const std::vector<ConfirmedCasesDataEntry>& case_data, const int region, Date date,
    std::vector<double>& num_Exposed, std::vector<double>& num_InfectedNoSymptoms,
    std::vector<double>& num_InfectedSymptoms, std::vector<double>& num_InfectedSevere,
    std::vector<double>& num_icu, std::vector<double>& num_death,
    std::vector<double>& num_rec, const std::vector<int>& t_Exposed,
    const std::vector<int>& t_InfectedNoSymptoms,
    const std::vector<int>& t_InfectedSymptoms, const std::vector<int>& t_InfectedSevere,
    const std::vector<int>& t_InfectedCritical, const std::vector<double>& mu_C_R,
    const std::vector<double>& mu_I_H, const std::vector<double>& mu_H_U,
    const std::vector<double>& scaling_factor_inf)
{
    auto max_date_entry = std::max_element(case_data.begin(), case_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == case_data.end()) {
        log_error("RKI data file is empty.");
        return failure(StatusCode::InvalidValue, "RKI data is empty.");
    }
    auto max_date = max_date_entry->date;
    if (max_date < date) {
        log_error("Specified date does not exist in RKI data");
        return failure(StatusCode::OutOfRange, "RKI data does not contain specified date.");
    }

    // shifts the initilization to the recent past if simulation starts
    // around current day and data of the future would be required.
    // Only needed for preinfection compartments, exposed and InfectedNoSymptoms.
    auto days_surplus = get_offset_in_days(max_date, date) - 6; // 6 > T_E + T_C
    if (days_surplus > 0) {
        days_surplus = 0;
    }

    for (auto&& entry : case_data) {
    
        auto age = (size_t)entry.age_group;
        if (entry.date == offset_date_by_days(date, 0)) {
            num_InfectedSymptoms[age] += scaling_factor_inf[age] * entry.num_confirmed;
            num_rec[age] += entry.num_confirmed;
        }
        if (entry.date == offset_date_by_days(date, t_InfectedNoSymptoms[age] + days_surplus)) {
            num_InfectedNoSymptoms[age] += 1 / (1 - mu_C_R[age]) * scaling_factor_inf[age] * entry.num_confirmed;
            num_Exposed[age] -= 1 / (1 - mu_C_R[age]) * scaling_factor_inf[age] * entry.num_confirmed;
        }
        if (entry.date == offset_date_by_days(date, days_surplus)) {
            num_InfectedNoSymptoms[age] -= 1 / (1 - mu_C_R[age]) * scaling_factor_inf[age] * entry.num_confirmed;
        }
        if (entry.date == offset_date_by_days(date, t_Exposed[age] + t_InfectedNoSymptoms[age] + days_surplus)) {
            num_Exposed[age] += 1 / (1 - mu_C_R[age]) * scaling_factor_inf[age] * entry.num_confirmed;
        }
        if (entry.date == offset_date_by_days(date, -t_InfectedSymptoms[age])) {
            num_InfectedSymptoms[age] -= scaling_factor_inf[age] * entry.num_confirmed;
            num_InfectedSevere[age] += mu_I_H[age] * scaling_factor_inf[age] * entry.num_confirmed;
        }
        if (entry.date == offset_date_by_days(date, -t_InfectedSymptoms[age] - t_InfectedSevere[age])) {
            num_InfectedSevere[age] -= mu_I_H[age] * scaling_factor_inf[age] * entry.num_confirmed;
            num_icu[age] += mu_I_H[age] * mu_H_U[age] * scaling_factor_inf[age] * entry.num_confirmed;
        }
        if (entry.date ==
            offset_date_by_days(date, -t_InfectedSymptoms[age] - t_InfectedSevere[age] - t_InfectedCritical[age])) {
            num_death[age] += entry.num_deaths;
            num_icu[age] -= mu_I_H[age] * mu_H_U[age] * scaling_factor_inf[age] * entry.num_confirmed;
        }
    }

    for (size_t i = 0; i < ConfirmedCasesDataEntry::age_group_names.size(); i++) {
        // subtract infected confirmed cases which are not yet recovered
        // and remove dark figure scaling factor
        num_rec[i] -= num_InfectedSymptoms[i] / scaling_factor_inf[i];
        num_rec[i] -= num_InfectedSevere[i] / scaling_factor_inf[i];
        num_rec[i] -=
            num_icu[i] /
            scaling_factor_inf[i]; // TODO: this has to be adapted for scaling_factor_inf != 1 or != ***_icu
        num_rec[i] -= num_death[i] / scaling_factor_inf[i];
        auto try_fix_constraints = [region, i](double& value, double error, auto str) {
            if (value < error) {
                // this should probably return a failure
                // but the algorithm is not robust enough to avoid large negative
                // values and there are tests that rely on it
                log_error("{:s} for age group {:s} is {:.4f} for region {:d}, "
                            "exceeds expected negative value.",
                            str, ConfirmedCasesDataEntry::age_group_names[i], value, region);
                value = 0.0;
            }
            else if (value < 0) {
                log_info("{:s} for age group {:s} is {:.4f} for region {:d}, "
                            "automatically corrected",
                            str, ConfirmedCasesDataEntry::age_group_names[i], value, region);
                value = 0.0;
            }
        };

        try_fix_constraints(num_InfectedSymptoms[i], -5, "InfectedSymptoms");
        try_fix_constraints(num_InfectedNoSymptoms[i], -5, "InfectedNoSymptoms");
        try_fix_constraints(num_Exposed[i], -5, "Exposed");
        try_fix_constraints(num_InfectedSevere[i], -5, "InfectedSevere");
        try_fix_constraints(num_death[i], -5, "Dead");
        try_fix_constraints(num_icu[i], -5, "InfectedCritical");
        try_fix_constraints(num_rec[i], -20, "Recovered or vaccinated");
    }

    return success();
}

IOResult<std::vector<double>> compute_confirmed_cases_data_fix_recovered(
    const std::vector<ConfirmedCasesDataEntry>& case_data, const int region, Date date, double delay)
{
    auto max_date_entry = std::max_element(case_data.begin(), case_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == case_data.end()) {
        log_error("RKI data is empty.");
        return failure(StatusCode::InvalidValue, "RKI data is empty.");
    }
    auto max_date = max_date_entry->date;
    if (max_date < date) {
        log_error("Specified date does not exist in RKI data");
        return failure(StatusCode::OutOfRange, "RKI data does not contain specified date.");
    }

    // shifts the initilization to the recent past if simulation starts
    // around current day and data of the future would be required.
    // Only needed for preinfection compartments, exposed and InfectedNoSymptoms.
    auto days_surplus = get_offset_in_days(max_date, date) - 6; // 6 > T_E^C + T_C^I
    if (days_surplus > 0) {
        days_surplus = 0;
    }

    auto num_groups = ConfirmedCasesDataEntry::age_group_names.size();
    std::vector<double> num_rec(num_groups, 0.0);

    for (auto&& rki_entry : case_data) {
        if (rki_entry.date == offset_date_by_days(date, int(-delay))) {
            num_rec[size_t(rki_entry.age_group)] = rki_entry.num_confirmed;
        }
    }

    for (size_t i = 0; i < num_groups; i++) {
        auto try_fix_constraints = [region, i](double& value, double error, auto str) {
            if (value < error) {
                // this should probably return a failure
                // but the algorithm is not robust enough to avoid large negative
                // values and there are tests that rely on it
                log_error("{:s} for age group {:s} is {:.4f} for region {:d}, "
                            "exceeds expected negative value.",
                            str, ConfirmedCasesDataEntry::age_group_names[i], value, region);
                value = 0.0;
            }
            else if (value < 0) {
                log_info("{:s} for age group {:s} is {:.4f} for region {:d}, "
                            "automatically corrected",
                            str, ConfirmedCasesDataEntry::age_group_names[i], value, region);
                value = 0.0;
            }
        };
        try_fix_constraints(num_rec[i], 0, "Recovered");
    }

    return success(num_rec);
}

} // namespace details
} // namespace osecirvvs
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP
