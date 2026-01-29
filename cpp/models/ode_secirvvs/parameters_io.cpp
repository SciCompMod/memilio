/*
* Copyright (C) 2020-2026 MEmilio
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
IOResult<void> read_confirmed_cases_data(
    std::string const& path, std::vector<int> const& vregion, Date date, std::vector<std::vector<double>>& vnum_Exposed,
    std::vector<std::vector<double>>& vnum_InfectedNoSymptoms, std::vector<std::vector<double>>& vnum_InfectedSymptoms,
    std::vector<std::vector<double>>& vnum_InfectedSevere, std::vector<std::vector<double>>& vnum_icu,
    std::vector<std::vector<double>>& vnum_death, std::vector<std::vector<double>>& vnum_rec,
    const std::vector<std::vector<int>>& vt_Exposed, const std::vector<std::vector<int>>& vt_InfectedNoSymptoms,
    const std::vector<std::vector<int>>& vt_InfectedSymptoms, const std::vector<std::vector<int>>& vt_InfectedSevere,
    const std::vector<std::vector<int>>& vt_InfectedCritical, const std::vector<std::vector<double>>& vmu_C_R,
    const std::vector<std::vector<double>>& vmu_I_H, const std::vector<std::vector<double>>& vmu_H_U,
    const std::vector<double>& scaling_factor_inf)
{
    BOOST_OUTCOME_TRY(auto&& rki_data, mio::read_confirmed_cases_data(path));
    return read_confirmed_cases_data(rki_data, vregion, date, vnum_Exposed, vnum_InfectedNoSymptoms,
                                     vnum_InfectedSymptoms, vnum_InfectedSevere, vnum_icu, vnum_death, vnum_rec,
                                     vt_Exposed, vt_InfectedNoSymptoms, vt_InfectedSymptoms, vt_InfectedSevere,
                                     vt_InfectedCritical, vmu_C_R, vmu_I_H, vmu_H_U, scaling_factor_inf);
}

IOResult<void> read_confirmed_cases_data(
    const std::vector<ConfirmedCasesDataEntry>& rki_data, std::vector<int> const& vregion, Date date,
    std::vector<std::vector<double>>& vnum_Exposed, std::vector<std::vector<double>>& vnum_InfectedNoSymptoms,
    std::vector<std::vector<double>>& vnum_InfectedSymptoms, std::vector<std::vector<double>>& vnum_InfectedSevere,
    std::vector<std::vector<double>>& vnum_icu, std::vector<std::vector<double>>& vnum_death,
    std::vector<std::vector<double>>& vnum_rec, const std::vector<std::vector<int>>& vt_Exposed,
    const std::vector<std::vector<int>>& vt_InfectedNoSymptoms,
    const std::vector<std::vector<int>>& vt_InfectedSymptoms, const std::vector<std::vector<int>>& vt_InfectedSevere,
    const std::vector<std::vector<int>>& vt_InfectedCritical, const std::vector<std::vector<double>>& vmu_C_R,
    const std::vector<std::vector<double>>& vmu_I_H, const std::vector<std::vector<double>>& vmu_H_U,
    const std::vector<double>& scaling_factor_inf)
{
    auto max_date_entry = std::max_element(rki_data.begin(), rki_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == rki_data.end()) {
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

    for (auto&& entry : rki_data) {
        auto it = std::find_if(vregion.begin(), vregion.end(), [&entry](auto r) {
            return r == 0 || get_region_id(entry) == r;
        });
        if (it != vregion.end()) {
            auto region_idx = size_t(it - vregion.begin());

            auto& t_Exposed            = vt_Exposed[region_idx];
            auto& t_InfectedNoSymptoms = vt_InfectedNoSymptoms[region_idx];
            auto& t_InfectedSymptoms   = vt_InfectedSymptoms[region_idx];
            auto& t_InfectedSevere     = vt_InfectedSevere[region_idx];
            auto& t_InfectedCritical   = vt_InfectedCritical[region_idx];

            auto& num_InfectedNoSymptoms = vnum_InfectedNoSymptoms[region_idx];
            auto& num_InfectedSymptoms   = vnum_InfectedSymptoms[region_idx];
            auto& num_rec                = vnum_rec[region_idx];
            auto& num_Exposed            = vnum_Exposed[region_idx];
            auto& num_InfectedSevere     = vnum_InfectedSevere[region_idx];
            auto& num_death              = vnum_death[region_idx];
            auto& num_icu                = vnum_icu[region_idx];

            auto& mu_C_R = vmu_C_R[region_idx];
            auto& mu_I_H = vmu_I_H[region_idx];
            auto& mu_H_U = vmu_H_U[region_idx];

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
    }

    for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
        auto region = vregion[region_idx];

        auto& num_InfectedNoSymptoms = vnum_InfectedNoSymptoms[region_idx];
        auto& num_InfectedSymptoms   = vnum_InfectedSymptoms[region_idx];
        auto& num_rec                = vnum_rec[region_idx];
        auto& num_Exposed            = vnum_Exposed[region_idx];
        auto& num_InfectedSevere     = vnum_InfectedSevere[region_idx];
        auto& num_death              = vnum_death[region_idx];
        auto& num_icu                = vnum_icu[region_idx];

        size_t num_groups = ConfirmedCasesDataEntry::age_group_names.size();
        for (size_t i = 0; i < num_groups; i++) {
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
    }

    return success();
}

IOResult<void> read_confirmed_cases_data_fix_recovered(std::string const& path, std::vector<int> const& vregion,
                                                       Date date, std::vector<std::vector<double>>& vnum_rec,
                                                       double delay)
{
    BOOST_OUTCOME_TRY(auto&& rki_data, mio::read_confirmed_cases_data(path));
    return read_confirmed_cases_data_fix_recovered(rki_data, vregion, date, vnum_rec, delay);
}

IOResult<void> read_confirmed_cases_data_fix_recovered(const std::vector<ConfirmedCasesDataEntry>& rki_data,
                                                       std::vector<int> const& vregion, Date date,
                                                       std::vector<std::vector<double>>& vnum_rec, double delay)
{
    auto max_date_entry = std::max_element(rki_data.begin(), rki_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == rki_data.end()) {
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

    for (auto&& rki_entry : rki_data) {
        auto it = std::find_if(vregion.begin(), vregion.end(), [&rki_entry](auto r) {
            return r == 0 || get_region_id(rki_entry) == r;
        });
        if (it != vregion.end()) {
            auto region_idx = size_t(it - vregion.begin());
            if (rki_entry.date == offset_date_by_days(date, int(-delay))) {
                vnum_rec[region_idx][size_t(rki_entry.age_group)] = rki_entry.num_confirmed;
            }
        }
    }

    for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
        auto region   = vregion[region_idx];
        auto& num_rec = vnum_rec[region_idx];

        size_t num_groups = ConfirmedCasesDataEntry::age_group_names.size();
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
    }

    return success();
}

} // namespace details
} // namespace osecirvvs
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP
