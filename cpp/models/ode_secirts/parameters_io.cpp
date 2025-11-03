/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker, Wadim Koslow, Daniel Abele, Martin J. Kühn
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

#include "ode_secirts/parameters_io.h"
#include "memilio/geography/regions.h"
#include "memilio/io/io.h"
#include "ode_secirts/parameters.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "memilio/io/epi_data.h"
#include "memilio/utils/memory.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/utils/stl_util.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/date.h"

#include <boost/filesystem.hpp>

#include <numeric>
#include <vector>
#include <iostream>
#include <string>
#include <random>
#include <fstream>

namespace mio
{
namespace osecirts
{
namespace details
{

IOResult<void> compute_confirmed_cases_data(
    const std::vector<ConfirmedCasesDataEntry>& case_data, const int region, Date date,
    std::vector<double>& num_Exposed, std::vector<double>& num_InfectedNoSymptoms,
    std::vector<double>& num_InfectedSymptoms, std::vector<double>& num_InfectedSevere,
    std::vector<double>& num_icu, std::vector<double>& num_death,
    std::vector<double>& num_imm, const std::vector<int>& t_Exposed,
    const std::vector<int>& t_InfectedNoSymptoms, const std::vector<int>& t_InfectedSymptoms, 
    const std::vector<int>& t_InfectedSevere, const std::vector<int>& t_InfectedCritical, 
    const std::vector<int>& t_imm_interval_i, const std::vector<double>& mu_C_R, 
    const std::vector<double>& mu_I_H, const std::vector<double>& mu_H_U,
    const std::vector<double>& reduc_t_Infected, const std::vector<double>& reduc_Exposed,
    const std::vector<double>& reduc_InfectedSymptoms, const std::vector<double>& reduc_icu_death,
    const std::vector<double>& scaling_factor_inf, const size_t layer)
{
    auto max_date_entry = std::max_element(case_data.begin(), case_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == case_data.end()) {
        log_error("Case data file is empty.");
        return failure(StatusCode::InvalidValue, "Case data is empty.");
    }
    auto max_date = max_date_entry->date;
    if (max_date < date) {
        log_error("Specified date does not exist in case data");
        return failure(StatusCode::OutOfRange, "Case data does not contain specified date.");
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

        // transition probabilities
        double recoveredPerInfectedNoSymptoms = mu_C_R[age];
        double severePerInfectedSymptoms = mu_I_H[age];
        double criticalPerSevere         = mu_H_U[age];

        // transition times
        int timeInfectedNoSymptoms = t_InfectedNoSymptoms[age];
        int timeInfectedSymptoms = t_InfectedSymptoms[age];

        // if we select a layer with better immunity (layer > 0), we need to adjust the times and transition rates
        if (layer > 0) {
            timeInfectedNoSymptoms = static_cast<int>(std::round(timeInfectedNoSymptoms * reduc_t_Infected[age]));
            timeInfectedSymptoms   = static_cast<int>(std::round(timeInfectedSymptoms * reduc_t_Infected[age]));

            const double red_fact_exp = reduc_Exposed[age];

            const double red_fact_inf = reduc_InfectedSymptoms[age];

            const double red_fact_sev = reduc_icu_death[age];

            recoveredPerInfectedNoSymptoms = 1 - red_fact_inf / red_fact_exp * (1 - recoveredPerInfectedNoSymptoms);
            severePerInfectedSymptoms      = red_fact_sev / red_fact_inf * severePerInfectedSymptoms;
        }

        if (entry.date == offset_date_by_days(date, 0)) {
            num_InfectedSymptoms[age] += scaling_factor_inf[age] * entry.num_confirmed;
            num_imm[age] += entry.num_confirmed;
        }
        if (entry.date == offset_date_by_days(date, timeInfectedNoSymptoms + days_surplus)) {
            num_InfectedNoSymptoms[age] +=
                1 / (1 - recoveredPerInfectedNoSymptoms) * scaling_factor_inf[age] * entry.num_confirmed;
            num_Exposed[age] -=
                1 / (1 - recoveredPerInfectedNoSymptoms) * scaling_factor_inf[age] * entry.num_confirmed;
        }
        if (entry.date == offset_date_by_days(date, days_surplus)) {
            num_InfectedNoSymptoms[age] -=
                1 / (1 - recoveredPerInfectedNoSymptoms) * scaling_factor_inf[age] * entry.num_confirmed;
        }
        if (entry.date == offset_date_by_days(date, t_Exposed[age] + timeInfectedNoSymptoms + days_surplus)) {
            num_Exposed[age] +=
                1 / (1 - recoveredPerInfectedNoSymptoms) * scaling_factor_inf[age] * entry.num_confirmed;
        }
        if (entry.date == offset_date_by_days(date, -timeInfectedSymptoms)) {
            num_InfectedSymptoms[age] -= scaling_factor_inf[age] * entry.num_confirmed;
            num_InfectedSevere[age] += severePerInfectedSymptoms * scaling_factor_inf[age] * entry.num_confirmed;
        }
        if (entry.date == offset_date_by_days(date, -timeInfectedSymptoms - t_InfectedSevere[age])) {
            num_InfectedSevere[age] -= severePerInfectedSymptoms * scaling_factor_inf[age] * entry.num_confirmed;
            num_icu[age] +=
                severePerInfectedSymptoms * criticalPerSevere * scaling_factor_inf[age] * entry.num_confirmed;
        }
        if (entry.date == offset_date_by_days(date, -timeInfectedSymptoms - t_InfectedSevere[age] - t_InfectedCritical[age])) {
            num_death[age] += entry.num_deaths;
            num_icu[age] -=
                severePerInfectedSymptoms * criticalPerSevere * scaling_factor_inf[age] * entry.num_confirmed;
        }
        if (entry.date == offset_date_by_days(date, 0 - t_imm_interval_i[age])) {
            num_imm[age] -= entry.num_confirmed;
        }
    }

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

        const double tol_error = -1e-8;
        try_fix_constraints(num_InfectedSymptoms[i], tol_error, "InfectedSymptoms");
        try_fix_constraints(num_InfectedNoSymptoms[i], tol_error, "InfectedNoSymptoms");
        try_fix_constraints(num_Exposed[i], tol_error, "Exposed");
        try_fix_constraints(num_InfectedSevere[i], tol_error, "InfectedSevere");
        try_fix_constraints(num_death[i], tol_error, "Dead");
        try_fix_constraints(num_icu[i], tol_error, "InfectedCritical");
        try_fix_constraints(num_timm_i[i], tol_error, "Recently Recovered or Vaccinated");
    }

    return success();
}

} // namespace details
} // namespace osecirts
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP
