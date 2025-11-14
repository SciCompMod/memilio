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
    std::vector<ScalarType>& num_Exposed, std::vector<ScalarType>& num_InfectedNoSymptoms,
    std::vector<ScalarType>& num_InfectedSymptoms, std::vector<ScalarType>& num_InfectedSevere,
    std::vector<ScalarType>& num_icu, std::vector<ScalarType>& num_death,
    std::vector<ScalarType>& num_imm, const std::vector<int>& t_Exposed,
    const std::vector<int>& t_InfectedNoSymptoms, const std::vector<int>& t_InfectedSymptoms, 
    const std::vector<int>& t_InfectedSevere, const std::vector<int>& t_InfectedCritical, 
    const std::vector<int>& t_imm_interval_i, const std::vector<ScalarType>& mu_C_R, 
    const std::vector<ScalarType>& mu_I_H, const std::vector<ScalarType>& mu_H_U,
    const std::vector<ScalarType>& reduc_t_Infected, const std::vector<ScalarType>& reduc_Exposed,
    const std::vector<ScalarType>& reduc_InfectedSymptoms, const std::vector<ScalarType>& reduc_icu_death,
    const std::vector<ScalarType>& scaling_factor_inf, const size_t layer)
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
        ScalarType recoveredPerInfectedNoSymptoms = mu_C_R[age];
        ScalarType severePerInfectedSymptoms = mu_I_H[age];
        ScalarType criticalPerSevere         = mu_H_U[age];

        // transition times
        int timeInfectedNoSymptoms = t_InfectedNoSymptoms[age];
        int timeInfectedSymptoms = t_InfectedSymptoms[age];

        // if we select a layer with better immunity (layer > 0), we need to adjust the times and transition rates
        if (layer > 0) {
            timeInfectedNoSymptoms = static_cast<int>(std::round(timeInfectedNoSymptoms * reduc_t_Infected[age]));
            timeInfectedSymptoms   = static_cast<int>(std::round(timeInfectedSymptoms * reduc_t_Infected[age]));

            const ScalarType red_fact_exp = reduc_Exposed[age];

            const ScalarType red_fact_inf = reduc_InfectedSymptoms[age];

            const ScalarType red_fact_sev = reduc_icu_death[age];

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
        auto try_fix_constraints = [region, i](ScalarType& value, ScalarType error, auto str) {
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

        const ScalarType tol_error = -1e-8;
        try_fix_constraints(num_InfectedSymptoms[i], tol_error, "InfectedSymptoms");
        try_fix_constraints(num_InfectedNoSymptoms[i], tol_error, "InfectedNoSymptoms");
        try_fix_constraints(num_Exposed[i], tol_error, "Exposed");
        try_fix_constraints(num_InfectedSevere[i], tol_error, "InfectedSevere");
        try_fix_constraints(num_death[i], tol_error, "Dead");
        try_fix_constraints(num_icu[i], tol_error, "InfectedCritical");
        try_fix_constraints(num_imm[i], tol_error, "Recently Recovered or Vaccinated");
    }

    return success();
}

IOResult<void>
set_confirmed_cases_data(Model<ScalarType>& model, const std::vector<ConfirmedCasesDataEntry>& case_data,
                         const int region, Date date, const std::vector<ScalarType>& scaling_factor_inf,
                         const std::vector<std::vector<ScalarType>>& immunity_population)
{
    auto num_age_groups = (size_t)model.parameters.get_num_groups();
    assert(scaling_factor_inf.size() == num_age_groups); //TODO: allow vector or scalar valued scaling factors
    assert(ConfirmedCasesDataEntry::age_group_names.size() == num_age_groups);

    std::vector<int> t_Exposed;
    std::vector<int> t_InfectedNoSymptoms;
    std::vector<int> t_InfectedSymptoms;
    std::vector<int> t_InfectedSevere;
    std::vector<int> t_InfectedCritical;
    std::vector<int> t_imm_interval_i;

    std::vector<ScalarType> mu_C_R;
    std::vector<ScalarType> mu_I_H;
    std::vector<ScalarType> mu_H_U;

    std::vector<ScalarType> reduc_t_Infected;
    std::vector<ScalarType> reduc_Exposed;
    std::vector<ScalarType> reduc_InfectedSymptoms;
    std::vector<ScalarType> reduc_icu_death;

    std::vector<ScalarType> num_InfectedSymptoms(num_age_groups, 0.0);
    std::vector<ScalarType> num_death(num_age_groups, 0.0);
    std::vector<ScalarType> num_Exposed(num_age_groups, 0.0);
    std::vector<ScalarType> num_InfectedNoSymptoms(num_age_groups, 0.0);
    std::vector<ScalarType> num_InfectedSevere(num_age_groups, 0.0);
    std::vector<ScalarType> num_icu(num_age_groups, 0.0);
    std::vector<ScalarType> num_timm1(num_age_groups, 0.0);
    std::vector<ScalarType> num_timm2(num_age_groups, 0.0);

    std::vector<ScalarType> denom_E(num_age_groups, 0.0);
    std::vector<ScalarType> denom_I_NS(num_age_groups, 0.0);
    std::vector<ScalarType> denom_I_Sy(num_age_groups, 0.0);
    std::vector<ScalarType> denom_I_Sev_Cr(num_age_groups, 0.0);

    // calculate the denominators to split the reported case numbers to the different immunity layers.
    for (size_t group = 0; group < num_age_groups; group++) {
        denom_E[group] =
            1 / (immunity_population[0][group] +
                    immunity_population[1][group] *
                        model.parameters.template get<ReducExposedPartialImmunity<ScalarType>>()[(AgeGroup)group] +
                    immunity_population[2][group] *
                        model.parameters.template get<ReducExposedImprovedImmunity<ScalarType>>()[(AgeGroup)group]);

        denom_I_NS[group] =
            1 / (immunity_population[0][group] +
                    immunity_population[1][group] *
                        model.parameters.template get<ReducExposedPartialImmunity<ScalarType>>()[(AgeGroup)group] +
                    immunity_population[2][group] *
                        model.parameters.template get<ReducExposedImprovedImmunity<ScalarType>>()[(AgeGroup)group]);

        denom_I_Sy[group] =
            1 / (immunity_population[0][group] +
                    immunity_population[1][group] *
                        model
                            .parameters.template get<ReducInfectedSymptomsPartialImmunity<ScalarType>>()[(AgeGroup)group] +
                    immunity_population[2][group] *
                        model
                            .parameters.template get<ReducInfectedSymptomsImprovedImmunity<ScalarType>>()[(AgeGroup)group]);

        denom_I_Sev_Cr[group] =
            1 / (immunity_population[0][group] +
                    immunity_population[1][group] *
                        model.parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<ScalarType>>()[(
                            AgeGroup)group] +
                    immunity_population[2][group] *
                        model.parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<ScalarType>>()[(
                            AgeGroup)group]);
    }

    /*----------- Naive immunity -----------*/
    for (size_t group = 0; group < num_age_groups; group++) {
        t_Exposed.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeExposed<ScalarType>>()[(AgeGroup)group])));
        t_InfectedNoSymptoms.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedNoSymptoms<ScalarType>>()[(AgeGroup)group])));
        t_InfectedSymptoms.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedSymptoms<ScalarType>>()[(AgeGroup)group])));
        t_InfectedSevere.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedSevere<ScalarType>>()[(AgeGroup)group])));
        t_InfectedCritical.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedCritical<ScalarType>>()[(AgeGroup)group])));
        t_imm_interval_i.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeTemporaryImmunityPI<ScalarType>>()[(AgeGroup)group])));

        mu_C_R.push_back(
            model.parameters.template get<RecoveredPerInfectedNoSymptoms<ScalarType>>()[(AgeGroup)group]);
        mu_I_H.push_back(
            model.parameters.template get<SeverePerInfectedSymptoms<ScalarType>>()[(AgeGroup)group]);
        mu_H_U.push_back(
            model.parameters.template get<CriticalPerSevere<ScalarType>>()[(AgeGroup)group]);

        reduc_t_Infected.push_back(
            model.parameters.template get<ReducTimeInfectedMild<ScalarType>>()[(AgeGroup)group]);
        reduc_Exposed.push_back(
            model.parameters.template get<ReducExposedPartialImmunity<ScalarType>>()[(AgeGroup)group]);
        reduc_InfectedSymptoms.push_back(
            model.parameters.template get<ReducInfectedSymptomsPartialImmunity<ScalarType>>()[(AgeGroup)group]);
        reduc_icu_death.push_back(
            model.parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<ScalarType>>()[(AgeGroup)group]);
    }

    BOOST_OUTCOME_TRY(compute_confirmed_cases_data(case_data, region, date, num_Exposed, num_InfectedNoSymptoms,
                                                   num_InfectedSymptoms, num_InfectedSevere, num_icu, num_death, num_timm1,
                                                   t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
                                                   t_InfectedCritical, t_imm_interval_i, reduc_t_Infected, reduc_Exposed,
                                                   reduc_InfectedSymptoms, reduc_icu_death, mu_C_R, mu_I_H, mu_H_U, 
                                                   scaling_factor_inf, 0));

    size_t num_groups = (size_t)model.parameters.get_num_groups();
    for (size_t i = 0; i < num_groups; i++) {
        model.populations[{AgeGroup(i), InfectionState::ExposedNaive}] =
            immunity_population[0][i] * denom_E[i] * num_Exposed[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsNaive}] =
            immunity_population[0][i] * denom_I_NS[i] * num_InfectedNoSymptoms[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsNaiveConfirmed}] = 0;
        model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsNaive}] =
            immunity_population[0][i] * denom_I_Sy[i] * num_InfectedSymptoms[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsNaiveConfirmed}] = 0;
        model.populations[{AgeGroup(i), InfectionState::InfectedSevereNaive}] =
            immunity_population[0][i] * denom_I_Sev_Cr[i] * num_InfectedSevere[i];
        // Only set the number of ICU patients here, if the date is not available in the data.
        if (!is_divi_data_available(date)) {
            model.populations[{AgeGroup(i), InfectionState::InfectedCriticalNaive}] =
                immunity_population[0][i] * denom_I_Sev_Cr[i] * num_icu[i];
        }
    }
    if (std::accumulate(num_InfectedSymptoms.begin(), num_InfectedSymptoms.end(), ScalarType(0.0),
                        [](const ScalarType& a, const ScalarType& b) {
                            return evaluate_intermediate<ScalarType>(a + b);
                        }) == 0) {
        log_warning(
            "No infections for unvaccinated reported on date {} for region {}. Population data has not been set.",
            date, region);
    }

    /*----------- PARTIAL Immunity -----------*/
    reduc_Exposed.clear();
    reduc_InfectedSymptoms.clear();
    reduc_icu_death.clear();

    num_InfectedSymptoms   = std::vector<ScalarType>(num_age_groups, 0.0);
    num_death              = std::vector<ScalarType>(num_age_groups, 0.0);
    num_Exposed            = std::vector<ScalarType>(num_age_groups, 0.0);
    num_InfectedNoSymptoms = std::vector<ScalarType>(num_age_groups, 0.0);
    num_InfectedSevere     = std::vector<ScalarType>(num_age_groups, 0.0);
    num_icu                = std::vector<ScalarType>(num_age_groups, 0.0);

    for (size_t group = 0; group < num_age_groups; group++) {
        reduc_Exposed.push_back(
            model.parameters.template get<ReducExposedPartialImmunity<ScalarType>>()[(AgeGroup)group]);
        reduc_InfectedSymptoms.push_back(
            model.parameters.template get<ReducInfectedSymptomsPartialImmunity<ScalarType>>()[(AgeGroup)group]);
        reduc_icu_death.push_back(
            model.parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<ScalarType>>()[(AgeGroup)group]);
    }

    BOOST_OUTCOME_TRY(compute_confirmed_cases_data(case_data, region, date, num_Exposed, num_InfectedNoSymptoms,
                                                   num_InfectedSymptoms, num_InfectedSevere, num_icu, num_death, num_timm1,
                                                   t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
                                                   t_InfectedCritical, t_imm_interval_i, reduc_t_Infected, reduc_Exposed,
                                                   reduc_InfectedSymptoms, reduc_icu_death, mu_C_R, mu_I_H, mu_H_U, 
                                                   scaling_factor_inf, 1));
    for (size_t i = 0; i < num_groups; i++) {
        model.populations[{AgeGroup(i), InfectionState::ExposedPartialImmunity}] =
            immunity_population[1][i] *
            model.parameters.template get<ReducExposedPartialImmunity<ScalarType>>()[(AgeGroup)i] * denom_E[i] *
            num_Exposed[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsPartialImmunity}] =
            immunity_population[1][i] *
            model.parameters.template get<ReducExposedPartialImmunity<ScalarType>>()[(AgeGroup)i] * denom_I_NS[i] *
            num_InfectedNoSymptoms[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}] = 0;
        model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsPartialImmunity}] =
            immunity_population[1][i] *
            model.parameters.template get<ReducInfectedSymptomsPartialImmunity<ScalarType>>()[(AgeGroup)i] *
            denom_I_Sy[i] * num_InfectedSymptoms[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsPartialImmunityConfirmed}] = 0;
        model.populations[{AgeGroup(i), InfectionState::InfectedSeverePartialImmunity}] =
            immunity_population[1][i] *
            model
                .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<ScalarType>>()[(AgeGroup)i] *
            denom_I_Sev_Cr[i] * num_InfectedSevere[i];
        // Only set the number of ICU patients here, if the date is not available in the data.
        if (!is_divi_data_available(date)) {
            model.populations[{AgeGroup(i), InfectionState::InfectedCriticalPartialImmunity}] =
                immunity_population[1][i] *
                model
                    .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<ScalarType>>()[(AgeGroup)i] *
                denom_I_Sev_Cr[i] * num_icu[i];
        }
        // the += is necessary because we already set the previous vaccinated individuals
        model.populations[{AgeGroup(i), InfectionState::TemporaryImmunePartialImmunity}] +=
            immunity_population[1][i] *
            model.parameters.template get<ReducExposedPartialImmunity<ScalarType>>()[(AgeGroup)i] * denom_E[i] *
            num_timm1[i];
    }
    if (std::accumulate(num_InfectedSymptoms.begin(), num_InfectedSymptoms.end(), ScalarType(0.0),
                        [](const ScalarType& a, const ScalarType& b) {
                            return evaluate_intermediate<ScalarType>(a + b);
                        }) == 0) {

        log_warning("No infections for partially vaccinated reported on date {} for region {}. "
                    "Population data has not been set.",
                    date, region);
    }

    /*----------- Improved Immunity -----------*/
    reduc_Exposed.clear();
    reduc_InfectedSymptoms.clear();
    reduc_icu_death.clear();

    num_InfectedSymptoms   = std::vector<ScalarType>(num_age_groups, 0.0);
    num_death              = std::vector<ScalarType>(num_age_groups, 0.0);
    num_Exposed            = std::vector<ScalarType>(num_age_groups, 0.0);
    num_InfectedNoSymptoms = std::vector<ScalarType>(num_age_groups, 0.0);
    num_InfectedSevere     = std::vector<ScalarType>(num_age_groups, 0.0);
    num_icu                = std::vector<ScalarType>(num_age_groups, 0.0);

    for (size_t group = 0; group < num_age_groups; group++) {
        reduc_Exposed.push_back(
            model.parameters.template get<ReducExposedImprovedImmunity<ScalarType>>()[(AgeGroup)group]);
        reduc_InfectedSymptoms.push_back(
            model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<ScalarType>>()[(AgeGroup)group]);
        reduc_icu_death.push_back(
            model.parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<ScalarType>>()[(AgeGroup)group]);
    }

    BOOST_OUTCOME_TRY(compute_confirmed_cases_data(case_data, region, date, num_Exposed, num_InfectedNoSymptoms,
                                                   num_InfectedSymptoms, num_InfectedSevere, num_icu, num_death, num_timm2,
                                                   t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
                                                   t_InfectedCritical, t_imm_interval_i, reduc_t_Infected, reduc_Exposed,
                                                   reduc_InfectedSymptoms, reduc_icu_death, mu_C_R, mu_I_H, mu_H_U, 
                                                   scaling_factor_inf, 2));

    for (size_t i = 0; i < num_groups; i++) {
        model.populations[{AgeGroup(i), InfectionState::ExposedImprovedImmunity}] =
            immunity_population[2][i] *
            model.parameters.template get<ReducExposedImprovedImmunity<ScalarType>>()[(AgeGroup)i] * denom_E[i] *
            num_Exposed[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsImprovedImmunity}] =
            immunity_population[2][i] *
            model.parameters.template get<ReducExposedImprovedImmunity<ScalarType>>()[(AgeGroup)i] * denom_I_NS[i] *
            num_InfectedNoSymptoms[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] = 0;
        model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsImprovedImmunity}] =
            immunity_population[2][i] *
            model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<ScalarType>>()[(AgeGroup)i] *
            denom_I_Sy[i] * num_InfectedSymptoms[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] = 0;
        model.populations[{AgeGroup(i), InfectionState::InfectedSevereImprovedImmunity}] =
            immunity_population[2][i] *
            model
                .parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<ScalarType>>()[(AgeGroup)i] *
            denom_I_Sev_Cr[i] * num_InfectedSevere[i];
        // Only set the number of ICU patients here, if the date is not available in the data.
        if (!is_divi_data_available(date)) {
            model.populations[{AgeGroup(i), InfectionState::InfectedCriticalImprovedImmunity}] =
                immunity_population[2][i] *
                model
                    .parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<ScalarType>>()[(AgeGroup)i] *
                denom_I_Sev_Cr[i] * num_icu[i];
        }

        // the += is necessary because we already set the previous vaccinated individuals
        model.populations[{AgeGroup(i), InfectionState::TemporaryImmuneImprovedImmunity}] +=
            immunity_population[2][i] *
            model.parameters.template get<ReducExposedImprovedImmunity<ScalarType>>()[(AgeGroup)i] * denom_E[i] *
            num_timm2[i];
    }
    if (std::accumulate(num_InfectedSymptoms.begin(), num_InfectedSymptoms.end(), ScalarType(0.0),
                        [](const ScalarType& a, const ScalarType& b) {
                            return evaluate_intermediate<ScalarType>(a + b);
                        }) == 0) {
        log_warning("No infections for vaccinated reported on date {} for region {}. "
                    "Population data has not been set.",
                    date, region);
    }
    return success();
}

IOResult<void> set_confirmed_cases_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, const std::string& path,
                                        Date date, const std::vector<ScalarType>& scaling_factor_inf,
                                        const std::vector<std::vector<ScalarType>>& immunity_population)
{
    BOOST_OUTCOME_TRY(auto&& case_data, mio::read_confirmed_cases_data(path));

    // sort case_data into regions and ignore once with no region associated
    std::vector<std::vector<ConfirmedCasesDataEntry>> vcase_data{model.size()};
    for (auto&& entry : case_data) {
        auto it = std::find_if(model.begin(), model.end(), [&entry](auto m) {
            return m.id == 0 || get_region_id(entry) == m.id;
        });
        if (it != model.end()) {
            auto region_idx = size_t(it - model.begin());
            vcase_data[region_idx].push_back(entry);
        }
    }

    for (size_t region_idx = 0; region_idx < model.size(); ++region_idx) {
        BOOST_OUTCOME_TRY(
            set_confirmed_cases_data(model[region_idx].property, vcase_data[region_idx], model[region_idx].id, date, scaling_factor_inf, immunity_population));
    }

    return success();
}

IOResult<void> set_population_data(Model<ScalarType>& model, const std::vector<ScalarType>& num_population,
                                   const int region, const std::vector<std::vector<ScalarType>>& immunity_population)
{
    if (std::accumulate(num_population.begin(), num_population.end(), ScalarType(0.0),
                        [](const ScalarType& a, const ScalarType& b) {
                            return evaluate_intermediate<ScalarType>(a + b);
                        }) <= 0)
    {    
        log_warning("No population data available for region " + std::to_string(region) +
                    ". Population data has not been set.");
        return success();
    }

    auto num_groups = model.parameters.get_num_groups();
    for (auto i = AgeGroup(0); i < num_groups; i++) {

        ScalarType SN  = num_population[size_t(i)] * immunity_population[0][size_t(i)];
        ScalarType SPI = num_population[size_t(i)] * immunity_population[1][size_t(i)];
        ScalarType SII = num_population[size_t(i)] - SN - SPI;

        model.populations[{i, InfectionState::SusceptibleImprovedImmunity}] = std::max(
            0.0,
            SII -
                (model.populations[{i, InfectionState::ExposedImprovedImmunity}] +
                model.populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunity}] +
                model.populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] +
                model.populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] +
                model.populations[{i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] +
                model.populations[{i, InfectionState::InfectedSevereImprovedImmunity}] +
                model.populations[{i, InfectionState::InfectedCriticalImprovedImmunity}] +
                model.populations[{i, InfectionState::DeadImprovedImmunity}] +
                model.populations[{i, InfectionState::TemporaryImmuneImprovedImmunity}]));

        model.populations[{i, InfectionState::SusceptiblePartialImmunity}] = std::max(
            0.0,
            SPI - model.populations[{i, InfectionState::ExposedPartialImmunity}] -
                model.populations[{i, InfectionState::InfectedNoSymptomsPartialImmunity}] -
                model.populations[{i, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}] -
                model.populations[{i, InfectionState::InfectedSymptomsPartialImmunity}] -
                model.populations[{i, InfectionState::InfectedSymptomsPartialImmunityConfirmed}] -
                model.populations[{i, InfectionState::InfectedSeverePartialImmunity}] -
                model.populations[{i, InfectionState::InfectedCriticalPartialImmunity}] -
                model.populations[{i, InfectionState::DeadPartialImmunity}] -
                model.populations[{i, InfectionState::TemporaryImmunePartialImmunity}]);

        model.populations.template set_difference_from_group_total<AgeGroup>(
            {i, InfectionState::SusceptibleNaive}, num_population[size_t(i)]);
    }

    for (auto i = AgeGroup(0); i < AgeGroup(6); i++) {
        for (auto j = Index<InfectionState>(0); j < InfectionState::Count; ++j) {
            if (model.populations[{i, j}] < 0) {
                log_warning("Compartment at age group {}, infection state {}, is negative: {}", size_t(i),
                            size_t(j), model.populations[{i, j}]);
            }
        }
    }

    return success();
}

IOResult<void> set_population_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, const std::string& path,
                                   const std::vector<std::vector<ScalarType>>& immunity_population)
{
    std::vector<int> vregion; 
    std::transform(model.begin(), model.end(), std::back_inserter(vregion), [](const auto& m) { return m.id; });
    BOOST_OUTCOME_TRY(auto&& num_population, mio::read_population_data(path, vregion));

    for (size_t region_idx = 0; region_idx < model.size(); ++region_idx) {
        BOOST_OUTCOME_TRY(set_population_data(model[region_idx].property, num_population[region_idx], model[region_idx].id, immunity_population));
    }
    return success();
}

IOResult<void> set_divi_data(Model<ScalarType>& model, const ScalarType num_icu, ScalarType scaling_factor_icu)
{
    ScalarType sum_mu_I_U = 0;
    std::vector<ScalarType> mu_I_U;
    auto num_groups = model.parameters.get_num_groups();
    for (auto i = AgeGroup(0); i < num_groups; i++) {
        sum_mu_I_U += model.parameters.template get<CriticalPerSevere<ScalarType>>()[i] *
                                model.parameters.template get<SeverePerInfectedSymptoms<ScalarType>>()[i];
        mu_I_U.push_back(model.parameters.template get<CriticalPerSevere<ScalarType>>()[i] *
                                    model.parameters.template get<SeverePerInfectedSymptoms<ScalarType>>()[i]);
    }

    for (auto i = AgeGroup(0); i < num_groups; i++) {
        model.populations[{i, InfectionState::InfectedCriticalNaive}] =
            scaling_factor_icu * num_icu * mu_I_U[(size_t)i] / sum_mu_I_U;
    }

    return success();
}

IOResult<void> set_divi_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, const std::string& path, Date date, 
                             ScalarType scaling_factor_icu)
{
    // DIVI dataset will no longer be updated from CW29 2024 on.
    if (!is_divi_data_available(date)) {
        log_warning("No DIVI data available for date: {}. "
                    "ICU compartment will be set based on Case data.",
                    date);
        return success();
    }

    std::vector<int> vregion; 
    std::transform(model.begin(), model.end(), std::back_inserter(vregion), [](const auto& m) { return m.id; });
    BOOST_OUTCOME_TRY(auto&& num_icu, read_divi_data(path, vregion, date));

    for (size_t region_idx = 0; region_idx < model.size(); ++region_idx) {
        BOOST_OUTCOME_TRY(set_divi_data(model[region_idx].property, num_icu[region_idx], scaling_factor_icu));
    }

    return success();
}

} // namespace details

IOResult<void> read_input_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, Date date,
                               const std::vector<ScalarType>& scaling_factor_inf, ScalarType scaling_factor_icu,
                               int num_days, const std::vector<std::vector<ScalarType>>& immunity_population, 
                               const mio::regions::de::EpidataFilenames& epidata_filenames)
{

    BOOST_OUTCOME_TRY(
        details::set_vaccination_data<ScalarType>(model, epidata_filenames.vaccination_data_path, date, num_days));

    // TODO: Reuse more code, e.g., set_divi_data (in secir) and a set_divi_data (here) only need a different ModelType.
    // TODO: add option to set ICU data from confirmed cases if DIVI or other data is not available.
    BOOST_OUTCOME_TRY(details::set_divi_data(model, epidata_filenames.divi_data_path, date, scaling_factor_icu));

    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(model, epidata_filenames.case_data_path, date, 
                                                        scaling_factor_inf, immunity_population));
    BOOST_OUTCOME_TRY(details::set_population_data(model, epidata_filenames.population_data_path, immunity_population));
    return success();
}

#ifdef MEMILIO_HAS_HDF5

IOResult<void> export_input_data_timeseries(
    const mio::VectorRange<Node<Model<ScalarType>>> model, const std::string& results_dir, Date date,
    const std::vector<ScalarType>& scaling_factor_inf, const ScalarType scaling_factor_icu, const int num_days,
    const std::vector<std::vector<ScalarType>>& immunity_population, const mio::regions::de::EpidataFilenames& epidata_filenames)
{
    const auto num_age_groups = (size_t)model[0].property.parameters.get_num_groups();
    assert(scaling_factor_inf.size() == num_age_groups);
    assert(num_age_groups == ConfirmedCasesDataEntry::age_group_names.size());
    std::vector<TimeSeries<ScalarType>> extrapolated_data(
        model.size(), TimeSeries<ScalarType>::zero(num_days + 1, (size_t)InfectionState::Count * num_age_groups));

    for (int t = 0; t <= num_days; ++t) {
        auto offset_day = offset_date_by_days(date, t);

        // TODO: empty vaccination data path guard
        BOOST_OUTCOME_TRY(read_input_data(model, date, scaling_factor_inf, scaling_factor_icu,
                                          num_days, immunity_population, epidata_filenames));

        for (size_t r = 0; r < model.size(); r++) {
            extrapolated_data[r][t] = model[r].property.get_initial_values();
            // in set_population_data the number of death individuals is subtracted from the SusceptibleImprovedImmunity compartment.
            // Since we should be independent whether we consider them or not, we add them back here before we save the data.
            for (size_t age = 0; age < num_age_groups; age++) {
                extrapolated_data[r][t][(size_t)InfectionState::SusceptibleImprovedImmunity +
                                        age * (size_t)InfectionState::Count] +=
                    extrapolated_data[r][t][(size_t)InfectionState::DeadNaive + age * (size_t)InfectionState::Count];
            }
        }
    }

    std::vector<int> vregion; 
    std::transform(model.begin(), model.end(), std::back_inserter(vregion), [](const auto& m) { return m.id; });
    BOOST_OUTCOME_TRY(save_result(extrapolated_data, vregion, static_cast<int>(num_age_groups),
                                  path_join(results_dir, "Results_rki.h5")));

    auto extrapolated_rki_data_sum = sum_nodes(std::vector<std::vector<TimeSeries<ScalarType>>>{extrapolated_data});
    BOOST_OUTCOME_TRY(save_result({extrapolated_rki_data_sum[0][0]}, {0}, static_cast<int>(num_age_groups),
                                  path_join(results_dir, "Results_rki_sum.h5")));

    return success();
}

#else
IOResult<void> export_input_data_timeseries(const mio::VectorRange<Node<Model<ScalarType>>>, const std::string&, Date, const std::vector<int>&,
                                            const std::vector<ScalarType>&, const ScalarType, const int,
                                            const std::vector<std::vector<ScalarType>>, const mio::regions::de::EpidataFilenames&)
{
    mio::log_warning("HDF5 not available. Cannot export time series of extrapolated real data.");
    return success();
}

#endif //MEMILIO_HAS_HDF5

} // namespace osecirts
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP
