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
    std::vector<ScalarType>& num_Exposed, std::vector<ScalarType>& num_InfectedNoSymptoms,
    std::vector<ScalarType>& num_InfectedSymptoms, std::vector<ScalarType>& num_InfectedSevere,
    std::vector<ScalarType>& num_icu, std::vector<ScalarType>& num_death,
    std::vector<ScalarType>& num_rec, const std::vector<int>& t_Exposed,
    const std::vector<int>& t_InfectedNoSymptoms,
    const std::vector<int>& t_InfectedSymptoms, const std::vector<int>& t_InfectedSevere,
    const std::vector<int>& t_InfectedCritical, const std::vector<ScalarType>& mu_C_R,
    const std::vector<ScalarType>& mu_I_H, const std::vector<ScalarType>& mu_H_U,
    const std::vector<ScalarType>& scaling_factor_inf)
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

IOResult<std::vector<ScalarType>> compute_confirmed_cases_data_fix_recovered(
    const std::vector<ConfirmedCasesDataEntry>& case_data, const int region, Date date, ScalarType delay)
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
    std::vector<ScalarType> num_rec(num_groups, 0.0);

    for (auto&& rki_entry : case_data) {
        if (rki_entry.date == offset_date_by_days(date, int(-delay))) {
            num_rec[size_t(rki_entry.age_group)] = rki_entry.num_confirmed;
        }
    }

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
        try_fix_constraints(num_rec[i], 0, "Recovered");
    }

    return success(num_rec);
}

IOResult<void> set_confirmed_cases_data(Model<ScalarType>& model, const std::vector<ConfirmedCasesDataEntry>& case_data,
                                        const int region, Date date, const std::vector<ScalarType>& scaling_factor_inf, 
                                        bool set_death = false)
{
    auto num_age_groups = (size_t)model.parameters.get_num_groups();
    assert(scaling_factor_inf.size() == num_age_groups); //TODO: allow vector or scalar valued scaling factors
    assert(ConfirmedCasesDataEntry::age_group_names.size() == num_age_groups);

    std::vector<int> t_Exposed;
    std::vector<int> t_InfectedNoSymptoms;
    std::vector<int> t_InfectedSymptoms;
    std::vector<int> t_InfectedSevere;
    std::vector<int> t_InfectedCritical;

    std::vector<ScalarType> mu_C_R;
    std::vector<ScalarType> mu_I_H;
    std::vector<ScalarType> mu_H_U;

    std::vector<ScalarType> num_InfectedSymptoms(num_age_groups, 0.0);
    std::vector<ScalarType> num_death(num_age_groups, 0.0);
    std::vector<ScalarType> num_rec(num_age_groups, 0.0);
    std::vector<ScalarType> num_Exposed(num_age_groups, 0.0);
    std::vector<ScalarType> num_InfectedNoSymptoms(num_age_groups, 0.0);
    std::vector<ScalarType> num_InfectedSevere(num_age_groups, 0.0);
    std::vector<ScalarType> num_icu(num_age_groups, 0.0);

    /*----------- UNVACCINATED -----------*/
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

        mu_C_R.push_back(
            model.parameters.template get<RecoveredPerInfectedNoSymptoms<ScalarType>>()[(AgeGroup)group]);
        mu_I_H.push_back(
            model.parameters.template get<SeverePerInfectedSymptoms<ScalarType>>()[(AgeGroup)group]);
        mu_H_U.push_back(
            model.parameters.template get<CriticalPerSevere<ScalarType>>()[(AgeGroup)group]);
    }

    BOOST_OUTCOME_TRY(compute_confirmed_cases_data(case_data, region, date, num_Exposed, num_InfectedNoSymptoms,
                                                   num_InfectedSymptoms, num_InfectedSevere, num_icu, num_death, num_rec,
                                                   t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
                                                   t_InfectedCritical, mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf));


    size_t num_groups = (size_t)model.parameters.get_num_groups();
    for (size_t i = 0; i < num_groups; i++) {
        model.populations[{AgeGroup(i), InfectionState::ExposedNaive}] = num_Exposed[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsNaive}] =
            num_InfectedNoSymptoms[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsNaiveConfirmed}] = 0;
        model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsNaive}] =
            num_InfectedSymptoms[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsNaiveConfirmed}] = 0;
        model.populations[{AgeGroup(i), InfectionState::InfectedSevereNaive}] =
            num_InfectedSevere[i];
        // Only set the number of ICU patients here, if the date is not available in the data.
        if (!is_divi_data_available(date)) {
            model.populations[{AgeGroup(i), InfectionState::InfectedCriticalNaive}] = num_icu[i];
        }
        model.populations[{AgeGroup(i), InfectionState::SusceptibleImprovedImmunity}] = num_rec[i];
        if (set_death) {
            // in set_confirmed_cases_data initilization, deaths are now set to 0. In order to visualize
            // the extrapolated real number of deaths, they have to be set here. In the comparison of data
            // it has to be paid attention to the fact, the the simulation starts with deaths=0
            // while this method starts with deaths=number of reported deaths so far...
            // Additionally, we set the number of reported deaths to DeadNaive since no information on that is
            // available here.
            // Do only add deaths after substraction.
            model.populations[{AgeGroup(i), InfectionState::DeadNaive}] = num_death[i];
        }
    }
    if (std::accumulate(num_InfectedSymptoms.begin(), num_InfectedSymptoms.end(), ScalarType(0.0),
                        [](const ScalarType& a, const ScalarType& b) {
                            return evaluate_intermediate<ScalarType>(a + b);
                        }) == 0.0) {
        log_warning(
            "No infections for unvaccinated reported on date {} for region {}. Population data has not been set.",
            date, region);
    }

    /*----------- PARTIALLY VACCINATED -----------*/
    t_InfectedNoSymptoms.clear();
    t_Exposed.clear();
    t_InfectedSymptoms.clear();
    t_InfectedSevere.clear();
    t_InfectedCritical.clear();

    mu_C_R.clear();
    mu_I_H.clear();
    mu_H_U.clear();

    num_InfectedSymptoms   = std::vector<ScalarType>(num_age_groups, 0.0);
    num_death              = std::vector<ScalarType>(num_age_groups, 0.0);
    num_rec                = std::vector<ScalarType>(num_age_groups, 0.0);
    num_Exposed            = std::vector<ScalarType>(num_age_groups, 0.0);
    num_InfectedNoSymptoms = std::vector<ScalarType>(num_age_groups, 0.0);
    num_InfectedSevere     = std::vector<ScalarType>(num_age_groups, 0.0);
    num_icu                = std::vector<ScalarType>(num_age_groups, 0.0);

    for (size_t group = 0; group < num_age_groups; group++) {
        ScalarType reduc_t = model.parameters.template get<ReducTimeInfectedMild<ScalarType>>()[(AgeGroup)group];
        t_Exposed.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeExposed<ScalarType>>()[(AgeGroup)group])));
        t_InfectedNoSymptoms.push_back(static_cast<int>(std::round(
            model.parameters.template get<TimeInfectedNoSymptoms<ScalarType>>()[(AgeGroup)group] * reduc_t)));
        t_InfectedSymptoms.push_back(static_cast<int>(std::round(
            model.parameters.template get<TimeInfectedSymptoms<ScalarType>>()[(AgeGroup)group] * reduc_t)));
        t_InfectedSevere.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedSevere<ScalarType>>()[(AgeGroup)group])));
        t_InfectedCritical.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedCritical<ScalarType>>()[(AgeGroup)group])));

        ScalarType exp_fac_part_immune =
            model.parameters.template get<ReducExposedPartialImmunity<ScalarType>>()[(AgeGroup)group];
        ScalarType inf_fac_part_immune =
            model.parameters.template get<ReducInfectedSymptomsPartialImmunity<ScalarType>>()[(AgeGroup)group];
        ScalarType hosp_fac_part_immune =
            model
                .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<ScalarType>>()[(AgeGroup)group];
        ScalarType icu_fac_part_immune =
            model
                .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<ScalarType>>()[(AgeGroup)group];
        mu_C_R.push_back((
            1 - inf_fac_part_immune / exp_fac_part_immune *
                    (1 - model
                                .parameters.template get<RecoveredPerInfectedNoSymptoms<ScalarType>>()[(AgeGroup)group])));
        mu_I_H.push_back(
            hosp_fac_part_immune / inf_fac_part_immune *
            model.parameters.template get<SeverePerInfectedSymptoms<ScalarType>>()[(AgeGroup)group]);
        // transfer from H to U, D unchanged.
        mu_H_U.push_back(
            icu_fac_part_immune / hosp_fac_part_immune *
            model.parameters.template get<CriticalPerSevere<ScalarType>>()[(AgeGroup)group]);
    }

    BOOST_OUTCOME_TRY(compute_confirmed_cases_data(case_data, region, date, num_Exposed, num_InfectedNoSymptoms,
                                                   num_InfectedSymptoms, num_InfectedSevere, num_icu, num_death, num_rec,
                                                   t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
                                                   t_InfectedCritical, mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf));


    for (size_t i = 0; i < num_groups; i++) {
        model.populations[{AgeGroup(i), InfectionState::ExposedPartialImmunity}] = num_Exposed[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsPartialImmunity}] =
            num_InfectedNoSymptoms[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}] = 0;
        model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsPartialImmunity}] =
            num_InfectedSymptoms[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsPartialImmunityConfirmed}] = 0;
        model.populations[{AgeGroup(i), InfectionState::InfectedSeverePartialImmunity}] =
            num_InfectedSevere[i];
        // Only set the number of ICU patients here, if the date is not available in the data.
        if (!is_divi_data_available(date)) {
            model.populations[{AgeGroup(i), InfectionState::InfectedCriticalPartialImmunity}] =
                num_icu[i];
        }
    }
    // }
    if (std::accumulate(num_InfectedSymptoms.begin(), num_InfectedSymptoms.end(), ScalarType(0.0),
                        [](const ScalarType& a, const ScalarType& b) {
                            return evaluate_intermediate<ScalarType>(a + b);
                        }) == 0.0) {
        log_warning("No infections for partially vaccinated reported on date {} for region {}. "
                    "Population data has not been set.",
                    date, region);
    }

    /*----------- FULLY VACCINATED -----------*/
    t_InfectedNoSymptoms.clear();
    t_Exposed.clear();
    t_InfectedSymptoms.clear();
    t_InfectedSevere.clear();
    t_InfectedCritical.clear();

    mu_C_R.clear();
    mu_I_H.clear();
    mu_H_U.clear();

    num_InfectedSymptoms   = std::vector<ScalarType>(num_age_groups, 0.0);
    num_death              = std::vector<ScalarType>(num_age_groups, 0.0);
    num_rec                = std::vector<ScalarType>(num_age_groups, 0.0);
    num_Exposed            = std::vector<ScalarType>(num_age_groups, 0.0);
    num_InfectedNoSymptoms = std::vector<ScalarType>(num_age_groups, 0.0);
    num_InfectedSevere     = std::vector<ScalarType>(num_age_groups, 0.0);
    num_icu                = std::vector<ScalarType>(num_age_groups, 0.0);

    for (size_t group = 0; group < num_age_groups; group++) {
        ScalarType reduc_t = model.parameters.template get<ReducTimeInfectedMild<ScalarType>>()[(AgeGroup)group];
        t_Exposed.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeExposed<ScalarType>>()[(AgeGroup)group])));
        t_InfectedNoSymptoms.push_back(static_cast<int>(std::round(
            model.parameters.template get<TimeInfectedNoSymptoms<ScalarType>>()[(AgeGroup)group] * reduc_t)));
        t_InfectedSymptoms.push_back(static_cast<int>(std::round(
            model.parameters.template get<TimeInfectedSymptoms<ScalarType>>()[(AgeGroup)group] * reduc_t)));
        t_InfectedSevere.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedSevere<ScalarType>>()[(AgeGroup)group])));
        t_InfectedCritical.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedCritical<ScalarType>>()[(AgeGroup)group])));

        ScalarType reduc_immune_exp =
            model.parameters.template get<ReducExposedImprovedImmunity<ScalarType>>()[(AgeGroup)group];
        ScalarType reduc_immune_inf =
            model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<ScalarType>>()[(AgeGroup)group];
        ScalarType reduc_immune_hosp =
            model.parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<ScalarType>>()[(
                AgeGroup)group];
        ScalarType reduc_immune_icu =
            model.parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<ScalarType>>()[(
                AgeGroup)group];
        mu_C_R.push_back((
            1 - reduc_immune_inf / reduc_immune_exp *
                    (1 - model
                                .parameters.template get<RecoveredPerInfectedNoSymptoms<ScalarType>>()[(AgeGroup)group])));
        mu_I_H.push_back(
            reduc_immune_hosp / reduc_immune_inf *
            model.parameters.template get<SeverePerInfectedSymptoms<ScalarType>>()[(AgeGroup)group]);
        // transfer from H to U, D unchanged.
        mu_H_U.push_back(
            reduc_immune_icu / reduc_immune_hosp *
            model.parameters.template get<CriticalPerSevere<ScalarType>>()[(AgeGroup)group]);
    }

    BOOST_OUTCOME_TRY(compute_confirmed_cases_data(case_data, region, date, num_Exposed, num_InfectedNoSymptoms,
                                                   num_InfectedSymptoms, num_InfectedSevere, num_icu, num_death, num_rec,
                                                   t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
                                                   t_InfectedCritical, mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf));

    for (size_t i = 0; i < num_groups; i++) {
        model.populations[{AgeGroup(i), InfectionState::ExposedImprovedImmunity}] = num_Exposed[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsImprovedImmunity}] =
            num_InfectedNoSymptoms[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] = 0;
        model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsImprovedImmunity}] =
            num_InfectedSymptoms[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] = 0;
        model.populations[{AgeGroup(i), InfectionState::InfectedSevereImprovedImmunity}] =
            num_InfectedSevere[i];
        // Only set the number of ICU patients here, if the date is not available in the data.
        if (!is_divi_data_available(date)) {
            model.populations[{AgeGroup(i), InfectionState::InfectedCriticalImprovedImmunity}] =
                num_icu[i];
        }
    }

    if (std::accumulate(num_InfectedSymptoms.begin(), num_InfectedSymptoms.end(), ScalarType(0.0),
                        [](const ScalarType& a, const ScalarType& b) {
                            return evaluate_intermediate<ScalarType>(a + b);
                        }) == 0.0) {
        log_warning("No infections for vaccinated reported on date {} for region {}. "
                    "Population data has not been set.",
                    date, region);
    }

    return success();
}

IOResult<void> set_confirmed_cases_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, const std::string& path, Date date,
                                        const std::vector<ScalarType>& scaling_factor_inf, bool set_death = false)
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
        BOOST_OUTCOME_TRY(set_confirmed_cases_data(model[region_idx].property, vcase_data[region_idx], model[region_idx].id, date, scaling_factor_inf, set_death));
    }
    return success();
}

IOResult<void> set_population_data(Model<ScalarType>& model, const std::vector<ScalarType>& num_population,
                                   const int region, const std::vector<ScalarType>& num_rec)
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

        ScalarType S_v = std::min(
            model.parameters.template get<DailyFullVaccinations<ScalarType>>()[{i, SimulationDay(0)}] +
                num_rec[size_t(i)],
            num_population[size_t(i)]);
        ScalarType S_pv = std::max(
            model.parameters.template get<DailyPartialVaccinations<ScalarType>>()[{i, SimulationDay(0)}] -
                model.parameters.template get<DailyFullVaccinations<ScalarType>>()[{i, SimulationDay(0)}],
            0.0); // use std::max with 0
        ScalarType S;
        if (num_population[size_t(i)] - S_pv - S_v < 0.0) {
            log_warning("Number of vaccinated persons greater than population in county {}, age group {}.",
                        region, size_t(i));
            S   = 0.0;
            S_v = num_population[size_t(i)] - S_pv;
        }
        else {
            S = num_population[size_t(i)] - S_pv - S_v;
        }

        ScalarType denom_E =
            1 / (S + S_pv * model.parameters.template get<ReducExposedPartialImmunity<ScalarType>>()[i] +
                    S_v * model.parameters.template get<ReducExposedImprovedImmunity<ScalarType>>()[i]);
        ScalarType denom_C =
            1 / (S + S_pv * model.parameters.template get<ReducExposedPartialImmunity<ScalarType>>()[i] +
                    S_v * model.parameters.template get<ReducExposedImprovedImmunity<ScalarType>>()[i]);
        ScalarType denom_I =
            1 /
            (S +
                S_pv * model.parameters.template get<ReducInfectedSymptomsPartialImmunity<ScalarType>>()[i] +
                S_v * model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<ScalarType>>()[i]);
        ScalarType denom_HU =
            1 /
            (S +
                S_pv * model
                        .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<ScalarType>>()[i] +
                S_v * model
                        .parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<ScalarType>>()[i]);

        model.populations[{i, InfectionState::ExposedNaive}] =
            S * model.populations[{i, InfectionState::ExposedNaive}] * denom_E;
        model.populations[{i, InfectionState::ExposedPartialImmunity}] =
            S_pv * model.parameters.template get<ReducExposedPartialImmunity<ScalarType>>()[i] *
            model.populations[{i, InfectionState::ExposedPartialImmunity}] * denom_E;
        model.populations[{i, InfectionState::ExposedImprovedImmunity}] =
            S_v * model.parameters.template get<ReducExposedImprovedImmunity<ScalarType>>()[i] *
            model.populations[{i, InfectionState::ExposedImprovedImmunity}] * denom_E;

        model.populations[{i, InfectionState::InfectedNoSymptomsNaive}] =
            S * model.populations[{i, InfectionState::InfectedNoSymptomsNaive}] * denom_C;
        model.populations[{i, InfectionState::InfectedNoSymptomsPartialImmunity}] =
            S_pv * model.populations[{i, InfectionState::InfectedNoSymptomsPartialImmunity}] * denom_C;
        model.populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunity}] =
            S_v * model.populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunity}] * denom_C;

        model.populations[{i, InfectionState::InfectedNoSymptomsNaiveConfirmed}] =
            S * model.populations[{i, InfectionState::InfectedNoSymptomsNaiveConfirmed}] * denom_C;
        model.populations[{i, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}] =
            S_pv * model.populations[{i, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}] *
            denom_C;
        model.populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] =
            S_v * model.populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] *
            denom_C;

        model.populations[{i, InfectionState::InfectedSymptomsNaive}] =
            S * model.populations[{i, InfectionState::InfectedSymptomsNaive}] * denom_I;
        model.populations[{i, InfectionState::InfectedSymptomsPartialImmunity}] =
            S_pv * model.parameters.template get<ReducInfectedSymptomsPartialImmunity<ScalarType>>()[i] *
            model.populations[{i, InfectionState::InfectedSymptomsPartialImmunity}] * denom_I;
        model.populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] =
            S_v * model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<ScalarType>>()[i] *
            model.populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] * denom_I;

        model.populations[{i, InfectionState::InfectedSymptomsNaiveConfirmed}] =
            S * model.populations[{i, InfectionState::InfectedSymptomsNaiveConfirmed}] * denom_I;
        model.populations[{i, InfectionState::InfectedSymptomsPartialImmunityConfirmed}] =
            S_pv * model.parameters.template get<ReducInfectedSymptomsPartialImmunity<ScalarType>>()[i] *
            model.populations[{i, InfectionState::InfectedSymptomsPartialImmunityConfirmed}] * denom_I;
        model.populations[{i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] =
            S_v * model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<ScalarType>>()[i] *
            model.populations[{i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] * denom_I;

        model.populations[{i, InfectionState::InfectedSevereNaive}] =
            S * model.populations[{i, InfectionState::InfectedSevereNaive}] * denom_HU;
        model.populations[{i, InfectionState::InfectedSeverePartialImmunity}] =
            S_pv *
            model.parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<ScalarType>>()[i] *
            model.populations[{i, InfectionState::InfectedSeverePartialImmunity}] * denom_HU;
        model.populations[{i, InfectionState::InfectedSevereImprovedImmunity}] =
            S_v *
            model
                .parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<ScalarType>>()[i] *
            model.populations[{i, InfectionState::InfectedSevereImprovedImmunity}] * denom_HU;

        model.populations[{i, InfectionState::InfectedCriticalPartialImmunity}] =
            S_pv *
            model.parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<ScalarType>>()[i] *
            model.populations[{i, InfectionState::InfectedCriticalNaive}] * denom_HU;
        model.populations[{i, InfectionState::InfectedCriticalImprovedImmunity}] =
            S_v *
            model
                .parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<ScalarType>>()[i] *
            model.populations[{i, InfectionState::InfectedCriticalNaive}] * denom_HU;
        model.populations[{i, InfectionState::InfectedCriticalNaive}] =
            S * model.populations[{i, InfectionState::InfectedCriticalNaive}] * denom_HU;

        model.populations[{i, InfectionState::SusceptibleImprovedImmunity}] =
            model.parameters.template get<DailyFullVaccinations<ScalarType>>()[{i, SimulationDay(0)}] +
            model.populations[{i, InfectionState::SusceptibleImprovedImmunity}] -
            (model.populations[{i, InfectionState::InfectedSymptomsNaive}] +
                model.populations[{i, InfectionState::InfectedSymptomsPartialImmunity}] +
                model.populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] +
                model.populations[{i, InfectionState::InfectedSymptomsNaiveConfirmed}] +
                model.populations[{i, InfectionState::InfectedSymptomsPartialImmunityConfirmed}] +
                model.populations[{i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] +
                model.populations[{i, InfectionState::InfectedSevereNaive}] +
                model.populations[{i, InfectionState::InfectedSeverePartialImmunity}] +
                model.populations[{i, InfectionState::InfectedSevereImprovedImmunity}] +
                model.populations[{i, InfectionState::InfectedCriticalNaive}] +
                model.populations[{i, InfectionState::InfectedCriticalPartialImmunity}] +
                model.populations[{i, InfectionState::InfectedCriticalImprovedImmunity}] +
                model.populations[{i, InfectionState::DeadNaive}] +
                model.populations[{i, InfectionState::DeadPartialImmunity}] +
                model.populations[{i, InfectionState::DeadImprovedImmunity}]);

        model.populations[{i, InfectionState::SusceptibleImprovedImmunity}] = std::min(
            S_v - model.populations[{i, InfectionState::ExposedImprovedImmunity}] -
                model.populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunity}] -
                model.populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] -
                model.populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] -
                model.populations[{i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] -
                model.populations[{i, InfectionState::InfectedSevereImprovedImmunity}] -
                model.populations[{i, InfectionState::InfectedCriticalImprovedImmunity}] -
                model.populations[{i, InfectionState::DeadImprovedImmunity}],
            std::max(0.0, ScalarType(model.populations[{i, InfectionState::SusceptibleImprovedImmunity}])));

        model.populations[{i, InfectionState::SusceptiblePartialImmunity}] = std::max(
            0.0,
            S_pv - model.populations[{i, InfectionState::ExposedPartialImmunity}] -
                model.populations[{i, InfectionState::InfectedNoSymptomsPartialImmunity}] -
                model.populations[{i, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}] -
                model.populations[{i, InfectionState::InfectedSymptomsPartialImmunity}] -
                model.populations[{i, InfectionState::InfectedSymptomsPartialImmunityConfirmed}] -
                model.populations[{i, InfectionState::InfectedSeverePartialImmunity}] -
                model.populations[{i, InfectionState::InfectedCriticalPartialImmunity}] -
                model.populations[{i, InfectionState::DeadPartialImmunity}]);

        model.populations.template set_difference_from_group_total<AgeGroup>(
            {i, InfectionState::SusceptibleNaive}, num_population[size_t(i)]);
    }

    for (auto i = AgeGroup(0); i < AgeGroup(6); i++) {
        for (auto j = Index<InfectionState>(0); j < InfectionState::Count; ++j) {
            if (model.populations[{i, j}] < 0) {
                log_warning("Compartment at age group {}, infection state {}, is negative: {}", size_t(i),
                            size_t(j), model.populations[{i, j}] / num_population[size_t(i)]);
            }
        }
    }

    return success();
}

IOResult<void> set_population_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, const std::string& path, const std::string& path_rki,
                                   Date date)
{
    std::vector<int> vregion; 
    std::transform(model.begin(), model.end(), std::back_inserter(vregion), [](const auto& m) { return m.id; });
    BOOST_OUTCOME_TRY(auto&& vnum_population, read_population_data(path, vregion));
    BOOST_OUTCOME_TRY(auto&& case_data, mio::read_confirmed_cases_data(path_rki));

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
        BOOST_OUTCOME_TRY(auto&& num_rec, compute_confirmed_cases_data_fix_recovered(vcase_data[region_idx], model[region_idx].id, date, 14.));
        BOOST_OUTCOME_TRY(set_population_data(model[region_idx].property, vnum_population[region_idx], model[region_idx].id, num_rec));
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
                               int num_days, const mio::regions::de::EpidataFilenames& epidata_filenames)
{

    BOOST_OUTCOME_TRY(
        details::set_vaccination_data<ScalarType>(model, epidata_filenames.vaccination_data_path, date, num_days));

    // TODO: Reuse more code, e.g., set_divi_data (in secir) and a set_divi_data (here) only need a different ModelType.
    // TODO: add option to set ICU data from confirmed cases if DIVI or other data is not available.
    BOOST_OUTCOME_TRY(details::set_divi_data(model, epidata_filenames.divi_data_path, date,
                                             scaling_factor_icu));

    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(model, epidata_filenames.case_data_path,
                                                        date, scaling_factor_inf));
    BOOST_OUTCOME_TRY(details::set_population_data(model, epidata_filenames.population_data_path, 
                                                   epidata_filenames.case_data_path, date));
    return success();
}

#ifdef MEMILIO_HAS_HDF5

IOResult<void> export_input_data_timeseries(
    const mio::VectorRange<Node<Model<ScalarType>>> model, const std::string& results_dir, Date date,
    const std::vector<ScalarType>& scaling_factor_inf, const ScalarType scaling_factor_icu, const int num_days,
    const mio::regions::de::EpidataFilenames& epidata_filenames)
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
                                      num_days, epidata_filenames));

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
IOResult<void> export_input_data_county_timeseries(const mio::VectorRange<Node<Model<ScalarType>>>, const std::string&, Date, const std::vector<int>&,
                                                   const std::vector<ScalarType>&, const ScalarType, const int,
                                                   const mio::regions::de::EpidataFilenames&)
{
    mio::log_warning("HDF5 not available. Cannot export time series of extrapolated real data.");
    return success();
}

#endif //MEMILIO_HAS_HDF5

} // namespace osecirvvs
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP
