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
#ifndef MIO_ODE_SECIRTS_PARAMETERS_IO_H
#define MIO_ODE_SECIRTS_PARAMETERS_IO_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "ode_secirts/model.h"
#include "ode_secirts/analyze_result.h"
#include "memilio/math/eigen_util.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"
#include "memilio/io/json_serializer.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/date.h"

namespace mio
{
namespace osecirts
{

namespace details
{

/**
 * @brief Gets the region ID (county, state, or district) of an EpiDataEntry.
 * 
 * If none are available, it defaults to 0 which is representing the whole country.
 * 
 * @tparam EpiDataEntry The type of the data entry.
 * @param data_entry The (RKI) data entry to extract the region ID from.
 * @return The region ID as integer, or 0 if no specific region information is available.
 */
template <class EpiDataEntry>
int get_region_id(const EpiDataEntry& data_entry)
{
    if (data_entry.county_id) {
        return data_entry.county_id->get();
    }
    if (data_entry.state_id) {
        return data_entry.state_id->get();
    }
    if (data_entry.district_id) {
        return data_entry.district_id->get();
    }
    return 0;
}

/**
 * @brief Computes the distribution of confirmed cases across infection states based on Case (RKI) data.
 *
 * This function processes case data for given regions and distributes the cases across different 
 * infection states, considering the corresponding transition times and probabilities defined in the model. 
 *
 * @tparam Model The type of the model used.
 * @tparam FP Floating point type (default: double).
 *
 * @param[in] case_data Vector of confirmed case data entries (defined in epi_data.h).
 * @param[out] vnum_Exposed Output vector for the number of exposed individuals per age group and region.
 * @param[out] vnum_InfectedNoSymptoms Output vector for the number of infected individuals without symptoms.
 * @param[out] vnum_InfectedSymptoms Output vector for the number of infected individuals with symptoms.
 * @param[out] vnum_InfectedSevere Output vector for the number of severely infected individuals.
 * @param[out] vnum_icu Output vector for the number of individuals in critical condition (ICU).
 * @param[out] vnum_death Output vector for the number of deaths.
 * @param[out] vnum_timm_i Output vector for the number of individuals in temporary immunity state.
 * @param[in] vregion Vector of region IDs representing the regions in the model vector.
 * @param[in] date Date for which the simulation starts.
 * @param[in] model Vector of models, each representing a region and containing the parameters.
 * @param[in] scaling_factor_inf Vector of scaling factors for confirmed cases for 
 * @param[in] layer Specifies the immunity layer: 0 (Naive), 1 (Partial Immunity), 2 (Improved Immunity).
 *
 * @return An IOResult showing success or failure.
 */
template <class Model, typename FP = double>
IOResult<void> compute_confirmed_cases_data(
    const std::vector<ConfirmedCasesDataEntry>& case_data, std::vector<std::vector<FP>>& vnum_Exposed,
    std::vector<std::vector<FP>>& vnum_InfectedNoSymptoms, std::vector<std::vector<FP>>& vnum_InfectedSymptoms,
    std::vector<std::vector<FP>>& vnum_InfectedSevere, std::vector<std::vector<FP>>& vnum_icu,
    std::vector<std::vector<FP>>& vnum_death, std::vector<std::vector<FP>>& vnum_timm_i,
    std::vector<int> const& vregion, Date date, const std::vector<Model>& model,
    const std::vector<FP>& scaling_factor_inf, const size_t layer)
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
        auto it = std::find_if(vregion.begin(), vregion.end(), [&entry](auto r) {
            return r == 0 || get_region_id(entry) == r;
        });
        if (it != vregion.end()) {
            auto region_idx = size_t(it - vregion.begin());

            auto params_region           = model[region_idx].parameters;
            auto& num_InfectedNoSymptoms = vnum_InfectedNoSymptoms[region_idx];
            auto& num_InfectedSymptoms   = vnum_InfectedSymptoms[region_idx];
            auto& num_Exposed            = vnum_Exposed[region_idx];
            auto& num_InfectedSevere     = vnum_InfectedSevere[region_idx];
            auto& num_death              = vnum_death[region_idx];
            auto& num_icu                = vnum_icu[region_idx];
            auto& num_imm                = vnum_timm_i[region_idx];

            auto age = (size_t)entry.age_group;
            // (rounded) transition times
            const int t_exposed =
                static_cast<int>(std::round(params_region.template get<TimeExposed<FP>>()[entry.age_group]));
            int t_InfectedNoSymptoms =
                static_cast<int>(std::round(params_region.template get<TimeInfectedNoSymptoms<FP>>()[entry.age_group]));
            int t_InfectedSymptoms =
                static_cast<int>(std::round(params_region.template get<TimeInfectedSymptoms<FP>>()[entry.age_group]));
            const int t_InfectedSevere =
                static_cast<int>(std::round(params_region.template get<TimeInfectedSevere<FP>>()[entry.age_group]));
            const int t_InfectedCritical =
                static_cast<int>(std::round(params_region.template get<TimeInfectedCritical<FP>>()[entry.age_group]));
            const int t_imm_interval_i = static_cast<int>(
                std::round(params_region.template get<TimeTemporaryImmunityPI<FP>>()[entry.age_group]));

            // transition probabilities
            FP recoveredPerInfectedNoSymptoms =
                params_region.template get<RecoveredPerInfectedNoSymptoms<FP>>()[entry.age_group];
            FP severePerInfectedSymptoms = params_region.template get<SeverePerInfectedSymptoms<FP>>()[entry.age_group];
            FP criticalPerSevere         = params_region.template get<CriticalPerSevere<FP>>()[entry.age_group];

            // if we select a layer with better immunity (layer > 0), we need to adjust the times and transition rates
            if (layer > 0) {
                t_InfectedNoSymptoms = static_cast<int>(std::round(
                    t_InfectedNoSymptoms * params_region.template get<ReducTimeInfectedMild<FP>>()[entry.age_group]));
                t_InfectedSymptoms   = static_cast<int>(std::round(
                    t_InfectedSymptoms * params_region.template get<ReducTimeInfectedMild<FP>>()[entry.age_group]));

                const FP red_fact_exp =
                    layer == 1 ? params_region.template get<ReducExposedPartialImmunity<FP>>()[entry.age_group]
                               : params_region.template get<ReducExposedImprovedImmunity<FP>>()[entry.age_group];

                const FP red_fact_inf =
                    layer == 1
                        ? params_region.template get<ReducInfectedSymptomsPartialImmunity<FP>>()[entry.age_group]
                        : params_region.template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[entry.age_group];

                const FP red_fact_sev =
                    layer == 1
                        ? params_region
                              .template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[entry.age_group]
                        : params_region
                              .template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[entry.age_group];

                recoveredPerInfectedNoSymptoms = 1 - red_fact_inf / red_fact_exp * (1 - recoveredPerInfectedNoSymptoms);
                severePerInfectedSymptoms      = red_fact_sev / red_fact_inf * severePerInfectedSymptoms;
            }

            if (entry.date == offset_date_by_days(date, 0)) {
                num_InfectedSymptoms[age] += scaling_factor_inf[age] * entry.num_confirmed;
                num_imm[age] += entry.num_confirmed;
            }
            if (entry.date == offset_date_by_days(date, t_InfectedNoSymptoms + days_surplus)) {
                num_InfectedNoSymptoms[age] +=
                    1 / (1 - recoveredPerInfectedNoSymptoms) * scaling_factor_inf[age] * entry.num_confirmed;
                num_Exposed[age] -=
                    1 / (1 - recoveredPerInfectedNoSymptoms) * scaling_factor_inf[age] * entry.num_confirmed;
            }
            if (entry.date == offset_date_by_days(date, days_surplus)) {
                num_InfectedNoSymptoms[age] -=
                    1 / (1 - recoveredPerInfectedNoSymptoms) * scaling_factor_inf[age] * entry.num_confirmed;
            }
            if (entry.date == offset_date_by_days(date, t_exposed + t_InfectedNoSymptoms + days_surplus)) {
                num_Exposed[age] +=
                    1 / (1 - recoveredPerInfectedNoSymptoms) * scaling_factor_inf[age] * entry.num_confirmed;
            }
            if (entry.date == offset_date_by_days(date, -t_InfectedSymptoms)) {
                num_InfectedSymptoms[age] -= scaling_factor_inf[age] * entry.num_confirmed;
                num_InfectedSevere[age] += severePerInfectedSymptoms * scaling_factor_inf[age] * entry.num_confirmed;
            }
            if (entry.date == offset_date_by_days(date, -t_InfectedSymptoms - t_InfectedSevere)) {
                num_InfectedSevere[age] -= severePerInfectedSymptoms * scaling_factor_inf[age] * entry.num_confirmed;
                num_icu[age] +=
                    severePerInfectedSymptoms * criticalPerSevere * scaling_factor_inf[age] * entry.num_confirmed;
            }
            if (entry.date == offset_date_by_days(date, -t_InfectedSymptoms - t_InfectedSevere - t_InfectedCritical)) {
                num_death[age] += entry.num_deaths;
                num_icu[age] -=
                    severePerInfectedSymptoms * criticalPerSevere * scaling_factor_inf[age] * entry.num_confirmed;
            }
            if (entry.date == offset_date_by_days(date, 0 - t_imm_interval_i)) {
                num_imm[age] -= entry.num_confirmed;
            }
        }
    }

    for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
        auto region = vregion[region_idx];

        auto& num_InfectedNoSymptoms = vnum_InfectedNoSymptoms[region_idx];
        auto& num_InfectedSymptoms   = vnum_InfectedSymptoms[region_idx];
        auto& num_Exposed            = vnum_Exposed[region_idx];
        auto& num_InfectedSevere     = vnum_InfectedSevere[region_idx];
        auto& num_death              = vnum_death[region_idx];
        auto& num_icu                = vnum_icu[region_idx];
        auto& num_timm_i             = vnum_timm_i[region_idx];

        size_t num_groups = ConfirmedCasesDataEntry::age_group_names.size();
        for (size_t i = 0; i < num_groups; i++) {
            auto try_fix_constraints = [region, i](FP& value, FP error, auto str) {
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

            const FP tol_error = -1e-8;
            try_fix_constraints(num_InfectedSymptoms[i], tol_error, "InfectedSymptoms");
            try_fix_constraints(num_InfectedNoSymptoms[i], tol_error, "InfectedNoSymptoms");
            try_fix_constraints(num_Exposed[i], tol_error, "Exposed");
            try_fix_constraints(num_InfectedSevere[i], tol_error, "InfectedSevere");
            try_fix_constraints(num_death[i], tol_error, "Dead");
            try_fix_constraints(num_icu[i], tol_error, "InfectedCritical");
            try_fix_constraints(num_timm_i[i], tol_error, "Recently Recovered or Vaccinated");
        }
    }

    return success();
}

/**
 * @brief Reads confirmed case data from a file and computes the distribution of cases across infection states.
 *
 * This function reads transformed RKI data from a specified file and processes the confirmed cases 
 * to distribute them across different infection states and age groups.
 *
 * @tparam Model The type of the model used.
 * @tparam FP Floating point type (default: double).
 *
 * @param[in] path Path to the file containing transformed case (rki) data.
 * @param[in] vregion Vector of region IDs to process.
 * @param[in] date Date for which the simulation starts.
 * @param[out] vnum_Exposed Output vector for the number of exposed individuals per age group and region.
 * @param[out] vnum_InfectedNoSymptoms Output vector for the number of infected individuals without symptoms.
 * @param[out] vnum_InfectedSymptoms Output vector for the number of infected individuals with symptoms.
 * @param[out] vnum_InfectedSevere Output vector for the number of severely infected individuals (Hospitalized).
 * @param[out] vnum_icu Output vector for the number of individuals in critical condition (ICU).
 * @param[out] vnum_death Output vector for the number of deaths.
 * @param[out] vnum_timm_i Output vector for the number of individuals in a temporary immunity state.
 * @param[in] model Vector of models, each representing a region and containing the parameters.
 * @param[in] scaling_factor_inf Vector of scaling factors for confirmed cases.
 * @param[in] layer Specifies the immunity layer: 0 (Naive), 1 (Partial Immunity), 2 (Improved Immunity).
 *
 * @return An IOResult indicating success or failure.
 */
template <class Model, typename FP = double>
IOResult<void> read_confirmed_cases_data(
    std::string const& path, std::vector<int> const& vregion, Date date, std::vector<std::vector<FP>>& vnum_Exposed,
    std::vector<std::vector<FP>>& vnum_InfectedNoSymptoms, std::vector<std::vector<FP>>& vnum_InfectedSymptoms,
    std::vector<std::vector<FP>>& vnum_InfectedSevere, std::vector<std::vector<FP>>& vnum_icu,
    std::vector<std::vector<FP>>& vnum_death, std::vector<std::vector<FP>>& vnum_timm_i,
    const std::vector<Model>& model, const std::vector<FP>& scaling_factor_inf, const size_t layer)
{
    BOOST_OUTCOME_TRY(auto&& case_data, mio::read_confirmed_cases_data(path));
    return compute_confirmed_cases_data(case_data, vnum_Exposed, vnum_InfectedNoSymptoms, vnum_InfectedSymptoms,
                                        vnum_InfectedSevere, vnum_icu, vnum_death, vnum_timm_i, vregion, date, model,
                                        scaling_factor_inf, layer);
}

/**
 * @brief Sets the confirmed cases data in the model considering different immunity layers.
 *
 * This function distributes confirmed case data across infection states for regions and age groups 
 * in the model. It considers different levels of immunity (naive, partial, and improved).
 *
 * @tparam Model The type of the model used.
 * @tparam FP Floating point type (default: double).
 *
 * @param[in,out] model Vector of models, each representing a region, where the compartments are updated.
 * @param[in] case_data Vector of confirmed case data entries.
 * @param[in] region Vector of region IDs for which the data is processed.
 * @param[in] date Date for which the confirmed cases are set in the model.
 * @param[in] scaling_factor_inf Vector of scaling factors for confirmed cases.
 * @param[in] immunity_population Vector containing the immunity distribution for naive, partial, and improved immunity layers.
 *
 * @return An IOResult indicating success or failure.
 */
template <class Model, typename FP = double>
IOResult<void>
set_confirmed_cases_data(std::vector<Model>& model, const std::vector<ConfirmedCasesDataEntry>& case_data,
                         std::vector<int> const& region, Date date, const std::vector<FP>& scaling_factor_inf,
                         const std::vector<std::vector<FP>> immunity_population)
{
    auto num_age_groups = (size_t)model[0].parameters.get_num_groups();
    assert(scaling_factor_inf.size() == num_age_groups); //TODO: allow vector or scalar valued scaling factors
    assert(ConfirmedCasesDataEntry::age_group_names.size() == num_age_groups);

    std::vector<std::vector<FP>> num_InfectedSymptoms(model.size());
    std::vector<std::vector<FP>> num_death(model.size());
    std::vector<std::vector<FP>> num_Exposed(model.size());
    std::vector<std::vector<FP>> num_InfectedNoSymptoms(model.size());
    std::vector<std::vector<FP>> num_InfectedSevere(model.size());
    std::vector<std::vector<FP>> num_icu(model.size());
    std::vector<std::vector<FP>> num_timm1(model.size());
    std::vector<std::vector<FP>> num_timm2(model.size());

    std::vector<FP> denom_E(num_age_groups, 0.0);
    std::vector<FP> denom_I_NS(num_age_groups, 0.0);
    std::vector<FP> denom_I_Sy(num_age_groups, 0.0);
    std::vector<FP> denom_I_Sev_Cr(num_age_groups, 0.0);

    /*----------- Naive immunity -----------*/
    for (size_t county = 0; county < model.size(); county++) {
        num_InfectedSymptoms[county]   = std::vector<FP>(num_age_groups, 0.0);
        num_death[county]              = std::vector<FP>(num_age_groups, 0.0);
        num_Exposed[county]            = std::vector<FP>(num_age_groups, 0.0);
        num_InfectedNoSymptoms[county] = std::vector<FP>(num_age_groups, 0.0);
        num_InfectedSevere[county]     = std::vector<FP>(num_age_groups, 0.0);
        num_icu[county]                = std::vector<FP>(num_age_groups, 0.0);
        num_timm1[county]              = std::vector<FP>(num_age_groups, 0.0);
        num_timm2[county]              = std::vector<FP>(num_age_groups, 0.0);
        for (size_t group = 0; group < num_age_groups; group++) {
            // calculate the denominators to split the reported case numbers to the different immunity layers.
            denom_E[group] =
                1 / (immunity_population[0][group] +
                     immunity_population[1][group] *
                         model[county].parameters.template get<ReducExposedPartialImmunity<FP>>()[(AgeGroup)group] +
                     immunity_population[2][group] *
                         model[county].parameters.template get<ReducExposedImprovedImmunity<FP>>()[(AgeGroup)group]);

            denom_I_NS[group] =
                1 / (immunity_population[0][group] +
                     immunity_population[1][group] *
                         model[county].parameters.template get<ReducExposedPartialImmunity<FP>>()[(AgeGroup)group] +
                     immunity_population[2][group] *
                         model[county].parameters.template get<ReducExposedImprovedImmunity<FP>>()[(AgeGroup)group]);

            denom_I_Sy[group] =
                1 / (immunity_population[0][group] +
                     immunity_population[1][group] *
                         model[county]
                             .parameters.template get<ReducInfectedSymptomsPartialImmunity<FP>>()[(AgeGroup)group] +
                     immunity_population[2][group] *
                         model[county]
                             .parameters.template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[(AgeGroup)group]);

            denom_I_Sev_Cr[group] =
                1 / (immunity_population[0][group] +
                     immunity_population[1][group] *
                         model[county].parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[(
                             AgeGroup)group] +
                     immunity_population[2][group] *
                         model[county].parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[(
                             AgeGroup)group]);
        }
    }

    BOOST_OUTCOME_TRY(compute_confirmed_cases_data(case_data, num_Exposed, num_InfectedNoSymptoms, num_InfectedSymptoms,
                                                   num_InfectedSevere, num_icu, num_death, num_timm1, region, date,
                                                   model, scaling_factor_inf, 0));

    for (size_t county = 0; county < model.size(); county++) {
        size_t num_groups = (size_t)model[county].parameters.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            model[county].populations[{AgeGroup(i), InfectionState::ExposedNaive}] =
                immunity_population[0][i] * denom_E[i] * num_Exposed[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsNaive}] =
                immunity_population[0][i] * denom_I_NS[i] * num_InfectedNoSymptoms[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsNaiveConfirmed}] = 0;
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsNaive}] =
                immunity_population[0][i] * denom_I_Sy[i] * num_InfectedSymptoms[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsNaiveConfirmed}] = 0;
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSevereNaive}] =
                immunity_population[0][i] * denom_I_Sev_Cr[i] * num_InfectedSevere[county][i];
            // Only set the number of ICU patients here, if the date is not available in the data.
            if (date <= Date(2020, 4, 23) || date >= Date(2024, 7, 21)) {
                model[county].populations[{AgeGroup(i), InfectionState::InfectedCriticalNaive}] =
                    immunity_population[0][i] * denom_I_Sev_Cr[i] * num_icu[county][i];
            }
        }
        if (std::accumulate(num_InfectedSymptoms[county].begin(), num_InfectedSymptoms[county].end(), 0.0) == 0) {
            log_warning("No infections for unvaccinated reported on date " + std::to_string(date.year) + "-" +
                        std::to_string(date.month) + "-" + std::to_string(date.day) + " for region " +
                        std::to_string(region[county]) + ". Population data has not been set.");
        }
    }

    /*----------- PARTIAL Immunity -----------*/
    for (size_t county = 0; county < model.size(); county++) {
        num_InfectedSymptoms[county]   = std::vector<FP>(num_age_groups, 0.0);
        num_death[county]              = std::vector<FP>(num_age_groups, 0.0);
        num_Exposed[county]            = std::vector<FP>(num_age_groups, 0.0);
        num_InfectedNoSymptoms[county] = std::vector<FP>(num_age_groups, 0.0);
        num_InfectedSevere[county]     = std::vector<FP>(num_age_groups, 0.0);
        num_icu[county]                = std::vector<FP>(num_age_groups, 0.0);
    }

    BOOST_OUTCOME_TRY(compute_confirmed_cases_data(case_data, num_Exposed, num_InfectedNoSymptoms, num_InfectedSymptoms,
                                                   num_InfectedSevere, num_icu, num_death, num_timm1, region, date,
                                                   model, scaling_factor_inf, 1));

    for (size_t county = 0; county < model.size(); county++) {
        size_t num_groups = (size_t)model[county].parameters.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            model[county].populations[{AgeGroup(i), InfectionState::ExposedPartialImmunity}] =
                immunity_population[1][i] *
                model[county].parameters.template get<ReducExposedPartialImmunity<FP>>()[(AgeGroup)i] * denom_E[i] *
                num_Exposed[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsPartialImmunity}] =
                immunity_population[1][i] *
                model[county].parameters.template get<ReducExposedPartialImmunity<FP>>()[(AgeGroup)i] * denom_I_NS[i] *
                num_InfectedNoSymptoms[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}] = 0;
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsPartialImmunity}] =
                immunity_population[1][i] *
                model[county].parameters.template get<ReducInfectedSymptomsPartialImmunity<FP>>()[(AgeGroup)i] *
                denom_I_Sy[i] * num_InfectedSymptoms[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsPartialImmunityConfirmed}] = 0;
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSeverePartialImmunity}] =
                immunity_population[1][i] *
                model[county]
                    .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[(AgeGroup)i] *
                denom_I_Sev_Cr[i] * num_InfectedSevere[county][i];
            // Only set the number of ICU patients here, if the date is not available in the data.
            if (date <= Date(2020, 4, 23) || date >= Date(2024, 7, 21)) {
                model[county].populations[{AgeGroup(i), InfectionState::InfectedCriticalPartialImmunity}] =
                    immunity_population[1][i] *
                    model[county]
                        .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[(AgeGroup)i] *
                    denom_I_Sev_Cr[i] * num_icu[county][i];
            }
            // the += is necessary because we already set the previous vaccinated individuals
            model[county].populations[{AgeGroup(i), InfectionState::TemporaryImmunePartialImmunity}] +=
                immunity_population[1][i] *
                model[county].parameters.template get<ReducExposedPartialImmunity<FP>>()[(AgeGroup)i] * denom_E[i] *
                num_timm1[county][i];
        }
        if (std::accumulate(num_InfectedSymptoms[county].begin(), num_InfectedSymptoms[county].end(), 0.0) == 0) {
            log_warning("No infections for partially vaccinated reported on date " + std::to_string(date.year) + "-" +
                        std::to_string(date.month) + "-" + std::to_string(date.day) + " for region " +
                        std::to_string(region[county]) + ". Population data has not been set.");
        }
    }

    /*----------- Improved Immunity -----------*/
    for (size_t county = 0; county < model.size(); county++) {
        num_InfectedSymptoms[county]   = std::vector<FP>(num_age_groups, 0.0);
        num_death[county]              = std::vector<FP>(num_age_groups, 0.0);
        num_Exposed[county]            = std::vector<FP>(num_age_groups, 0.0);
        num_InfectedNoSymptoms[county] = std::vector<FP>(num_age_groups, 0.0);
        num_InfectedSevere[county]     = std::vector<FP>(num_age_groups, 0.0);
        num_icu[county]                = std::vector<FP>(num_age_groups, 0.0);
    }

    BOOST_OUTCOME_TRY(compute_confirmed_cases_data(case_data, num_Exposed, num_InfectedNoSymptoms, num_InfectedSymptoms,
                                                   num_InfectedSevere, num_icu, num_death, num_timm2, region, date,
                                                   model, scaling_factor_inf, 2));

    for (size_t county = 0; county < model.size(); county++) {
        size_t num_groups = (size_t)model[county].parameters.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            model[county].populations[{AgeGroup(i), InfectionState::ExposedImprovedImmunity}] =
                immunity_population[2][i] *
                model[county].parameters.template get<ReducExposedImprovedImmunity<FP>>()[(AgeGroup)i] * denom_E[i] *
                num_Exposed[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsImprovedImmunity}] =
                immunity_population[2][i] *
                model[county].parameters.template get<ReducExposedImprovedImmunity<FP>>()[(AgeGroup)i] * denom_I_NS[i] *
                num_InfectedNoSymptoms[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] = 0;
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsImprovedImmunity}] =
                immunity_population[2][i] *
                model[county].parameters.template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[(AgeGroup)i] *
                denom_I_Sy[i] * num_InfectedSymptoms[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] = 0;
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSevereImprovedImmunity}] =
                immunity_population[2][i] *
                model[county]
                    .parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[(AgeGroup)i] *
                denom_I_Sev_Cr[i] * num_InfectedSevere[county][i];
            // Only set the number of ICU patients here, if the date is not available in the data.
            if (date <= Date(2020, 4, 23) || date >= Date(2024, 7, 21)) {
                model[county].populations[{AgeGroup(i), InfectionState::InfectedCriticalImprovedImmunity}] =
                    immunity_population[2][i] *
                    model[county]
                        .parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[(AgeGroup)i] *
                    denom_I_Sev_Cr[i] * num_icu[county][i];
            }

            // the += is necessary because we already set the previous vaccinated individuals
            model[county].populations[{AgeGroup(i), InfectionState::TemporaryImmuneImprovedImmunity}] +=
                immunity_population[2][i] *
                model[county].parameters.template get<ReducExposedImprovedImmunity<FP>>()[(AgeGroup)i] * denom_E[i] *
                num_timm2[county][i];
        }
        if (std::accumulate(num_InfectedSymptoms[county].begin(), num_InfectedSymptoms[county].end(), 0.0) == 0) {
            log_warning("No infections for vaccinated reported on date " + std::to_string(date.year) + "-" +
                        std::to_string(date.month) + "-" + std::to_string(date.day) + " for region " +
                        std::to_string(region[county]) + ". Population data has not been set.");
        }
    }
    return success();
}

/**
 * @brief Reads confirmed case data from a file and sets it in the model.
 *
 * This function reads transformed RKI data from the specified file and distributes the confirmed case data 
 * across different infection states for regions and age groups in the model. It considers naive, partial, 
 * and improved immunity layers.
 *
 * @tparam Model The type of the model used.
 * @tparam FP Floating point type (default: double).
 *
 * @param[in,out] model Vector of models, each representing a region, where the compartments are updated.
 * @param[in] path Path to the file containing case (RKI) data.
 * @param[in] region Vector of region IDs for which the data is processed.
 * @param[in] date Date for which the confirmed cases are set in the model.
 * @param[in] scaling_factor_inf Vector of scaling factors for confirmed cases.
 * @param[in] immunity_population Vector containing the immunity distribution for naive, partial, and improved immunity layers.
 *
 * @return An IOResult indicating success or failure.
 */
template <class Model, typename FP = double>
IOResult<void> set_confirmed_cases_data(std::vector<Model>& model, const std::string& path,
                                        std::vector<int> const& region, Date date,
                                        const std::vector<FP>& scaling_factor_inf,
                                        const std::vector<std::vector<FP>> immunity_population)
{
    BOOST_OUTCOME_TRY(auto&& case_data, mio::read_confirmed_cases_data(path));
    BOOST_OUTCOME_TRY(
        set_confirmed_cases_data(model, case_data, region, date, scaling_factor_inf, immunity_population));
    return success();
}

/**
 * @brief Extracts the number of individuals in critical condition (ICU) for each region 
 * on a specified date from the provided DIVI data.
 *
 * @tparam FP Floating point type (default: double).
 *
 * @param[in] divi_data Vector of DIVI data entries containing date, region, and ICU information.
 * @param[in] vregion Vector of region IDs for which the data is computed.
 * @param[in] date Date for which the ICU data is computed.
 * @param[out] vnum_icu Output vector containing the number of ICU cases for each region.
 *
 * @return An IOResult indicating success or failure.
 */
template <typename FP = double>
IOResult<void> compute_divi_data(const std::vector<DiviEntry>& divi_data, const std::vector<int>& vregion, Date date,
                                 std::vector<FP>& vnum_icu)
{
    auto max_date_entry = std::max_element(divi_data.begin(), divi_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == divi_data.end()) {
        log_error("DIVI data is empty.");
        return failure(StatusCode::InvalidValue, "DIVI data is empty.");
    }
    auto max_date = max_date_entry->date;
    if (max_date < date) {
        log_error("DIVI data does not contain the specified date.");
        return failure(StatusCode::OutOfRange, "DIVI data does not contain the specified date.");
    }

    for (auto&& entry : divi_data) {
        auto it      = std::find_if(vregion.begin(), vregion.end(), [&entry](auto r) {
            return r == 0 || r == get_region_id(entry);
        });
        auto date_df = entry.date;
        if (it != vregion.end() && date_df == date) {
            auto region_idx      = size_t(it - vregion.begin());
            vnum_icu[region_idx] = entry.num_icu;
        }
    }

    return success();
}

/**
 * @brief Reads DIVI data from a file and computes the ICU data for specified regions and date.
 *
 * @tparam FP Floating point type (default: double).
 *
 * @param[in] path Path to the file containing DIVI data.
 * @param[in] vregion Vector of region IDs for which the data is computed.
 * @param[in] date Date for which the ICU data is computed.
 * @param[out] vnum_icu Output vector containing the number of ICU cases for each region.
 *
 * @return An IOResult indicating success or failure.
 */
template <typename FP = double>
IOResult<void> read_divi_data(const std::string& path, const std::vector<int>& vregion, Date date,
                              std::vector<FP>& vnum_icu)
{
    BOOST_OUTCOME_TRY(auto&& divi_data, mio::read_divi_data(path));
    return compute_divi_data(divi_data, vregion, date, vnum_icu);
}

/**
 * @brief Sets ICU data from DIVI data into the a vector of models, distributed across age groups.
 *
 * This function reads DIVI data from a file, computes the number of individuals in critical condition (ICU) 
 * for each region, and sets these values in the model. The ICU cases are distributed across age groups 
 * using the transition probabilities from severe to critical.
 *
 * @tparam Model The type of the model used.
 * @tparam FP Floating point type (default: double).
 *
 * @param[in,out] model Vector of models, each representing a region, where the ICU population is updated.
 * @param[in] path Path to the file containing DIVI data.
 * @param[in] vregion Vector of region IDs for which the data is computed.
 * @param[in] date Date for which the ICU data is computed.
 * @param[in] scaling_factor_icu Scaling factor for reported ICU cases.
 *
 * @return An IOResult indicating success or failure.
 */
template <class Model, typename FP = double>
IOResult<void> set_divi_data(std::vector<Model>& model, const std::string& path, const std::vector<int>& vregion,
                             Date date, FP scaling_factor_icu)
{
    // DIVI dataset will no longer be updated from CW29 2024 on.
    if (date <= Date(2020, 4, 23) || date >= Date(2024, 7, 21)) {
        log_warning("No DIVI data available for date: {}-{}-{}", date.day, date.month, date.year,
                    ". ICU compartment will be set based on Case data.");
        return success();
    }
    std::vector<FP> sum_mu_I_U(vregion.size(), 0);
    std::vector<std::vector<FP>> mu_I_U{model.size()};
    for (size_t region = 0; region < vregion.size(); region++) {
        auto num_groups = model[region].parameters.get_num_groups();
        for (auto i = AgeGroup(0); i < num_groups; i++) {
            sum_mu_I_U[region] += model[region].parameters.template get<CriticalPerSevere<FP>>()[i] *
                                  model[region].parameters.template get<SeverePerInfectedSymptoms<FP>>()[i];
            mu_I_U[region].push_back(model[region].parameters.template get<CriticalPerSevere<FP>>()[i] *
                                     model[region].parameters.template get<SeverePerInfectedSymptoms<FP>>()[i]);
        }
    }
    std::vector<FP> num_icu(model.size(), 0.0);
    BOOST_OUTCOME_TRY(read_divi_data(path, vregion, date, num_icu));

    for (size_t region = 0; region < vregion.size(); region++) {
        auto num_groups = model[region].parameters.get_num_groups();
        for (auto i = AgeGroup(0); i < num_groups; i++) {
            model[region].populations[{i, InfectionState::InfectedCriticalNaive}] =
                scaling_factor_icu * num_icu[region] * mu_I_U[region][(size_t)i] / sum_mu_I_U[region];
        }
    }

    return success();
}

/**
 * @brief Reads population data from census data.
 * 
 * @param[in] path Path to the population data file.
 * @param[in] vregion Vector of keys representing the regions of interest.
 * @return An IOResult containing a vector of vectors, where each inner vector represents the population
 *         distribution across age groups for a specific region, or an error if the function fails.
 * @see mio::read_population_data
 */
IOResult<std::vector<std::vector<double>>> read_population_data(const std::string& path,
                                                                const std::vector<int>& vregion);

/**
 * @brief Reads population data from a vector of population data entries.
 * 
 * @param[in] population_data Vector of population data entries.
 * @param[in] vregion Vector of keys representing the regions of interest.
 * @return An IOResult containing a vector of vectors, where each inner vector represents the population
 *         distribution across age groups for a specific region, or an error if the function fails.
 * @see mio::read_population_data
 */
IOResult<std::vector<std::vector<double>>> read_population_data(const std::vector<PopulationDataEntry>& population_data,
                                                                const std::vector<int>& vregion);

/**
 * @brief Sets the population data for the given models based on the provided population distribution and immunity levels.
 *
 * @tparam Model The type of the model used.
 * @tparam FP Floating point type (default: double).
 * 
 * @param[in,out] model A vector of models for which population data will be set.
 * @param[in] num_population A 2D vector where each row represents the age group population distribution for a specific region.
 * @param[in] vregion A vector of region identifiers corresponding to the population data.
 * @param[in] immunity_population A 2D vector where each row represents the immunity distribution for a specific region 
 *                                 across different levels of immunity (e.g., naive, partial, improved immunity).
 * 
 * @return An IOResult indicating success or failure.
 */
template <class Model, typename FP = double>
IOResult<void> set_population_data(std::vector<Model>& model, const std::vector<std::vector<FP>>& num_population,
                                   const std::vector<int>& vregion,
                                   const std::vector<std::vector<FP>> immunity_population)
{
    for (size_t region = 0; region < vregion.size(); region++) {
        if (std::accumulate(num_population[region].begin(), num_population[region].end(), 0.0) > 0) {
            auto num_groups = model[region].parameters.get_num_groups();
            for (auto i = AgeGroup(0); i < num_groups; i++) {

                FP SN  = num_population[region][size_t(i)] * immunity_population[0][size_t(i)];
                FP SPI = num_population[region][size_t(i)] * immunity_population[1][size_t(i)];
                FP SII = num_population[region][size_t(i)] - SN - SPI;

                model[region].populations[{i, InfectionState::SusceptibleImprovedImmunity}] = std::max(
                    0.0,
                    FP(SII -
                       (model[region].populations[{i, InfectionState::ExposedImprovedImmunity}] +
                        model[region].populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunity}] +
                        model[region].populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] +
                        model[region].populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] +
                        model[region].populations[{i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] +
                        model[region].populations[{i, InfectionState::InfectedSevereImprovedImmunity}] +
                        model[region].populations[{i, InfectionState::InfectedCriticalImprovedImmunity}] +
                        model[region].populations[{i, InfectionState::DeadImprovedImmunity}] +
                        model[region].populations[{i, InfectionState::TemporaryImmuneImprovedImmunity}])));

                model[region].populations[{i, InfectionState::SusceptiblePartialImmunity}] = std::max(
                    0.0,
                    SPI - model[region].populations[{i, InfectionState::ExposedPartialImmunity}] -
                        model[region].populations[{i, InfectionState::InfectedNoSymptomsPartialImmunity}] -
                        model[region].populations[{i, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}] -
                        model[region].populations[{i, InfectionState::InfectedSymptomsPartialImmunity}] -
                        model[region].populations[{i, InfectionState::InfectedSymptomsPartialImmunityConfirmed}] -
                        model[region].populations[{i, InfectionState::InfectedSeverePartialImmunity}] -
                        model[region].populations[{i, InfectionState::InfectedCriticalPartialImmunity}] -
                        model[region].populations[{i, InfectionState::DeadPartialImmunity}] -
                        model[region].populations[{i, InfectionState::TemporaryImmunePartialImmunity}]);

                model[region].populations.template set_difference_from_group_total<AgeGroup>(
                    {i, InfectionState::SusceptibleNaive}, num_population[region][size_t(i)]);
            }

            for (auto i = AgeGroup(0); i < AgeGroup(6); i++) {
                for (auto j = Index<InfectionState>(0); j < InfectionState::Count; ++j) {
                    if (model[region].populations[{i, j}] < 0) {
                        log_warning("Compartment at age group {}, infection state {}, is negative: {}", size_t(i),
                                    size_t(j), model[region].populations[{i, j}]);
                    }
                }
            }
        }
        else {
            log_warning("No population data available for region " + std::to_string(region) +
                        ". Population data has not been set.");
        }
    }

    return success();
}

/**
 * @brief Reads population data from a file and sets it for the each given model.
 *
 * @tparam Model The type of the model used.
 * 
 * @param[in,out] model A vector of models for which population data will be set.
 * @param[in] path The file path to the population data.
 * @param[in] vregion A vector of region identifiers corresponding to the population data.
 * @param[in] immunity_population A 2D vector where each row represents the immunity distribution for a specific region 
 *                                 across different levels of immunity (e.g., naive, partial, improved).
 * 
 * @return An IOResult indicating success or failure.
 */
template <class Model>
IOResult<void> set_population_data(std::vector<Model>& model, const std::string& path, const std::vector<int>& vregion,
                                   const std::vector<std::vector<double>> immunity_population)
{
    BOOST_OUTCOME_TRY(auto&& num_population, details::read_population_data(path, vregion));
    BOOST_OUTCOME_TRY(set_population_data(model, num_population, vregion, immunity_population));
    return success();
}

/**
 * @brief Sets vaccination data for the given models using provided vaccination (partial, full, and booster) data.
 *
 *
 * @tparam FP Floating point type (default: double).
 *
 * @param[in,out] model A vector of models for which vaccination data will be set.
 * @param[in] vacc_data A vector of VaccinationDataEntry objects containing the vaccination data.
 * @param[in] date The starting date for the simulation.
 * @param[in] vregion A vector of region identifiers corresponding to the vaccination data.
 * @param[in] num_days The number of days for which the simulation runs.
 *
 * @return An IOResult indicating success or failure.
 */
template <typename FP = double>
IOResult<void> set_vaccination_data(std::vector<Model<FP>>& model, const std::vector<VaccinationDataEntry>& vacc_data,
                                    Date date, const std::vector<int>& vregion, int num_days)
{
    auto num_groups = model[0].parameters.get_num_groups();

    auto days_until_effective_n =
        (int)(double)model[0].parameters.template get<DaysUntilEffectivePartialVaccination<FP>>()[AgeGroup(0)];
    auto days_until_effective_pi =
        (int)(double)model[0].parameters.template get<DaysUntilEffectiveImprovedVaccination<FP>>()[AgeGroup(0)];
    auto days_until_effective_ii =
        (int)(double)model[0].parameters.template get<DaysUntilEffectiveBoosterImmunity<FP>>()[AgeGroup(0)];
    // iterate over regions (e.g., counties)
    for (size_t i = 0; i < model.size(); ++i) {
        // iterate over age groups in region
        for (auto g = AgeGroup(0); g < num_groups; ++g) {

            model[i].parameters.template get<DailyPartialVaccinations<FP>>().resize(SimulationDay(num_days + 1));
            model[i].parameters.template get<DailyFullVaccinations<FP>>().resize(SimulationDay(num_days + 1));
            model[i].parameters.template get<DailyBoosterVaccinations<FP>>().resize(SimulationDay(num_days + 1));
            for (auto d = SimulationDay(0); d < SimulationDay(num_days + 1); ++d) {
                model[i].parameters.template get<DailyPartialVaccinations<FP>>()[{g, d}] = 0.0;
                model[i].parameters.template get<DailyFullVaccinations<FP>>()[{g, d}]    = 0.0;
                model[i].parameters.template get<DailyBoosterVaccinations<FP>>()[{g, d}] = 0.0;
            }
        }
    }

    auto max_date_entry = std::max_element(vacc_data.begin(), vacc_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == vacc_data.end()) {
        return failure(StatusCode::InvalidFileFormat, "Vaccination data file is empty.");
    }
    auto max_date = max_date_entry->date;
    if (max_date < offset_date_by_days(date, num_days)) {
        log_error("Vaccination data does not contain all dates. After the last day the data, vaccinations per day are "
                  "set to 0.");
    }

    for (auto&& vacc_data_entry : vacc_data) {
        auto it      = std::find_if(vregion.begin(), vregion.end(), [&vacc_data_entry](auto&& r) {
            return r == 0 || (vacc_data_entry.county_id && vacc_data_entry.county_id == regions::CountyId(r)) ||
                   (vacc_data_entry.state_id && vacc_data_entry.state_id == regions::StateId(r)) ||
                   (vacc_data_entry.district_id && vacc_data_entry.district_id == regions::DistrictId(r));
        });
        auto date_df = vacc_data_entry.date;
        if (it != vregion.end()) {
            auto region_idx = size_t(it - vregion.begin());
            AgeGroup age    = vacc_data_entry.age_group;

            // get daily vaccinations for each layer
            for (size_t d = 0; d < (size_t)num_days + 1; ++d) {
                auto offset_first_date = offset_date_by_days(date, (int)d - days_until_effective_n);
                if (max_date >= offset_first_date) {
                    if (date_df == offset_first_date) {
                        model[region_idx]
                            .parameters.template get<DailyPartialVaccinations<FP>>()[{age, SimulationDay(d)}] =
                            vacc_data_entry.num_vaccinations_partial;
                    }
                }
                else {
                    if (date_df == offset_first_date) {
                        model[region_idx]
                            .parameters.template get<DailyPartialVaccinations<FP>>()[{age, SimulationDay(d)}] = 0;
                    }
                }

                auto offset_full_date = offset_date_by_days(date, (int)d - days_until_effective_pi);
                if (max_date >= offset_full_date) {
                    if (date_df == offset_full_date) {
                        model[region_idx]
                            .parameters.template get<DailyFullVaccinations<FP>>()[{age, SimulationDay(d)}] =
                            vacc_data_entry.num_vaccinations_completed;
                    }
                }
                else {
                    if (date_df == offset_first_date) {
                        model[region_idx]
                            .parameters.template get<DailyFullVaccinations<FP>>()[{age, SimulationDay(d)}] = 0;
                    }
                }

                auto offset_booster_date = offset_date_by_days(date, (int)d - days_until_effective_ii);
                if (max_date >= offset_booster_date) {
                    if (date_df == offset_booster_date) {
                        model[region_idx]
                            .parameters.template get<DailyBoosterVaccinations<FP>>()[{age, SimulationDay(d)}] =
                            vacc_data_entry.num_vaccinations_refreshed_first +
                            vacc_data_entry.num_vaccinations_refreshed_additional;
                    }
                }
                else {
                    if (date_df == offset_first_date) {
                        model[region_idx]
                            .parameters.template get<DailyBoosterVaccinations<FP>>()[{age, SimulationDay(d)}] = 0;
                    }
                }
            }
        }
    }
    return success();
}

/**
 * @brief Sets vaccination data for the given models using vaccination data from a file.
 *
 * This function reads vaccination data from a specified file, and assigns daily vaccination numbers
 * (partial, full, and booster) to each region and age group in the models.
 *
 * @tparam FP Floating point type (default: double).
 *
 * @param[in,out] model A vector of models for which vaccination data will be set.
 * @param[in] path The file path to the vaccination data.
 * @param[in] date The starting date for the simulation.
 * @param[in] vregion A vector of region identifiers corresponding to the vaccination data.
 * @param[in] num_days The number of days for which the simulation runs.
 *
 * @return An IOResult indicating success or failure.
 */
template <typename FP = double>
IOResult<void> set_vaccination_data(std::vector<Model<FP>>& model, const std::string& path, Date date,
                                    const std::vector<int>& vregion, int num_days)
{
    BOOST_OUTCOME_TRY(auto&& vacc_data, read_vaccination_data(path));
    BOOST_OUTCOME_TRY(set_vaccination_data(model, vacc_data, date, vregion, num_days));
    return success();
}

} // namespace details

#ifdef MEMILIO_HAS_HDF5

/**
 * @brief Exports extrapolated real-world data time series for specified regions.
 *
 * This function generates and exports time series data based on the extrapolation and approximation methods
 * used to initialize the model with real-world data. The resulting data represents the initialized states of
 * the model over the specified time range.
 *
 * @tparam Model Type of the model used.
 *
 * @param[in] models A vector of models for which the extrapolated data is set.
 * @param[in] results_dir Path to the directory where the extrapolated results will be saved in a h5 file.
 * @param[in] counties A vector of region identifiers for which the time series will be exported.
 * @param[in] date The starting date of the time series.
 * @param[in] scaling_factor_inf A vector of scaling factors applied to confirmed cases.
 * @param[in] scaling_factor_icu A scaling factor applied to ICU cases.
 * @param[in] num_days The number of days for which will be extrapolated.
 * @param[in] divi_data_path Path to the DIVI ICU data file.
 * @param[in] confirmed_cases_path Path to the confirmed cases data file.
 * @param[in] population_data_path Path to the population data file.
 * @param[in] immunity_population A vector of vectors specifying immunity for each age group and immunity layer.
 * @param[in] vaccination_data_path Path to the vaccination data file (optional).
 *
 * @return An IOResult indicating success or failure.
 */
template <class Model>
IOResult<void> export_input_data_county_timeseries(
    std::vector<Model> models, const std::string& results_dir, const std::vector<int>& counties, Date date,
    const std::vector<double>& scaling_factor_inf, const double scaling_factor_icu, const int num_days,
    const std::string& divi_data_path, const std::string& confirmed_cases_path, const std::string& population_data_path,
    const std::vector<std::vector<double>> immunity_population, const std::string& vaccination_data_path = "")
{
    const auto num_groups = (size_t)models[0].parameters.get_num_groups();
    assert(scaling_factor_inf.size() == num_groups);
    assert(num_groups == ConfirmedCasesDataEntry::age_group_names.size());
    assert(models.size() == counties.size());
    std::vector<TimeSeries<double>> extrapolated_data(
        models.size(), TimeSeries<double>::zero(num_days + 1, (size_t)InfectionState::Count * num_groups));

    BOOST_OUTCOME_TRY(auto&& case_data, read_confirmed_cases_data(confirmed_cases_path));
    BOOST_OUTCOME_TRY(auto&& population_data, details::read_population_data(population_data_path, counties));

    // empty vector if set_vaccination_data is not set
    std::vector<VaccinationDataEntry> vacc_data;
    if (!vaccination_data_path.empty()) {
        BOOST_OUTCOME_TRY(vacc_data, read_vaccination_data(vaccination_data_path));
    }

    for (int t = 0; t <= num_days; ++t) {
        auto offset_day = offset_date_by_days(date, t);

        if (!vaccination_data_path.empty()) {
            BOOST_OUTCOME_TRY(details::set_vaccination_data(models, vacc_data, offset_day, counties, num_days));
        }

        // TODO: Reuse more code, e.g., set_divi_data (in secir) and a set_divi_data (here) only need a different ModelType.
        // TODO: add option to set ICU data from confirmed cases if DIVI or other data is not available.
        BOOST_OUTCOME_TRY(details::set_divi_data(models, divi_data_path, counties, offset_day, scaling_factor_icu));

        BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(models, case_data, counties, offset_day, scaling_factor_inf,
                                                            immunity_population));

        BOOST_OUTCOME_TRY(details::set_population_data(models, population_data, counties, immunity_population));

        for (size_t r = 0; r < counties.size(); r++) {
            extrapolated_data[r][t] = models[r].get_initial_values();
            // in set_population_data the number of death individuals is subtracted from the SusceptibleImprovedImmunity compartment.
            // Since we should be independent whether we consider them or not, we add them back here before we save the data.
            for (size_t age = 0; age < num_groups; age++) {
                extrapolated_data[r][t][(size_t)InfectionState::SusceptibleImprovedImmunity +
                                        age * (size_t)InfectionState::Count] +=
                    extrapolated_data[r][t][(size_t)InfectionState::DeadNaive + age * (size_t)InfectionState::Count];
            }
        }
    }
    BOOST_OUTCOME_TRY(save_result(extrapolated_data, counties, static_cast<int>(num_groups),
                                  path_join(results_dir, "Results_rki.h5")));

    auto extrapolated_rki_data_sum = sum_nodes(std::vector<std::vector<TimeSeries<double>>>{extrapolated_data});
    BOOST_OUTCOME_TRY(save_result({extrapolated_rki_data_sum[0][0]}, {0}, static_cast<int>(num_groups),
                                  path_join(results_dir, "Results_rki_sum.h5")));

    return success();
}

#else
template <class Model>
IOResult<void> export_input_data_county_timeseries(std::vector<Model>, const std::string&, const std::vector<int>&,
                                                   Date, const std::vector<double>&, const double, const int,
                                                   const std::string&, const std::string&, const std::string&,
                                                   const std::vector<std::vector<double>>, const std::string&)
{
    mio::log_warning("HDF5 not available. Cannot export time series of extrapolated real data.");
    return success();
}

#endif //MEMILIO_HAS_HDF5

/**
 * @brief Reads input data for specified counties and initializes the model accordingly.
 *
 * This function loads real-world data for specified counties and initializes
 * the corresponding model compartments. Optionally, it can extrapolate real-world data
 * using the export_input_data_county_timeseries function.
 *
 * @tparam Model The model type used.
 *
 * @param[in,out] model Vector of model objects that will be initialized with the input data.
 * @param[in] date The starting date for the simulation.
 * @param[in] county Vector of county IDs for which the data is read.
 * @param[in] scaling_factor_inf Vector of scaling factors for confirmed cases.
 * @param[in] scaling_factor_icu Scaling factor for ICU data.
 * @param[in] dir Path to the directory containing input data files.
 * @param[in] num_days The number of days for which the simulation runs.
 * @param[in] immunity_population Vector of vectors representing immunity proportions for each age group and immunity layer.
 * @param[in] export_time_series Boolean flag indicating whether to export time series of extrapolated data (default: false).
 *
 * @return An IOResult indicating success or failure.
 */
template <class Model>
IOResult<void> read_input_data_county(std::vector<Model>& model, Date date, const std::vector<int>& county,
                                      const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                      const std::string& dir, int num_days,
                                      const std::vector<std::vector<double>> immunity_population,
                                      bool export_time_series = false)
{
    BOOST_OUTCOME_TRY(details::set_vaccination_data(
        model, path_join(dir, "pydata/Germany", "all_county_ageinf_vacc_ma7.json"), date, county, num_days));

    // TODO: Reuse more code, e.g., set_divi_data (in secir) and a set_divi_data (here) only need a different ModelType.
    // TODO: add option to set ICU data from confirmed cases if DIVI or other data is not available.
    BOOST_OUTCOME_TRY(details::set_divi_data(model, path_join(dir, "pydata/Germany", "county_divi_ma7.json"), county,
                                             date, scaling_factor_icu));

    BOOST_OUTCOME_TRY(
        details::set_confirmed_cases_data(model, path_join(dir, "pydata/Germany", "cases_all_county_age_ma7.json"),
                                          county, date, scaling_factor_inf, immunity_population));
    BOOST_OUTCOME_TRY(details::set_population_data(
        model, path_join(dir, "pydata/Germany", "county_current_population.json"), county, immunity_population));

    if (export_time_series) {
        // Use only if extrapolated real data is needed for comparison. EXPENSIVE !
        // Run time equals run time of the previous functions times the num_days !
        // (This only represents the vectorization of the previous function over all simulation days...)
        log_info("Exporting time series of extrapolated real data. This may take some minutes. "
                 "For simulation runs over the same time period, deactivate it.");
        BOOST_OUTCOME_TRY(export_input_data_county_timeseries(
            model, dir, county, date, scaling_factor_inf, scaling_factor_icu, num_days,
            path_join(dir, "pydata/Germany", "county_divi_ma7.json"),
            path_join(dir, "pydata/Germany", "cases_all_county_age_ma7.json"),
            path_join(dir, "pydata/Germany", "county_current_population.json"), immunity_population,
            path_join(dir, "pydata/Germany", "all_county_ageinf_vacc_ma7.json")));
    }

    return success();
}

/**
 * @brief Reads compartments for geographic units at a specified date from data files.
 *
 * This function estimates all compartments from available data using the provided model parameters.
 *
 * @tparam Model The type of tmodel used.
 *
 * @param[in,out] model Vector of models, one per county, to be initialized with data.
 * @param[in] date Date for which the data should be read.
 * @param[in] node_ids Vector of IDs of the units for which data is read.
 * @param[in] scaling_factor_inf Vector of scaling factors for confirmed cases.
 * @param[in] scaling_factor_icu Scaling factor for ICU cases.
 * @param[in] data_dir Directory containing the input data files.
 * @param[in] num_days Number of days to simulate.
 * @param[in] immunity_population Matrix containing immunity proportions for each age group and immunity layer.
 * @param[in] export_time_series Boolean flag indicating whether to export time series of extrapolated data (default: false).
 *
 * @return An IOResult indicating success or failure.
    */
template <class Model>
IOResult<void> read_input_data(std::vector<Model>& model, Date date, const std::vector<int>& node_ids,
                               const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                               const std::string& data_dir, int num_days,
                               const std::vector<std::vector<double>> immunity_population,
                               bool export_time_series = false)
{

    BOOST_OUTCOME_TRY(
        details::set_vaccination_data(model, path_join(data_dir, "vaccination_data.json"), date, node_ids, num_days));

    // TODO: Reuse more code, e.g., set_divi_data (in secir) and a set_divi_data (here) only need a different ModelType.
    // TODO: add option to set ICU data from confirmed cases if DIVI or other data is not available.
    BOOST_OUTCOME_TRY(
        details::set_divi_data(model, path_join(data_dir, "critical_cases.json"), node_ids, date, scaling_factor_icu));

    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(model, path_join(data_dir, "confirmed_cases.json"), node_ids,
                                                        date, scaling_factor_inf, immunity_population));
    BOOST_OUTCOME_TRY(details::set_population_data(model, path_join(data_dir, "population_data.json"), node_ids,
                                                   immunity_population));

    if (export_time_series) {
        // Use only if extrapolated real data is needed for comparison. EXPENSIVE !
        // Run time equals run time of the previous functions times the num_days !
        // (This only represents the vectorization of the previous function over all simulation days...)
        log_info("Exporting time series of extrapolated real data. This may take some minutes. "
                 "For simulation runs over the same time period, deactivate it.");
        BOOST_OUTCOME_TRY(export_input_data_county_timeseries(
            model, data_dir, node_ids, date, scaling_factor_inf, scaling_factor_icu, num_days,
            path_join(data_dir, "critical_cases.json"), path_join(data_dir, "confirmed_cases.json"),
            path_join(data_dir, "population_data.json"), immunity_population,
            path_join(data_dir, "vaccination_data.json")));
    }

    return success();
}

} // namespace osecirts
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // MIO_ODE_SECIRTS_PARAMETERS_IO_H
