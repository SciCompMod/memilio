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
#include "memilio/geography/regions.h"
#include "memilio/math/eigen_util.h"
#include "memilio/math/math_utils.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/parameters_io.h"
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
    const std::vector<double>& scaling_factor_inf, const size_t layer);
/**@}*/

/**
 * @brief Sets the confirmed cases data in the model considering different immunity layers.
 *
 * This function distributes confirmed case data across infection states for regions and age groups
 * in the model. It considers different levels of immunity (naive, partial, and improved).
 *
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
template <typename FP>
IOResult<void>
set_confirmed_cases_data(Model<FP>& model, const std::vector<ConfirmedCasesDataEntry>& case_data,
                         const int region, Date date, const std::vector<FP>& scaling_factor_inf,
                         const std::vector<std::vector<double>>& immunity_population)
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

    std::vector<double> mu_C_R;
    std::vector<double> mu_I_H;
    std::vector<double> mu_H_U;

    std::vector<double> reduc_t_Infected;
    std::vector<double> reduc_Exposed;
    std::vector<double> reduc_InfectedSymptoms;
    std::vector<double> reduc_icu_death;

    std::vector<double> num_InfectedSymptoms(num_age_groups, 0.0);
    std::vector<double> num_death(num_age_groups, 0.0);
    std::vector<double> num_Exposed(num_age_groups, 0.0);
    std::vector<double> num_InfectedNoSymptoms(num_age_groups, 0.0);
    std::vector<double> num_InfectedSevere(num_age_groups, 0.0);
    std::vector<double> num_icu(num_age_groups, 0.0);
    std::vector<double> num_timm1(num_age_groups, 0.0);
    std::vector<double> num_timm2(num_age_groups, 0.0);

    std::vector<double> denom_E(num_age_groups, 0.0);
    std::vector<double> denom_I_NS(num_age_groups, 0.0);
    std::vector<double> denom_I_Sy(num_age_groups, 0.0);
    std::vector<double> denom_I_Sev_Cr(num_age_groups, 0.0);

    // calculate the denominators to split the reported case numbers to the different immunity layers.
    for (size_t group = 0; group < num_age_groups; group++) {
        denom_E[group] =
            1 / (immunity_population[0][group] +
                    immunity_population[1][group] *
                        model.parameters.template get<ReducExposedPartialImmunity<FP>>()[(AgeGroup)group] +
                    immunity_population[2][group] *
                        model.parameters.template get<ReducExposedImprovedImmunity<FP>>()[(AgeGroup)group]);

        denom_I_NS[group] =
            1 / (immunity_population[0][group] +
                    immunity_population[1][group] *
                        model.parameters.template get<ReducExposedPartialImmunity<FP>>()[(AgeGroup)group] +
                    immunity_population[2][group] *
                        model.parameters.template get<ReducExposedImprovedImmunity<FP>>()[(AgeGroup)group]);

        denom_I_Sy[group] =
            1 / (immunity_population[0][group] +
                    immunity_population[1][group] *
                        model
                            .parameters.template get<ReducInfectedSymptomsPartialImmunity<FP>>()[(AgeGroup)group] +
                    immunity_population[2][group] *
                        model
                            .parameters.template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[(AgeGroup)group]);

        denom_I_Sev_Cr[group] =
            1 / (immunity_population[0][group] +
                    immunity_population[1][group] *
                        model.parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[(
                            AgeGroup)group] +
                    immunity_population[2][group] *
                        model.parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[(
                            AgeGroup)group]);
    }

    /*----------- Naive immunity -----------*/
    for (size_t group = 0; group < num_age_groups; group++) {
        t_Exposed.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeExposed<FP>>()[(AgeGroup)group])));
        t_InfectedNoSymptoms.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedNoSymptoms<FP>>()[(AgeGroup)group])));
        t_InfectedSymptoms.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedSymptoms<FP>>()[(AgeGroup)group])));
        t_InfectedSevere.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedSevere<FP>>()[(AgeGroup)group])));
        t_InfectedCritical.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedCritical<FP>>()[(AgeGroup)group])));
        t_imm_interval_i.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeTemporaryImmunityPI<FP>>()[(AgeGroup)group])));

        mu_C_R.push_back(
            model.parameters.template get<RecoveredPerInfectedNoSymptoms<FP>>()[(AgeGroup)group]);
        mu_I_H.push_back(
            model.parameters.template get<SeverePerInfectedSymptoms<FP>>()[(AgeGroup)group]);
        mu_H_U.push_back(
            model.parameters.template get<CriticalPerSevere<FP>>()[(AgeGroup)group]);

        reduc_t_Infected.push_back(
            model.parameters.template get<ReducTimeInfectedMild<FP>>()[(AgeGroup)group]);
        reduc_Exposed.push_back(
            model.parameters.template get<ReducExposedPartialImmunity<FP>>()[(AgeGroup)group]);
        reduc_InfectedSymptoms.push_back(
            model.parameters.template get<ReducInfectedSymptomsPartialImmunity<FP>>()[(AgeGroup)group]);
        reduc_icu_death.push_back(
            model.parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[(AgeGroup)group]);
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
    if (std::accumulate(num_InfectedSymptoms.begin(), num_InfectedSymptoms.end(), FP(0.0),
                        [](const FP& a, const FP& b) {
                            return evaluate_intermediate<FP>(a + b);
                        }) == 0) {
        log_warning(
            "No infections for unvaccinated reported on date {} for region {}. Population data has not been set.",
            date, region);
    }

    /*----------- PARTIAL Immunity -----------*/
    reduc_Exposed.clear();
    reduc_InfectedSymptoms.clear();
    reduc_icu_death.clear();

    num_InfectedSymptoms   = std::vector<double>(num_age_groups, 0.0);
    num_death              = std::vector<double>(num_age_groups, 0.0);
    num_Exposed            = std::vector<double>(num_age_groups, 0.0);
    num_InfectedNoSymptoms = std::vector<double>(num_age_groups, 0.0);
    num_InfectedSevere     = std::vector<double>(num_age_groups, 0.0);
    num_icu                = std::vector<double>(num_age_groups, 0.0);

    for (size_t group = 0; group < num_age_groups; group++) {
        reduc_Exposed.push_back(
            model.parameters.template get<ReducExposedPartialImmunity<FP>>()[(AgeGroup)group]);
        reduc_InfectedSymptoms.push_back(
            model.parameters.template get<ReducInfectedSymptomsPartialImmunity<FP>>()[(AgeGroup)group]);
        reduc_icu_death.push_back(
            model.parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[(AgeGroup)group]);
    }

    BOOST_OUTCOME_TRY(compute_confirmed_cases_data(case_data, region, date, num_Exposed, num_InfectedNoSymptoms,
                                                   num_InfectedSymptoms, num_InfectedSevere, num_icu, num_death, num_timm1,
                                                   t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
                                                   t_InfectedCritical, t_imm_interval_i, reduc_t_Infected, reduc_Exposed,
                                                   reduc_InfectedSymptoms, reduc_icu_death, mu_C_R, mu_I_H, mu_H_U, 
                                                   scaling_factor_inf, 1));

    size_t num_groups = (size_t)model.parameters.get_num_groups();
    for (size_t i = 0; i < num_groups; i++) {
        model.populations[{AgeGroup(i), InfectionState::ExposedPartialImmunity}] =
            immunity_population[1][i] *
            model.parameters.template get<ReducExposedPartialImmunity<FP>>()[(AgeGroup)i] * denom_E[i] *
            num_Exposed[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsPartialImmunity}] =
            immunity_population[1][i] *
            model.parameters.template get<ReducExposedPartialImmunity<FP>>()[(AgeGroup)i] * denom_I_NS[i] *
            num_InfectedNoSymptoms[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}] = 0;
        model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsPartialImmunity}] =
            immunity_population[1][i] *
            model.parameters.template get<ReducInfectedSymptomsPartialImmunity<FP>>()[(AgeGroup)i] *
            denom_I_Sy[i] * num_InfectedSymptoms[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsPartialImmunityConfirmed}] = 0;
        model.populations[{AgeGroup(i), InfectionState::InfectedSeverePartialImmunity}] =
            immunity_population[1][i] *
            model
                .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[(AgeGroup)i] *
            denom_I_Sev_Cr[i] * num_InfectedSevere[i];
        // Only set the number of ICU patients here, if the date is not available in the data.
        if (!is_divi_data_available(date)) {
            model.populations[{AgeGroup(i), InfectionState::InfectedCriticalPartialImmunity}] =
                immunity_population[1][i] *
                model
                    .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[(AgeGroup)i] *
                denom_I_Sev_Cr[i] * num_icu[i];
        }
        // the += is necessary because we already set the previous vaccinated individuals
        model.populations[{AgeGroup(i), InfectionState::TemporaryImmunePartialImmunity}] +=
            immunity_population[1][i] *
            model.parameters.template get<ReducExposedPartialImmunity<FP>>()[(AgeGroup)i] * denom_E[i] *
            num_timm1[i];
    }
    if (std::accumulate(num_InfectedSymptoms.begin(), num_InfectedSymptoms.end(), FP(0.0),
                        [](const FP& a, const FP& b) {
                            return evaluate_intermediate<FP>(a + b);
                        }) == 0) {

        log_warning("No infections for partially vaccinated reported on date {} for region {}. "
                    "Population data has not been set.",
                    date, region);
    }

    /*----------- Improved Immunity -----------*/
    reduc_Exposed.clear();
    reduc_InfectedSymptoms.clear();
    reduc_icu_death.clear();

    num_InfectedSymptoms   = std::vector<double>(num_age_groups, 0.0);
    num_death              = std::vector<double>(num_age_groups, 0.0);
    num_Exposed            = std::vector<double>(num_age_groups, 0.0);
    num_InfectedNoSymptoms = std::vector<double>(num_age_groups, 0.0);
    num_InfectedSevere     = std::vector<double>(num_age_groups, 0.0);
    num_icu                = std::vector<double>(num_age_groups, 0.0);

    for (size_t group = 0; group < num_age_groups; group++) {
        reduc_Exposed.push_back(
            model.parameters.template get<ReducExposedImprovedImmunity<FP>>()[(AgeGroup)group]);
        reduc_InfectedSymptoms.push_back(
            model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[(AgeGroup)group]);
        reduc_icu_death.push_back(
            model.parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[(AgeGroup)group]);
    }

    BOOST_OUTCOME_TRY(compute_confirmed_cases_data(case_data, region, date, num_Exposed, num_InfectedNoSymptoms,
                                                   num_InfectedSymptoms, num_InfectedSevere, num_icu, num_death, num_timm2,
                                                   t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
                                                   t_InfectedCritical, t_imm_interval_i, reduc_t_Infected, reduc_Exposed,
                                                   reduc_InfectedSymptoms, reduc_icu_death, mu_C_R, mu_I_H, mu_H_U, 
                                                   scaling_factor_inf, 2));

    size_t num_groups = (size_t)model.parameters.get_num_groups();
    for (size_t i = 0; i < num_groups; i++) {
        model.populations[{AgeGroup(i), InfectionState::ExposedImprovedImmunity}] =
            immunity_population[2][i] *
            model.parameters.template get<ReducExposedImprovedImmunity<FP>>()[(AgeGroup)i] * denom_E[i] *
            num_Exposed[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsImprovedImmunity}] =
            immunity_population[2][i] *
            model.parameters.template get<ReducExposedImprovedImmunity<FP>>()[(AgeGroup)i] * denom_I_NS[i] *
            num_InfectedNoSymptoms[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] = 0;
        model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsImprovedImmunity}] =
            immunity_population[2][i] *
            model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[(AgeGroup)i] *
            denom_I_Sy[i] * num_InfectedSymptoms[i];
        model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] = 0;
        model.populations[{AgeGroup(i), InfectionState::InfectedSevereImprovedImmunity}] =
            immunity_population[2][i] *
            model
                .parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[(AgeGroup)i] *
            denom_I_Sev_Cr[i] * num_InfectedSevere[i];
        // Only set the number of ICU patients here, if the date is not available in the data.
        if (!is_divi_data_available(date)) {
            model.populations[{AgeGroup(i), InfectionState::InfectedCriticalImprovedImmunity}] =
                immunity_population[2][i] *
                model
                    .parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[(AgeGroup)i] *
                denom_I_Sev_Cr[i] * num_icu[i];
        }

        // the += is necessary because we already set the previous vaccinated individuals
        model.populations[{AgeGroup(i), InfectionState::TemporaryImmuneImprovedImmunity}] +=
            immunity_population[2][i] *
            model.parameters.template get<ReducExposedImprovedImmunity<FP>>()[(AgeGroup)i] * denom_E[i] *
            num_timm2[i];
    }
    if (std::accumulate(num_InfectedSymptoms.begin(), num_InfectedSymptoms.end(), FP(0.0),
                        [](const FP& a, const FP& b) {
                            return evaluate_intermediate<FP>(a + b);
                        }) == 0) {
        log_warning("No infections for vaccinated reported on date {} for region {}. "
                    "Population data has not been set.",
                    date, region);
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
 * @tparam FP Floating point type (default: double).
 *
 * @param[in,out] model Vector of models, each representing a region, where the compartments are updated.
 * @param[in] path Path to the file containing case (RKI) data.
 * @param[in] vregion Vector of region IDs for which the data is processed.
 * @param[in] date Date for which the confirmed cases are set in the model.
 * @param[in] scaling_factor_inf Vector of scaling factors for confirmed cases.
 * @param[in] immunity_population Vector containing the immunity distribution for naive, partial, and improved immunity layers.
 *
 * @return An IOResult indicating success or failure.
 */
template <typename FP>
IOResult<void> set_confirmed_cases_data(mio::VectorRange<Model<FP>>& model, const std::string& path,
                                        std::vector<int> const& vregion, Date date,
                                        const std::vector<FP>& scaling_factor_inf,
                                        const std::vector<std::vector<double>>& immunity_population)
{
    BOOST_OUTCOME_TRY(auto&& case_data, mio::read_confirmed_cases_data(path));

    // sort case_data into regions and ignore once with no region associated
    std::vector<std::vector<ConfirmedCasesDataEntry>> vcase_data{model.size()};
    for (auto&& entry : case_data) {
        auto it = std::find_if(vregion.begin(), vregion.end(), [&entry](auto r) {
            return r == 0 || get_region_id(entry) == r;
        });
        if (it != vregion.end()) {
            auto region_idx = size_t(it - vregion.begin());
            vcase_data[region_idx].pushback(entry);
        }
    }

    for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
        BOOST_OUTCOME_TRY(
            set_confirmed_cases_data(model[region_idx], vcase_data[region_idx], vregion[region_idx], date, scaling_factor_inf, immunity_population));
    }

    return success();
}

/**
 * @brief Sets the population data for the given models based on the provided population distribution and immunity levels.
 *
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
template <typename FP>
IOResult<void> set_population_data(Model<FP>& model, const std::vector<double>& num_population,
                                   const int region, const std::vector<std::vector<double>>& immunity_population)
{
    if (std::accumulate(num_population.begin(), num_population.end(), double(0.0),
                        [](const double& a, const double& b) {
                            return evaluate_intermediate<double>(a + b);
                        }) <= 0)
    {    
        log_warning("No population data available for region " + std::to_string(region) +
                    ". Population data has not been set.");
        return success();
    }

    auto num_groups = model.parameters.get_num_groups();
    for (auto i = AgeGroup(0); i < num_groups; i++) {

        FP SN  = num_population[size_t(i)] * immunity_population[0][size_t(i)];
        FP SPI = num_population[size_t(i)] * immunity_population[1][size_t(i)];
        FP SII = num_population[size_t(i)] - SN - SPI;

        model.populations[{i, InfectionState::SusceptibleImprovedImmunity}] = std::max(
            0.0,
            FP(SII -
                (model.populations[{i, InfectionState::ExposedImprovedImmunity}] +
                model.populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunity}] +
                model.populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] +
                model.populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] +
                model.populations[{i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] +
                model.populations[{i, InfectionState::InfectedSevereImprovedImmunity}] +
                model.populations[{i, InfectionState::InfectedCriticalImprovedImmunity}] +
                model.populations[{i, InfectionState::DeadImprovedImmunity}] +
                model.populations[{i, InfectionState::TemporaryImmuneImprovedImmunity}])));

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

/**
 * @brief Reads population data from a file and sets it for the each given model.
 *
 * @param[in,out] model A vector of models for which population data will be set.
 * @param[in] path The file path to the population data.
 * @param[in] vregion A vector of region identifiers corresponding to the population data.
 * @param[in] immunity_population A 2D vector where each row represents the immunity distribution for a specific region
 *                                 across different levels of immunity (e.g., naive, partial, improved).
 *
 * @return An IOResult indicating success or failure.
 */
template <typename FP>
IOResult<void> set_population_data(mio::VectorRange<Model<FP>>& model, const std::string& path, const std::vector<int>& vregion,
                                   const std::vector<std::vector<double>>& immunity_population)
{
    BOOST_OUTCOME_TRY(auto&& num_population, mio::read_population_data(path, vregion));

    for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
        BOOST_OUTCOME_TRY(set_population_data(model[region_idx], num_population[region_idx], vregion[region_idx], immunity_population));
    }
    return success();
}

template <typename FP>
IOResult<void> set_vaccination_data(Model<FP>& model, const VaccinationDataEntry& vacc_data_entry,
                                    Date date, int num_days, Date max_date, const int days_until_effective_n, 
                                    const int days_until_effective_pi, const int days_until_effective_ii)
{
    auto num_groups = model.parameters.get_num_groups();

    auto date_df = vacc_data_entry.date;
    AgeGroup age    = vacc_data_entry.age_group;

    // get daily vaccinations for each layer
    for (size_t d = 0; d < (size_t)num_days + 1; ++d) {
        auto offset_first_date = offset_date_by_days(date, (int)d - days_until_effective_n);
        if (max_date >= offset_first_date) {
            if (date_df == offset_first_date) {
                model.parameters.template get<DailyPartialVaccinations<FP>>()[{age, SimulationDay(d)}] =
                    vacc_data_entry.num_vaccinations_partial;
            }
        }
        else {
            if (date_df == offset_first_date) {
                model.parameters.template get<DailyPartialVaccinations<FP>>()[{age, SimulationDay(d)}] = 0;
            }
        }

        auto offset_full_date = offset_date_by_days(date, (int)d - days_until_effective_pi);
        if (max_date >= offset_full_date) {
            if (date_df == offset_full_date) {
                model.parameters.template get<DailyFullVaccinations<FP>>()[{age, SimulationDay(d)}] =
                    vacc_data_entry.num_vaccinations_completed;
            }
        }
        else {
            if (date_df == offset_first_date) {
                model.parameters.template get<DailyFullVaccinations<FP>>()[{age, SimulationDay(d)}] = 0;
            }
        }

        auto offset_booster_date = offset_date_by_days(date, (int)d - days_until_effective_ii);
        if (max_date >= offset_booster_date) {
            if (date_df == offset_booster_date) {
                model.parameters.template get<DailyBoosterVaccinations<FP>>()[{age, SimulationDay(d)}] =
                    vacc_data_entry.num_vaccinations_refreshed_first +
                    vacc_data_entry.num_vaccinations_refreshed_additional;
            }
        }
        else {
            if (date_df == offset_first_date) {
                model.parameters.template get<DailyBoosterVaccinations<FP>>()[{age, SimulationDay(d)}] = 0;
            }
        }
    }
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
 * @param[in] num_days The number of days for which the simulation runs.
 *
 * @return An IOResult indicating success or failure.
 */
template <typename FP>
IOResult<void> set_vaccination_data(Model<FP>& model, const std::vector<VaccinationDataEntry>& vacc_data,
                                    Date date, int num_days)
{

    auto max_date_entry = std::max_element(vacc_data.begin(), vacc_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == vacc_data.end()) {
        return failure(StatusCode::InvalidFileFormat, "Vaccination data file is empty.");
    }
    auto max_date       = max_date_entry->date;
    auto min_date_entry = std::min_element(vacc_data.begin(), vacc_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    auto min_date       = min_date_entry->date;
    if (min_date > date || max_date < offset_date_by_days(date, num_days)) {
        log_warning("Vaccination data only available from {} to {}. "
                    "For days before and after, vaccinations will be set to 0.",
                    min_date, max_date);
    }

    auto days_until_effective_n =
        (int)(double)model.parameters.template get<DaysUntilEffectivePartialVaccination<FP>>()[AgeGroup(0)];
    auto days_until_effective_pi =
        (int)(double)model.parameters.template get<DaysUntilEffectiveImprovedVaccination<FP>>()[AgeGroup(0)];
    auto days_until_effective_ii =
        (int)(double)model.parameters.template get<DaysUntilEffectiveBoosterImmunity<FP>>()[AgeGroup(0)];

    for (auto&& vacc_data_entry : vacc_data) {
        BOOST_OUTCOME_TRY(set_vaccination_data(model, vacc_data_entry, date, num_days, max_date,
                                                days_until_effective_n, days_until_effective_pi, days_until_effective_ii));
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
template <typename FP>
IOResult<void> set_vaccination_data(mio::VectorRange<Model<FP>>& model, const std::string& path, Date date,
                                    const std::vector<int>& vregion, int num_days)
{
    // Set vaccination data to 0 for all models
    for (auto& m : model) {
        m.parameters.template get<DailyPartialVaccinations<FP>>().resize(SimulationDay(num_days + 1));
        m.parameters.template get<DailyFullVaccinations<FP>>().resize(SimulationDay(num_days + 1));
        m.parameters.template get<DailyBoosterVaccinations<FP>>().resize(SimulationDay(num_days + 1));
        for (auto d = SimulationDay(0); d < SimulationDay(num_days + 1); ++d) {
            for (auto a = AgeGroup(0); a < m.parameters.get_num_groups(); ++a) {
                m.parameters.template get<DailyPartialVaccinations<FP>>()[{a, d}] = 0.0;
                m.parameters.template get<DailyFullVaccinations<FP>>()[{a, d}]    = 0.0;
                m.parameters.template get<DailyBoosterVaccinations<FP>>()[{a, d}] = 0.0;
            }
        }
    }

    // Check if vaccination data is available for the given date range
    auto end_date = offset_date_by_days(date, num_days);
    if (!is_vaccination_data_available(date, end_date)) {
        log_warning("No vaccination data available in range from {} to {}. "
                    "Vaccination data will be set to 0.",
                    date, end_date);
        return success();
    }

    BOOST_OUTCOME_TRY(auto&& vacc_data, read_vaccination_data(path));
    
    // Sort case_data into regions and ignore once with no region associated
    std::vector<std::vector<VaccinationDataEntry>> vvacc_data{model.size()};
    for (auto&& vacc_data_entry : vacc_data) {
        auto it      = std::find_if(vregion.begin(), vregion.end(), [&vacc_data_entry](auto&& r) {
            return r == 0 || (vacc_data_entry.county_id && vacc_data_entry.county_id == regions::CountyId(r)) ||
                   (vacc_data_entry.state_id && vacc_data_entry.state_id == regions::StateId(r)) ||
                   (vacc_data_entry.district_id && vacc_data_entry.district_id == regions::DistrictId(r));
        });
        if (it != vregion.end()) {
            auto region_idx = size_t(it - vregion.begin());
            vvacc_data[region_idx].pushback(vacc_data_entry);
        }
    }

    for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
        BOOST_OUTCOME_TRY(set_vaccination_data(model[region_idx], vacc_data[region_idx], date, num_days));
    }
    
    return success();
}

/**
 * @brief Sets ICU data from DIVI data into the a vector of models, distributed across age groups.
 *
 * This function reads DIVI data from a file, computes the number of individuals in critical condition (ICU)
 * for each region, and sets these values in the model. The ICU cases are distributed across age groups
 * using the transition probabilities from severe to critical.
 * @tparam FP Floating point type (default: double).
 *
 * @param[in,out] model Vector of models, each representing a region, where the ICU population is updated.
 * @param[in] num_icu icu data
 * @param[in] scaling_factor_icu Scaling factor for reported ICU cases.
 *
 * @return An IOResult indicating success or failure.
 */
template <typename FP>
IOResult<void> set_divi_data(Model<FP>& model, const double num_icu, double scaling_factor_icu)
{
    FP sum_mu_I_U = 0;
    std::vector<FP> mu_I_U;
    auto num_groups = model.parameters.get_num_groups();
    for (auto i = AgeGroup(0); i < num_groups; i++) {
        sum_mu_I_U += model.parameters.template get<CriticalPerSevere<FP>>()[i] *
                                model.parameters.template get<SeverePerInfectedSymptoms<FP>>()[i];
        mu_I_U.push_back(model.parameters.template get<CriticalPerSevere<FP>>()[i] *
                                    model.parameters.template get<SeverePerInfectedSymptoms<FP>>()[i]);
    }

    for (auto i = AgeGroup(0); i < num_groups; i++) {
        model.populations[{i, InfectionState::InfectedCriticalNaive}] =
            scaling_factor_icu * num_icu * mu_I_U[(size_t)i] / sum_mu_I_U;
    }

    return success();
}

/**
 * @brief sets populations data from DIVI register into Model
 * @param[in, out] model vector of objects in which the data is set
 * @param[in] path Path to transformed DIVI file
 * @param[in] vregion vector of keys of the regions of interest
 * @param[in] date Date for which the arrays are initialized
 * @param[in] scaling_factor_icu factor by which to scale the icu cases of divi data
 */
template <class FP>
IOResult<void> set_divi_data(mio::VectorRange<Model<FP>>& model, const std::string& path, const std::vector<int>& vregion,
                             Date date, double scaling_factor_icu)
{
    // DIVI dataset will no longer be updated from CW29 2024 on.
    if (!is_divi_data_available(date)) {
        log_warning("No DIVI data available for date: {}. "
                    "ICU compartment will be set based on Case data.",
                    date);
        return success();
    }
    BOOST_OUTCOME_TRY(auto&& num_icu, read_divi_data(path, vregion, date));

    for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
        BOOST_OUTCOME_TRY(set_divi_data(model[region_idx], num_icu[region_idx], vregion[region_idx], scaling_factor_icu));
    }

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
 * @param[in] models A vector of models for which the extrapolated data is set.
 * @param[in] results_dir Path to the directory where the extrapolated results will be saved in a h5 file.
 * @param[in] date Date for which the data should be read.
 * @param[in] node_ids Vector of keys of the node_ids of interest.
 * @param[in] scaling_factor_inf A vector of scaling factors applied to confirmed cases.
 * @param[in] scaling_factor_icu A scaling factor applied to ICU cases.
 * @param[in] num_days The number of days for which will be extrapolated.
 * @param[in] pydata_dir Directory that contains the data files.
 * @param[in] immunity_population A vector of vectors specifying immunity for each age group and immunity layer.
 *
 * @return An IOResult indicating success or failure.
 */
template <typename FP>
IOResult<void> export_input_data_timeseries(
    mio::VectorRange<Model<FP>> models, const std::string& results_dir, Date date, const std::vector<int>& node_ids,
    const std::vector<double>& scaling_factor_inf, const double scaling_factor_icu, const int num_days,
    const std::vector<std::vector<double>>& immunity_population, const mio::regions::de::EpidataFilenames& epidata_filenames)
{
    const auto num_age_groups = (size_t)models[0].parameters.get_num_groups();
    assert(scaling_factor_inf.size() == num_age_groups);
    assert(num_age_groups == ConfirmedCasesDataEntry::age_group_names.size());
    assert(models.size() == node_ids.size());
    std::vector<TimeSeries<double>> extrapolated_data(
        models.size(), TimeSeries<double>::zero(num_days + 1, (size_t)InfectionState::Count * num_age_groups));

    for (int t = 0; t <= num_days; ++t) {
        auto offset_day = offset_date_by_days(date, t);

        // TODO: empty vaccination data path guard
        BOOST_OUTCOME_TRY(read_input_data(model, date, county, scaling_factor_inf, scaling_factor_icu,
                                      num_days, immunity_population, epidata_filenames));

        for (size_t r = 0; r < node_ids.size(); r++) {
            extrapolated_data[r][t] = models[r].get_initial_values();
            // in set_population_data the number of death individuals is subtracted from the SusceptibleImprovedImmunity compartment.
            // Since we should be independent whether we consider them or not, we add them back here before we save the data.
            for (size_t age = 0; age < num_age_groups; age++) {
                extrapolated_data[r][t][(size_t)InfectionState::SusceptibleImprovedImmunity +
                                        age * (size_t)InfectionState::Count] +=
                    extrapolated_data[r][t][(size_t)InfectionState::DeadNaive + age * (size_t)InfectionState::Count];
            }
        }
    }
    BOOST_OUTCOME_TRY(save_result(extrapolated_data, node_ids, static_cast<int>(num_age_groups),
                                  path_join(results_dir, "Results_rki.h5")));

    auto extrapolated_rki_data_sum = sum_nodes(std::vector<std::vector<TimeSeries<double>>>{extrapolated_data});
    BOOST_OUTCOME_TRY(save_result({extrapolated_rki_data_sum[0][0]}, {0}, static_cast<int>(num_age_groups),
                                  path_join(results_dir, "Results_rki_sum.h5")));

    return success();
}

#else
template <typename FP>
IOResult<void> export_input_data_timeseries(mio::VectorRange<Model<FP>>, const std::string&, Date, const std::vector<int>&,
                                            const std::vector<double>&, const double, const int,
                                            const std::vector<std::vector<double>>, const mio::regions::de::EpidataFilenames&)
{
    mio::log_warning("HDF5 not available. Cannot export time series of extrapolated real data.");
    return success();
}

#endif //MEMILIO_HAS_HDF5

/**
 * @brief Reads compartments for geographic units at a specified date from data files.
 *
 * This function estimates all compartments from available data using the provided model parameters.
 *
 * @param[in,out] model Vector of models, one per county, to be initialized with data.
 * @param[in] date Date for which the data should be read.
 * @param[in] node_ids Vector of IDs of the units for which data is read.
 * @param[in] scaling_factor_inf Vector of scaling factors for confirmed cases.
 * @param[in] scaling_factor_icu Scaling factor for ICU cases.
 * @param[in] pydata_dir Directory containing the input data files.
 * @param[in] num_days Number of days to simulate.
 * @param[in] immunity_population Matrix containing immunity proportions for each age group and immunity layer.
 *
 * @return An IOResult indicating success or failure.
    */
template <typename FP>
IOResult<void> read_input_data(mio::VectorRange<Model<FP>>& model, Date date, const std::vector<int>& node_ids,
                               const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                               int num_days, const std::vector<std::vector<double>>& immunity_population, 
                               const mio::regions::de::EpidataFilenames& epidata_filenames)
{

    BOOST_OUTCOME_TRY(
        details::set_vaccination_data(model, epidata_filenames.vaccination_data_path, date, node_ids, num_days));

    // TODO: Reuse more code, e.g., set_divi_data (in secir) and a set_divi_data (here) only need a different ModelType.
    // TODO: add option to set ICU data from confirmed cases if DIVI or other data is not available.
    BOOST_OUTCOME_TRY(details::set_divi_data(model, epidata_filenames.divi_data_path, node_ids, date,
                                             scaling_factor_icu));

    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(model, epidata_filenames.case_data_path, node_ids,
                                                        date, scaling_factor_inf, immunity_population));
    BOOST_OUTCOME_TRY(details::set_population_data(model, epidata_filenames.population_data_path, node_ids,
                                                   immunity_population));
    return success();
}

} // namespace osecirts
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // MIO_ODE_SECIRTS_PARAMETERS_IO_H
