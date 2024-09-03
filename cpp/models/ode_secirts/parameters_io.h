/* 
* Copyright (C) 2020-2024 MEmilio
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
#ifndef ODESECIRTS_PARAMETERS_IO_H
#define ODESECIRTS_PARAMETERS_IO_H

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
        * @brief Reads subpopulations of infection states from transformed RKI cases file.
        * @param path Path to transformed RKI cases file.
        * @param vregion vector of keys of the region of interest     
        * @param date Date for which the arrays are initialized
        * @param num_* output vector for number of people in the corresponding compartement
        * @param t_* average time it takes to get from one compartement to another (vector with one element per age group)
        * @param mu_* probabilities to get from one compartement to another (vector with one element per age group)
        * @param scaling_factor_inf Factor for scaling the confirmed cases to account for estimated undetected cases.
        * @see mio::read_confirmed_cases_data
        * @{
        */
IOResult<void> read_confirmed_cases_data(
    std::string const& path, std::vector<int> const& vregion, Date date, std::vector<std::vector<double>>& num_Exposed,
    std::vector<std::vector<double>>& num_InfectedNoSymptoms, std::vector<std::vector<double>>& num_InfectedSymptoms,
    std::vector<std::vector<double>>& num_InfectedSevere, std::vector<std::vector<double>>& num_icu,
    std::vector<std::vector<double>>& num_death, std::vector<std::vector<double>>& vnum_timm_i,
    const std::vector<std::vector<int>>& t_Exposed, const std::vector<std::vector<int>>& t_InfectedNoSymptoms,
    const std::vector<std::vector<int>>& t_InfectedSymptoms, const std::vector<std::vector<int>>& t_InfectedSevere,
    const std::vector<std::vector<int>>& t_InfectedCritical, const std::vector<std::vector<int>>& vt_imm_interval_i,
    const std::vector<std::vector<double>>& mu_C_R, const std::vector<std::vector<double>>& mu_I_H,
    const std::vector<std::vector<double>>& mu_H_U, const std::vector<double>& scaling_factor_inf);

IOResult<void> read_confirmed_cases_data(
    const std::vector<ConfirmedCasesDataEntry>& rki_data, std::vector<int> const& vregion, Date date,
    std::vector<std::vector<double>>& num_Exposed, std::vector<std::vector<double>>& num_InfectedNoSymptoms,
    std::vector<std::vector<double>>& num_InfectedSymptoms, std::vector<std::vector<double>>& num_InfectedSevere,
    std::vector<std::vector<double>>& num_icu, std::vector<std::vector<double>>& num_death,
    std::vector<std::vector<double>>& vnum_timm_i, const std::vector<std::vector<int>>& t_Exposed,
    const std::vector<std::vector<int>>& t_InfectedNoSymptoms, const std::vector<std::vector<int>>& t_InfectedSymptoms,
    const std::vector<std::vector<int>>& t_InfectedSevere, const std::vector<std::vector<int>>& t_InfectedCritical,
    const std::vector<std::vector<int>>& vt_imm_interval_i, const std::vector<std::vector<double>>& mu_C_R,
    const std::vector<std::vector<double>>& mu_I_H, const std::vector<std::vector<double>>& mu_H_U,
    const std::vector<double>& scaling_factor_inf);
/**@}*/

IOResult<std::vector<std::vector<double>>> read_immunity_population(const std::string& path,
                                                                    const size_t& num_age_groups);
/**@}*/

/**
        * @brief sets infected and immun compartment(s) from a transformed RKI cases file into a Model.
        * @param model vector of objects in which the data is set
        * @param path Path to transformed RKI cases file
        * @param region vector of keys of the region of interest
        * @param date Date for which the arrays are initialized
        * @param scaling_factor_inf factors by which to scale the confirmed cases of
        * @param immunity_data_path Path to immunity data file
        */
template <class Model>
IOResult<void>
set_confirmed_cases_data(std::vector<Model>& model, const std::string& path, std::vector<int> const& region, Date date,
                         const std::vector<double>& scaling_factor_inf, const std::string& immunity_data_path)
{
    auto num_age_groups = (size_t)model[0].parameters.get_num_groups();
    assert(scaling_factor_inf.size() == num_age_groups); //TODO: allow vector or scalar valued scaling factors
    assert(ConfirmedCasesDataEntry::age_group_names.size() == num_age_groups);

    BOOST_OUTCOME_TRY(auto&& immunity_population, read_immunity_population(immunity_data_path, num_age_groups));

    BOOST_OUTCOME_TRY(auto&& rki_data, mio::read_confirmed_cases_data(path));

    std::vector<std::vector<int>> t_Exposed{model.size()};
    std::vector<std::vector<int>> t_InfectedNoSymptoms{model.size()};
    std::vector<std::vector<int>> t_InfectedSymptoms{model.size()};
    std::vector<std::vector<int>> t_InfectedSevere{model.size()};
    std::vector<std::vector<int>> t_InfectedCritical{model.size()};
    std::vector<std::vector<int>> t_imm_interval1{model.size()};
    std::vector<std::vector<int>> t_imm_interval2{model.size()};

    std::vector<std::vector<double>> mu_C_R{model.size()};
    std::vector<std::vector<double>> mu_I_H{model.size()};
    std::vector<std::vector<double>> mu_H_U{model.size()};

    std::vector<std::vector<double>> num_InfectedSymptoms(model.size());
    std::vector<std::vector<double>> num_death(model.size());
    std::vector<std::vector<double>> num_Exposed(model.size());
    std::vector<std::vector<double>> num_InfectedNoSymptoms(model.size());
    std::vector<std::vector<double>> num_InfectedSevere(model.size());
    std::vector<std::vector<double>> num_icu(model.size());
    std::vector<std::vector<double>> num_timm1(model.size());
    std::vector<std::vector<double>> num_timm2(model.size());

    std::vector<double> denom_E(num_age_groups, 0.0);
    std::vector<double> denom_I_NS(num_age_groups, 0.0);
    std::vector<double> denom_I_Sy(num_age_groups, 0.0);
    std::vector<double> denom_I_Sev(num_age_groups, 0.0);

    /*----------- UNVACCINATED -----------*/
    for (size_t county = 0; county < model.size(); county++) {
        num_InfectedSymptoms[county]   = std::vector<double>(num_age_groups, 0.0);
        num_death[county]              = std::vector<double>(num_age_groups, 0.0);
        num_Exposed[county]            = std::vector<double>(num_age_groups, 0.0);
        num_InfectedNoSymptoms[county] = std::vector<double>(num_age_groups, 0.0);
        num_InfectedSevere[county]     = std::vector<double>(num_age_groups, 0.0);
        num_icu[county]                = std::vector<double>(num_age_groups, 0.0);
        num_timm1[county]              = std::vector<double>(num_age_groups, 0.0);
        num_timm2[county]              = std::vector<double>(num_age_groups, 0.0);
        for (size_t group = 0; group < num_age_groups; group++) {

            t_Exposed[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeExposed<double>>()[(AgeGroup)group])));
            t_InfectedNoSymptoms[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeInfectedNoSymptoms<double>>()[(AgeGroup)group])));
            t_InfectedSymptoms[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeInfectedSymptoms<double>>()[(AgeGroup)group])));
            t_InfectedSevere[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeInfectedSevere<double>>()[(AgeGroup)group])));
            t_InfectedCritical[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeInfectedCritical<double>>()[(AgeGroup)group])));
            t_imm_interval1[county].push_back(static_cast<int>(
                model[county].parameters.template get<TimeTemporaryImmunityPI<double>>()[(AgeGroup)group]));
            t_imm_interval2[county].push_back(static_cast<int>(
                model[county].parameters.template get<TimeTemporaryImmunityII<double>>()[(AgeGroup)group]));

            mu_C_R[county].push_back(
                model[county].parameters.template get<RecoveredPerInfectedNoSymptoms<double>>()[(AgeGroup)group]);
            mu_I_H[county].push_back(
                model[county].parameters.template get<SeverePerInfectedSymptoms<double>>()[(AgeGroup)group]);
            mu_H_U[county].push_back(
                model[county].parameters.template get<CriticalPerSevere<double>>()[(AgeGroup)group]);

            // calculate the denominators to split the reported case numbers to the different immunity layers.
            denom_E[group] =
                1 /
                (immunity_population[0][group] +
                 immunity_population[1][group] *
                     model[county].parameters.template get<ReducExposedPartialImmunity<double>>()[(AgeGroup)group] +
                 immunity_population[2][group] *
                     model[county].parameters.template get<ReducExposedImprovedImmunity<double>>()[(AgeGroup)group]);

            denom_I_NS[group] =
                1 /
                (immunity_population[0][group] +
                 immunity_population[1][group] *
                     model[county].parameters.template get<ReducExposedPartialImmunity<double>>()[(AgeGroup)group] +
                 immunity_population[2][group] *
                     model[county].parameters.template get<ReducExposedImprovedImmunity<double>>()[(AgeGroup)group]);

            denom_I_Sy[group] =
                1 /
                (immunity_population[0][group] +
                 immunity_population[1][group] *
                     model[county]
                         .parameters.template get<ReducInfectedSymptomsPartialImmunity<double>>()[(AgeGroup)group] +
                 immunity_population[2][group] *
                     model[county]
                         .parameters.template get<ReducInfectedSymptomsImprovedImmunity<double>>()[(AgeGroup)group]);

            denom_I_Sev[group] =
                1 /
                (immunity_population[0][group] +
                 immunity_population[1][group] *
                     model[county].parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[(
                         AgeGroup)group] +
                 immunity_population[2][group] *
                     model[county].parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[(
                         AgeGroup)group]);
        }
    }

    BOOST_OUTCOME_TRY(read_confirmed_cases_data(
        rki_data, region, date, num_Exposed, num_InfectedNoSymptoms, num_InfectedSymptoms, num_InfectedSevere, num_icu,
        num_death, num_timm1, t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere, t_InfectedCritical,
        t_imm_interval1, mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf));

    for (size_t county = 0; county < model.size(); county++) {
        // if (std::accumulate(num_InfectedSymptoms[county].begin(), num_InfectedSymptoms[county].end(), 0.0) > 0) {
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
                immunity_population[0][i] * denom_I_Sev[i] * num_InfectedSevere[county][i];
        }

        // }
        if (std::accumulate(num_InfectedSymptoms[county].begin(), num_InfectedSymptoms[county].end(), 0.0) == 0) {
            log_warning("No infections for unvaccinated reported on date " + std::to_string(date.year) + "-" +
                        std::to_string(date.month) + "-" + std::to_string(date.day) + " for region " +
                        std::to_string(region[county]) + ". Population data has not been set.");
        }
    }

    /*----------- PARTIALLY VACCINATED -----------*/
    for (size_t county = 0; county < model.size(); county++) {
        t_InfectedNoSymptoms[county].clear();
        t_Exposed[county].clear();
        t_InfectedSymptoms[county].clear();
        t_InfectedSevere[county].clear();
        t_InfectedCritical[county].clear();

        mu_C_R[county].clear();
        mu_I_H[county].clear();
        mu_H_U[county].clear();

        num_InfectedSymptoms[county]   = std::vector<double>(num_age_groups, 0.0);
        num_death[county]              = std::vector<double>(num_age_groups, 0.0);
        num_Exposed[county]            = std::vector<double>(num_age_groups, 0.0);
        num_InfectedNoSymptoms[county] = std::vector<double>(num_age_groups, 0.0);
        num_InfectedSevere[county]     = std::vector<double>(num_age_groups, 0.0);
        num_icu[county]                = std::vector<double>(num_age_groups, 0.0);
        for (size_t group = 0; group < num_age_groups; group++) {

            double reduc_t = model[0].parameters.template get<ReducTimeInfectedMild<double>>()[(AgeGroup)group];
            t_Exposed[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeExposed<double>>()[(AgeGroup)group])));
            t_InfectedNoSymptoms[county].push_back(static_cast<int>(std::round(
                model[county].parameters.template get<TimeInfectedNoSymptoms<double>>()[(AgeGroup)group] * reduc_t)));
            t_InfectedSymptoms[county].push_back(static_cast<int>(std::round(
                model[county].parameters.template get<TimeInfectedSymptoms<double>>()[(AgeGroup)group] * reduc_t)));
            t_InfectedSevere[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeInfectedSevere<double>>()[(AgeGroup)group])));
            t_InfectedCritical[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeInfectedCritical<double>>()[(AgeGroup)group])));

            double exp_fac_part_immune =
                model[county].parameters.template get<ReducExposedPartialImmunity<double>>()[(AgeGroup)group];
            double inf_fac_part_immune =
                model[county].parameters.template get<ReducInfectedSymptomsPartialImmunity<double>>()[(AgeGroup)group];
            double hosp_fac_part_immune =
                model[county]
                    .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[(AgeGroup)group];
            double icu_fac_part_immune =
                model[county]
                    .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[(AgeGroup)group];
            mu_C_R[county].push_back((
                1 - inf_fac_part_immune / exp_fac_part_immune *
                        (1 - model[county]
                                 .parameters.template get<RecoveredPerInfectedNoSymptoms<double>>()[(AgeGroup)group])));
            mu_I_H[county].push_back(
                hosp_fac_part_immune / inf_fac_part_immune *
                model[county].parameters.template get<SeverePerInfectedSymptoms<double>>()[(AgeGroup)group]);
            // transfer from H to U, D unchanged.
            mu_H_U[county].push_back(
                icu_fac_part_immune / hosp_fac_part_immune *
                model[county].parameters.template get<CriticalPerSevere<double>>()[(AgeGroup)group]);
        }
    }

    BOOST_OUTCOME_TRY(read_confirmed_cases_data(
        rki_data, region, date, num_Exposed, num_InfectedNoSymptoms, num_InfectedSymptoms, num_InfectedSevere, num_icu,
        num_death, num_timm1, t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere, t_InfectedCritical,
        t_imm_interval1, mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf));

    for (size_t county = 0; county < model.size(); county++) {
        // if (std::accumulate(num_InfectedSymptoms[county].begin(), num_InfectedSymptoms[county].end(), 0.0) > 0) {
        size_t num_groups = (size_t)model[county].parameters.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            model[county].populations[{AgeGroup(i), InfectionState::ExposedPartialImmunity}] =
                immunity_population[1][i] *
                model[county].parameters.template get<ReducExposedPartialImmunity<double>>()[(AgeGroup)i] * denom_E[i] *
                num_Exposed[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsPartialImmunity}] =
                immunity_population[1][i] *
                model[county].parameters.template get<ReducExposedPartialImmunity<double>>()[(AgeGroup)i] *
                denom_I_NS[i] * num_InfectedNoSymptoms[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}] = 0;
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsPartialImmunity}] =
                immunity_population[1][i] *
                model[county].parameters.template get<ReducInfectedSymptomsPartialImmunity<double>>()[(AgeGroup)i] *
                denom_I_Sy[i] * num_InfectedSymptoms[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsPartialImmunityConfirmed}] = 0;
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSeverePartialImmunity}] =
                immunity_population[1][i] *
                model[county]
                    .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[(AgeGroup)i] *
                denom_I_Sev[i] * num_InfectedSevere[county][i];
            // the += is necessary because we already set the previous vaccinated individuals
            model[county].populations[{AgeGroup(i), InfectionState::TemporaryImmunPartialImmunity}] +=
                immunity_population[1][i] *
                model[county].parameters.template get<ReducExposedPartialImmunity<double>>()[(AgeGroup)i] * denom_E[i] *
                num_timm1[county][i];
        }
        // }
        if (std::accumulate(num_InfectedSymptoms[county].begin(), num_InfectedSymptoms[county].end(), 0.0) == 0) {
            log_warning("No infections for partially vaccinated reported on date " + std::to_string(date.year) + "-" +
                        std::to_string(date.month) + "-" + std::to_string(date.day) + " for region " +
                        std::to_string(region[county]) + ". Population data has not been set.");
        }
    }

    /*----------- FULLY VACCINATED -----------*/
    for (size_t county = 0; county < model.size(); county++) {
        t_InfectedNoSymptoms[county].clear();
        t_Exposed[county].clear();
        t_InfectedSymptoms[county].clear();
        t_InfectedSevere[county].clear();
        t_InfectedCritical[county].clear();

        mu_C_R[county].clear();
        mu_I_H[county].clear();
        mu_H_U[county].clear();

        num_InfectedSymptoms[county]   = std::vector<double>(num_age_groups, 0.0);
        num_death[county]              = std::vector<double>(num_age_groups, 0.0);
        num_Exposed[county]            = std::vector<double>(num_age_groups, 0.0);
        num_InfectedNoSymptoms[county] = std::vector<double>(num_age_groups, 0.0);
        num_InfectedSevere[county]     = std::vector<double>(num_age_groups, 0.0);
        num_icu[county]                = std::vector<double>(num_age_groups, 0.0);
        for (size_t group = 0; group < num_age_groups; group++) {

            double reduc_t = model[0].parameters.template get<ReducTimeInfectedMild<double>>()[(AgeGroup)group];
            t_Exposed[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeExposed<double>>()[(AgeGroup)group])));
            t_InfectedNoSymptoms[county].push_back(static_cast<int>(std::round(
                model[county].parameters.template get<TimeInfectedNoSymptoms<double>>()[(AgeGroup)group] * reduc_t)));
            t_InfectedSymptoms[county].push_back(static_cast<int>(std::round(
                model[county].parameters.template get<TimeInfectedSymptoms<double>>()[(AgeGroup)group] * reduc_t)));
            t_InfectedSevere[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeInfectedSevere<double>>()[(AgeGroup)group])));
            t_InfectedCritical[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeInfectedCritical<double>>()[(AgeGroup)group])));

            double reduc_immune_exp =
                model[county].parameters.template get<ReducExposedImprovedImmunity<double>>()[(AgeGroup)group];
            double reduc_immune_inf =
                model[county].parameters.template get<ReducInfectedSymptomsImprovedImmunity<double>>()[(AgeGroup)group];
            double reduc_immune_hosp =
                model[county].parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[(
                    AgeGroup)group];
            double reduc_immune_icu =
                model[county].parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[(
                    AgeGroup)group];
            mu_C_R[county].push_back((
                1 - reduc_immune_inf / reduc_immune_exp *
                        (1 - model[county]
                                 .parameters.template get<RecoveredPerInfectedNoSymptoms<double>>()[(AgeGroup)group])));
            mu_I_H[county].push_back(
                reduc_immune_hosp / reduc_immune_inf *
                model[county].parameters.template get<SeverePerInfectedSymptoms<double>>()[(AgeGroup)group]);
            // transfer from H to U, D unchanged.
            mu_H_U[county].push_back(
                reduc_immune_icu / reduc_immune_hosp *
                model[county].parameters.template get<CriticalPerSevere<double>>()[(AgeGroup)group]);
        }
    }

    BOOST_OUTCOME_TRY(read_confirmed_cases_data(
        rki_data, region, date, num_Exposed, num_InfectedNoSymptoms, num_InfectedSymptoms, num_InfectedSevere, num_icu,
        num_death, num_timm2, t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere, t_InfectedCritical,
        t_imm_interval2, mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf));

    for (size_t county = 0; county < model.size(); county++) {
        // if (std::accumulate(num_InfectedSymptoms[county].begin(), num_InfectedSymptoms[county].end(), 0.0) > 0) {
        size_t num_groups = (size_t)model[county].parameters.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            model[county].populations[{AgeGroup(i), InfectionState::ExposedImprovedImmunity}] =
                immunity_population[2][i] *
                model[county].parameters.template get<ReducExposedImprovedImmunity<double>>()[(AgeGroup)i] *
                denom_E[i] * num_Exposed[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsImprovedImmunity}] =
                immunity_population[2][i] *
                model[county].parameters.template get<ReducExposedImprovedImmunity<double>>()[(AgeGroup)i] *
                denom_I_NS[i] * num_InfectedNoSymptoms[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] = 0;
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsImprovedImmunity}] =
                immunity_population[2][i] *
                model[county].parameters.template get<ReducInfectedSymptomsImprovedImmunity<double>>()[(AgeGroup)i] *
                denom_I_Sy[i] * num_InfectedSymptoms[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] = 0;
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSevereImprovedImmunity}] =
                immunity_population[2][i] *
                model[county]
                    .parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[(AgeGroup)i] *
                denom_I_Sev[i] * num_InfectedSevere[county][i];
            // the += is necessary because we already set the previous vaccinated individuals
            model[county].populations[{AgeGroup(i), InfectionState::TemporaryImmunImprovedImmunity}] +=
                immunity_population[2][i] *
                model[county].parameters.template get<ReducExposedImprovedImmunity<double>>()[(AgeGroup)i] *
                denom_E[i] * num_timm2[county][i];
        }
        // }
        if (std::accumulate(num_InfectedSymptoms[county].begin(), num_InfectedSymptoms[county].end(), 0.0) == 0) {
            log_warning("No infections for vaccinated reported on date " + std::to_string(date.year) + "-" +
                        std::to_string(date.month) + "-" + std::to_string(date.day) + " for region " +
                        std::to_string(region[county]) + ". Population data has not been set.");
        }
    }
    return success();
}

/**
        * @brief reads number of ICU patients from DIVI register into Parameters
        * @param path Path to transformed DIVI file
        * @param vregion Keys of the region of interest
        * @param date Date for which the arrays are initialized
        * @param vnum_icu number of ICU patients
        * @see mio::read_divi_data
        * @{
        */
IOResult<void> read_divi_data(const std::string& path, const std::vector<int>& vregion, Date date,
                              std::vector<double>& vnum_icu);
IOResult<void> read_divi_data(const std::vector<DiviEntry>& divi_data, const std::vector<int>& vregion, Date date,
                              std::vector<double>& vnum_icu);
/**@}*/

/**
        * @brief sets populations data from DIVI register into Model
        * @param model vector of objects in which the data is set
        * @param path Path to transformed DIVI file
        * @param vregion vector of keys of the regions of interest
        * @param date Date for which the arrays are initialized
        * @param scaling_factor_icu factor by which to scale the icu cases of divi data
        */
template <class Model>
IOResult<void> set_divi_data(std::vector<Model>& model, const std::string& path, const std::vector<int>& vregion,
                             Date date, double scaling_factor_icu)
{
    std::vector<double> sum_mu_I_U(vregion.size(), 0);
    std::vector<std::vector<double>> mu_I_U{model.size()};
    for (size_t region = 0; region < vregion.size(); region++) {
        auto num_groups = model[region].parameters.get_num_groups();
        for (auto i = AgeGroup(0); i < num_groups; i++) {
            sum_mu_I_U[region] += model[region].parameters.template get<CriticalPerSevere<double>>()[i] *
                                  model[region].parameters.template get<SeverePerInfectedSymptoms<double>>()[i];
            mu_I_U[region].push_back(model[region].parameters.template get<CriticalPerSevere<double>>()[i] *
                                     model[region].parameters.template get<SeverePerInfectedSymptoms<double>>()[i]);
        }
    }
    std::vector<double> num_icu(model.size(), 0.0);
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
        * @brief reads population data from census data.
        * @param path Path to population data file.
        * @param vregion vector of keys of the regions of interest
        * @see mio::read_population_data
        * @{
        */
IOResult<std::vector<std::vector<double>>> read_population_data(const std::string& path,
                                                                const std::vector<int>& vregion);
IOResult<std::vector<std::vector<double>>> read_population_data(const std::vector<PopulationDataEntry>& population_data,
                                                                const std::vector<int>& vregion);
/**@}*/

template <class Model>
IOResult<void> set_population_data(std::vector<Model>& model, const std::string& path, const std::vector<int>& vregion,
                                   const std::string& immunity_data_path)
{
    BOOST_OUTCOME_TRY(auto&& num_population, read_population_data(path, vregion));

    auto num_age_groups = ConfirmedCasesDataEntry::age_group_names.size();
    BOOST_OUTCOME_TRY(auto&& immunity_population, read_immunity_population(immunity_data_path, num_age_groups));

    for (size_t region = 0; region < vregion.size(); region++) {
        if (std::accumulate(num_population[region].begin(), num_population[region].end(), 0.0) > 0) {
            auto num_groups = model[region].parameters.get_num_groups();
            for (auto i = AgeGroup(0); i < num_groups; i++) {

                double SN  = num_population[region][size_t(i)] * immunity_population[0][size_t(i)];
                double SPI = num_population[region][size_t(i)] * immunity_population[1][size_t(i)];
                double SII = num_population[region][size_t(i)] - SN - SPI;

                model[region].populations[{i, InfectionState::SusceptibleImprovedImmunity}] = std::max(
                    0.0,
                    double(
                        SII -
                        (model[region].populations[{i, InfectionState::ExposedImprovedImmunity}] +
                         model[region].populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunity}] +
                         model[region].populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] +
                         model[region].populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] +
                         model[region].populations[{i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] +
                         model[region].populations[{i, InfectionState::InfectedSevereImprovedImmunity}] +
                         model[region].populations[{i, InfectionState::InfectedCriticalImprovedImmunity}] +
                         model[region].populations[{i, InfectionState::DeadImprovedImmunity}] +
                         model[region].populations[{i, InfectionState::TemporaryImmunImprovedImmunity}])));

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
                        model[region].populations[{i, InfectionState::TemporaryImmunPartialImmunity}]);

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

template <typename FP = double>
IOResult<void> set_vaccination_data(std::vector<Model<FP>>& model, const std::string& path, Date date,
                                    const std::vector<int>& vregion, int num_days)
{
    BOOST_OUTCOME_TRY(auto&& vacc_data, read_vaccination_data(path));

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
} // namespace details

#ifdef MEMILIO_HAS_HDF5

/**
    * @brief Exports the time series of extrapolated real data according to
    *   the extrapolation / approximation method used to initialize the model from
    *   real world data.
        (This is the vector-valued functionality of set_confirmed_cases_data())
    * @param model vector of objects in which the data is set
    * @param data_dir Path to transformed RKI cases files
    * @param results_dir Path to result files
    * @param start_date Start date of the time series to be exported.
    * @param region vector of keys of the region of interest
    * @param scaling_factor_inf Factor for scaling the confirmed cases to account for an estimated number of undetected cases.
    * @param scaling_factor_icu Factor for scaling the reported ICU cases to account for possibly unreported ICU cases.
    * @param num_days Number of days for which the time series is exported.
    * @param divi_data_path path to divi data file
    * @param confirmed_cases_path path to confirmed cases file
    * @param population_data_path path to population data file
    * @param immunity_data_path path to immunity data file
    */
template <class Model>
IOResult<void> export_input_data_county_timeseries(
    const std::vector<Model>& model, const std::string& dir, std::vector<int> const& region, Date start_date,
    const std::vector<double>& scaling_factor_inf, double scaling_factor_icu, int num_days,
    const std::string& divi_data_path, const std::string& confirmed_cases_path, const std::string& population_data_path,
    const std::string& immunity_data_path)
{
    auto num_age_groups = (size_t)model[0].parameters.get_num_groups();
    BOOST_OUTCOME_TRY(auto&& immunity_population,
                      details::read_immunity_population(immunity_data_path, num_age_groups));

    assert(scaling_factor_inf.size() == num_age_groups);
    assert(num_age_groups == ConfirmedCasesDataEntry::age_group_names.size());
    assert(model.size() == region.size());

    BOOST_OUTCOME_TRY(auto&& rki_data, read_confirmed_cases_data(confirmed_cases_path));
    BOOST_OUTCOME_TRY(auto&& population_data, read_population_data(population_data_path));
    BOOST_OUTCOME_TRY(auto&& divi_data, read_divi_data(divi_data_path));

    /* functionality copy from set_confirmed_cases_data() here splitted in params */
    /* which do not need to be reset for each day and compartments sizes that are */
    /* set later for each day */
    /*----------- UNVACCINATED -----------*/
    // data needs to be int, because access to data-specific confirmed cases
    // is done with these parameters. TODO: Rounding instead
    // of casting to int should be introduced in the future.
    std::vector<std::vector<int>> t_InfectedNoSymptoms_n{model.size()};
    std::vector<std::vector<int>> t_Exposed_n{model.size()};
    std::vector<std::vector<int>> t_InfectedSymptoms_n{model.size()};
    std::vector<std::vector<int>> t_InfectedSevere_n{model.size()};
    std::vector<std::vector<int>> t_InfectedCritical_n{model.size()};
    std::vector<std::vector<int>> t_imm_interval1{model.size()};
    std::vector<std::vector<int>> t_imm_interval2{model.size()};

    std::vector<std::vector<double>> mu_C_R_n{model.size()};
    std::vector<std::vector<double>> mu_I_H_n{model.size()};
    std::vector<std::vector<double>> mu_H_U_n{model.size()};
    // ICU data is not age-resolved. Use a partition of unity defined by
    // the age-dependent probability I->H->U divided by the sum over all
    // age groups of all of these probabilities.
    std::vector<double> sum_mu_I_U_n(model.size(), 0);
    std::vector<std::vector<double>> mu_I_U_n{model.size()};

    for (size_t county = 0; county < model.size(); county++) {
        for (size_t group = 0; group < num_age_groups; group++) {

            t_Exposed_n[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeExposed<double>>()[(AgeGroup)group])));
            t_InfectedNoSymptoms_n[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeInfectedNoSymptoms<double>>()[(AgeGroup)group])));
            t_InfectedSymptoms_n[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeInfectedSymptoms<double>>()[(AgeGroup)group])));
            t_InfectedSevere_n[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeInfectedSevere<double>>()[(AgeGroup)group])));
            t_InfectedCritical_n[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeInfectedCritical<double>>()[(AgeGroup)group])));
            t_imm_interval1[county].push_back(static_cast<int>(
                model[county].parameters.template get<TimeTemporaryImmunityPI<double>>()[(AgeGroup)group]));
            t_imm_interval2[county].push_back(static_cast<int>(
                model[county].parameters.template get<TimeTemporaryImmunityII<double>>()[(AgeGroup)group]));

            mu_C_R_n[county].push_back(
                model[county].parameters.template get<RecoveredPerInfectedNoSymptoms<double>>()[(AgeGroup)group]);
            mu_I_H_n[county].push_back(
                model[county].parameters.template get<SeverePerInfectedSymptoms<double>>()[(AgeGroup)group]);
            mu_H_U_n[county].push_back(
                model[county].parameters.template get<CriticalPerSevere<double>>()[(AgeGroup)group]);

            /* begin: NOT in set_confirmed_cases_data() */
            sum_mu_I_U_n[county] +=
                model[county].parameters.template get<CriticalPerSevere<double>>()[AgeGroup(group)] *
                model[county].parameters.template get<SeverePerInfectedSymptoms<double>>()[AgeGroup(group)];
            mu_I_U_n[county].push_back(
                model[county].parameters.template get<CriticalPerSevere<double>>()[AgeGroup(group)] *
                model[county].parameters.template get<SeverePerInfectedSymptoms<double>>()[AgeGroup(group)]);
            /* end: NOT in set_confirmed_cases_data() */
        }
    }

    /*----------- PARTIALLY VACCINATED -----------*/
    // data needs to be int, because access to data-specific confirmed cases
    // is done with these parameters. TODO: Rounding instead
    // of casting to int should be introduced in the future.
    std::vector<std::vector<int>> t_InfectedNoSymptoms_pi{model.size()};
    std::vector<std::vector<int>> t_Exposed_pi{model.size()};
    std::vector<std::vector<int>> t_InfectedSymptoms_pi{model.size()};
    std::vector<std::vector<int>> t_InfectedSevere_pi{model.size()};
    std::vector<std::vector<int>> t_InfectedCritical_pi{model.size()};

    std::vector<std::vector<double>> mu_C_R_pi{model.size()};
    std::vector<std::vector<double>> mu_I_H_pi{model.size()};
    std::vector<std::vector<double>> mu_H_U_pi{model.size()};

    // ICU data is not age-resolved. Use a partition of unity defined by
    // the age-dependent probability I->H->U divided by the sum over all
    // age groups of all of these probabilities.
    std::vector<double> sum_mu_I_U_pi(model.size(), 0);
    std::vector<std::vector<double>> mu_I_U_pi{model.size()};
    for (size_t county = 0; county < model.size(); county++) {
        for (size_t group = 0; group < num_age_groups; group++) {

            double reduc_t = model[0].parameters.template get<ReducTimeInfectedMild<double>>()[(AgeGroup)group];
            t_Exposed_pi[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeExposed<double>>()[(AgeGroup)group])));
            t_InfectedNoSymptoms_pi[county].push_back(static_cast<int>(std::round(
                model[county].parameters.template get<TimeInfectedNoSymptoms<double>>()[(AgeGroup)group] * reduc_t)));
            t_InfectedSymptoms_pi[county].push_back(static_cast<int>(std::round(
                model[county].parameters.template get<TimeInfectedSymptoms<double>>()[(AgeGroup)group] * reduc_t)));
            t_InfectedSevere_pi[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeInfectedSevere<double>>()[(AgeGroup)group])));
            t_InfectedCritical_pi[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeInfectedCritical<double>>()[(AgeGroup)group])));

            double exp_fact_part_immune =
                model[county].parameters.template get<ReducExposedPartialImmunity<double>>()[(AgeGroup)group];
            double inf_fact_part_immune =
                model[county].parameters.template get<ReducInfectedSymptomsPartialImmunity<double>>()[(AgeGroup)group];
            double hosp_fact_part_immune =
                model[county]
                    .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[(AgeGroup)group];
            double icu_fact_part_immune =
                model[county]
                    .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[(AgeGroup)group];
            mu_C_R_pi[county].push_back((
                1 - inf_fact_part_immune / exp_fact_part_immune *
                        (1 - model[county]
                                 .parameters.template get<RecoveredPerInfectedNoSymptoms<double>>()[(AgeGroup)group])));
            mu_I_H_pi[county].push_back(
                hosp_fact_part_immune / inf_fact_part_immune *
                model[county].parameters.template get<SeverePerInfectedSymptoms<double>>()[(AgeGroup)group]);
            // transfer from H to U, D unchanged.
            mu_H_U_pi[county].push_back(
                icu_fact_part_immune / hosp_fact_part_immune *
                model[county].parameters.template get<CriticalPerSevere<double>>()[(AgeGroup)group]);

            sum_mu_I_U_pi[county] +=
                icu_fact_part_immune / hosp_fact_part_immune *
                model[county].parameters.template get<CriticalPerSevere<double>>()[AgeGroup(group)] *
                hosp_fact_part_immune / inf_fact_part_immune *
                model[county].parameters.template get<SeverePerInfectedSymptoms<double>>()[AgeGroup(group)];
            mu_I_U_pi[county].push_back(
                icu_fact_part_immune / hosp_fact_part_immune *
                model[county].parameters.template get<CriticalPerSevere<double>>()[AgeGroup(group)] *
                hosp_fact_part_immune / inf_fact_part_immune *
                model[county].parameters.template get<SeverePerInfectedSymptoms<double>>()[AgeGroup(group)]);
        }
    }

    /*----------- FULLY VACCINATED -----------*/
    // data needs to be int, because access to data-specific confirmed cases
    // is done with these parameters. TODO: Rounding instead
    // of casting to int should be introduced in the future.
    std::vector<std::vector<int>> t_InfectedNoSymptoms_ii{model.size()};
    std::vector<std::vector<int>> t_Exposed_ii{model.size()};
    std::vector<std::vector<int>> t_InfectedSymptoms_ii{model.size()};
    std::vector<std::vector<int>> t_InfectedSevere_ii{model.size()};
    std::vector<std::vector<int>> t_InfectedCritical_ii{model.size()};

    std::vector<std::vector<double>> mu_C_R_ii{model.size()};
    std::vector<std::vector<double>> mu_I_H_ii{model.size()};
    std::vector<std::vector<double>> mu_H_U_ii{model.size()};
    // ICU data is not age-resolved. Use a partition of unity defined by
    // the age-dependent probability I->H->U divided by the sum over all
    // age groups of all of these probabilities.
    std::vector<double> sum_mu_I_U_ii(model.size(), 0);
    std::vector<std::vector<double>> mu_I_U_ii{model.size()};
    for (size_t county = 0; county < model.size(); county++) {
        for (size_t group = 0; group < num_age_groups; group++) {

            double reduc_t = model[0].parameters.template get<ReducTimeInfectedMild<double>>()[(AgeGroup)group];
            t_Exposed_ii[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeExposed<double>>()[(AgeGroup)group])));
            t_InfectedNoSymptoms_ii[county].push_back(static_cast<int>(std::round(
                model[county].parameters.template get<TimeInfectedNoSymptoms<double>>()[(AgeGroup)group] * reduc_t)));
            t_InfectedSymptoms_ii[county].push_back(static_cast<int>(std::round(
                model[county].parameters.template get<TimeInfectedSymptoms<double>>()[(AgeGroup)group] * reduc_t)));
            t_InfectedSevere_ii[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeInfectedSevere<double>>()[(AgeGroup)group])));
            t_InfectedCritical_ii[county].push_back(static_cast<int>(
                std::round(model[county].parameters.template get<TimeInfectedCritical<double>>()[(AgeGroup)group])));

            double reduc_immune_exp =
                model[county].parameters.template get<ReducExposedImprovedImmunity<double>>()[(AgeGroup)group];
            double reduc_immune_inf =
                model[county].parameters.template get<ReducInfectedSymptomsImprovedImmunity<double>>()[(AgeGroup)group];
            double reduc_immune_hosp =
                model[county].parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[(
                    AgeGroup)group];
            double reduc_immune_icu =
                model[county].parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[(
                    AgeGroup)group];
            mu_C_R_ii[county].push_back((
                1 - reduc_immune_inf / reduc_immune_exp *
                        (1 - model[county]
                                 .parameters.template get<RecoveredPerInfectedNoSymptoms<double>>()[(AgeGroup)group])));
            mu_I_H_ii[county].push_back(
                reduc_immune_hosp / reduc_immune_inf *
                model[county].parameters.template get<SeverePerInfectedSymptoms<double>>()[(AgeGroup)group]);
            // transfer from H to U, D unchanged.
            mu_H_U_ii[county].push_back(
                reduc_immune_icu / reduc_immune_hosp *
                model[county].parameters.template get<CriticalPerSevere<double>>()[(AgeGroup)group]);

            sum_mu_I_U_ii[county] +=
                reduc_immune_icu / reduc_immune_hosp *
                model[county].parameters.template get<CriticalPerSevere<double>>()[AgeGroup(group)] *
                reduc_immune_hosp / reduc_immune_inf *
                model[county].parameters.template get<SeverePerInfectedSymptoms<double>>()[AgeGroup(group)];
            mu_I_U_ii[county].push_back(
                reduc_immune_icu / reduc_immune_hosp *
                model[county].parameters.template get<CriticalPerSevere<double>>()[AgeGroup(group)] *
                reduc_immune_hosp / reduc_immune_inf *
                model[county].parameters.template get<SeverePerInfectedSymptoms<double>>()[AgeGroup(group)]);
        }
    }
    std::vector<TimeSeries<double>> extrapolated_rki(
        model.size(), TimeSeries<double>::zero(num_days + 1, (size_t)InfectionState::Count * num_age_groups));

    for (size_t day = 0; day <= static_cast<size_t>(num_days); day++) {
        auto date = offset_date_by_days(start_date, int(day));

        // unvaccinated
        std::vector<std::vector<double>> num_Exposed_n(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_InfectedNoSymptoms_n(model.size(),
                                                                  std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_InfectedSymptoms_n(model.size(), std::vector<double>(num_age_groups, 0.0));
        // potential TODO: these confirmed are only confirmed by commuting, set to zero here. Adapt if generalized!
        std::vector<std::vector<double>> num_InfectedNoSymptomsConfirmed_n(model.size(),
                                                                           std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_InfectedSymptomsConfirmed_n(model.size(),
                                                                         std::vector<double>(num_age_groups, 0.0));
        // end TODO
        std::vector<std::vector<double>> num_InfectedSevere_n(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_death_n(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> dummy_icu(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> dummy_timm(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_timm1(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_timm2(model.size(), std::vector<double>(num_age_groups, 0.0));
        BOOST_OUTCOME_TRY(details::read_confirmed_cases_data(
            rki_data, region, date, num_Exposed_n, num_InfectedNoSymptoms_n, num_InfectedSymptoms_n,
            num_InfectedSevere_n, dummy_icu, num_death_n, dummy_timm, t_Exposed_n, t_InfectedNoSymptoms_n,
            t_InfectedSymptoms_n, t_InfectedSevere_n, t_InfectedCritical_n, t_imm_interval1, mu_C_R_n, mu_I_H_n,
            mu_H_U_n, scaling_factor_inf));

        // partially vaccinated
        std::vector<std::vector<double>> num_Exposed_pi(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_InfectedNoSymptoms_pi(model.size(),
                                                                   std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_InfectedSymptoms_pi(model.size(),
                                                                 std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_InfectedSevere_pi(model.size(), std::vector<double>(num_age_groups, 0.0));
        // potential TODO: these confirmed are only confirmed by commuting, set to zero here. Adapt if generalized!
        std::vector<std::vector<double>> num_InfectedNoSymptomsConfirmed_pi(model.size(),
                                                                            std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_InfectedSymptomsConfirmed_pi(model.size(),
                                                                          std::vector<double>(num_age_groups, 0.0));
        // end TODO
        std::vector<std::vector<double>> dummy_death(model.size(), std::vector<double>(num_age_groups, 0.0));
        for (size_t county = 0; county < model.size(); county++) {
            dummy_death[county] = std::vector<double>(num_age_groups, 0.0);
            dummy_icu[county]   = std::vector<double>(num_age_groups, 0.0);
        }
        BOOST_OUTCOME_TRY(details::read_confirmed_cases_data(
            rki_data, region, date, num_Exposed_pi, num_InfectedNoSymptoms_pi, num_InfectedSymptoms_pi,
            num_InfectedSevere_pi, dummy_icu, dummy_death, num_timm1, t_Exposed_pi, t_InfectedNoSymptoms_pi,
            t_InfectedSymptoms_pi, t_InfectedSevere_pi, t_InfectedCritical_pi, t_imm_interval1, mu_C_R_pi, mu_I_H_pi,
            mu_H_U_pi, scaling_factor_inf));

        // fully vaccinated
        std::vector<std::vector<double>> num_Exposed_ii(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_InfectedNoSymptoms_ii(model.size(),
                                                                   std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_InfectedSymptoms_ii(model.size(),
                                                                 std::vector<double>(num_age_groups, 0.0));
        // potential TODO: these confirmed are only confirmed by commuting, set to zero here. Adapt if generalized!
        std::vector<std::vector<double>> num_InfectedNoSymptomsConfirmed_ii(model.size(),
                                                                            std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_InfectedSymptomsConfirmed_ii(model.size(),
                                                                          std::vector<double>(num_age_groups, 0.0));
        // end TODO
        std::vector<std::vector<double>> num_InfectedSevere_ii(model.size(), std::vector<double>(num_age_groups, 0.0));
        for (size_t county = 0; county < model.size(); county++) {
            dummy_death[county] = std::vector<double>(num_age_groups, 0.0);
            dummy_icu[county]   = std::vector<double>(num_age_groups, 0.0);
        }
        BOOST_OUTCOME_TRY(details::read_confirmed_cases_data(
            rki_data, region, date, num_Exposed_ii, num_InfectedNoSymptoms_ii, num_InfectedSymptoms_ii,
            num_InfectedSevere_ii, dummy_icu, dummy_death, num_timm2, t_Exposed_ii, t_InfectedNoSymptoms_ii,
            t_InfectedSymptoms_ii, t_InfectedSevere_ii, t_InfectedCritical_ii, t_imm_interval2, mu_C_R_ii, mu_I_H_ii,
            mu_H_U_ii, scaling_factor_inf));

        // ICU only read for compartment InfectionState::InfectedCritical and then distributed later
        std::vector<double> dummy_icu2(model.size(), 0.0);
        BOOST_OUTCOME_TRY(details::read_divi_data(divi_data, region, date, dummy_icu2));

        std::vector<std::vector<double>> num_icu(model.size(), std::vector<double>(num_age_groups, 0.0));
        for (size_t county = 0; county < region.size(); county++) {
            for (size_t age = 0; age < num_age_groups; age++) {
                num_icu[county][age] =
                    scaling_factor_icu * dummy_icu2[county] * mu_I_U_n[county][age] / sum_mu_I_U_n[county];
            }
        }

        // read population basics
        BOOST_OUTCOME_TRY(auto&& num_population, details::read_population_data(population_data, region));

        for (size_t county = 0; county < region.size(); county++) {
            if (std::accumulate(num_population[county].begin(), num_population[county].end(), 0.0) > 0) {
                for (size_t age = 0; age < num_age_groups; age++) {
                    auto age_group_offset = age * (size_t)InfectionState::Count;
                    double SN             = num_population[county][size_t(age)] * immunity_population[0][size_t(age)];
                    double SPI            = num_population[county][size_t(age)] * immunity_population[1][size_t(age)];
                    double SII            = num_population[county][size_t(age)] - SN - SPI;

                    double denom_E =
                        1 / (SN +
                             SPI * model[county]
                                       .parameters.template get<ReducExposedPartialImmunity<double>>()[AgeGroup(age)] +
                             SII * model[county]
                                       .parameters.template get<ReducExposedImprovedImmunity<double>>()[AgeGroup(age)]);
                    double denom_C =
                        1 / (SN +
                             SPI * model[county]
                                       .parameters.template get<ReducExposedPartialImmunity<double>>()[AgeGroup(age)] +
                             SII * model[county]
                                       .parameters.template get<ReducExposedImprovedImmunity<double>>()[AgeGroup(age)]);
                    double denom_I =
                        1 / (SN +
                             SPI * model[county]
                                       .parameters
                                       .template get<ReducInfectedSymptomsPartialImmunity<double>>()[AgeGroup(age)] +
                             SII * model[county]
                                       .parameters
                                       .template get<ReducInfectedSymptomsImprovedImmunity<double>>()[AgeGroup(age)]);
                    double denom_HU =
                        1 / (SN +
                             SPI * model[county]
                                       .parameters.template get<
                                           ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[AgeGroup(age)] +
                             SII * model[county]
                                       .parameters.template get<
                                           ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[AgeGroup(age)]);

                    extrapolated_rki[county][day]((size_t)InfectionState::ExposedNaive + age_group_offset) =
                        SN * denom_E * num_Exposed_n[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::ExposedPartialImmunity + age_group_offset) =
                        SPI *
                        model[county].parameters.template get<ReducExposedPartialImmunity<double>>()[AgeGroup(age)] *
                        denom_E * num_Exposed_pi[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::ExposedImprovedImmunity + age_group_offset) =
                        SII *
                        model[county].parameters.template get<ReducExposedImprovedImmunity<double>>()[AgeGroup(age)] *
                        denom_E * num_Exposed_ii[county][age];

                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedNoSymptomsNaive + age_group_offset) =
                        SN * denom_C * num_InfectedNoSymptoms_n[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedNoSymptomsPartialImmunity +
                                                  age_group_offset) =
                        SPI * denom_C * num_InfectedNoSymptoms_pi[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedNoSymptomsImprovedImmunity +
                                                  age_group_offset) =
                        SII * denom_C * num_InfectedNoSymptoms_ii[county][age];

                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedNoSymptomsNaiveConfirmed +
                                                  age_group_offset) =
                        SN * denom_C * num_InfectedNoSymptomsConfirmed_n[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedNoSymptomsPartialImmunityConfirmed +
                                                  age_group_offset) =
                        SPI * denom_C * num_InfectedNoSymptomsConfirmed_pi[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed +
                                                  age_group_offset) =
                        SII * denom_C * num_InfectedNoSymptomsConfirmed_ii[county][age];

                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsNaive + age_group_offset) =
                        SN * denom_I * num_InfectedSymptoms_n[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsPartialImmunity +
                                                  age_group_offset) =
                        SPI *
                        model[county]
                            .parameters.template get<ReducInfectedSymptomsPartialImmunity<double>>()[AgeGroup(age)] *
                        denom_I * num_InfectedSymptoms_pi[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsImprovedImmunity +
                                                  age_group_offset) =
                        SII *
                        model[county]
                            .parameters.template get<ReducInfectedSymptomsImprovedImmunity<double>>()[AgeGroup(age)] *
                        denom_I * num_InfectedSymptoms_ii[county][age];

                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsNaiveConfirmed +
                                                  age_group_offset) =
                        SN * denom_I * num_InfectedSymptomsConfirmed_n[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsPartialImmunityConfirmed +
                                                  age_group_offset) =
                        SPI *
                        model[county]
                            .parameters.template get<ReducInfectedSymptomsPartialImmunity<double>>()[AgeGroup(age)] *
                        denom_I * num_InfectedSymptomsConfirmed_pi[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsImprovedImmunityConfirmed +
                                                  age_group_offset) =
                        SII *
                        model[county]
                            .parameters.template get<ReducInfectedSymptomsImprovedImmunity<double>>()[AgeGroup(age)] *
                        denom_I * num_InfectedSymptomsConfirmed_ii[county][age];

                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedSevereNaive + age_group_offset) =
                        SN * denom_HU * num_InfectedSevere_n[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedSeverePartialImmunity +
                                                  age_group_offset) =
                        SPI *
                        model[county]
                            .parameters
                            .template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[AgeGroup(age)] *
                        denom_HU * num_InfectedSevere_pi[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedSevereImprovedImmunity +
                                                  age_group_offset) =
                        SII *
                        model[county]
                            .parameters
                            .template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[AgeGroup(age)] *
                        denom_HU * num_InfectedSevere_ii[county][age];

                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedCriticalNaive + age_group_offset) =
                        SN * denom_HU * num_icu[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedCriticalPartialImmunity +
                                                  age_group_offset) =
                        SPI *
                        model[county]
                            .parameters
                            .template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[AgeGroup(age)] *
                        denom_HU * num_icu[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedCriticalImprovedImmunity +
                                                  age_group_offset) =
                        SII *
                        model[county]
                            .parameters
                            .template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[AgeGroup(age)] *
                        denom_HU * num_icu[county][age];

                    extrapolated_rki[county][day]((size_t)InfectionState::TemporaryImmunPartialImmunity +
                                                  age_group_offset) =
                        num_timm1[county][age] * immunity_population[0][size_t(age)];
                    extrapolated_rki[county][day]((size_t)InfectionState::TemporaryImmunImprovedImmunity +
                                                  age_group_offset) =
                        num_timm2[county][age] * immunity_population[1][size_t(age)] +
                        immunity_population[2][size_t(age)] * immunity_population[2][size_t(age)];

                    extrapolated_rki[county][day]((size_t)InfectionState::SusceptibleImprovedImmunity +
                                                  age_group_offset) =
                        SII -
                        (extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsImprovedImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day](
                             (size_t)InfectionState::InfectedSymptomsImprovedImmunityConfirmed + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedSevereImprovedImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedCriticalImprovedImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::DeadImprovedImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::TemporaryImmunImprovedImmunity +
                                                       age_group_offset));

                    extrapolated_rki[county][day]((size_t)InfectionState::SusceptibleImprovedImmunity +
                                                  age_group_offset) =
                        std::min(
                            SII - (extrapolated_rki[county][day](
                                       (size_t)InfectionState::InfectedSymptomsImprovedImmunity + age_group_offset) +
                                   extrapolated_rki[county][day](
                                       (size_t)InfectionState::InfectedSymptomsImprovedImmunityConfirmed +
                                       age_group_offset) +
                                   extrapolated_rki[county][day](
                                       (size_t)InfectionState::InfectedSevereImprovedImmunity + age_group_offset) +
                                   extrapolated_rki[county][day](
                                       (size_t)InfectionState::InfectedCriticalImprovedImmunity + age_group_offset) +
                                   extrapolated_rki[county][day]((size_t)InfectionState::DeadImprovedImmunity +
                                                                 age_group_offset) +
                                   extrapolated_rki[county][day](
                                       (size_t)InfectionState::TemporaryImmunImprovedImmunity + age_group_offset)),
                            std::max(0.0,
                                     double(extrapolated_rki[county][day](
                                         (size_t)InfectionState::SusceptibleImprovedImmunity + age_group_offset))));

                    /////////////////
                    extrapolated_rki[county][day]((size_t)InfectionState::SusceptiblePartialImmunity +
                                                  age_group_offset) =
                        SPI -
                        (extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsPartialImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day](
                             (size_t)InfectionState::InfectedSymptomsPartialImmunityConfirmed + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedSeverePartialImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedCriticalPartialImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::DeadPartialImmunity + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::TemporaryImmunImprovedImmunity +
                                                       age_group_offset));

                    extrapolated_rki[county][day]((size_t)InfectionState::SusceptiblePartialImmunity +
                                                  age_group_offset) =
                        std::min(SPI, std::max(0.0, double(extrapolated_rki[county][day](
                                                        (size_t)InfectionState::SusceptiblePartialImmunity +
                                                        age_group_offset))));

                    /////////////////
                    extrapolated_rki[county][day]((size_t)InfectionState::SusceptibleNaive + age_group_offset) =
                        SPI -
                        (extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsNaive +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsNaiveConfirmed +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedSevereNaive + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedCriticalNaive +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::DeadNaive + age_group_offset));

                    extrapolated_rki[county][day]((size_t)InfectionState::SusceptibleNaive + age_group_offset) =
                        std::min(SN, std::max(0.0, double(extrapolated_rki[county][day](
                                                       (size_t)InfectionState::SusceptibleNaive + age_group_offset))));

                    // in set_confirmed_cases_data initilization, deaths are now set to 0. In order to visualize
                    // the extrapolated real number of deaths, they have to be set here. In the comparison of data
                    // it has to be paid attention to the fact, the the simulation starts with deaths=0
                    // while this method starts with deaths=number of reported deaths so far...
                    // Additionally, we set the number of reported deaths to DeadNaive since no information on that is
                    // available here.
                    // Do only add deaths after substraction.
                    extrapolated_rki[county][day]((size_t)InfectionState::DeadNaive + age_group_offset) =
                        num_death_n[county][age];
                }
            }
            else {
                log_warning("No population data available for region " + std::to_string(county) +
                            ". Population data has not been set.");
            }
        }
        log_info("extrapolated real data for date: {}-{}-{}", date.day, date.month, date.year);
    }
    /* end: similar functionality in set_confirmed_cases_data(), here only for vector of TimeSeries */
    auto num_groups = (int)(size_t)model[0].parameters.get_num_groups();
    BOOST_OUTCOME_TRY(save_result(extrapolated_rki, region, num_groups, path_join(dir, "Results_rki.h5")));

    auto extrapolated_rki_data_sum = sum_nodes(std::vector<std::vector<TimeSeries<double>>>{extrapolated_rki});
    BOOST_OUTCOME_TRY(
        save_result({extrapolated_rki_data_sum[0][0]}, {0}, num_groups, path_join(dir, "Results_rki_sum.h5")));

    return success();
}

template <class Model>
IOResult<void> export_input_data_county_timeseries(
    std::vector<Model>&& model, const std::string& dir, std::vector<int> const& region, Date date,
    const std::vector<double>& scaling_factor_inf, double scaling_factor_icu, int num_days,
    const std::string& divi_data_path, const std::string& confirmed_cases_path, const std::string& population_data_path,
    bool set_vaccination_data, const std::string& vaccination_data_path, const std::string& immunity_data_path)
{
    if (set_vaccination_data) {
        BOOST_OUTCOME_TRY(details::set_vaccination_data(model, vaccination_data_path, date, region, num_days));
    }

    BOOST_OUTCOME_TRY(export_input_data_county_timeseries(
        model, dir, region, date, scaling_factor_inf, scaling_factor_icu, num_days, divi_data_path,
        confirmed_cases_path, population_data_path, immunity_data_path));

    return success();
}

#endif //MEMILIO_HAS_HDF5

/**
    * Reads compartments for German counties at a specified date from data files.
    * Estimates all compartments from available data using the model parameters, so the 
    * model parameters must be set before calling this function.
    * Uses data files that contain centered 7-day moving average.
    * @param model Vector of SECIRS-type models, one per county.
    * @param date Date for which the data should be read.
    * @param county Ids of the counties.
    * @param scaling_factor_inf Factor of confirmed cases to account for undetected cases in each county.
    * @param scaling_factor_icu Factor of ICU cases to account for underreporting.
    * @param dir Directory that contains the data files.
    * @param num_days Number of days to be simulated; required to load data for vaccinations during the simulation.
    * @param export_time_series If true, reads data for each day of simulation and writes it in the same directory as the input files.
    */
template <class Model>
IOResult<void> read_input_data_county(std::vector<Model>& model, Date date, const std::vector<int>& county,
                                      const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                      const std::string& dir, int num_days, bool export_time_series = false)
{
    BOOST_OUTCOME_TRY(details::set_vaccination_data(
        model, path_join(dir, "pydata/Germany", "all_county_ageinf_vacc_ma7.json"), date, county, num_days));

    // TODO: Reuse more code, e.g., set_divi_data (in secir) and a set_divi_data (here) only need a different ModelType.
    // TODO: add option to set ICU data from confirmed cases if DIVI or other data is not available.
    if (date > Date(2020, 4, 23)) {
        BOOST_OUTCOME_TRY(details::set_divi_data(model, path_join(dir, "pydata/Germany", "county_divi_ma7.json"),
                                                 county, date, scaling_factor_icu));
    }
    else {
        log_warning("No DIVI data available for this date");
    }

    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(
        model, path_join(dir, "pydata/Germany", "cases_all_county_age_ma7.json"), county, date, scaling_factor_inf,
        path_join(dir, "pydata/Germany", "immunity_population.txt")));
    BOOST_OUTCOME_TRY(
        details::set_population_data(model, path_join(dir, "pydata/Germany", "county_current_population.json"), county,
                                     path_join(dir, "pydata/Germany", "immunity_population.txt")));

    if (export_time_series) {
        // Use only if extrapolated real data is needed for comparison. EXPENSIVE !
        // Run time equals run time of the previous functions times the num_days !
        // (This only represents the vectorization of the previous function over all simulation days...)
        log_warning("Exporting time series of extrapolated real data. This may take some minutes. "
                    "For simulation runs over the same time period, deactivate it.");
        BOOST_OUTCOME_TRY(
            export_input_data_county_timeseries(model, dir, county, date, scaling_factor_inf, scaling_factor_icu,
                                                num_days, path_join(dir, "pydata/Germany", "county_divi_ma7.json"),
                                                path_join(dir, "pydata/Germany", "cases_all_county_age_ma7.json"),
                                                path_join(dir, "pydata/Germany", "county_current_population.json"),
                                                path_join(dir, "pydata/Germany", "immunity_population.txt")));
    }

    return success();
}

/**
    * Reads compartments for German counties at a specified date from data files.
    * Estimates all compartments from available data using the model parameters, so the 
    * model parameters must be set before calling this function.
    * Uses data files that contain centered 7-day moving average.
    * @param model Vector of SECIRS-type models, one per county.
    * @param date Date for which the data should be read.
    * @param county Ids of the counties.
    * @param scaling_factor_inf Factor of confirmed cases to account for undetected cases in each county.
    * @param scaling_factor_icu Factor of ICU cases to account for underreporting.
    * @param dir Directory that contains the data files.
    * @param num_days Number of days to be simulated; required to load data for vaccinations during the simulation.
    * @param export_time_series If true, reads data for each day of simulation and writes it in the same directory as the input files.
    */
template <class Model>
IOResult<void> read_input_data(std::vector<Model>& model, Date date, const std::vector<int>& node_ids,
                               const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                               const std::string& data_dir, int num_days, bool export_time_series = false)
{

    BOOST_OUTCOME_TRY(
        details::set_vaccination_data(model, path_join(data_dir, "vaccination_data.json"), date, node_ids, num_days));

    // TODO: Reuse more code, e.g., set_divi_data (in secir) and a set_divi_data (here) only need a different ModelType.
    // TODO: add option to set ICU data from confirmed cases if DIVI or other data is not available.
    if (date > Date(2020, 4, 23)) {
        BOOST_OUTCOME_TRY(details::set_divi_data(model, path_join(data_dir, "critical_cases.json"), node_ids, date,
                                                 scaling_factor_icu));
    }
    else {
        log_warning("No DIVI data available for this date");
    }

    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(model, path_join(data_dir, "confirmed_cases.json"), node_ids,
                                                        date, scaling_factor_inf,
                                                        path_join(data_dir, "immunity_population.txt")));
    BOOST_OUTCOME_TRY(details::set_population_data(model, path_join(data_dir, "population_data.json"), node_ids,
                                                   path_join(data_dir, "immunity_population.txt")));

    if (export_time_series) {
        // Use only if extrapolated real data is needed for comparison. EXPENSIVE !
        // Run time equals run time of the previous functions times the num_days !
        // (This only represents the vectorization of the previous function over all simulation days...)
        log_warning("Exporting time series of extrapolated real data. This may take some minutes. "
                    "For simulation runs over the same time period, deactivate it.");
        BOOST_OUTCOME_TRY(export_input_data_county_timeseries(
            model, data_dir, node_ids, date, scaling_factor_inf, scaling_factor_icu, num_days,
            path_join(data_dir, "critical_cases.json"), path_join(data_dir, "confirmed_cases.json"),
            path_join(data_dir, "population_data.json"), path_join(data_dir, "immunity_population.txt")));
    }

    return success();
}

} // namespace osecirts
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // ODESECIRTS_PARAMETERS_IO_H
