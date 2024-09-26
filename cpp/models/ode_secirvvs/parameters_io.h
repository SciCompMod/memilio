/* 
* Copyright (C) 2020-2024 MEmilio
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
#ifndef ODESECIRVVS_PARAMETERS_IO_H
#define ODESECIRVVS_PARAMETERS_IO_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "ode_secirvvs/model.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/date.h"
#include "memilio/utils/logging.h"

namespace mio
{
namespace osecirvvs
{

namespace details
{
/**
        * @brief Reads subpopulations of infection states from transformed RKI cases file.
        * @param[in] path Path to transformed RKI cases file.
        * @param[in] vregion Vector of keys of the region of interest.
        * @param[in] date Date for which the arrays are initialized.
        * @param[in, out] num_* Output vector for number of people in the corresponding compartement.
        * @param[in] vt_* Average time it takes to get from one compartement to another (vector with one element per age group).
        * @param[in] vmu_* Probabilities to get from one compartement to another (vector with one element per age group).
        * @param[in] scaling_factor_inf Factor for scaling the confirmed cases to account for estimated undetected cases.
        * @see mio::read_confirmed_cases_data
        * @{
        */
IOResult<void> read_confirmed_cases_data(
    std::string const& path, std::vector<int> const& vregion, Date date, std::vector<std::vector<double>>& num_Exposed,
    std::vector<std::vector<double>>& num_InfectedNoSymptoms, std::vector<std::vector<double>>& num_InfectedSymptoms,
    std::vector<std::vector<double>>& num_InfectedSevere, std::vector<std::vector<double>>& num_icu,
    std::vector<std::vector<double>>& num_death, std::vector<std::vector<double>>& num_rec,
    const std::vector<std::vector<int>>& t_Exposed, const std::vector<std::vector<int>>& t_InfectedNoSymptoms,
    const std::vector<std::vector<int>>& t_InfectedSymptoms, const std::vector<std::vector<int>>& t_InfectedSevere,
    const std::vector<std::vector<int>>& t_InfectedCritical, const std::vector<std::vector<double>>& mu_C_R,
    const std::vector<std::vector<double>>& mu_I_H, const std::vector<std::vector<double>>& mu_H_U,
    const std::vector<double>& scaling_factor_inf);

IOResult<void> read_confirmed_cases_data(
    const std::vector<ConfirmedCasesDataEntry>& rki_data, std::vector<int> const& vregion, Date date,
    std::vector<std::vector<double>>& num_Exposed, std::vector<std::vector<double>>& num_InfectedNoSymptoms,
    std::vector<std::vector<double>>& num_InfectedSymptoms, std::vector<std::vector<double>>& num_InfectedSevere,
    std::vector<std::vector<double>>& num_icu, std::vector<std::vector<double>>& num_death,
    std::vector<std::vector<double>>& num_rec, const std::vector<std::vector<int>>& t_Exposed,
    const std::vector<std::vector<int>>& t_InfectedNoSymptoms, const std::vector<std::vector<int>>& t_InfectedSymptoms,
    const std::vector<std::vector<int>>& t_InfectedSevere, const std::vector<std::vector<int>>& t_InfectedCritical,
    const std::vector<std::vector<double>>& mu_C_R, const std::vector<std::vector<double>>& mu_I_H,
    const std::vector<std::vector<double>>& mu_H_U, const std::vector<double>& scaling_factor_inf);
/**@}*/

/**
        * @brief Reads confirmed cases data and translates data of day t0-delay to recovered compartment.
        * @param[in] path Path to RKI confirmed cases file.
        * @param[in] vregion Vector of keys of the region of interest.
        * @param[in] date Date for which the arrays are initialized.
        * @param[in, out] num_rec Output vector for number of people in the compartement recovered.
        * @param[in] delay Number of days in the past the are used to set recovered compartment.
        * @see mio::read_confirmed_cases_data
        * @{
        */
IOResult<void> read_confirmed_cases_data_fix_recovered(const std::vector<ConfirmedCasesDataEntry>& rki_data,
                                                       std::vector<int> const& vregion, Date date,
                                                       std::vector<std::vector<double>>& vnum_rec, double delay = 14.);
IOResult<void> read_confirmed_cases_data_fix_recovered(std::string const& path, std::vector<int> const& vregion,
                                                       Date date, std::vector<std::vector<double>>& vnum_rec,
                                                       double delay = 14.);
/**@}*/

/**
        * @brief Sets the confirmed cases data for a vector of models based on input data.
        * @param[in, out] model Vector of objects in which the data is set.
        * @param[in] case_data Vector of case data. Each inner vector represents a different region.
        * @param[in] region Vector of keys of the region of interest.
        * @param[in] date Date for which the arrays are initialized.
        * @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of RKI data.
        * @param[in] set_death If true, set the number of deaths.
        */
template <class Model>
IOResult<void> set_confirmed_cases_data(std::vector<Model>& model,
                                        const std::vector<ConfirmedCasesDataEntry>& case_data,
                                        std::vector<int> const& region, Date date,
                                        const std::vector<double>& scaling_factor_inf, bool set_death = false)
{
    auto num_age_groups = (size_t)model[0].parameters.get_num_groups();
    assert(scaling_factor_inf.size() == num_age_groups); //TODO: allow vector or scalar valued scaling factors
    assert(ConfirmedCasesDataEntry::age_group_names.size() == num_age_groups);

    std::vector<std::vector<int>> t_Exposed{model.size()};
    std::vector<std::vector<int>> t_InfectedNoSymptoms{model.size()};
    std::vector<std::vector<int>> t_InfectedSymptoms{model.size()};
    std::vector<std::vector<int>> t_InfectedSevere{model.size()};
    std::vector<std::vector<int>> t_InfectedCritical{model.size()};

    std::vector<std::vector<double>> mu_C_R{model.size()};
    std::vector<std::vector<double>> mu_I_H{model.size()};
    std::vector<std::vector<double>> mu_H_U{model.size()};

    std::vector<std::vector<double>> num_InfectedSymptoms(model.size());
    std::vector<std::vector<double>> num_death(model.size());
    std::vector<std::vector<double>> num_rec(model.size());
    std::vector<std::vector<double>> num_Exposed(model.size());
    std::vector<std::vector<double>> num_InfectedNoSymptoms(model.size());
    std::vector<std::vector<double>> num_InfectedSevere(model.size());
    std::vector<std::vector<double>> num_icu(model.size());

    /*----------- UNVACCINATED -----------*/
    for (size_t county = 0; county < model.size(); county++) {
        num_InfectedSymptoms[county]   = std::vector<double>(num_age_groups, 0.0);
        num_death[county]              = std::vector<double>(num_age_groups, 0.0);
        num_rec[county]                = std::vector<double>(num_age_groups, 0.0);
        num_Exposed[county]            = std::vector<double>(num_age_groups, 0.0);
        num_InfectedNoSymptoms[county] = std::vector<double>(num_age_groups, 0.0);
        num_InfectedSevere[county]     = std::vector<double>(num_age_groups, 0.0);
        num_icu[county]                = std::vector<double>(num_age_groups, 0.0);
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

            mu_C_R[county].push_back(
                model[county].parameters.template get<RecoveredPerInfectedNoSymptoms<double>>()[(AgeGroup)group]);
            mu_I_H[county].push_back(
                model[county].parameters.template get<SeverePerInfectedSymptoms<double>>()[(AgeGroup)group]);
            mu_H_U[county].push_back(
                model[county].parameters.template get<CriticalPerSevere<double>>()[(AgeGroup)group]);
        }
    }

    BOOST_OUTCOME_TRY(read_confirmed_cases_data(case_data, region, date, num_Exposed, num_InfectedNoSymptoms,
                                                num_InfectedSymptoms, num_InfectedSevere, num_icu, num_death, num_rec,
                                                t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
                                                t_InfectedCritical, mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf));

    for (size_t county = 0; county < model.size(); county++) {
        // if (std::accumulate(num_InfectedSymptoms[county].begin(), num_InfectedSymptoms[county].end(), 0.0) > 0) {
        size_t num_groups = (size_t)model[county].parameters.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            model[county].populations[{AgeGroup(i), InfectionState::ExposedNaive}] = num_Exposed[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsNaive}] =
                num_InfectedNoSymptoms[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsNaiveConfirmed}] = 0;
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsNaive}] =
                num_InfectedSymptoms[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsNaiveConfirmed}] = 0;
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSevereNaive}] =
                num_InfectedSevere[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::SusceptibleImprovedImmunity}] = num_rec[county][i];
            if (set_death) {
                // in set_confirmed_cases_data initilization, deaths are now set to 0. In order to visualize
                // the extrapolated real number of deaths, they have to be set here. In the comparison of data
                // it has to be paid attention to the fact, the the simulation starts with deaths=0
                // while this method starts with deaths=number of reported deaths so far...
                // Additionally, we set the number of reported deaths to DeadNaive since no information on that is
                // available here.
                // Do only add deaths after substraction.
                model[county].populations[{AgeGroup(i), InfectionState::DeadNaive}] = num_death[county][i];
            }
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
        num_rec[county]                = std::vector<double>(num_age_groups, 0.0);
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

    BOOST_OUTCOME_TRY(read_confirmed_cases_data(case_data, region, date, num_Exposed, num_InfectedNoSymptoms,
                                                num_InfectedSymptoms, num_InfectedSevere, num_icu, num_death, num_rec,
                                                t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
                                                t_InfectedCritical, mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf));

    for (size_t county = 0; county < model.size(); county++) {
        // if (std::accumulate(num_InfectedSymptoms[county].begin(), num_InfectedSymptoms[county].end(), 0.0) > 0) {
        size_t num_groups = (size_t)model[county].parameters.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            model[county].populations[{AgeGroup(i), InfectionState::ExposedPartialImmunity}] = num_Exposed[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsPartialImmunity}] =
                num_InfectedNoSymptoms[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}] = 0;
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsPartialImmunity}] =
                num_InfectedSymptoms[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsPartialImmunityConfirmed}] = 0;
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSeverePartialImmunity}] =
                num_InfectedSevere[county][i];
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
        num_rec[county]                = std::vector<double>(num_age_groups, 0.0);
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

    BOOST_OUTCOME_TRY(read_confirmed_cases_data(case_data, region, date, num_Exposed, num_InfectedNoSymptoms,
                                                num_InfectedSymptoms, num_InfectedSevere, num_icu, num_death, num_rec,
                                                t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
                                                t_InfectedCritical, mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf));

    for (size_t county = 0; county < model.size(); county++) {
        // if (std::accumulate(num_InfectedSymptoms[county].begin(), num_InfectedSymptoms[county].end(), 0.0) > 0) {
        size_t num_groups = (size_t)model[county].parameters.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            model[county].populations[{AgeGroup(i), InfectionState::ExposedImprovedImmunity}] = num_Exposed[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsImprovedImmunity}] =
                num_InfectedNoSymptoms[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] = 0;
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsImprovedImmunity}] =
                num_InfectedSymptoms[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] = 0;
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSevereImprovedImmunity}] =
                num_InfectedSevere[county][i];
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
        * @brief sets populations data from a transformed RKI cases file into a Model.
        * @param[in, out] model vector of objects in which the data is set
        * @param[in] path Path to transformed RKI cases file
        * @param[in] region vector of keys of the region of interest
        * @param[in] date Date for which the arrays are initialized
        * @param[in] scaling_factor_inf factors by which to scale the confirmed cases of
        * rki data
        * @param set_death[in] If true, set the number of deaths.
        */
template <class Model>
IOResult<void> set_confirmed_cases_data(std::vector<Model>& model, const std::string& path,
                                        std::vector<int> const& region, Date date,
                                        const std::vector<double>& scaling_factor_inf, bool set_death = false)
{
    BOOST_OUTCOME_TRY(auto&& case_data, mio::read_confirmed_cases_data(path));
    BOOST_OUTCOME_TRY(set_confirmed_cases_data(model, case_data, region, date, scaling_factor_inf, set_death));
    return success();
}

/**
        * @brief reads number of ICU patients from DIVI register into Parameters
        * @param[in] path Path to transformed DIVI file
        * @param[in] vregion Keys of the region of interest
        * @param[in] date Date for which the arrays are initialized
        * @param[in, out] vnum_icu number of ICU patients
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
        * @param[in, out] model vector of objects in which the data is set
        * @param[in] path Path to transformed DIVI file
        * @param[in] vregion vector of keys of the regions of interest
        * @param[in] date Date for which the arrays are initialized
        * @param[in] scaling_factor_icu factor by which to scale the icu cases of divi data
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
        * @param[in] path Path to population data file.
        * @param[in] vregion vector of keys of the regions of interest
        * @see mio::read_population_data
        * @{
        */
IOResult<std::vector<std::vector<double>>> read_population_data(const std::string& path,
                                                                const std::vector<int>& vregion);
IOResult<std::vector<std::vector<double>>> read_population_data(const std::vector<PopulationDataEntry>& population_data,
                                                                const std::vector<int>& vregion);
/**@}*/

/**
* @brief sets population data from census data which has been read into num_population
* @param[in, out] model vector of objects in which the data is set
* @param[in] num_population vector of population data
* @param[in] case_data vector of case data. Each inner vector represents a different region
* @param[in] vregion vector of keys of the regions of interest
* @param[in] date Date for which the arrays are initialized
*/
template <class Model>
IOResult<void> set_population_data(std::vector<Model>& model, const std::vector<std::vector<double>>& num_population,
                                   const std::vector<ConfirmedCasesDataEntry>& case_data,
                                   const std::vector<int>& vregion, Date date)
{
    auto num_age_groups = ConfirmedCasesDataEntry::age_group_names.size();
    std::vector<std::vector<double>> vnum_rec(model.size(), std::vector<double>(num_age_groups, 0.0));

    BOOST_OUTCOME_TRY(read_confirmed_cases_data_fix_recovered(case_data, vregion, date, vnum_rec, 14.));

    for (size_t region = 0; region < vregion.size(); region++) {
        if (std::accumulate(num_population[region].begin(), num_population[region].end(), 0.0) > 0) {
            auto num_groups = model[region].parameters.get_num_groups();
            for (auto i = AgeGroup(0); i < num_groups; i++) {

                double S_v = std::min(
                    model[region].parameters.template get<DailyFullVaccination<double>>()[{i, SimulationDay(0)}] +
                        vnum_rec[region][size_t(i)],
                    num_population[region][size_t(i)]);
                double S_pv = std::max(
                    model[region].parameters.template get<DailyFirstVaccination<double>>()[{i, SimulationDay(0)}] -
                        model[region].parameters.template get<DailyFullVaccination<double>>()[{i, SimulationDay(0)}],
                    0.0); // use std::max with 0
                double S;
                if (num_population[region][size_t(i)] - S_pv - S_v < 0.0) {
                    log_warning("Number of vaccinated persons greater than population in county {}, age group {}.",
                                region, size_t(i));
                    S   = 0.0;
                    S_v = num_population[region][size_t(i)] - S_pv;
                }
                else {
                    S = num_population[region][size_t(i)] - S_pv - S_v;
                }

                double denom_E =
                    1 / (S + S_pv * model[region].parameters.template get<ReducExposedPartialImmunity<double>>()[i] +
                         S_v * model[region].parameters.template get<ReducExposedImprovedImmunity<double>>()[i]);
                double denom_C =
                    1 / (S + S_pv * model[region].parameters.template get<ReducExposedPartialImmunity<double>>()[i] +
                         S_v * model[region].parameters.template get<ReducExposedImprovedImmunity<double>>()[i]);
                double denom_I =
                    1 /
                    (S +
                     S_pv * model[region].parameters.template get<ReducInfectedSymptomsPartialImmunity<double>>()[i] +
                     S_v * model[region].parameters.template get<ReducInfectedSymptomsImprovedImmunity<double>>()[i]);
                double denom_HU =
                    1 /
                    (S +
                     S_pv * model[region]
                                .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[i] +
                     S_v * model[region]
                               .parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[i]);

                model[region].populations[{i, InfectionState::ExposedNaive}] =
                    S * model[region].populations[{i, InfectionState::ExposedNaive}] * denom_E;
                model[region].populations[{i, InfectionState::ExposedPartialImmunity}] =
                    S_pv * model[region].parameters.template get<ReducExposedPartialImmunity<double>>()[i] *
                    model[region].populations[{i, InfectionState::ExposedPartialImmunity}] * denom_E;
                model[region].populations[{i, InfectionState::ExposedImprovedImmunity}] =
                    S_v * model[region].parameters.template get<ReducExposedImprovedImmunity<double>>()[i] *
                    model[region].populations[{i, InfectionState::ExposedImprovedImmunity}] * denom_E;

                model[region].populations[{i, InfectionState::InfectedNoSymptomsNaive}] =
                    S * model[region].populations[{i, InfectionState::InfectedNoSymptomsNaive}] * denom_C;
                model[region].populations[{i, InfectionState::InfectedNoSymptomsPartialImmunity}] =
                    S_pv * model[region].populations[{i, InfectionState::InfectedNoSymptomsPartialImmunity}] * denom_C;
                model[region].populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunity}] =
                    S_v * model[region].populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunity}] * denom_C;

                model[region].populations[{i, InfectionState::InfectedNoSymptomsNaiveConfirmed}] =
                    S * model[region].populations[{i, InfectionState::InfectedNoSymptomsNaiveConfirmed}] * denom_C;
                model[region].populations[{i, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}] =
                    S_pv * model[region].populations[{i, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}] *
                    denom_C;
                model[region].populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] =
                    S_v * model[region].populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] *
                    denom_C;

                model[region].populations[{i, InfectionState::InfectedSymptomsNaive}] =
                    S * model[region].populations[{i, InfectionState::InfectedSymptomsNaive}] * denom_I;
                model[region].populations[{i, InfectionState::InfectedSymptomsPartialImmunity}] =
                    S_pv * model[region].parameters.template get<ReducInfectedSymptomsPartialImmunity<double>>()[i] *
                    model[region].populations[{i, InfectionState::InfectedSymptomsPartialImmunity}] * denom_I;
                model[region].populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] =
                    S_v * model[region].parameters.template get<ReducInfectedSymptomsImprovedImmunity<double>>()[i] *
                    model[region].populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] * denom_I;

                model[region].populations[{i, InfectionState::InfectedSymptomsNaiveConfirmed}] =
                    S * model[region].populations[{i, InfectionState::InfectedSymptomsNaiveConfirmed}] * denom_I;
                model[region].populations[{i, InfectionState::InfectedSymptomsPartialImmunityConfirmed}] =
                    S_pv * model[region].parameters.template get<ReducInfectedSymptomsPartialImmunity<double>>()[i] *
                    model[region].populations[{i, InfectionState::InfectedSymptomsPartialImmunityConfirmed}] * denom_I;
                model[region].populations[{i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] =
                    S_v * model[region].parameters.template get<ReducInfectedSymptomsImprovedImmunity<double>>()[i] *
                    model[region].populations[{i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] * denom_I;

                model[region].populations[{i, InfectionState::InfectedSevereNaive}] =
                    S * model[region].populations[{i, InfectionState::InfectedSevereNaive}] * denom_HU;
                model[region].populations[{i, InfectionState::InfectedSeverePartialImmunity}] =
                    S_pv *
                    model[region].parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[i] *
                    model[region].populations[{i, InfectionState::InfectedSeverePartialImmunity}] * denom_HU;
                model[region].populations[{i, InfectionState::InfectedSevereImprovedImmunity}] =
                    S_v *
                    model[region]
                        .parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[i] *
                    model[region].populations[{i, InfectionState::InfectedSevereImprovedImmunity}] * denom_HU;

                model[region].populations[{i, InfectionState::InfectedCriticalPartialImmunity}] =
                    S_pv *
                    model[region].parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[i] *
                    model[region].populations[{i, InfectionState::InfectedCriticalNaive}] * denom_HU;
                model[region].populations[{i, InfectionState::InfectedCriticalImprovedImmunity}] =
                    S_v *
                    model[region]
                        .parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[i] *
                    model[region].populations[{i, InfectionState::InfectedCriticalNaive}] * denom_HU;
                model[region].populations[{i, InfectionState::InfectedCriticalNaive}] =
                    S * model[region].populations[{i, InfectionState::InfectedCriticalNaive}] * denom_HU;

                model[region].populations[{i, InfectionState::SusceptibleImprovedImmunity}] =
                    model[region].parameters.template get<DailyFullVaccination<double>>()[{i, SimulationDay(0)}] +
                    model[region].populations[{i, InfectionState::SusceptibleImprovedImmunity}] -
                    (model[region].populations[{i, InfectionState::InfectedSymptomsNaive}] +
                     model[region].populations[{i, InfectionState::InfectedSymptomsPartialImmunity}] +
                     model[region].populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] +
                     model[region].populations[{i, InfectionState::InfectedSymptomsNaiveConfirmed}] +
                     model[region].populations[{i, InfectionState::InfectedSymptomsPartialImmunityConfirmed}] +
                     model[region].populations[{i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] +
                     model[region].populations[{i, InfectionState::InfectedSevereNaive}] +
                     model[region].populations[{i, InfectionState::InfectedSeverePartialImmunity}] +
                     model[region].populations[{i, InfectionState::InfectedSevereImprovedImmunity}] +
                     model[region].populations[{i, InfectionState::InfectedCriticalNaive}] +
                     model[region].populations[{i, InfectionState::InfectedCriticalPartialImmunity}] +
                     model[region].populations[{i, InfectionState::InfectedCriticalImprovedImmunity}] +
                     model[region].populations[{i, InfectionState::DeadNaive}] +
                     model[region].populations[{i, InfectionState::DeadPartialImmunity}] +
                     model[region].populations[{i, InfectionState::DeadImprovedImmunity}]);

                model[region].populations[{i, InfectionState::SusceptibleImprovedImmunity}] = std::min(
                    S_v - model[region].populations[{i, InfectionState::ExposedImprovedImmunity}] -
                        model[region].populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunity}] -
                        model[region].populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] -
                        model[region].populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] -
                        model[region].populations[{i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] -
                        model[region].populations[{i, InfectionState::InfectedSevereImprovedImmunity}] -
                        model[region].populations[{i, InfectionState::InfectedCriticalImprovedImmunity}] -
                        model[region].populations[{i, InfectionState::DeadImprovedImmunity}],
                    std::max(0.0, double(model[region].populations[{i, InfectionState::SusceptibleImprovedImmunity}])));

                model[region].populations[{i, InfectionState::SusceptiblePartialImmunity}] = std::max(
                    0.0,
                    S_pv - model[region].populations[{i, InfectionState::ExposedPartialImmunity}] -
                        model[region].populations[{i, InfectionState::InfectedNoSymptomsPartialImmunity}] -
                        model[region].populations[{i, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}] -
                        model[region].populations[{i, InfectionState::InfectedSymptomsPartialImmunity}] -
                        model[region].populations[{i, InfectionState::InfectedSymptomsPartialImmunityConfirmed}] -
                        model[region].populations[{i, InfectionState::InfectedSeverePartialImmunity}] -
                        model[region].populations[{i, InfectionState::InfectedCriticalPartialImmunity}] -
                        model[region].populations[{i, InfectionState::DeadPartialImmunity}]);

                model[region].populations.template set_difference_from_group_total<AgeGroup>(
                    {i, InfectionState::SusceptibleNaive}, num_population[region][size_t(i)]);
            }

            for (auto i = AgeGroup(0); i < AgeGroup(6); i++) {
                for (auto j = Index<InfectionState>(0); j < InfectionState::Count; ++j) {
                    if (model[region].populations[{i, j}] < 0) {
                        log_warning("Compartment at age group {}, infection state {}, is negative: {}", size_t(i),
                                    size_t(j), model[region].populations[{i, j}] / num_population[region][size_t(i)]);
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
* @brief Sets population data from census data which has been read into num_population.
* @param[in, out] model Vector of objects in which the data is set.
* @param[in] path Path to population data file.
* @param[in] path_rki Path to RKI cases data file.
* @param[in] vregion Vector of keys of the regions of interest.
* @param[in] date Date for which the arrays are initialized.
*/
template <class Model>
IOResult<void> set_population_data(std::vector<Model>& model, const std::string& path, const std::string& path_rki,
                                   const std::vector<int>& vregion, Date date)
{
    BOOST_OUTCOME_TRY(auto&& num_population, details::read_population_data(path, vregion));
    BOOST_OUTCOME_TRY(auto&& rki_data, mio::read_confirmed_cases_data(path_rki));

    BOOST_OUTCOME_TRY(set_population_data(model, num_population, rki_data, vregion, date));
    return success();
}

/**
 * @brief Sets vaccination data for models stored in a vector.
 * 
 * @tparam FP Floating point type used in the Model objects.
 * @param[in, out] model Vector of Model objects in which the vaccination data is set.
 * @param[in] vacc_data Vector of VaccinationDataEntry objects containing the vaccination data.
 * @param[in] date Start date for the simulation.
 * @param[in] vregion Vector of region identifiers.
 * @param[in] num_days Number of days for which the simulation is run.
 */
template <typename FP = double>
IOResult<void> set_vaccination_data(std::vector<Model<FP>>& model, const std::vector<VaccinationDataEntry>& vacc_data,
                                    Date date, const std::vector<int>& vregion, int num_days)
{
    auto num_groups = model[0].parameters.get_num_groups();

    // type conversion from UncertainValue -> FP -> int
    auto days_until_effective1 = static_cast<int>(
        static_cast<FP>(model[0].parameters.template get<DaysUntilEffectivePartialImmunity<FP>>()[AgeGroup(0)]));
    auto days_until_effective2 = static_cast<int>(
        static_cast<FP>(model[0].parameters.template get<DaysUntilEffectiveImprovedImmunity<FP>>()[AgeGroup(0)]));
    auto vaccination_distance =
        static_cast<int>(static_cast<FP>(model[0].parameters.template get<VaccinationGap<FP>>()[AgeGroup(0)]));

    // iterate over regions (e.g., counties)
    for (size_t i = 0; i < model.size(); ++i) {
        // iterate over age groups in region
        for (auto g = AgeGroup(0); g < num_groups; ++g) {

            model[i].parameters.template get<DailyFirstVaccination<FP>>().resize(SimulationDay(num_days + 1));
            model[i].parameters.template get<DailyFullVaccination<FP>>().resize(SimulationDay(num_days + 1));
            for (auto d = SimulationDay(0); d < SimulationDay(num_days + 1); ++d) {
                model[i].parameters.template get<DailyFirstVaccination<FP>>()[{g, d}] = 0.0;
                model[i].parameters.template get<DailyFullVaccination<FP>>()[{g, d}]  = 0.0;
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

    for (auto&& vacc_data_entry : vacc_data) {
        auto it      = std::find_if(vregion.begin(), vregion.end(), [&vacc_data_entry](auto&& r) {
            return r == 0 || (vacc_data_entry.county_id && vacc_data_entry.county_id == regions::CountyId(r)) ||
                   (vacc_data_entry.state_id && vacc_data_entry.state_id == regions::StateId(r)) ||
                   (vacc_data_entry.district_id && vacc_data_entry.district_id == regions::DistrictId(r));
        });
        auto date_df = vacc_data_entry.date;
        if (it != vregion.end()) {
            auto region_idx = size_t(it - vregion.begin());
            auto age        = vacc_data_entry.age_group;

            for (size_t d = 0; d < (size_t)num_days + 1; ++d) {
                int days_plus;
                // In the following, second dose means previous 'full immunization', now 'Grundimmunisierung'.
                // ---
                // date: start_date of the simulation (Input from IO call read_input_data_county_vaccmodel())
                // d: day of simulation, counted from 0 to num_days (for which we need (approximated) vaccination numbers)
                // root[i]["Vacc_completed"]: accumulated number of total second doses up to day date_df;
                //                               taken from input dataframe, single value, per county and age group
                // ----
                // An averaged distance between first and second doses (vaccination_distance) is assumed in the following
                // and the first doses are computed based on the second doses given 'vaccination_distance' days later.
                // ----
                // a person whose second dose is reported at start_date + simulation_day - days_until_effective1 + vaccination_distance
                // had the first dose on start_date + simulation_day - days_until_effective1. Furthermore, he/she has the full protection
                // of the first dose at day X = start_date + simulation_day
                // Storing its value in get<DailyFirstVaccination>() will eventually (in the simulation)
                // transfer the difference (between get<DailyFirstVaccination>() at d and d-1) of
                // N susceptible individuals to 'Susceptible Partially Vaccinated' state at day d; see secir_vaccinated.h
                auto offset_first_date =
                    offset_date_by_days(date, (int)d - days_until_effective1 + vaccination_distance);
                if (max_date >= offset_first_date) {
                    // Option 1: considered offset_first_date is available in input data frame
                    if (date_df == offset_first_date) {
                        model[region_idx]
                            .parameters.template get<DailyFirstVaccination<FP>>()[{age, SimulationDay(d)}] =
                            vacc_data_entry.num_vaccinations_completed;
                    }
                }
                else { // offset_first_date > max_date
                    // Option 2: considered offset_first_date is NOT available in input data frame
                    // Here, a constant number of first and second doses is assumed, i.e.,
                    // the the number of vaccinationes at day d (N days after max_date) will be:
                    // total number of vaccinations up to day max_date + N * number of vaccinations ON max_date
                    // (where the latter is computed as the difference between the total number at max_date and max_date-1)
                    days_plus = get_offset_in_days(offset_first_date, max_date);
                    if (date_df == offset_date_by_days(max_date, -1)) {
                        model[region_idx]
                            .parameters.template get<DailyFirstVaccination<FP>>()[{age, SimulationDay(d)}] -=
                            days_plus * vacc_data_entry.num_vaccinations_completed;
                    }
                    else if (date_df == max_date) {
                        model[region_idx]
                            .parameters.template get<DailyFirstVaccination<FP>>()[{age, SimulationDay(d)}] +=
                            (days_plus + 1) * vacc_data_entry.num_vaccinations_completed;
                    }
                }

                // a person whose second dose is reported at start_date + simulation_day - days_until_effective2
                // has the full protection of the second dose at day X = start_date + simulation_day
                // Storing its value in get<DailyFullVaccination>() will eventually (in the simulation)
                // transfer the difference (between get<DailyFullVaccination>() at d and d-1) of
                // N susceptible, partially vaccinated individuals to 'SusceptibleImprovedImmunity' state at day d; see secir_vaccinated.h
                auto offset_full_date = offset_date_by_days(date, (int)d - days_until_effective2);
                if (max_date >= offset_full_date) {
                    // Option 1: considered offset_full_date is available in input data frame
                    if (date_df == offset_full_date) {
                        model[region_idx].parameters.template get<DailyFullVaccination<FP>>()[{age, SimulationDay(d)}] =
                            vacc_data_entry.num_vaccinations_completed;
                    }
                }
                else { // offset_full_date > max_full_date
                    // Option 2: considered offset_full_date is NOT available in input data frame
                    days_plus = get_offset_in_days(offset_full_date, max_date);
                    if (date_df == offset_date_by_days(max_date, -1)) {
                        model[region_idx]
                            .parameters.template get<DailyFullVaccination<FP>>()[{age, SimulationDay(d)}] -=
                            days_plus * vacc_data_entry.num_vaccinations_completed;
                    }
                    else if (date_df == max_date) {
                        model[region_idx]
                            .parameters.template get<DailyFullVaccination<FP>>()[{age, SimulationDay(d)}] +=
                            (days_plus + 1) * vacc_data_entry.num_vaccinations_completed;
                    }
                }
            }
        }
    }
    return success();
}

/**
 * @brief Reads vaccination data from a file and sets it for each model.
 * 
 * @tparam FP Floating point type used in the Model objects.
 * @param[in, out] model Vector of Model objects in which the vaccination data is set.
 * @param[in] path Path to vaccination data file.
 * @param[in] date Start date for the simulation.
 * @param[in] vregion Vector of region identifiers.
 * @param[in] num_days Number of days for which the simulation is run.
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
* @brief Uses the initialisation method, which uses the reported data to set the initial conditions for the model for a given day. 
* The initialisation is applied for a predefined number of days and finally saved in a timeseries for each region. In the end,
* we save the files "Results_rki.h5" and "Results_rki_sum.h5" in the results_dir.
* Results_rki.h5 contains a time series for each region and Results_rki_sum.h5 contains the sum of all regions.
* @param[in] models Vector of models in which the data is set. Copy is made to avoid changing the original model.
* @param[in] results_dir Path to result files.
* @param[in] counties Vector of keys of the counties of interest.
* @param[in] date Date for which the data should be read.
* @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
* @param[in] scaling_factor_icu Factor by which to scale the icu cases of divi data.
* @param[in] num_days Number of days to be simulated/initialized.
* @param[in] divi_data_path Path to DIVI file.
* @param[in] confirmed_cases_path Path to confirmed cases file.
* @param[in] population_data_path Path to population data file.
* @param[in] vaccination_data_path Path to vaccination data file.
*/
template <class Model>
IOResult<void> export_input_data_county_timeseries(
    std::vector<Model> models, const std::string& results_dir, const std::vector<int>& counties, Date date,
    const std::vector<double>& scaling_factor_inf, const double scaling_factor_icu, const int num_days,
    const std::string& divi_data_path, const std::string& confirmed_cases_path, const std::string& population_data_path,
    const std::string& vaccination_data_path = "")
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
        if (offset_day > Date(2020, 4, 23)) {
            BOOST_OUTCOME_TRY(details::set_divi_data(models, divi_data_path, counties, offset_day, scaling_factor_icu));
        }
        else {
            log_warning("No DIVI data available for for date: {}-{}-{}", offset_day.day, offset_day.month,
                        offset_day.year);
        }

        BOOST_OUTCOME_TRY(
            details::set_confirmed_cases_data(models, case_data, counties, offset_day, scaling_factor_inf, true));
        BOOST_OUTCOME_TRY(details::set_population_data(models, population_data, case_data, counties, offset_day));

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
                                                   const std::string&)
{
    mio::log_warning("HDF5 not available. Cannot export time series of extrapolated real data.");
    return success();
}

#endif //MEMILIO_HAS_HDF5

/**
    * Reads compartments for German counties at a specified date from data files.
    * Estimates all compartments from available data using the model parameters, so the 
    * model parameters must be set before calling this function.
    * Uses data files that contain centered 7-day moving average.
    * @param[in, out] model Vector of SECIRVVS models, one per county.
    * @param[in] date Date for which the data should be read.
    * @param[in] county Ids of the counties.
    * @param[in] scaling_factor_inf Factor of confirmed cases to account for undetected cases in each county.
    * @param[in] scaling_factor_icu Factor of ICU cases to account for underreporting.
    * @param[in] dir Directory that contains the data files.
    * @param[in] num_days Number of days to be simulated; required to load data for vaccinations during the simulation.
    * @param[in] export_time_series If true, reads data for each day of simulation and writes it in the same directory as the input files.
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
        model, path_join(dir, "pydata/Germany", "cases_all_county_age_ma7.json"), county, date, scaling_factor_inf));
    BOOST_OUTCOME_TRY(
        details::set_population_data(model, path_join(dir, "pydata/Germany", "county_current_population.json"),
                                     path_join(dir, "pydata/Germany", "cases_all_county_age_ma7.json"), county, date));

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
                                                path_join(dir, "pydata/Germany", "all_county_ageinf_vacc_ma7.json")));
    }

    return success();
}

/**
    * Reads compartments for German counties at a specified date from data files.
    * Estimates all compartments from available data using the model parameters, so the 
    * model parameters must be set before calling this function.
    * Uses data files that contain centered 7-day moving average.
    * @param[in, out] model Vector of SECIRVVS models, one per county.
    * @param[in] date Date for which the data should be read.
    * @param[in] node_ids Ids of the nodes.
    * @param[in] scaling_factor_inf Factor of confirmed cases to account for undetected cases in each county.
    * @param[in] scaling_factor_icu Factor of ICU cases to account for underreporting.
    * @param[in] dir Directory that contains the data files.
    * @param[in] num_days Number of days to be simulated; required to load data for vaccinations during the simulation.
    * @param[in] export_time_series If true, reads data for each day of simulation and writes it in the same directory as the input files.
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
                                                        date, scaling_factor_inf));
    BOOST_OUTCOME_TRY(details::set_population_data(model, path_join(data_dir, "population_data.json"),
                                                   path_join(data_dir, "confirmed_cases.json"), node_ids, date));

    if (export_time_series) {
        // Use only if extrapolated real data is needed for comparison. EXPENSIVE !
        // Run time equals run time of the previous functions times the num_days !
        // (This only represents the vectorization of the previous function over all simulation days...)
        log_warning("Exporting time series of extrapolated real data. This may take some minutes. "
                    "For simulation runs over the same time period, deactivate it.");
        BOOST_OUTCOME_TRY(export_input_data_county_timeseries(
            model, data_dir, node_ids, date, scaling_factor_inf, scaling_factor_icu, num_days,
            path_join(data_dir, "divi_data.json"), path_join(data_dir, "confirmed_cases.json"),
            path_join(data_dir, "population_data.json"), path_join(data_dir, "vaccination_data.json")));
    }

    return success();
}

} // namespace osecirvvs
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // ODESECIRVVS_PARAMETERS_IO_H
