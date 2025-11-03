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
#ifndef MIO_ODE_SECIRVVS_PARAMETERS_IO_H
#define MIO_ODE_SECIRVVS_PARAMETERS_IO_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "ode_secirvvs/model.h"
#include "memilio/geography/regions.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/parameters_io.h"
#include "memilio/io/io.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/date.h"
#include "memilio/utils/logging.h"
#include "memilio/math/math_utils.h"

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
    const std::vector<double>& scaling_factor_inf);
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
IOResult<void> compute_confirmed_cases_data_fix_recovered(const std::vector<ConfirmedCasesDataEntry>& rki_data,
                                                          std::vector<int> const& vregion, Date date,
                                                          std::vector<std::vector<double>>& vnum_rec, double delay = 14.);
IOResult<void> compute_confirmed_cases_data_fix_recovered(std::string const& path, std::vector<int> const& vregion,
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
template <class FP>
IOResult<void> set_confirmed_cases_data(Model<FP>& model, const std::vector<ConfirmedCasesDataEntry>& case_data,
                                        const int region, Date date, const std::vector<double>& scaling_factor_inf, 
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

    std::vector<double> mu_C_R;
    std::vector<double> mu_I_H;
    std::vector<double> mu_H_U;

    std::vector<double> num_InfectedSymptoms(num_age_groups, 0.0);
    std::vector<double> num_death(num_age_groups, 0.0);
    std::vector<double> num_rec(num_age_groups, 0.0);
    std::vector<double> num_Exposed(num_age_groups, 0.0);
    std::vector<double> num_InfectedNoSymptoms(num_age_groups, 0.0);
    std::vector<double> num_InfectedSevere(num_age_groups, 0.0);
    std::vector<double> num_icu(num_age_groups, 0.0);

    /*----------- UNVACCINATED -----------*/
    for (size_t group = 0; group < num_age_groups; group++) {

        t_Exposed.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeExposed<double>>()[(AgeGroup)group])));
        t_InfectedNoSymptoms.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedNoSymptoms<double>>()[(AgeGroup)group])));
        t_InfectedSymptoms.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedSymptoms<double>>()[(AgeGroup)group])));
        t_InfectedSevere.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedSevere<double>>()[(AgeGroup)group])));
        t_InfectedCritical.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedCritical<double>>()[(AgeGroup)group])));

        mu_C_R.push_back(
            model.parameters.template get<RecoveredPerInfectedNoSymptoms<double>>()[(AgeGroup)group]);
        mu_I_H.push_back(
            model.parameters.template get<SeverePerInfectedSymptoms<double>>()[(AgeGroup)group]);
        mu_H_U.push_back(
            model.parameters.template get<CriticalPerSevere<double>>()[(AgeGroup)group]);
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
    if (std::accumulate(num_InfectedSymptoms.begin(), num_InfectedSymptoms.end(), double(0.0),
                        [](const double& a, const double& b) {
                            return evaluate_intermediate<double>(a + b);
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

    num_InfectedSymptoms   = std::vector<double>(num_age_groups, 0.0);
    num_death              = std::vector<double>(num_age_groups, 0.0);
    num_rec                = std::vector<double>(num_age_groups, 0.0);
    num_Exposed            = std::vector<double>(num_age_groups, 0.0);
    num_InfectedNoSymptoms = std::vector<double>(num_age_groups, 0.0);
    num_InfectedSevere     = std::vector<double>(num_age_groups, 0.0);
    num_icu                = std::vector<double>(num_age_groups, 0.0);
    for (size_t group = 0; group < num_age_groups; group++) {
        double reduc_t = model[0].parameters.template get<ReducTimeInfectedMild<double>>()[(AgeGroup)group];
        t_Exposed.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeExposed<double>>()[(AgeGroup)group])));
        t_InfectedNoSymptoms.push_back(static_cast<int>(std::round(
            model.parameters.template get<TimeInfectedNoSymptoms<double>>()[(AgeGroup)group] * reduc_t)));
        t_InfectedSymptoms.push_back(static_cast<int>(std::round(
            model.parameters.template get<TimeInfectedSymptoms<double>>()[(AgeGroup)group] * reduc_t)));
        t_InfectedSevere.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedSevere<double>>()[(AgeGroup)group])));
        t_InfectedCritical.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedCritical<double>>()[(AgeGroup)group])));

        double exp_fac_part_immune =
            model.parameters.template get<ReducExposedPartialImmunity<double>>()[(AgeGroup)group];
        double inf_fac_part_immune =
            model.parameters.template get<ReducInfectedSymptomsPartialImmunity<double>>()[(AgeGroup)group];
        double hosp_fac_part_immune =
            model
                .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[(AgeGroup)group];
        double icu_fac_part_immune =
            model
                .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[(AgeGroup)group];
        mu_C_R.push_back((
            1 - inf_fac_part_immune / exp_fac_part_immune *
                    (1 - model
                                .parameters.template get<RecoveredPerInfectedNoSymptoms<double>>()[(AgeGroup)group])));
        mu_I_H.push_back(
            hosp_fac_part_immune / inf_fac_part_immune *
            model.parameters.template get<SeverePerInfectedSymptoms<double>>()[(AgeGroup)group]);
        // transfer from H to U, D unchanged.
        mu_H_U.push_back(
            icu_fac_part_immune / hosp_fac_part_immune *
            model.parameters.template get<CriticalPerSevere<double>>()[(AgeGroup)group]);
    }

    BOOST_OUTCOME_TRY(compute_confirmed_cases_data(case_data, region, date, num_Exposed, num_InfectedNoSymptoms,
                                                   num_InfectedSymptoms, num_InfectedSevere, num_icu, num_death, num_rec,
                                                   t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
                                                   t_InfectedCritical, mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf));


    size_t num_groups = (size_t)model.parameters.get_num_groups();
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
    if (std::accumulate(num_InfectedSymptoms.begin(), num_InfectedSymptoms.end(), double(0.0),
                        [](const double& a, const double& b) {
                            return evaluate_intermediate<double>(a + b);
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

    num_InfectedSymptoms   = std::vector<double>(num_age_groups, 0.0);
    num_death              = std::vector<double>(num_age_groups, 0.0);
    num_rec                = std::vector<double>(num_age_groups, 0.0);
    num_Exposed            = std::vector<double>(num_age_groups, 0.0);
    num_InfectedNoSymptoms = std::vector<double>(num_age_groups, 0.0);
    num_InfectedSevere     = std::vector<double>(num_age_groups, 0.0);
    num_icu                = std::vector<double>(num_age_groups, 0.0);
    for (size_t group = 0; group < num_age_groups; group++) {
        double reduc_t = model[0].parameters.template get<ReducTimeInfectedMild<double>>()[(AgeGroup)group];
        t_Exposed.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeExposed<double>>()[(AgeGroup)group])));
        t_InfectedNoSymptoms.push_back(static_cast<int>(std::round(
            model.parameters.template get<TimeInfectedNoSymptoms<double>>()[(AgeGroup)group] * reduc_t)));
        t_InfectedSymptoms.push_back(static_cast<int>(std::round(
            model.parameters.template get<TimeInfectedSymptoms<double>>()[(AgeGroup)group] * reduc_t)));
        t_InfectedSevere.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedSevere<double>>()[(AgeGroup)group])));
        t_InfectedCritical.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedCritical<double>>()[(AgeGroup)group])));

        double reduc_immune_exp =
            model.parameters.template get<ReducExposedImprovedImmunity<double>>()[(AgeGroup)group];
        double reduc_immune_inf =
            model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<double>>()[(AgeGroup)group];
        double reduc_immune_hosp =
            model.parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[(
                AgeGroup)group];
        double reduc_immune_icu =
            model.parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[(
                AgeGroup)group];
        mu_C_R.push_back((
            1 - reduc_immune_inf / reduc_immune_exp *
                    (1 - model
                                .parameters.template get<RecoveredPerInfectedNoSymptoms<double>>()[(AgeGroup)group])));
        mu_I_H.push_back(
            reduc_immune_hosp / reduc_immune_inf *
            model.parameters.template get<SeverePerInfectedSymptoms<double>>()[(AgeGroup)group]);
        // transfer from H to U, D unchanged.
        mu_H_U.push_back(
            reduc_immune_icu / reduc_immune_hosp *
            model.parameters.template get<CriticalPerSevere<double>>()[(AgeGroup)group]);
    }

    BOOST_OUTCOME_TRY(compute_confirmed_cases_data(case_data, region, date, num_Exposed, num_InfectedNoSymptoms,
                                                   num_InfectedSymptoms, num_InfectedSevere, num_icu, num_death, num_rec,
                                                   t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
                                                   t_InfectedCritical, mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf));

    for (size_t county = 0; county < model.size(); county++) {
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
            // Only set the number of ICU patients here, if the date is not available in the data.
            if (!is_divi_data_available(date)) {
                model[county].populations[{AgeGroup(i), InfectionState::InfectedCriticalImprovedImmunity}] =
                    num_icu[county][i];
            }
        }

        if (std::accumulate(num_InfectedSymptoms[county].begin(), num_InfectedSymptoms[county].end(), double(0.0),
                            [](const double& a, const double& b) {
                                return evaluate_intermediate<double>(a + b);
                            }) == 0.0) {
            log_warning("No infections for vaccinated reported on date {} for region {}. "
                        "Population data has not been set.",
                        date, region[county]);
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
template <class FP>
IOResult<void> set_confirmed_cases_data(mio::VectorRange<Node<Model<FP>>>& model, const std::string& path,
                                        std::vector<int> const& vregion, Date date,
                                        const std::vector<double>& scaling_factor_inf, bool set_death = false)
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
        BOOST_OUTCOME_TRY(set_confirmed_cases_data(model[region_idx].property, vcase_data[region_idx], vregion[region_idx], date, scaling_factor_inf, set_death));
    }
    return success();
}

/**
* @brief sets population data from census data which has been read into num_population
* @param[in, out] model vector of objects in which the data is set
* @param[in] num_population vector of population data
* @param[in] vregion vector of keys of the regions of interest
* @param[in] date Date for which the arrays are initialized
*/
template <class FP>
IOResult<void> set_population_data(Model<FP>& model, const std::vector<double>& num_population,
                                   const int region, const std::vector<double>& num_rec)
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

        double S_v = std::min(
            model.parameters.template get<DailyFullVaccinations<double>>()[{i, SimulationDay(0)}] +
                num_rec[size_t(i)],
            num_population[size_t(i)]);
        double S_pv = std::max(
            model.parameters.template get<DailyPartialVaccinations<double>>()[{i, SimulationDay(0)}] -
                model.parameters.template get<DailyFullVaccinations<double>>()[{i, SimulationDay(0)}],
            0.0); // use std::max with 0
        double S;
        if (num_population[size_t(i)] - S_pv - S_v < 0.0) {
            log_warning("Number of vaccinated persons greater than population in county {}, age group {}.",
                        region, size_t(i));
            S   = 0.0;
            S_v = num_population[size_t(i)] - S_pv;
        }
        else {
            S = num_population[size_t(i)] - S_pv - S_v;
        }

        double denom_E =
            1 / (S + S_pv * model.parameters.template get<ReducExposedPartialImmunity<double>>()[i] +
                    S_v * model.parameters.template get<ReducExposedImprovedImmunity<double>>()[i]);
        double denom_C =
            1 / (S + S_pv * model.parameters.template get<ReducExposedPartialImmunity<double>>()[i] +
                    S_v * model.parameters.template get<ReducExposedImprovedImmunity<double>>()[i]);
        double denom_I =
            1 /
            (S +
                S_pv * model.parameters.template get<ReducInfectedSymptomsPartialImmunity<double>>()[i] +
                S_v * model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<double>>()[i]);
        double denom_HU =
            1 /
            (S +
                S_pv * model
                        .parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[i] +
                S_v * model
                        .parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[i]);

        model.populations[{i, InfectionState::ExposedNaive}] =
            S * model.populations[{i, InfectionState::ExposedNaive}] * denom_E;
        model.populations[{i, InfectionState::ExposedPartialImmunity}] =
            S_pv * model.parameters.template get<ReducExposedPartialImmunity<double>>()[i] *
            model.populations[{i, InfectionState::ExposedPartialImmunity}] * denom_E;
        model.populations[{i, InfectionState::ExposedImprovedImmunity}] =
            S_v * model.parameters.template get<ReducExposedImprovedImmunity<double>>()[i] *
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
            S_pv * model.parameters.template get<ReducInfectedSymptomsPartialImmunity<double>>()[i] *
            model.populations[{i, InfectionState::InfectedSymptomsPartialImmunity}] * denom_I;
        model.populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] =
            S_v * model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<double>>()[i] *
            model.populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] * denom_I;

        model.populations[{i, InfectionState::InfectedSymptomsNaiveConfirmed}] =
            S * model.populations[{i, InfectionState::InfectedSymptomsNaiveConfirmed}] * denom_I;
        model.populations[{i, InfectionState::InfectedSymptomsPartialImmunityConfirmed}] =
            S_pv * model.parameters.template get<ReducInfectedSymptomsPartialImmunity<double>>()[i] *
            model.populations[{i, InfectionState::InfectedSymptomsPartialImmunityConfirmed}] * denom_I;
        model.populations[{i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] =
            S_v * model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<double>>()[i] *
            model.populations[{i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] * denom_I;

        model.populations[{i, InfectionState::InfectedSevereNaive}] =
            S * model.populations[{i, InfectionState::InfectedSevereNaive}] * denom_HU;
        model.populations[{i, InfectionState::InfectedSeverePartialImmunity}] =
            S_pv *
            model.parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[i] *
            model.populations[{i, InfectionState::InfectedSeverePartialImmunity}] * denom_HU;
        model.populations[{i, InfectionState::InfectedSevereImprovedImmunity}] =
            S_v *
            model
                .parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[i] *
            model.populations[{i, InfectionState::InfectedSevereImprovedImmunity}] * denom_HU;

        model.populations[{i, InfectionState::InfectedCriticalPartialImmunity}] =
            S_pv *
            model.parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[i] *
            model.populations[{i, InfectionState::InfectedCriticalNaive}] * denom_HU;
        model.populations[{i, InfectionState::InfectedCriticalImprovedImmunity}] =
            S_v *
            model
                .parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[i] *
            model.populations[{i, InfectionState::InfectedCriticalNaive}] * denom_HU;
        model.populations[{i, InfectionState::InfectedCriticalNaive}] =
            S * model.populations[{i, InfectionState::InfectedCriticalNaive}] * denom_HU;

        model.populations[{i, InfectionState::SusceptibleImprovedImmunity}] =
            model.parameters.template get<DailyFullVaccinations<double>>()[{i, SimulationDay(0)}] +
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
            std::max(0.0, double(model.populations[{i, InfectionState::SusceptibleImprovedImmunity}])));

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

/**
* @brief Sets population data from census data which has been read into num_population.
* @param[in, out] model Vector of objects in which the data is set.
* @param[in] path Path to population data file.
* @param[in] path_rki Path to RKI cases data file.
* @param[in] vregion Vector of keys of the regions of interest.
* @param[in] date Date for which the arrays are initialized.
*/
template <class FP>
IOResult<void> set_population_data(mio::VectorRange<Node<Model<FP>>>& model, const std::string& path, const std::string& path_rki,
                                   const std::vector<int>& vregion, Date date)
{
    BOOST_OUTCOME_TRY(auto&& vnum_population, read_population_data(path, vregion));
    BOOST_OUTCOME_TRY(auto&& case_data, mio::read_confirmed_cases_data(path_rki));

    auto num_age_groups = ConfirmedCasesDataEntry::age_group_names.size();
    std::vector<std::vector<double>> vnum_rec(model.size(), std::vector<double>(num_age_groups, 0.0));

    BOOST_OUTCOME_TRY(compute_confirmed_cases_data_fix_recovered(case_data, vregion, date, vnum_rec, 14.));

    for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
        BOOST_OUTCOME_TRY(set_population_data(model[region_idx].property, vnum_population[region_idx], vregion[region_idx], vnum_rec[region_idx]));
    }
    return success();
}

/**
 * @brief Sets vaccination data for models stored in a vector.
 *
 * @tparam FP Floating point type used in the Model objects.
 * @param[in, out] model Vector of Model objects in which the vaccination data is set.
 * @param[in] vacc_data Vector of VaccinationDataEntry objects containing the vaccination data.
 * @param[in] date Start date for the simulation.
 * @param[in] num_days Number of days for which the simulation is run.
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
    auto max_date = max_date_entry->date;

    auto min_date_entry = std::min_element(vacc_data.begin(), vacc_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    auto min_date       = min_date_entry->date;
    if (min_date > date || max_date < offset_date_by_days(date, num_days)) {
        log_warning("Vaccination data only available from {} to {}. "
                    "For days before, vaccination data will be set to 0. For days after, "
                    "vaccination data will be set to the last available date.",
                    min_date, max_date);
    }

    auto num_groups = model.parameters.get_num_groups();

    // type conversion from UncertainValue -> FP -> int
    auto days_until_effective1 = static_cast<int>(
        static_cast<FP>(model.parameters.template get<DaysUntilEffectivePartialImmunity<FP>>()[AgeGroup(0)]));
    auto days_until_effective2 = static_cast<int>(
        static_cast<FP>(model.parameters.template get<DaysUntilEffectiveImprovedImmunity<FP>>()[AgeGroup(0)]));
    auto vaccination_distance =
        static_cast<int>(static_cast<FP>(model.parameters.template get<VaccinationGap<FP>>()[AgeGroup(0)]));

    for (auto&& vacc_data_entry : vacc_data) {
        auto date_df = vacc_data_entry.date;
        AgeGroup age        = vacc_data_entry.age_group;

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
            // Storing its value in get<DailyPartialVaccinations>() will eventually (in the simulation)
            // transfer the difference (between get<DailyPartialVaccinations>() at d and d-1) of
            // N susceptible individuals to 'Susceptible Partially Vaccinated' state at day d; see secir_vaccinated.h
            auto offset_first_date =
                offset_date_by_days(date, (int)d - days_until_effective1 + vaccination_distance);
            if (max_date >= offset_first_date) {
                // Option 1: considered offset_first_date is available in input data frame
                if (date_df == offset_first_date) {
                    model.parameters.template get<DailyPartialVaccinations<FP>>()[{age, SimulationDay(d)}] =
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
                    model.parameters.template get<DailyPartialVaccinations<FP>>()[{age, SimulationDay(d)}] -=
                        days_plus * vacc_data_entry.num_vaccinations_completed;
                }
                else if (date_df == max_date) {
                    model.parameters.template get<DailyPartialVaccinations<FP>>()[{age, SimulationDay(d)}] +=
                        (days_plus + 1) * vacc_data_entry.num_vaccinations_completed;
                }
            }

            // a person whose second dose is reported at start_date + simulation_day - days_until_effective2
            // has the full protection of the second dose at day X = start_date + simulation_day
            // Storing its value in get<DailyFullVaccinations>() will eventually (in the simulation)
            // transfer the difference (between get<DailyFullVaccinations>() at d and d-1) of
            // N susceptible, partially vaccinated individuals to 'SusceptibleImprovedImmunity' state at day d; see secir_vaccinated.h
            auto offset_full_date = offset_date_by_days(date, (int)d - days_until_effective2);
            if (max_date >= offset_full_date) {
                // Option 1: considered offset_full_date is available in input data frame
                if (date_df == offset_full_date) {
                    model.parameters.template get<DailyFullVaccinations<FP>>()[{age, SimulationDay(d)}] =
                        vacc_data_entry.num_vaccinations_completed;
                }
            }
            else { // offset_full_date > max_full_date
                // Option 2: considered offset_full_date is NOT available in input data frame
                days_plus = get_offset_in_days(offset_full_date, max_date);
                if (date_df == offset_date_by_days(max_date, -1)) {
                    model.parameters.template get<DailyFullVaccinations<FP>>()[{age, SimulationDay(d)}] -=
                        days_plus * vacc_data_entry.num_vaccinations_completed;
                }
                else if (date_df == max_date) {
                    model.parameters.template get<DailyFullVaccinations<FP>>()[{age, SimulationDay(d)}] +=
                        (days_plus + 1) * vacc_data_entry.num_vaccinations_completed;
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
template <typename FP>
IOResult<void> set_vaccination_data(mio::VectorRange<Node<Model<FP>>>& model, const std::string& path, Date date,
                                    const std::vector<int>& vregion, int num_days)
{
    // Set vaccination data to 0 for all models
    for (auto& m : model) {
        m.property.parameters.template get<DailyPartialVaccinations<FP>>().resize(SimulationDay(num_days + 1));
        m.property.parameters.template get<DailyFullVaccinations<FP>>().resize(SimulationDay(num_days + 1));
        for (auto d = SimulationDay(0); d < SimulationDay(num_days + 1); ++d) {
            for (auto a = AgeGroup(0); a < m.property.parameters.get_num_groups(); ++a) {
                m.property.parameters.template get<DailyPartialVaccinations<FP>>()[{a, d}] = 0.0;
                m.property.parameters.template get<DailyFullVaccinations<FP>>()[{a, d}]    = 0.0;
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
            return r == 0 || (vacc_data_entry.county_id && vacc_data_entry.county_id == regions::de::CountyId(r)) ||
                   (vacc_data_entry.state_id && vacc_data_entry.state_id == regions::de::StateId(r)) ||
                   (vacc_data_entry.district_id && vacc_data_entry.district_id == regions::de::DistrictId(r));
        });
        if (it != vregion.end()) {
            auto region_idx = size_t(it - vregion.begin());
            vvacc_data[region_idx].pushback(vacc_data_entry);
        }
    }

    for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
        BOOST_OUTCOME_TRY(set_vaccination_data(model[region_idx].property, vacc_data[region_idx], date, num_days));
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
IOResult<void> set_divi_data(mio::VectorRange<Node<Model<FP>>>& model, const std::string& path, const std::vector<int>& vregion,
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
        BOOST_OUTCOME_TRY(set_divi_data(model[region_idx].property, num_icu[region_idx], vregion[region_idx], scaling_factor_icu));
    }

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
* @param[in] date Date for which the data should be read.
* @param[in] node_ids Vector of keys of the node_ids of interest.
* @param[in] scaling_factor_inf Factors by which to scale the confirmed cases of rki data.
* @param[in] scaling_factor_icu Factor by which to scale the icu cases of divi data.
* @param[in] num_days Number of days to be simulated/initialized.
* @param[in] pydata_dir Directory that contains the data files.
*/
template <class FP>
IOResult<void> export_input_data_timeseries(
    mio::VectorRange<Node<Model<FP>>> models, const std::string& results_dir, Date date, const std::vector<int>& node_ids,
    const std::vector<double>& scaling_factor_inf, const double scaling_factor_icu, const int num_days,
    const mio::regions::de::EpidataFilenames& epidata_filenames)
{
    const auto num_age_groups = (size_t)models[0].property.parameters.get_num_groups();
    assert(scaling_factor_inf.size() == num_age_groups);
    assert(num_age_groups == ConfirmedCasesDataEntry::age_group_names.size());
    assert(models.size() == node_ids.size());
    std::vector<TimeSeries<double>> extrapolated_data(
        models.size(), TimeSeries<double>::zero(num_days + 1, (size_t)InfectionState::Count * num_age_groups));

    for (int t = 0; t <= num_days; ++t) {
        auto offset_day = offset_date_by_days(date, t);

        // TODO: empty vaccination data path guard
        BOOST_OUTCOME_TRY(read_input_data(model, date, county, scaling_factor_inf, scaling_factor_icu,
                                      num_days, epidata_filenames));

        for (size_t r = 0; r < node_ids.size(); r++) {
            extrapolated_data[r][t] = models[r].property.get_initial_values();
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
template <class FP>
IOResult<void> export_input_data_county_timeseries(mio::VectorRange<Node<Model<FP>>>, const std::string&, Date, const std::vector<int>&,
                                                   const std::vector<double>&, const double, const int,
                                                   const mio::regions::de::EpidataFilenames&)
{
    mio::log_warning("HDF5 not available. Cannot export time series of extrapolated real data.");
    return success();
}

#endif //MEMILIO_HAS_HDF5

template <typename FP>
IOResult<void> read_input_data(mio::VectorRange<Node<Model<FP>>>& model, Date date, const std::vector<int>& node_ids,
                               const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                               int num_days, const mio::regions::de::EpidataFilenames& epidata_filenames)
{

    BOOST_OUTCOME_TRY(
        details::set_vaccination_data(model, epidata_filenames.vaccination_data_path, date, node_ids, num_days));

    // TODO: Reuse more code, e.g., set_divi_data (in secir) and a set_divi_data (here) only need a different ModelType.
    // TODO: add option to set ICU data from confirmed cases if DIVI or other data is not available.
    BOOST_OUTCOME_TRY(details::set_divi_data(model, epidata_filenames.divi_data_path, node_ids, date,
                                             scaling_factor_icu));

    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(model, epidata_filenames.case_data_path, node_ids,
                                                        date, scaling_factor_inf));
    BOOST_OUTCOME_TRY(details::set_population_data(model, epidata_filenames.population_data_path, 
                                                   epidata_filenames.case_data_path, node_ids, date));
    return success();
}


} // namespace osecirvvs
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // MIO_ODE_SECIRVVS_PARAMETERS_IO_H
