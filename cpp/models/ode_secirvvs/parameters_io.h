/* 
* Copyright (C) 2020-2022 German Aerospace Center (DLR-SC)
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
#include "ode_secirvvs/analyze_result.h"
#include "memilio/math/eigen_util.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/mobility.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"
#include "memilio/io/json_serializer.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/date.h"

namespace mio
{
namespace osecirvvs
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
    std::string const& path, std::vector<int> const& vregion, Date date, std::vector<std::vector<double>>& num_exp,
    std::vector<std::vector<double>>& num_car, std::vector<std::vector<double>>& num_inf,
    std::vector<std::vector<double>>& num_hosp, std::vector<std::vector<double>>& num_icu,
    std::vector<std::vector<double>>& num_death, std::vector<std::vector<double>>& num_rec,
    const std::vector<std::vector<int>>& t_car_to_rec, const std::vector<std::vector<int>>& t_car_to_inf,
    const std::vector<std::vector<int>>& t_exp_to_car, const std::vector<std::vector<int>>& t_inf_to_rec,
    const std::vector<std::vector<int>>& t_inf_to_hosp, const std::vector<std::vector<int>>& t_hosp_to_rec,
    const std::vector<std::vector<int>>& t_hosp_to_icu, const std::vector<std::vector<int>>& t_icu_to_dead,
    const std::vector<std::vector<int>>& t_icu_to_rec, const std::vector<std::vector<double>>& mu_C_R,
    const std::vector<std::vector<double>>& mu_I_H, const std::vector<std::vector<double>>& mu_H_U,
    const std::vector<std::vector<double>>& mu_U_D, const std::vector<double>& scaling_factor_inf);

IOResult<void> read_confirmed_cases_data(
    const std::vector<ConfirmedCasesDataEntry>& rki_data, std::vector<int> const& vregion, Date date,
    std::vector<std::vector<double>>& num_exp, std::vector<std::vector<double>>& num_car,
    std::vector<std::vector<double>>& num_inf, std::vector<std::vector<double>>& num_hosp,
    std::vector<std::vector<double>>& num_icu, std::vector<std::vector<double>>& num_death,
    std::vector<std::vector<double>>& num_rec, const std::vector<std::vector<int>>& t_car_to_rec,
    const std::vector<std::vector<int>>& t_car_to_inf, const std::vector<std::vector<int>>& t_exp_to_car,
    const std::vector<std::vector<int>>& t_inf_to_rec, const std::vector<std::vector<int>>& t_inf_to_hosp,
    const std::vector<std::vector<int>>& t_hosp_to_rec, const std::vector<std::vector<int>>& t_hosp_to_icu,
    const std::vector<std::vector<int>>& t_icu_to_dead, const std::vector<std::vector<int>>& t_icu_to_rec,
    const std::vector<std::vector<double>>& mu_C_R, const std::vector<std::vector<double>>& mu_I_H,
    const std::vector<std::vector<double>>& mu_H_U, const std::vector<std::vector<double>>& mu_U_D,
    const std::vector<double>& scaling_factor_inf);
/**@}*/

/**
        * @brief Reads confirmed cases data and translates data of day t0-delay to recovered compartment,
        * @param path Path to RKI confirmed cases file.
        * @param vregion vector of keys of the region of interest     
        * @param date Date for which the arrays are initialized
        * @param num_rec output vector for number of people in the compartement recovered
        * @param delay number of days in the past the are used to set recovered compartment.
        * @see mio::read_rki_data
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
        * @brief sets populations data from a transformed RKI cases file into a Model.
        * @param model vector of objects in which the data is set
        * @param path Path to transformed RKI cases file
        * @param region vector of keys of the region of interest
        * @param date Date for which the arrays are initialized
        * @param scaling_factor_inf factors by which to scale the confirmed cases of
        * rki data
        */
template <class Model>
IOResult<void> set_confirmed_cases_data(std::vector<Model>& model, const std::string& path,
                                        std::vector<int> const& region, Date date,
                                        const std::vector<double>& scaling_factor_inf)
{
    auto num_age_groups = (size_t)model[0].parameters.get_num_groups();
    assert(scaling_factor_inf.size() == num_age_groups); //TODO: allow vector or scalar valued scaling factors
    assert(ConfirmedCasesDataEntry::age_group_names.size() == num_age_groups);

    BOOST_OUTCOME_TRY(rki_data, mio::read_confirmed_cases_data(path));

    std::vector<std::vector<int>> t_car_to_rec{model.size()}; // R9
    std::vector<std::vector<int>> t_car_to_inf{model.size()}; // R3
    std::vector<std::vector<int>> t_exp_to_car{model.size()}; // R2
    std::vector<std::vector<int>> t_inf_to_rec{model.size()}; // R4
    std::vector<std::vector<int>> t_inf_to_hosp{model.size()}; // R6
    std::vector<std::vector<int>> t_hosp_to_rec{model.size()}; // R5
    std::vector<std::vector<int>> t_hosp_to_icu{model.size()}; // R7
    std::vector<std::vector<int>> t_icu_to_rec{model.size()};

    std::vector<std::vector<double>> mu_C_R{model.size()};
    std::vector<std::vector<double>> mu_I_H{model.size()};
    std::vector<std::vector<double>> mu_H_U{model.size()};
    std::vector<std::vector<double>> mu_U_D{model.size()};

    std::vector<std::vector<double>> num_inf(model.size());
    std::vector<std::vector<double>> num_death(model.size());
    std::vector<std::vector<double>> num_rec(model.size());
    std::vector<std::vector<double>> num_exp(model.size());
    std::vector<std::vector<double>> num_car(model.size());
    std::vector<std::vector<double>> num_hosp(model.size());
    std::vector<std::vector<double>> num_icu(model.size());

    /*----------- UNVACCINATED -----------*/
    for (size_t county = 0; county < model.size(); county++) {
        num_inf[county]   = std::vector<double>(num_age_groups, 0.0);
        num_death[county] = std::vector<double>(num_age_groups, 0.0);
        num_rec[county]   = std::vector<double>(num_age_groups, 0.0);
        num_exp[county]   = std::vector<double>(num_age_groups, 0.0);
        num_car[county]   = std::vector<double>(num_age_groups, 0.0);
        num_hosp[county]  = std::vector<double>(num_age_groups, 0.0);
        num_icu[county]   = std::vector<double>(num_age_groups, 0.0);
        for (size_t group = 0; group < num_age_groups; group++) {

            t_car_to_inf[county].push_back(
                static_cast<int>(2 * (model[county].parameters.template get<IncubationTime>()[(AgeGroup)group] -
                                      model[county].parameters.template get<SerialInterval>()[(AgeGroup)group])));
            t_car_to_rec[county].push_back(
                static_cast<int>(t_car_to_inf[county][group] +
                                 0.5 * model[county].parameters.template get<TimeInfectedSymptoms>()[(AgeGroup)group]));
            t_exp_to_car[county].push_back(
                static_cast<int>(2 * model[county].parameters.template get<SerialInterval>()[(AgeGroup)group] -
                                 model[county].parameters.template get<IncubationTime>()[(AgeGroup)group]));
            t_inf_to_rec[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedSymptoms>()[(AgeGroup)group]));
            t_hosp_to_rec[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedSevere>()[(AgeGroup)group]));
            t_hosp_to_icu[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedSevere>()[(AgeGroup)group]));
            t_icu_to_dead[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedCritical>()[(AgeGroup)group]));
            t_icu_to_rec[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedCritical>()[(AgeGroup)group]));

            mu_C_R[county].push_back(
                model[county].parameters.template get<AsymptoticCasesPerInfectious>()[(AgeGroup)group]);
            mu_I_H[county].push_back(
                model[county].parameters.template get<HospitalizedCasesPerInfectious>()[(AgeGroup)group]);
            mu_H_U[county].push_back(model[county].parameters.template get<ICUCasesPerHospitalized>()[(AgeGroup)group]);
            mu_U_D[county].push_back(model[county].parameters.template get<DeathsPerICU>()[(AgeGroup)group]);
        }
    }

    BOOST_OUTCOME_TRY(read_confirmed_cases_data(
        rki_data, region, date, num_exp, num_car, num_inf, num_hosp, num_icu, num_death, num_rec, t_car_to_rec,
        t_car_to_inf, t_exp_to_car, t_inf_to_rec, t_inf_to_hosp, t_hosp_to_rec, t_hosp_to_icu, t_icu_to_dead,
        t_icu_to_rec, mu_C_R, mu_I_H, mu_H_U, mu_U_D, scaling_factor_inf));

    for (size_t county = 0; county < model.size(); county++) {
        // if (std::accumulate(num_inf[county].begin(), num_inf[county].end(), 0.0) > 0) {
        size_t num_groups = (size_t)model[county].parameters.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            model[county].populations[{AgeGroup(i), InfectionState::ExposedNaive}]      = num_exp[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsNaive}]      = num_car[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsNaive}]     = num_inf[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSevereNaive}] = num_hosp[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::SusceptibleImprovedImmunity}]         = num_rec[county][i];
        }
        // }
        if (std::accumulate(num_inf[county].begin(), num_inf[county].end(), 0.0) == 0) {
            log_warning("No infections for unvaccinated reported on date " + std::to_string(date.year) + "-" +
                        std::to_string(date.month) + "-" + std::to_string(date.day) + " for region " +
                        std::to_string(region[county]) + ". Population data has not been set.");
        }
    }

    /*----------- PARTIALLY VACCINATED -----------*/
    for (size_t county = 0; county < model.size(); county++) {
        t_car_to_rec[county].clear();
        t_car_to_inf[county].clear();
        t_exp_to_car[county].clear();
        t_inf_to_rec[county].clear();
        t_inf_to_hosp[county].clear();
        t_hosp_to_rec[county].clear();
        t_hosp_to_icu[county].clear();
        t_icu_to_dead[county].clear();
        t_icu_to_rec[county].clear();

        mu_C_R[county].clear();
        mu_I_H[county].clear();
        mu_H_U[county].clear();
        mu_U_D[county].clear();

        num_inf[county]   = std::vector<double>(num_age_groups, 0.0);
        num_death[county] = std::vector<double>(num_age_groups, 0.0);
        num_rec[county]   = std::vector<double>(num_age_groups, 0.0);
        num_exp[county]   = std::vector<double>(num_age_groups, 0.0);
        num_car[county]   = std::vector<double>(num_age_groups, 0.0);
        num_hosp[county]  = std::vector<double>(num_age_groups, 0.0);
        num_icu[county]   = std::vector<double>(num_age_groups, 0.0);
        for (size_t group = 0; group < num_age_groups; group++) {

            double reduc_t = model[0].parameters.template get<InfectiousTimeFactorImmune>()[(AgeGroup)group];
            t_car_to_inf[county].push_back(
                static_cast<int>(2 * (model[county].parameters.template get<IncubationTime>()[(AgeGroup)group] -
                                      model[county].parameters.template get<SerialInterval>()[(AgeGroup)group])));
            t_car_to_rec[county].push_back(static_cast<int>(
                t_car_to_inf[county][group] +
                0.5 * model[county].parameters.template get<TimeInfectedSymptoms>()[(AgeGroup)group] * reduc_t));
            t_exp_to_car[county].push_back(
                static_cast<int>(2 * model[county].parameters.template get<SerialInterval>()[(AgeGroup)group] -
                                 model[county].parameters.template get<IncubationTime>()[(AgeGroup)group]));
            t_inf_to_rec[county].push_back(static_cast<int>(
                model[county].parameters.template get<TimeInfectedSymptoms>()[(AgeGroup)group] * reduc_t));
            t_hosp_to_rec[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedSevere>()[(AgeGroup)group]));
            t_hosp_to_icu[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedSevere>()[(AgeGroup)group]));
            t_icu_to_dead[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedCritical>()[(AgeGroup)group]));
            t_icu_to_rec[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedCritical>()[(AgeGroup)group]));

            double exp_fac_part_immune =
                model[county].parameters.template get<ExposedFactorPartialImmunity>()[(AgeGroup)group];
            double inf_fac_part_immune =
                model[county].parameters.template get<InfectedFactorPartialImmunity>()[(AgeGroup)group];
            double hosp_fac_part_immune =
                model[county].parameters.template get<HospitalizedFactorPartialImmunity>()[(AgeGroup)group];
            double icu_fac_part_immune =
                model[county].parameters.template get<HospitalizedFactorPartialImmunity>()[(AgeGroup)group];
            double death_fac_part_immune =
                model[county].parameters.template get<HospitalizedFactorPartialImmunity>()[(AgeGroup)group];
            mu_C_R[county].push_back(
                (1 - inf_fac_part_immune / exp_fac_part_immune *
                         (1 - model[county].parameters.template get<AsymptoticCasesPerInfectious>()[(AgeGroup)group])));
            mu_I_H[county].push_back(
                hosp_fac_part_immune / inf_fac_part_immune *
                model[county].parameters.template get<HospitalizedCasesPerInfectious>()[(AgeGroup)group]);
            // transfer from H to U, D unchanged.
            mu_H_U[county].push_back(icu_fac_part_immune / hosp_fac_part_immune *
                                     model[county].parameters.template get<ICUCasesPerHospitalized>()[(AgeGroup)group]);
            mu_U_D[county].push_back(death_fac_part_immune / icu_fac_part_immune *
                                     model[county].parameters.template get<DeathsPerICU>()[(AgeGroup)group]);
        }
    }

    BOOST_OUTCOME_TRY(read_confirmed_cases_data(
        rki_data, region, date, num_exp, num_car, num_inf, num_hosp, num_icu, num_death, num_rec, t_car_to_rec,
        t_car_to_inf, t_exp_to_car, t_inf_to_rec, t_inf_to_hosp, t_hosp_to_rec, t_hosp_to_icu, t_icu_to_dead,
        t_icu_to_rec, mu_C_R, mu_I_H, mu_H_U, mu_U_D, scaling_factor_inf));

    for (size_t county = 0; county < model.size(); county++) {
        // if (std::accumulate(num_inf[county].begin(), num_inf[county].end(), 0.0) > 0) {
        size_t num_groups = (size_t)model[county].parameters.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            model[county].populations[{AgeGroup(i), InfectionState::ExposedPartialImmunity}]      = num_exp[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsPartialImmunity}]      = num_car[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsPartialImmunity}]     = num_inf[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSeverePartialImmunity}] = num_hosp[county][i];
        }
        // }
        if (std::accumulate(num_inf[county].begin(), num_inf[county].end(), 0.0) == 0) {
            log_warning("No infections for partially vaccinated reported on date " + std::to_string(date.year) + "-" +
                        std::to_string(date.month) + "-" + std::to_string(date.day) + " for region " +
                        std::to_string(region[county]) + ". Population data has not been set.");
        }
    }

    /*----------- FULLY VACCINATED -----------*/
    for (size_t county = 0; county < model.size(); county++) {
        t_car_to_rec[county].clear();
        t_car_to_inf[county].clear();
        t_exp_to_car[county].clear();
        t_inf_to_rec[county].clear();
        t_inf_to_hosp[county].clear();
        t_hosp_to_rec[county].clear();
        t_hosp_to_icu[county].clear();
        t_icu_to_dead[county].clear();
        t_icu_to_rec[county].clear();

        mu_C_R[county].clear();
        mu_I_H[county].clear();
        mu_H_U[county].clear();
        mu_U_D[county].clear();

        num_inf[county]   = std::vector<double>(num_age_groups, 0.0);
        num_death[county] = std::vector<double>(num_age_groups, 0.0);
        num_rec[county]   = std::vector<double>(num_age_groups, 0.0);
        num_exp[county]   = std::vector<double>(num_age_groups, 0.0);
        num_car[county]   = std::vector<double>(num_age_groups, 0.0);
        num_hosp[county]  = std::vector<double>(num_age_groups, 0.0);
        num_icu[county]   = std::vector<double>(num_age_groups, 0.0);
        for (size_t group = 0; group < num_age_groups; group++) {

            double reduc_t = model[0].parameters.template get<InfectiousTimeFactorImmune>()[(AgeGroup)group];
            t_car_to_inf[county].push_back(
                static_cast<int>(2 * (model[county].parameters.template get<IncubationTime>()[(AgeGroup)group] -
                                      model[county].parameters.template get<SerialInterval>()[(AgeGroup)group])));
            t_car_to_rec[county].push_back(static_cast<int>(
                t_car_to_inf[county][group] +
                0.5 * model[county].parameters.template get<TimeInfectedSymptoms>()[(AgeGroup)group] * reduc_t));
            t_exp_to_car[county].push_back(
                static_cast<int>(2 * model[county].parameters.template get<SerialInterval>()[(AgeGroup)group] -
                                 model[county].parameters.template get<IncubationTime>()[(AgeGroup)group]));
            t_inf_to_rec[county].push_back(static_cast<int>(
                model[county].parameters.template get<TimeInfectedSymptoms>()[(AgeGroup)group] * reduc_t));
            t_hosp_to_rec[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedSevere>()[(AgeGroup)group]));
            t_hosp_to_icu[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedSevere>()[(AgeGroup)group]));
            t_icu_to_dead[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedCritical>()[(AgeGroup)group]));
            t_icu_to_rec[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedCritical>()[(AgeGroup)group]));

            double reduc_immune_exp =
                model[county].parameters.template get<ExposedFactorImprovedImmunity>()[(AgeGroup)group];
            double reduc_immune_inf =
                model[county].parameters.template get<InfectedFactorImprovedImmunity>()[(AgeGroup)group];
            double reduc_immune_hosp =
                model[county].parameters.template get<HospitalizedFactorImprovedImmunity>()[(AgeGroup)group];
            double reduc_immune_icu =
                model[county].parameters.template get<HospitalizedFactorImprovedImmunity>()[(AgeGroup)group];
            double reduc_immune_dead =
                model[county].parameters.template get<HospitalizedFactorImprovedImmunity>()[(AgeGroup)group];
            mu_C_R[county].push_back(
                (1 - reduc_immune_inf / reduc_immune_exp *
                         (1 - model[county].parameters.template get<AsymptoticCasesPerInfectious>()[(AgeGroup)group])));
            mu_I_H[county].push_back(
                reduc_immune_hosp / reduc_immune_inf *
                model[county].parameters.template get<HospitalizedCasesPerInfectious>()[(AgeGroup)group]);
            // transfer from H to U, D unchanged.
            mu_H_U[county].push_back(reduc_immune_icu / reduc_immune_hosp *
                                     model[county].parameters.template get<ICUCasesPerHospitalized>()[(AgeGroup)group]);
            mu_U_D[county].push_back(reduc_immune_dead / reduc_immune_icu *
                                     model[county].parameters.template get<DeathsPerICU>()[(AgeGroup)group]);
        }
    }

    BOOST_OUTCOME_TRY(read_confirmed_cases_data(
        rki_data, region, date, num_exp, num_car, num_inf, num_hosp, num_icu, num_death, num_rec, t_car_to_rec,
        t_car_to_inf, t_exp_to_car, t_inf_to_rec, t_inf_to_hosp, t_hosp_to_rec, t_hosp_to_icu, t_icu_to_dead,
        t_icu_to_rec, mu_C_R, mu_I_H, mu_H_U, mu_U_D, scaling_factor_inf));

    for (size_t county = 0; county < model.size(); county++) {
        // if (std::accumulate(num_inf[county].begin(), num_inf[county].end(), 0.0) > 0) {
        size_t num_groups = (size_t)model[county].parameters.get_num_groups();
        for (size_t i = 0; i < num_groups; i++) {
            model[county].populations[{AgeGroup(i), InfectionState::ExposedImprovedImmunity}]  = num_exp[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsImprovedImmunity}]  = num_car[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSymptomsImprovedImmunity}] = num_inf[county][i];
            model[county].populations[{AgeGroup(i), InfectionState::InfectedSevereImprovedImmunity}] =
                num_hosp[county][i];
        }
        // }
        if (std::accumulate(num_inf[county].begin(), num_inf[county].end(), 0.0) == 0) {
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
            sum_mu_I_U[region] += model[region].parameters.template get<ICUCasesPerHospitalized>()[i] *
                                  model[region].parameters.template get<HospitalizedCasesPerInfectious>()[i];
            mu_I_U[region].push_back(model[region].parameters.template get<ICUCasesPerHospitalized>()[i] *
                                     model[region].parameters.template get<HospitalizedCasesPerInfectious>()[i]);
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
IOResult<void> set_population_data(std::vector<Model>& model, const std::string& path, const std::string& path_rki,
                                   const std::vector<int>& vregion, Date date)
{
    BOOST_OUTCOME_TRY(num_population, read_population_data(path, vregion));

    auto num_age_groups = ConfirmedCasesDataEntry::age_group_names.size();
    std::vector<std::vector<double>> num_rec(model.size(), std::vector<double>(num_age_groups, 0.0));

    BOOST_OUTCOME_TRY(read_confirmed_cases_data_fix_recovered(path_rki, vregion, date, num_rec, 14.));

    for (size_t region = 0; region < vregion.size(); region++) {
        if (std::accumulate(num_population[region].begin(), num_population[region].end(), 0.0) > 0) {
            auto num_groups = model[region].parameters.get_num_groups();
            for (auto i = AgeGroup(0); i < num_groups; i++) {

                double S_v =
                    std::min(model[region].parameters.template get<DailyFullVaccination>()[{i, SimulationDay(0)}] +
                                 num_rec[region][size_t(i)],
                             num_population[region][size_t(i)]);
                double S_pv =
                    std::max(model[region].parameters.template get<DailyFirstVaccination>()[{i, SimulationDay(0)}] -
                                 model[region].parameters.template get<DailyFullVaccination>()[{i, SimulationDay(0)}],
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
                    1 / (S + S_pv * model[region].parameters.template get<ExposedFactorPartialImmunity>()[i] +
                         S_v * model[region].parameters.template get<ExposedFactorImprovedImmunity>()[i]);
                double denom_C = 1 / (S + S_pv + S_v);
                double denom_I =
                    1 / (S + S_pv * model[region].parameters.template get<InfectedFactorPartialImmunity>()[i] +
                         S_v * model[region].parameters.template get<InfectedFactorImprovedImmunity>()[i]);
                double denom_HU =
                    1 / (S + S_pv * model[region].parameters.template get<HospitalizedFactorPartialImmunity>()[i] +
                         S_v * model[region].parameters.template get<HospitalizedFactorImprovedImmunity>()[i]);

                model[region].populations[{i, InfectionState::ExposedNaive}] =
                    S * model[region].populations[{i, InfectionState::ExposedNaive}] * denom_E;
                model[region].populations[{i, InfectionState::ExposedPartialImmunity}] =
                    S_pv * model[region].parameters.template get<ExposedFactorPartialImmunity>()[i] *
                    model[region].populations[{i, InfectionState::ExposedPartialImmunity}] * denom_E;
                model[region].populations[{i, InfectionState::ExposedImprovedImmunity}] =
                    S_v * model[region].parameters.template get<ExposedFactorImprovedImmunity>()[i] *
                    model[region].populations[{i, InfectionState::ExposedImprovedImmunity}] * denom_E;

                model[region].populations[{i, InfectionState::InfectedNoSymptomsNaive}] =
                    S * model[region].populations[{i, InfectionState::InfectedNoSymptomsNaive}] * denom_C;
                model[region].populations[{i, InfectionState::InfectedNoSymptomsPartialImmunity}] =
                    S_pv * model[region].populations[{i, InfectionState::InfectedNoSymptomsPartialImmunity}] * denom_C;
                model[region].populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunity}] =
                    S_v * model[region].populations[{i, InfectionState::InfectedNoSymptomsImprovedImmunity}] * denom_C;

                model[region].populations[{i, InfectionState::InfectedSymptomsNaive}] =
                    S * model[region].populations[{i, InfectionState::InfectedSymptomsNaive}] * denom_I;
                model[region].populations[{i, InfectionState::InfectedSymptomsPartialImmunity}] =
                    S_pv * model[region].parameters.template get<InfectedFactorPartialImmunity>()[i] *
                    model[region].populations[{i, InfectionState::InfectedSymptomsPartialImmunity}] * denom_I;
                model[region].populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] =
                    S_v * model[region].parameters.template get<InfectedFactorImprovedImmunity>()[i] *
                    model[region].populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] * denom_I;

                model[region].populations[{i, InfectionState::InfectedSevereNaive}] =
                    S * model[region].populations[{i, InfectionState::InfectedSevereNaive}] * denom_HU;
                model[region].populations[{i, InfectionState::InfectedSeverePartialImmunity}] =
                    S_pv * model[region].parameters.template get<HospitalizedFactorPartialImmunity>()[i] *
                    model[region].populations[{i, InfectionState::InfectedSeverePartialImmunity}] * denom_HU;
                model[region].populations[{i, InfectionState::InfectedSevereImprovedImmunity}] =
                    S_v * model[region].parameters.template get<HospitalizedFactorImprovedImmunity>()[i] *
                    model[region].populations[{i, InfectionState::InfectedSevereImprovedImmunity}] * denom_HU;

                model[region].populations[{i, InfectionState::InfectedCriticalNaive}] =
                    S * model[region].populations[{i, InfectionState::InfectedCriticalNaive}] * denom_HU;
                model[region].populations[{i, InfectionState::InfectedCriticalPartialImmunity}] =
                    S_pv * model[region].parameters.template get<HospitalizedFactorPartialImmunity>()[i] *
                    model[region].populations[{i, InfectionState::InfectedCriticalNaive}] * denom_HU;
                model[region].populations[{i, InfectionState::InfectedCriticalImprovedImmunity}] =
                    S_v * model[region].parameters.template get<HospitalizedFactorImprovedImmunity>()[i] *
                    model[region].populations[{i, InfectionState::InfectedCriticalNaive}] * denom_HU;

                model[region].populations[{i, InfectionState::SusceptibleImprovedImmunity}] =
                    model[region].parameters.template get<DailyFullVaccination>()[{i, SimulationDay(0)}] +
                    model[region].populations[{i, InfectionState::SusceptibleImprovedImmunity}] -
                    (model[region].populations[{i, InfectionState::InfectedSymptomsNaive}] +
                     model[region].populations[{i, InfectionState::InfectedSymptomsPartialImmunity}] +
                     model[region].populations[{i, InfectionState::InfectedSymptomsImprovedImmunity}] +
                     model[region].populations[{i, InfectionState::InfectedSevereNaive}] +
                     model[region].populations[{i, InfectionState::InfectedSeverePartialImmunity}] +
                     model[region].populations[{i, InfectionState::InfectedSevereImprovedImmunity}] +
                     model[region].populations[{i, InfectionState::InfectedCriticalNaive}] +
                     model[region].populations[{i, InfectionState::InfectedCriticalPartialImmunity}] +
                     model[region].populations[{i, InfectionState::InfectedCriticalImprovedImmunity}] +
                     model[region].populations[{i, InfectionState::Dead}]);

                model[region].populations[{i, InfectionState::SusceptibleImprovedImmunity}] = std::min(
                    S + S_pv + S_v, std::max(0.0, double(model[region].populations[{i, InfectionState::SusceptibleImprovedImmunity}])));

                model[region].populations[{i, InfectionState::SusceptiblePartialImmunity}] =
                    std::max(0.0, S_pv - model[region].populations[{i, InfectionState::ExposedPartialImmunity}] -
                                      model[region].populations[{i, InfectionState::InfectedNoSymptomsPartialImmunity}] -
                                      model[region].populations[{i, InfectionState::InfectedSymptomsPartialImmunity}] -
                                      model[region].populations[{i, InfectionState::InfectedSeverePartialImmunity}] -
                                      model[region].populations[{i, InfectionState::InfectedCriticalPartialImmunity}]);

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

IOResult<void> set_vaccination_data(std::vector<Model>& model, const std::string& path, Date date,
                                    const std::vector<int>& vregion, int num_days);
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
    */
template <class Model>
IOResult<void> export_input_data_county_timeseries(const std::vector<Model>& model, const std::string& data_dir,
                                                   const std::string& results_dir, std::vector<int> const& region,
                                                   Date start_date, const std::vector<double>& scaling_factor_inf,
                                                   double scaling_factor_icu, int num_days)
{
    auto num_age_groups = (size_t)model[0].parameters.get_num_groups();
    assert(scaling_factor_inf.size() == num_age_groups);
    assert(num_age_groups == ConfirmedCasesDataEntry::age_group_names.size());
    assert(model.size() == region.size());

    BOOST_OUTCOME_TRY(rki_data, read_confirmed_cases_data(path_join(data_dir, "cases_all_county_age_ma7.json")));
    BOOST_OUTCOME_TRY(population_data, read_population_data(path_join(data_dir, "county_current_population.json")));
    BOOST_OUTCOME_TRY(divi_data, read_divi_data(path_join(data_dir, "county_divi_ma7.json")));

    /* functionality copy from set_rki_data() here splitted in params */
    /* which do not need to be reset for each day and compartments sizes that are */
    /* set later for each day */
    /*----------- UNVACCINATED -----------*/
    // data needs to be int, because access to data-specific confirmed cases
    // is done with these parameters. TODO: Rounding instead
    // of casting to int should be introduced in the future.
    std::vector<std::vector<int>> t_car_to_rec_uv{model.size()}; // R9
    std::vector<std::vector<int>> t_car_to_inf_uv{model.size()}; // R3
    std::vector<std::vector<int>> t_exp_to_car_uv{model.size()}; // R2
    std::vector<std::vector<int>> t_inf_to_rec_uv{model.size()}; // R4
    std::vector<std::vector<int>> t_inf_to_hosp_uv{model.size()}; // R6
    std::vector<std::vector<int>> t_hosp_to_rec_uv{model.size()}; // R5
    std::vector<std::vector<int>> t_icu_to_rec_uv{model.size()};

    std::vector<std::vector<double>> mu_C_R_uv{model.size()};
    std::vector<std::vector<double>> mu_I_H_uv{model.size()};
    std::vector<std::vector<double>> mu_H_U_uv{model.size()};
    std::vector<std::vector<double>> mu_U_D_uv{model.size()};
    // ICU data is not age-resolved. Use a partition of unity defined by
    // the age-dependent probability I->H->U divided by the sum over all
    // age groups of all of these probabilities.
    std::vector<double> sum_mu_I_U_uv(model.size(), 0);
    std::vector<std::vector<double>> mu_I_U_uv{model.size()};

    for (size_t county = 0; county < model.size(); county++) {
        for (size_t group = 0; group < num_age_groups; group++) {

            t_car_to_inf_uv[county].push_back(
                static_cast<int>(2 * (model[county].parameters.template get<IncubationTime>()[(AgeGroup)group] -
                                      model[county].parameters.template get<SerialInterval>()[(AgeGroup)group])));
            t_car_to_rec_uv[county].push_back(
                static_cast<int>(t_car_to_inf_uv[county][group] +
                                 0.5 * model[county].parameters.template get<TimeInfectedSymptoms>()[(AgeGroup)group]));
            t_exp_to_car_uv[county].push_back(
                static_cast<int>(2 * model[county].parameters.template get<SerialInterval>()[(AgeGroup)group] -
                                 model[county].parameters.template get<IncubationTime>()[(AgeGroup)group]));
            t_inf_to_rec_uv[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedSymptoms>()[(AgeGroup)group]));
            t_hosp_to_rec_uv[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedSevere>()[(AgeGroup)group]));
            t_hosp_to_icu_uv[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedSevere>()[(AgeGroup)group]));
            t_icu_to_dead_uv[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedCritical>()[(AgeGroup)group]));
            t_icu_to_rec_uv[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedCritical>()[(AgeGroup)group]));

            mu_C_R_uv[county].push_back(
                model[county].parameters.template get<AsymptoticCasesPerInfectious>()[(AgeGroup)group]);
            mu_I_H_uv[county].push_back(
                model[county].parameters.template get<HospitalizedCasesPerInfectious>()[(AgeGroup)group]);
            mu_H_U_uv[county].push_back(
                model[county].parameters.template get<ICUCasesPerHospitalized>()[(AgeGroup)group]);
            mu_U_D_uv[county].push_back(model[county].parameters.template get<DeathsPerICU>()[(AgeGroup)group]);

            /* begin: NOT in set_rki_data() */
            sum_mu_I_U_uv[county] +=
                model[county].parameters.template get<ICUCasesPerHospitalized>()[AgeGroup(group)] *
                model[county].parameters.template get<HospitalizedCasesPerInfectious>()[AgeGroup(group)];
            mu_I_U_uv[county].push_back(
                model[county].parameters.template get<ICUCasesPerHospitalized>()[AgeGroup(group)] *
                model[county].parameters.template get<HospitalizedCasesPerInfectious>()[AgeGroup(group)]);
            /* end: NOT in set_rki_data() */
        }
    }

    /*----------- PARTIALLY VACCINATED -----------*/
    // data needs to be int, because access to data-specific confirmed cases
    // is done with these parameters. TODO: Rounding instead
    // of casting to int should be introduced in the future.
    std::vector<std::vector<int>> t_car_to_rec_pv{model.size()}; // R9
    std::vector<std::vector<int>> t_car_to_inf_pv{model.size()}; // R3
    std::vector<std::vector<int>> t_exp_to_car_pv{model.size()}; // R2
    std::vector<std::vector<int>> t_inf_to_rec_pv{model.size()}; // R4
    std::vector<std::vector<int>> t_inf_to_hosp_pv{model.size()}; // R6
    std::vector<std::vector<int>> t_hosp_to_rec_pv{model.size()}; // R5
    std::vector<std::vector<int>> t_icu_to_rec_pv{model.size()};

    std::vector<std::vector<double>> mu_C_R_pv{model.size()};
    std::vector<std::vector<double>> mu_I_H_pv{model.size()};
    std::vector<std::vector<double>> mu_H_U_pv{model.size()};
    std::vector<std::vector<double>> mu_U_D_pv{model.size()};
    // ICU data is not age-resolved. Use a partition of unity defined by
    // the age-dependent probability I->H->U divided by the sum over all
    // age groups of all of these probabilities.
    std::vector<double> sum_mu_I_U_pv(model.size(), 0);
    std::vector<std::vector<double>> mu_I_U_pv{model.size()};
    for (size_t county = 0; county < model.size(); county++) {
        for (size_t group = 0; group < num_age_groups; group++) {

            double reduc_t = model[0].parameters.template get<InfectiousTimeFactorImmune>()[(AgeGroup)group];
            t_car_to_inf_pv[county].push_back(
                static_cast<int>(2 * (model[county].parameters.template get<IncubationTime>()[(AgeGroup)group] -
                                      model[county].parameters.template get<SerialInterval>()[(AgeGroup)group])));
            t_car_to_rec_pv[county].push_back(static_cast<int>(
                t_car_to_inf_pv[county][group] +
                0.5 * model[county].parameters.template get<TimeInfectedSymptoms>()[(AgeGroup)group] * reduc_t));
            t_exp_to_car_pv[county].push_back(
                static_cast<int>(2 * model[county].parameters.template get<SerialInterval>()[(AgeGroup)group] -
                                 model[county].parameters.template get<IncubationTime>()[(AgeGroup)group]));
            t_inf_to_rec_pv[county].push_back(static_cast<int>(
                model[county].parameters.template get<TimeInfectedSymptoms>()[(AgeGroup)group] * reduc_t));
            t_hosp_to_rec_pv[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedSevere>()[(AgeGroup)group]));
            t_hosp_to_icu_pv[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedSevere>()[(AgeGroup)group]));
            t_icu_to_dead_pv[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedCritical>()[(AgeGroup)group]));
            t_icu_to_rec_pv[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedCritical>()[(AgeGroup)group]));

            double exp_fact_part_immune =
                model[county].parameters.template get<ExposedFactorPartialImmunity>()[(AgeGroup)group];
            double inf_fact_part_immune =
                model[county].parameters.template get<InfectedFactorPartialImmunity>()[(AgeGroup)group];
            double hosp_fact_part_immune =
                model[county].parameters.template get<HospitalizedFactorPartialImmunity>()[(AgeGroup)group];
            double icu_fact_part_immune =
                model[county].parameters.template get<HospitalizedFactorPartialImmunity>()[(AgeGroup)group];
            double death_fact_part_immune =
                model[county].parameters.template get<HospitalizedFactorPartialImmunity>()[(AgeGroup)group];
            mu_C_R_pv[county].push_back(
                (1 - inf_fact_part_immune / exp_fact_part_immune *
                         (1 - model[county].parameters.template get<AsymptoticCasesPerInfectious>()[(AgeGroup)group])));
            mu_I_H_pv[county].push_back(
                hosp_fact_part_immune / inf_fact_part_immune *
                model[county].parameters.template get<HospitalizedCasesPerInfectious>()[(AgeGroup)group]);
            // transfer from H to U, D unchanged.
            mu_H_U_pv[county].push_back(
                icu_fact_part_immune / hosp_fact_part_immune *
                model[county].parameters.template get<ICUCasesPerHospitalized>()[(AgeGroup)group]);
            mu_U_D_pv[county].push_back(death_fact_part_immune / icu_fact_part_immune *
                                        model[county].parameters.template get<DeathsPerICU>()[(AgeGroup)group]);

            sum_mu_I_U_pv[county] +=
                icu_fact_part_immune / hosp_fact_part_immune *
                model[county].parameters.template get<ICUCasesPerHospitalized>()[AgeGroup(group)] *
                hosp_fact_part_immune / inf_fact_part_immune *
                model[county].parameters.template get<HospitalizedCasesPerInfectious>()[AgeGroup(group)];
            mu_I_U_pv[county].push_back(
                icu_fact_part_immune / hosp_fact_part_immune *
                model[county].parameters.template get<ICUCasesPerHospitalized>()[AgeGroup(group)] *
                hosp_fact_part_immune / inf_fact_part_immune *
                model[county].parameters.template get<HospitalizedCasesPerInfectious>()[AgeGroup(group)]);
        }
    }

    /*----------- FULLY VACCINATED -----------*/
    // data needs to be int, because access to data-specific confirmed cases
    // is done with these parameters. TODO: Rounding instead
    // of casting to int should be introduced in the future.
    std::vector<std::vector<int>> t_car_to_rec_fv{model.size()}; // R9
    std::vector<std::vector<int>> t_car_to_inf_fv{model.size()}; // R3
    std::vector<std::vector<int>> t_exp_to_car_fv{model.size()}; // R2
    std::vector<std::vector<int>> t_inf_to_rec_fv{model.size()}; // R4
    std::vector<std::vector<int>> t_inf_to_hosp_fv{model.size()}; // R6
    std::vector<std::vector<int>> t_hosp_to_rec_fv{model.size()}; // R5
    std::vector<std::vector<int>> t_icu_to_rec_fv{model.size()};

    std::vector<std::vector<double>> mu_C_R_fv{model.size()};
    std::vector<std::vector<double>> mu_I_H_fv{model.size()};
    std::vector<std::vector<double>> mu_H_U_fv{model.size()};
    std::vector<std::vector<double>> mu_U_D_fv{model.size()};
    // ICU data is not age-resolved. Use a partition of unity defined by
    // the age-dependent probability I->H->U divided by the sum over all
    // age groups of all of these probabilities.
    std::vector<double> sum_mu_I_U_fv(model.size(), 0);
    std::vector<std::vector<double>> mu_I_U_fv{model.size()};
    for (size_t county = 0; county < model.size(); county++) {
        for (size_t group = 0; group < num_age_groups; group++) {

            double reduc_t = model[0].parameters.template get<InfectiousTimeFactorImmune>()[(AgeGroup)group];
            t_car_to_inf_fv[county].push_back(
                static_cast<int>(2 * (model[county].parameters.template get<IncubationTime>()[(AgeGroup)group] -
                                      model[county].parameters.template get<SerialInterval>()[(AgeGroup)group])));
            t_car_to_rec_fv[county].push_back(static_cast<int>(
                t_car_to_inf_fv[county][group] +
                0.5 * model[county].parameters.template get<TimeInfectedSymptoms>()[(AgeGroup)group] * reduc_t));
            t_exp_to_car_fv[county].push_back(
                static_cast<int>(2 * model[county].parameters.template get<SerialInterval>()[(AgeGroup)group] -
                                 model[county].parameters.template get<IncubationTime>()[(AgeGroup)group]));
            t_inf_to_rec_fv[county].push_back(static_cast<int>(
                model[county].parameters.template get<TimeInfectedSymptoms>()[(AgeGroup)group] * reduc_t));
            t_hosp_to_rec_fv[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedSevere>()[(AgeGroup)group]));
            t_hosp_to_icu_fv[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedSevere>()[(AgeGroup)group]));
            t_icu_to_dead_fv[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedCritical>()[(AgeGroup)group]));
            t_icu_to_rec_fv[county].push_back(
                static_cast<int>(model[county].parameters.template get<TimeInfectedCritical>()[(AgeGroup)group]));

            double reduc_immune_exp =
                model[county].parameters.template get<ExposedFactorImprovedImmunity>()[(AgeGroup)group];
            double reduc_immune_inf =
                model[county].parameters.template get<InfectedFactorImprovedImmunity>()[(AgeGroup)group];
            double reduc_immune_hosp =
                model[county].parameters.template get<HospitalizedFactorImprovedImmunity>()[(AgeGroup)group];
            double reduc_immune_icu =
                model[county].parameters.template get<HospitalizedFactorImprovedImmunity>()[(AgeGroup)group];
            double reduc_immune_dead =
                model[county].parameters.template get<HospitalizedFactorImprovedImmunity>()[(AgeGroup)group];
            mu_C_R_fv[county].push_back(
                (1 - reduc_immune_inf / reduc_immune_exp *
                         (1 - model[county].parameters.template get<AsymptoticCasesPerInfectious>()[(AgeGroup)group])));
            mu_I_H_fv[county].push_back(
                reduc_immune_hosp / reduc_immune_inf *
                model[county].parameters.template get<HospitalizedCasesPerInfectious>()[(AgeGroup)group]);
            // transfer from H to U, D unchanged.
            mu_H_U_fv[county].push_back(
                reduc_immune_icu / reduc_immune_hosp *
                model[county].parameters.template get<ICUCasesPerHospitalized>()[(AgeGroup)group]);
            mu_U_D_fv[county].push_back(reduc_immune_dead / reduc_immune_icu *
                                        model[county].parameters.template get<DeathsPerICU>()[(AgeGroup)group]);

            sum_mu_I_U_fv[county] +=
                reduc_immune_icu / reduc_immune_hosp *
                model[county].parameters.template get<ICUCasesPerHospitalized>()[AgeGroup(group)] * reduc_immune_hosp /
                reduc_immune_inf *
                model[county].parameters.template get<HospitalizedCasesPerInfectious>()[AgeGroup(group)];
            mu_I_U_fv[county].push_back(
                reduc_immune_icu / reduc_immune_hosp *
                model[county].parameters.template get<ICUCasesPerHospitalized>()[AgeGroup(group)] * reduc_immune_hosp /
                reduc_immune_inf *
                model[county].parameters.template get<HospitalizedCasesPerInfectious>()[AgeGroup(group)]);
        }
    }
    std::vector<TimeSeries<double>> extrapolated_rki(
        model.size(), TimeSeries<double>::zero(num_days + 1, (size_t)InfectionState::Count * num_age_groups));

    for (size_t day = 0; day <= static_cast<size_t>(num_days); day++) {
        auto date = offset_date_by_days(start_date, int(day));

        // unvaccinated
        std::vector<std::vector<double>> num_exp_uv(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_car_uv(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_inf_uv(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_rec_uv(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_hosp_uv(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_death_uv(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> dummy_icu(model.size(), std::vector<double>(num_age_groups, 0.0));
        BOOST_OUTCOME_TRY(details::read_confirmed_cases_data(
            rki_data, region, date, num_exp_uv, num_car_uv, num_inf_uv, num_hosp_uv, dummy_icu, num_death_uv,
            num_rec_uv, t_car_to_rec_uv, t_car_to_inf_uv, t_exp_to_car_uv, t_inf_to_rec_uv, t_inf_to_hosp_uv,
            t_hosp_to_rec_uv, t_hosp_to_icu_uv, t_icu_to_dead_uv, t_icu_to_rec_uv, mu_C_R_uv, mu_I_H_uv, mu_H_U_uv,
            mu_U_D_uv, scaling_factor_inf));

        // partially vaccinated
        std::vector<std::vector<double>> num_exp_pv(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_car_pv(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_inf_pv(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_hosp_pv(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> dummy_death(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> dummy_rec(model.size(), std::vector<double>(num_age_groups, 0.0));
        for (size_t county = 0; county < model.size(); county++) {
            dummy_death[county] = std::vector<double>(num_age_groups, 0.0);
            dummy_icu[county]   = std::vector<double>(num_age_groups, 0.0);
        }
        BOOST_OUTCOME_TRY(details::read_confirmed_cases_data(
            rki_data, region, date, num_exp_pv, num_car_pv, num_inf_pv, num_hosp_pv, dummy_icu, dummy_death, dummy_rec,
            t_car_to_rec_pv, t_car_to_inf_pv, t_exp_to_car_pv, t_inf_to_rec_pv, t_inf_to_hosp_pv, t_hosp_to_rec_pv,
            t_hosp_to_icu_pv, t_icu_to_dead_pv, t_icu_to_rec_pv, mu_C_R_pv, mu_I_H_pv, mu_H_U_pv, mu_U_D_pv,
            scaling_factor_inf));

        // fully vaccinated
        std::vector<std::vector<double>> num_exp_fv(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_car_fv(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_inf_fv(model.size(), std::vector<double>(num_age_groups, 0.0));
        std::vector<std::vector<double>> num_hosp_fv(model.size(), std::vector<double>(num_age_groups, 0.0));
        for (size_t county = 0; county < model.size(); county++) {
            dummy_rec[county]   = std::vector<double>(num_age_groups, 0.0);
            dummy_death[county] = std::vector<double>(num_age_groups, 0.0);
            dummy_icu[county]   = std::vector<double>(num_age_groups, 0.0);
        }
        BOOST_OUTCOME_TRY(details::read_confirmed_cases_data(
            rki_data, region, date, num_exp_fv, num_car_fv, num_inf_fv, num_hosp_fv, dummy_icu, dummy_death, dummy_rec,
            t_car_to_rec_fv, t_car_to_inf_fv, t_exp_to_car_fv, t_inf_to_rec_fv, t_inf_to_hosp_fv, t_hosp_to_rec_fv,
            t_hosp_to_icu_fv, t_icu_to_dead_fv, t_icu_to_rec_fv, mu_C_R_fv, mu_I_H_fv, mu_H_U_fv, mu_U_D_fv,
            scaling_factor_inf));

        // ICU only read for compartment InfectionState::ICU and then distributed later
        std::vector<double> dummy_icu2(model.size(), 0.0);
        BOOST_OUTCOME_TRY(details::read_divi_data(divi_data, region, date, dummy_icu2));

        std::vector<std::vector<double>> num_icu(model.size(), std::vector<double>(num_age_groups, 0.0));
        for (size_t county = 0; county < region.size(); county++) {
            for (size_t age = 0; age < num_age_groups; age++) {
                num_icu[county][age] =
                    scaling_factor_icu * dummy_icu2[county] * mu_I_U_uv[county][age] / sum_mu_I_U_uv[county];
            }
        }

        // read population basics
        BOOST_OUTCOME_TRY(num_population, details::read_population_data(population_data, region));

        std::vector<std::vector<double>> num_rec(model.size(), std::vector<double>(num_age_groups, 0.0));
        BOOST_OUTCOME_TRY(details::read_confirmed_cases_data_fix_recovered(rki_data, region, date, num_rec, 14.));

        for (size_t county = 0; county < region.size(); county++) {
            if (std::accumulate(num_population[county].begin(), num_population[county].end(), 0.0) > 0) {
                for (size_t age = 0; age < num_age_groups; age++) {
                    auto age_group_offset = age * (size_t)InfectionState::Count;
                    double S_v            = std::min(
                        model[county]
                                .parameters.template get<DailyFullVaccination>()[{AgeGroup(age), SimulationDay(day)}] +
                            num_rec[county][age],
                        num_population[county][age]);
                    double S_pv = std::max(
                        model[county]
                                .parameters.template get<DailyFirstVaccination>()[{AgeGroup(age), SimulationDay(day)}] -
                            model[county]
                                .parameters.template get<DailyFullVaccination>()[{AgeGroup(age), SimulationDay(day)}],
                        0.0); // use std::max with 0
                    double S;
                    if (num_population[county][age] - S_pv - S_v < 0.0) {
                        log_warning("Number of vaccinated greater than population at county {}, age group {}: {} + "
                                    "{} > {}.",
                                    county, age, S_pv, S_v, num_population[county][age]);
                        S   = 0.0;
                        S_v = num_population[county][age] - S_pv;
                    }
                    else {
                        S = num_population[county][age] - S_pv - S_v;
                    }

                    double denom_E =
                        1 /
                        (S +
                         S_pv * model[county].parameters.template get<ExposedFactorPartialImmunity>()[AgeGroup(age)] +
                         S_v * model[county].parameters.template get<ExposedFactorImprovedImmunity>()[AgeGroup(age)]);
                    double denom_C = 1 / (S + S_pv + S_v);
                    double denom_I =
                        1 /
                        (S +
                         S_pv * model[county].parameters.template get<InfectedFactorPartialImmunity>()[AgeGroup(age)] +
                         S_v * model[county].parameters.template get<InfectedFactorImprovedImmunity>()[AgeGroup(age)]);
                    double denom_HU =
                        1 /
                        (S +
                         S_pv *
                             model[county].parameters.template get<HospitalizedFactorPartialImmunity>()[AgeGroup(age)] +
                         S_v * model[county]
                                   .parameters.template get<HospitalizedFactorImprovedImmunity>()[AgeGroup(age)]);

                    extrapolated_rki[county][day]((size_t)InfectionState::ExposedNaive + age_group_offset) =
                        S * denom_E * num_exp_uv[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::ExposedPartialImmunity + age_group_offset) =
                        S_pv * model[county].parameters.template get<ExposedFactorPartialImmunity>()[AgeGroup(age)] *
                        denom_E * num_exp_pv[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::ExposedImprovedImmunity + age_group_offset) =
                        S_v * model[county].parameters.template get<ExposedFactorImprovedImmunity>()[AgeGroup(age)] *
                        denom_E * num_exp_fv[county][age];

                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedNoSymptomsNaive + age_group_offset) =
                        S * denom_C * num_car_uv[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedNoSymptomsPartialImmunity + age_group_offset) =
                        S_pv * denom_C * num_car_pv[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedNoSymptomsImprovedImmunity + age_group_offset) =
                        S_v * denom_C * num_car_fv[county][age];

                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsNaive + age_group_offset) =
                        S * denom_I * num_inf_uv[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsPartialImmunity + age_group_offset) =
                        S_pv * model[county].parameters.template get<InfectedFactorPartialImmunity>()[AgeGroup(age)] *
                        denom_I * num_inf_pv[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsImprovedImmunity + age_group_offset) =
                        S_v * model[county].parameters.template get<InfectedFactorImprovedImmunity>()[AgeGroup(age)] *
                        denom_I * num_inf_fv[county][age];

                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedSevereNaive + age_group_offset) =
                        S * denom_HU * num_hosp_uv[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedSeverePartialImmunity +
                                                  age_group_offset) =
                        S_pv *
                        model[county].parameters.template get<HospitalizedFactorPartialImmunity>()[AgeGroup(age)] *
                        denom_HU * num_hosp_pv[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedSevereImprovedImmunity +
                                                  age_group_offset) =
                        S_v *
                        model[county].parameters.template get<HospitalizedFactorImprovedImmunity>()[AgeGroup(age)] *
                        denom_HU * num_hosp_fv[county][age];

                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedCriticalNaive + age_group_offset) =
                        S * denom_HU * num_icu[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedCriticalPartialImmunity + age_group_offset) =
                        S_pv *
                        model[county].parameters.template get<HospitalizedFactorPartialImmunity>()[AgeGroup(age)] *
                        denom_HU * num_icu[county][age];
                    extrapolated_rki[county][day]((size_t)InfectionState::InfectedCriticalImprovedImmunity + age_group_offset) =
                        S_v *
                        model[county].parameters.template get<HospitalizedFactorImprovedImmunity>()[AgeGroup(age)] *
                        denom_HU * num_icu[county][age];

                    // in set_rki_data initilization, deaths are now set to 0. In order to visualize
                    // the extrapolated real number of deaths, they have to be set here. In the comparison of data
                    // it has to be paid attention to the fact, the the simulation starts with deaths=0
                    // while this method starts with deaths=number of reported deaths so far...
                    extrapolated_rki[county][day]((size_t)InfectionState::Dead + age_group_offset) =
                        num_death_uv[county][age];

                    extrapolated_rki[county][day]((size_t)InfectionState::SusceptibleImprovedImmunity + age_group_offset) =
                        model[county]
                            .parameters.template get<DailyFullVaccination>()[{AgeGroup(age), SimulationDay(day)}] +
                        num_rec_uv[county][age] -
                        (extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsNaive + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsPartialImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsImprovedImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedSevereNaive + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedSeverePartialImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedSevereImprovedImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedCriticalNaive + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedCriticalPartialImmunity + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedCriticalImprovedImmunity + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::Dead + age_group_offset));

                    extrapolated_rki[county][day]((size_t)InfectionState::SusceptibleImprovedImmunity + age_group_offset) = std::min(
                        S + S_pv + S_v, std::max(0.0, double(extrapolated_rki[county][day](
                                                          (size_t)InfectionState::SusceptibleImprovedImmunity + age_group_offset))));

                    extrapolated_rki[county][day]((size_t)InfectionState::SusceptiblePartialImmunity +
                                                  age_group_offset) =
                        std::max(0.0, S_pv -
                                          extrapolated_rki[county][day]((size_t)InfectionState::ExposedPartialImmunity +
                                                                        age_group_offset) -
                                          extrapolated_rki[county][day]((size_t)InfectionState::InfectedNoSymptomsPartialImmunity +
                                                                        age_group_offset) -
                                          extrapolated_rki[county][day](
                                              (size_t)InfectionState::InfectedSymptomsPartialImmunity + age_group_offset) -
                                          extrapolated_rki[county][day](
                                              (size_t)InfectionState::InfectedSeverePartialImmunity + age_group_offset) -
                                          extrapolated_rki[county][day]((size_t)InfectionState::InfectedCriticalPartialImmunity +
                                                                        age_group_offset));

                    extrapolated_rki[county][day]((size_t)InfectionState::SusceptibleNaive + age_group_offset) =
                        num_population[county][age] -
                        (extrapolated_rki[county][day]((size_t)InfectionState::SusceptiblePartialImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::SusceptibleImprovedImmunity + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::ExposedNaive + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::ExposedPartialImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::ExposedImprovedImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedNoSymptomsNaive + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedNoSymptomsPartialImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedNoSymptomsImprovedImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsNaive + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsPartialImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedSymptomsImprovedImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedSevereNaive + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedSeverePartialImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedSevereImprovedImmunity +
                                                       age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedCriticalNaive + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedCriticalPartialImmunity + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedCriticalImprovedImmunity + age_group_offset) +
                         extrapolated_rki[county][day]((size_t)InfectionState::Dead + age_group_offset));
                }
            }
            else {
                log_warning("No population data available for region " + std::to_string(county) +
                            ". Population data has not been set.");
            }
        }
        log_info("extrapolated real data for date: {}-{}-{}", date.day, date.month, date.year);
    }
    /* end: similar functionality in set_rki_data(), here only for vector of TimeSeries */
    auto num_groups = (int)(size_t)model[0].parameters.get_num_groups();
    BOOST_OUTCOME_TRY(save_result(extrapolated_rki, region, num_groups, path_join(results_dir, "Results_rki.h5")));

    auto extrapolated_rki_data_sum = sum_nodes(std::vector<std::vector<TimeSeries<double>>>{extrapolated_rki});
    BOOST_OUTCOME_TRY(
        save_result({extrapolated_rki_data_sum[0][0]}, {0}, num_groups, path_join(results_dir, "Results_rki_sum.h5")));

    return success();
}

#endif //MEMILIO_HAS_HDF5

/**
    * Reads compartments for German counties at a specified date from data files.
    * Estimates all compartments from available data using the model parameters, so the 
    * model parameters must be set before calling this function.
    * Uses data files that contain centered 7-day moving average.
    * @param model Vector of SECIRVVS models, one per county.
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
    BOOST_OUTCOME_TRY(details::set_vaccination_data(model, path_join(dir, "all_county_ageinf_vacc_ma7.json"), date,
                                                    county, num_days));

    // TODO: Reuse more code, e.g., set_divi_data (in secir) and a set_divi_data (here) only need a different ModelType.
    // TODO: add option to set ICU data from confirmed cases if DIVI or other data is not available.
    if (date > Date(2020, 4, 23)) {
        BOOST_OUTCOME_TRY(
            details::set_divi_data(model, path_join(dir, "county_divi_ma7.json"), county, date, scaling_factor_icu));
    }
    else {
        log_warning("No DIVI data available for this date");
    }

    BOOST_OUTCOME_TRY(details::set_confirmed_cases_data(model, path_join(dir, "cases_all_county_age_ma7.json"), county,
                                                        date, scaling_factor_inf));
    BOOST_OUTCOME_TRY(details::set_population_data(model, path_join(dir, "county_current_population.json"),
                                                   path_join(dir, "cases_all_county_age_ma7.json"), county, date));

    if (export_time_series) {
        // Use only if extrapolated real data is needed for comparison. EXPENSIVE !
        // Run time equals run time of the previous functions times the num_days !
        // (This only represents the vectorization of the previous function over all simulation days...)
        BOOST_OUTCOME_TRY(export_input_data_county_timeseries(model, dir, dir, county, date, scaling_factor_inf,
                                                              scaling_factor_icu, num_days));
    }

    return success();
}

} // namespace osecirvvs
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // ODESECIRVVS_PARAMETERS_IO_H
