/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Wadim Koslow, Daniel Abele, Martin J. Kuehn, Lena Ploetzke
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
#ifndef SECIRV_PARAMETERS_IO_H
#define SECIRV_PARAMETERS_IO_H

#include "memilio/config.h"
#include "memilio/io/epi_data.h"
#include "memilio/utils/logging.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "secir_vaccine/model.h"
#include "secir_vaccine/analyze_result.h"
#include "memilio/math/eigen_util.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/mobility.h"
#include "memilio/io/io.h"
#include "memilio/io/json_serializer.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/date.h"

namespace mio
{
namespace secirv
{

    namespace details
    {
        /**
        * @brief reads populations data from RKI
        * @param path Path to RKI file
        * @param vregion vector of keys of the region of interest     
        * @param date Date for which the arrays are initialized
        * @param num_* output vector for number of people in the corresponding compartement
        * @param t_* vector average time it takes to get from one compartement to another for each age group
        * @param mu_* vector probabilities to get from one compartement to another for each age group
        * @see mio::read_rki_data
        * @{
        */
        IOResult<void> read_rki_data(
            std::string const& path, std::vector<int> const& vregion, Date date,
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

        IOResult<void> read_rki_data(
            const std::vector<RkiEntry>& rki_data, std::vector<int> const& vregion, Date date,
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
        * @brief reads confirmed cases data and translates data of day t0-delay to recovered compartment,
        * @param path Path to RKI file
        * @param vregion vector of keys of the region of interest     
        * @param date Date for which the arrays are initialized
        * @param num_rec output vector for number of people in the compartement recovered
        * @param delay number of days in the past the are used to set recovered compartment.
        * @see mio::read_rki_data
        * @{
        */
        IOResult<void> read_rki_data_confirmed_to_recovered(const std::vector<RkiEntry>& rki_data,
                                                            std::vector<int> const& vregion, Date date,
                                                            std::vector<std::vector<double>>& vnum_rec,
                                                            double delay = 14.);
        IOResult<void> read_rki_data_confirmed_to_recovered(std::string const& path, std::vector<int> const& vregion,
                                                            Date date, std::vector<std::vector<double>>& vnum_rec,
                                                            double delay = 14.);
        /**@}*/

        /**
        * @brief sets populations data from RKI into a SecirModel
        * @param model vector of objects in which the data is set
        * @param path Path to RKI file
        * @param region vector of keys of the region of interest
        * @param date Date for which the arrays are initialized
        * @param scaling_factor_inf factors by which to scale the confirmed cases of
        * rki data
        */
        template <class Model>
        IOResult<void> set_rki_data(std::vector<Model>& model, const std::string& path, std::vector<int> const& region,
                                    Date date, const std::vector<double>& scaling_factor_inf)
        {
            auto num_age_groups = (size_t)model[0].parameters.get_num_groups();
            assert(scaling_factor_inf.size() == num_age_groups);
            assert(StringRkiAgeGroup::age_group_names.size() == num_age_groups);

            BOOST_OUTCOME_TRY(rki_data, mio::read_rki_data(path));

            std::vector<std::vector<int>> t_car_to_rec{model.size()}; // R9
            std::vector<std::vector<int>> t_car_to_inf{model.size()}; // R3
            std::vector<std::vector<int>> t_exp_to_car{model.size()}; // R2
            std::vector<std::vector<int>> t_inf_to_rec{model.size()}; // R4
            std::vector<std::vector<int>> t_inf_to_hosp{model.size()}; // R6
            std::vector<std::vector<int>> t_hosp_to_rec{model.size()}; // R5
            std::vector<std::vector<int>> t_hosp_to_icu{model.size()}; // R7
            std::vector<std::vector<int>> t_icu_to_dead{model.size()}; // R10
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

                    t_car_to_inf[county].push_back(static_cast<int>(
                        2 * (model[county].parameters.template get<IncubationTime>()[(AgeGroup)group] -
                             model[county].parameters.template get<SerialInterval>()[(AgeGroup)group])));
                    t_car_to_rec[county].push_back(static_cast<int>(
                        t_car_to_inf[county][group] +
                        0.5 * model[county].parameters.template get<InfectiousTimeMild>()[(AgeGroup)group]));
                    t_exp_to_car[county].push_back(
                        static_cast<int>(2 * model[county].parameters.template get<SerialInterval>()[(AgeGroup)group] -
                                         model[county].parameters.template get<IncubationTime>()[(AgeGroup)group]));
                    t_inf_to_rec[county].push_back(
                        static_cast<int>(model[county].parameters.template get<InfectiousTimeMild>()[(AgeGroup)group]));
                    t_inf_to_hosp[county].push_back(static_cast<int>(
                        model[county].parameters.template get<HomeToHospitalizedTime>()[(AgeGroup)group]));
                    t_hosp_to_rec[county].push_back(static_cast<int>(
                        model[county].parameters.template get<HospitalizedToHomeTime>()[(AgeGroup)group]));
                    t_hosp_to_icu[county].push_back(static_cast<int>(
                        model[county].parameters.template get<HospitalizedToICUTime>()[(AgeGroup)group]));
                    t_icu_to_dead[county].push_back(
                        static_cast<int>(model[county].parameters.template get<ICUToDeathTime>()[(AgeGroup)group]));
                    t_icu_to_rec[county].push_back(
                        static_cast<int>(model[county].parameters.template get<ICUToHomeTime>()[(AgeGroup)group]));

                    mu_C_R[county].push_back(
                        model[county].parameters.template get<AsymptoticCasesPerInfectious>()[(AgeGroup)group]);
                    mu_I_H[county].push_back(
                        model[county].parameters.template get<HospitalizedCasesPerInfectious>()[(AgeGroup)group]);
                    mu_H_U[county].push_back(
                        model[county].parameters.template get<ICUCasesPerHospitalized>()[(AgeGroup)group]);
                    mu_U_D[county].push_back(model[county].parameters.template get<DeathsPerICU>()[(AgeGroup)group]);
                }
            }

            BOOST_OUTCOME_TRY(read_rki_data(rki_data, region, date, num_exp, num_car, num_inf, num_hosp, num_icu,
                                            num_death, num_rec, t_car_to_rec, t_car_to_inf, t_exp_to_car, t_inf_to_rec,
                                            t_inf_to_hosp, t_hosp_to_rec, t_hosp_to_icu, t_icu_to_dead, t_icu_to_rec,
                                            mu_C_R, mu_I_H, mu_H_U, mu_U_D, scaling_factor_inf));

            for (size_t county = 0; county < model.size(); county++) {
                // if (std::accumulate(num_inf[county].begin(), num_inf[county].end(), 0.0) > 0) {
                size_t num_groups = (size_t)model[county].parameters.get_num_groups();
                for (size_t i = 0; i < num_groups; i++) {
                    model[county].populations[{AgeGroup(i), InfectionState::Exposed}]      = num_exp[county][i];
                    model[county].populations[{AgeGroup(i), InfectionState::Carrier}]      = num_car[county][i];
                    model[county].populations[{AgeGroup(i), InfectionState::Infected}]     = num_inf[county][i];
                    model[county].populations[{AgeGroup(i), InfectionState::Hospitalized}] = num_hosp[county][i];
                    model[county].populations[{AgeGroup(i), InfectionState::Recovered}]    = num_rec[county][i];
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
                    t_car_to_inf[county].push_back(static_cast<int>(
                        2 * (model[county].parameters.template get<IncubationTime>()[(AgeGroup)group] -
                             model[county].parameters.template get<SerialInterval>()[(AgeGroup)group])));
                    t_car_to_rec[county].push_back(static_cast<int>(
                        t_car_to_inf[county][group] +
                        0.5 * model[county].parameters.template get<InfectiousTimeMild>()[(AgeGroup)group] * reduc_t));
                    t_exp_to_car[county].push_back(
                        static_cast<int>(2 * model[county].parameters.template get<SerialInterval>()[(AgeGroup)group] -
                                         model[county].parameters.template get<IncubationTime>()[(AgeGroup)group]));
                    t_inf_to_rec[county].push_back(static_cast<int>(
                        model[county].parameters.template get<InfectiousTimeMild>()[(AgeGroup)group] * reduc_t));
                    t_inf_to_hosp[county].push_back(static_cast<int>(
                        model[county].parameters.template get<HomeToHospitalizedTime>()[(AgeGroup)group]));
                    t_hosp_to_rec[county].push_back(static_cast<int>(
                        model[county].parameters.template get<HospitalizedToHomeTime>()[(AgeGroup)group]));
                    t_hosp_to_icu[county].push_back(static_cast<int>(
                        model[county].parameters.template get<HospitalizedToICUTime>()[(AgeGroup)group]));
                    t_icu_to_dead[county].push_back(
                        static_cast<int>(model[county].parameters.template get<ICUToDeathTime>()[(AgeGroup)group]));
                    t_icu_to_rec[county].push_back(
                        static_cast<int>(model[county].parameters.template get<ICUToHomeTime>()[(AgeGroup)group]));

                    double reduc_vacc_exp  = model[county].parameters.template get<ExposedFactorPartiallyImmune>()[(AgeGroup)group];
                    double reduc_vacc_inf  = model[county].parameters.template get<InfectedFactorPartiallyImmune>()[(AgeGroup)group];
                    double reduc_vacc_hosp = model[county].parameters.template get<HospitalizedFactorPartiallyImmune>()[(AgeGroup)group];
                    double reduc_vacc_icu  = model[county].parameters.template get<HospitalizedFactorPartiallyImmune>()[(AgeGroup)group];
                    double reduc_vacc_dead = model[county].parameters.template get<HospitalizedFactorPartiallyImmune>()[(AgeGroup)group];
                    mu_C_R[county].push_back(
                        (1 -
                         reduc_vacc_inf / reduc_vacc_exp *
                             (1 -
                              model[county].parameters.template get<AsymptoticCasesPerInfectious>()[(AgeGroup)group])));
                    mu_I_H[county].push_back(
                        reduc_vacc_hosp / reduc_vacc_inf *
                        model[county].parameters.template get<HospitalizedCasesPerInfectious>()[(AgeGroup)group]);
                    // transfer from H to U, D unchanged.
                    mu_H_U[county].push_back(
                        reduc_vacc_icu / reduc_vacc_hosp *
                        model[county].parameters.template get<ICUCasesPerHospitalized>()[(AgeGroup)group]);
                    mu_U_D[county].push_back(reduc_vacc_dead / reduc_vacc_icu *
                                             model[county].parameters.template get<DeathsPerICU>()[(AgeGroup)group]);
                }
            }

            BOOST_OUTCOME_TRY(read_rki_data(rki_data, region, date, num_exp, num_car, num_inf, num_hosp, num_icu,
                                            num_death, num_rec, t_car_to_rec, t_car_to_inf, t_exp_to_car, t_inf_to_rec,
                                            t_inf_to_hosp, t_hosp_to_rec, t_hosp_to_icu, t_icu_to_dead, t_icu_to_rec,
                                            mu_C_R, mu_I_H, mu_H_U, mu_U_D, scaling_factor_inf));

            for (size_t county = 0; county < model.size(); county++) {
                // if (std::accumulate(num_inf[county].begin(), num_inf[county].end(), 0.0) > 0) {
                size_t num_groups = (size_t)model[county].parameters.get_num_groups();
                for (size_t i = 0; i < num_groups; i++) {
                    model[county].populations[{AgeGroup(i), InfectionState::ExposedPartiallyImmune}] =
                        num_exp[county][i];
                    model[county].populations[{AgeGroup(i), InfectionState::CarrierPartiallyImmune}] =
                        num_car[county][i];
                    model[county].populations[{AgeGroup(i), InfectionState::InfectedPartiallyImmune}] =
                        num_inf[county][i];
                    model[county].populations[{AgeGroup(i), InfectionState::HospitalizedPartiallyImmune}] =
                        num_hosp[county][i];
                }
                // }
                if (std::accumulate(num_inf[county].begin(), num_inf[county].end(), 0.0) == 0) {
                    log_warning("No infections for partially vaccinated reported on date " + std::to_string(date.year) +
                                "-" + std::to_string(date.month) + "-" + std::to_string(date.day) + " for region " +
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
                    t_car_to_inf[county].push_back(static_cast<int>(
                        2 * (model[county].parameters.template get<IncubationTime>()[(AgeGroup)group] -
                             model[county].parameters.template get<SerialInterval>()[(AgeGroup)group])));
                    t_car_to_rec[county].push_back(static_cast<int>(
                        t_car_to_inf[county][group] +
                        0.5 * model[county].parameters.template get<InfectiousTimeMild>()[(AgeGroup)group] * reduc_t));
                    t_exp_to_car[county].push_back(
                        static_cast<int>(2 * model[county].parameters.template get<SerialInterval>()[(AgeGroup)group] -
                                         model[county].parameters.template get<IncubationTime>()[(AgeGroup)group]));
                    t_inf_to_rec[county].push_back(static_cast<int>(
                        model[county].parameters.template get<InfectiousTimeMild>()[(AgeGroup)group] * reduc_t));
                    t_inf_to_hosp[county].push_back(static_cast<int>(
                        model[county].parameters.template get<HomeToHospitalizedTime>()[(AgeGroup)group]));
                    t_hosp_to_rec[county].push_back(static_cast<int>(
                        model[county].parameters.template get<HospitalizedToHomeTime>()[(AgeGroup)group]));
                    t_hosp_to_icu[county].push_back(static_cast<int>(
                        model[county].parameters.template get<HospitalizedToICUTime>()[(AgeGroup)group]));
                    t_icu_to_dead[county].push_back(
                        static_cast<int>(model[county].parameters.template get<ICUToDeathTime>()[(AgeGroup)group]));
                    t_icu_to_rec[county].push_back(
                        static_cast<int>(model[county].parameters.template get<ICUToHomeTime>()[(AgeGroup)group]));

                    double reduc_immune_exp = model[county].parameters.template get<ExposedFactorFullyImmune>()[(AgeGroup)group];
                    double reduc_immune_inf = model[county].parameters.template get<InfectedFactorFullyImmune>()[(AgeGroup)group];
                    double reduc_immune_hosp =
                        model[county].parameters.template get<HospitalizedFactorFullyImmune>()[(AgeGroup)group];
                    double reduc_immune_icu = model[county].parameters.template get<HospitalizedFactorFullyImmune>()[(AgeGroup)group];
                    double reduc_immune_dead =
                        model[county].parameters.template get<HospitalizedFactorFullyImmune>()[(AgeGroup)group];
                    mu_C_R[county].push_back(
                        (1 -
                         reduc_immune_inf / reduc_immune_exp *
                             (1 -
                              model[county].parameters.template get<AsymptoticCasesPerInfectious>()[(AgeGroup)group])));
                    mu_I_H[county].push_back(
                        reduc_immune_hosp / reduc_immune_inf *
                        model[county].parameters.template get<HospitalizedCasesPerInfectious>()[(AgeGroup)group]);
                    // transfer from H to U, D unchanged.
                    mu_H_U[county].push_back(
                        reduc_immune_icu / reduc_immune_hosp *
                        model[county].parameters.template get<ICUCasesPerHospitalized>()[(AgeGroup)group]);
                    mu_U_D[county].push_back(reduc_immune_dead / reduc_immune_icu *
                                             model[county].parameters.template get<DeathsPerICU>()[(AgeGroup)group]);
                }
            }

            BOOST_OUTCOME_TRY(read_rki_data(rki_data, region, date, num_exp, num_car, num_inf, num_hosp, num_icu,
                                            num_death, num_rec, t_car_to_rec, t_car_to_inf, t_exp_to_car, t_inf_to_rec,
                                            t_inf_to_hosp, t_hosp_to_rec, t_hosp_to_icu, t_icu_to_dead, t_icu_to_rec,
                                            mu_C_R, mu_I_H, mu_H_U, mu_U_D, scaling_factor_inf));

            for (size_t county = 0; county < model.size(); county++) {
                // if (std::accumulate(num_inf[county].begin(), num_inf[county].end(), 0.0) > 0) {
                size_t num_groups = (size_t)model[county].parameters.get_num_groups();
                for (size_t i = 0; i < num_groups; i++) {
                    model[county].populations[{AgeGroup(i), InfectionState::ExposedFullyImmune}]  = num_exp[county][i];
                    model[county].populations[{AgeGroup(i), InfectionState::CarrierFullyImmune}]  = num_car[county][i];
                    model[county].populations[{AgeGroup(i), InfectionState::InfectedFullyImmune}] = num_inf[county][i];
                    model[county].populations[{AgeGroup(i), InfectionState::HospitalizedFullyImmune}] =
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
        * @brief reads number of ICU patients from DIVI register into SecirParams
        * @param path Path to DIVI file
        * @param vregion Keys of the region of interest
        * @param date Date for which the arrays are initialized
        * @param vnum_icu number of ICU patients
        * @see mio::read_divi_data
        * @{
        */
        IOResult<void> read_divi_data(const std::string& path, const std::vector<int>& vregion, Date date,
                                      std::vector<double>& vnum_icu);
        IOResult<void> read_divi_data(const std::vector<DiviEntry>& divi_data, const std::vector<int>& vregion,
                                      Date date, std::vector<double>& vnum_icu);
        /**@}*/

        /**
        * @brief sets populations data from DIVI register into Model
        * @param model vector of objects in which the data is set
        * @param path Path to DIVI file
        * @param vregion vector of keys of the regions of interest
        * @param date Date for which the arrays are initialized
        * @param scaling_factor_icu factor by which to scale the icu cases of divi data
        */
        template <class Model>
        IOResult<void> set_divi_data(std::vector<Model>& model, const std::string& path,
                                     const std::vector<int>& vregion, Date date, double scaling_factor_icu)
        {
            std::vector<double> sum_mu_I_U(vregion.size(), 0);
            std::vector<std::vector<double>> mu_I_U{model.size()};
            for (size_t region = 0; region < vregion.size(); region++) {
                auto num_groups = model[region].parameters.get_num_groups();
                for (auto i = AgeGroup(0); i < num_groups; i++) {
                    sum_mu_I_U[region] += model[region].parameters.template get<ICUCasesPerHospitalized>()[i] *
                                          model[region].parameters.template get<HospitalizedCasesPerInfectious>()[i];
                    mu_I_U[region].push_back(
                        model[region].parameters.template get<ICUCasesPerHospitalized>()[i] *
                        model[region].parameters.template get<HospitalizedCasesPerInfectious>()[i]);
                }
            }
            std::vector<double> num_icu(model.size(), 0.0);
            BOOST_OUTCOME_TRY(read_divi_data(path, vregion, date, num_icu));

            for (size_t region = 0; region < vregion.size(); region++) {
                auto num_groups = model[region].parameters.get_num_groups();
                for (auto i = AgeGroup(0); i < num_groups; i++) {
                    model[region].populations[{i, InfectionState::ICU}] =
                        scaling_factor_icu * num_icu[region] * mu_I_U[region][(size_t)i] / sum_mu_I_U[region];
                }
            }

            return success();
        }

        /**
        * @brief reads population data from census data
        * @param path Path to RKI file
        * @param vregion vector of keys of the regions of interest
        * @see mio::read_population_data
        * @{
        */
        IOResult<std::vector<std::vector<double>>> read_population_data(const std::string& path,
                                                                        const std::vector<int>& vregion);
        IOResult<std::vector<std::vector<double>>>
        read_population_data(const std::vector<PopulationDataEntry>& population_data, const std::vector<int>& vregion);
        /**@}*/

        template <class Model>
        IOResult<void> set_population_data(std::vector<Model>& model, const std::string& path,
                                           const std::string& path_rki, const std::vector<int>& vregion, Date date)
        {
            BOOST_OUTCOME_TRY(num_population, read_population_data(path, vregion));

            auto num_age_groups = StringRkiAgeGroup::age_group_names.size();
            std::vector<std::vector<double>> num_rec(model.size(), std::vector<double>(num_age_groups, 0.0));

            BOOST_OUTCOME_TRY(read_rki_data_confirmed_to_recovered(path_rki, vregion, date, num_rec, 14.));

            for (size_t region = 0; region < vregion.size(); region++) {
                if (std::accumulate(num_population[region].begin(), num_population[region].end(), 0.0) > 0) {
                    auto num_groups = model[region].parameters.get_num_groups();
                    for (auto i = AgeGroup(0); i < num_groups; i++) {

                        double S_v  = std::min(model[region].parameters.template get<DailyFullVaccination>()[{i, SimulationDay(0)}] +
                                                  num_rec[region][size_t(i)],
                                              num_population[region][size_t(i)]);
                        double S_pv = std::max(model[region].parameters.template get<DailyFirstVaccination>()[{i, SimulationDay(0)}] -
                                                   model[region].parameters.template get<DailyFullVaccination>()[{i, SimulationDay(0)}],
                                               0.0); // use std::max with 0
                        double S;
                        if (num_population[region][size_t(i)] - S_pv - S_v < 0.0) {
                            log_warning(
                                "Number of vaccinated persons greater than population in county {}, age group {}.",
                                region, size_t(i));
                            S   = 0.0;
                            S_v = num_population[region][size_t(i)] - S_pv;
                        }
                        else {
                            S = num_population[region][size_t(i)] - S_pv - S_v;
                        }

                        double denom_E  = 1 / (S + S_pv * model[region].parameters.template get<ExposedFactorPartiallyImmune>()[i] +
                                              S_v * model[region].parameters.template get<ExposedFactorFullyImmune>()[i]);
                        double denom_C  = 1 / (S + S_pv + S_v);
                        double denom_I  = 1 / (S + S_pv * model[region].parameters.template get<InfectedFactorPartiallyImmune>()[i] +
                                              S_v * model[region].parameters.template get<InfectedFactorFullyImmune>()[i]);
                        double denom_HU = 1 / (S + S_pv * model[region].parameters.template get<HospitalizedFactorPartiallyImmune>()[i] +
                                               S_v * model[region].parameters.template get<HospitalizedFactorFullyImmune>()[i]);

                        model[region].populations[{i, InfectionState::Exposed}] =
                            S * model[region].populations[{i, InfectionState::Exposed}] * denom_E;
                        model[region].populations[{i, InfectionState::ExposedPartiallyImmune}] =
                            S_pv * model[region].parameters.template get<ExposedFactorPartiallyImmune>()[i] *
                            model[region].populations[{i, InfectionState::ExposedPartiallyImmune}] * denom_E;
                        model[region].populations[{i, InfectionState::ExposedFullyImmune}] =
                            S_v * model[region].parameters.template get<ExposedFactorFullyImmune>()[i] *
                            model[region].populations[{i, InfectionState::ExposedFullyImmune}] * denom_E;

                        model[region].populations[{i, InfectionState::Carrier}] =
                            S * model[region].populations[{i, InfectionState::Carrier}] * denom_C;
                        model[region].populations[{i, InfectionState::CarrierPartiallyImmune}] =
                            S_pv * model[region].populations[{i, InfectionState::CarrierPartiallyImmune}] * denom_C;
                        model[region].populations[{i, InfectionState::CarrierFullyImmune}] =
                            S_v * model[region].populations[{i, InfectionState::CarrierFullyImmune}] * denom_C;

                        model[region].populations[{i, InfectionState::Infected}] =
                            S * model[region].populations[{i, InfectionState::Infected}] * denom_I;
                        model[region].populations[{i, InfectionState::InfectedPartiallyImmune}] =
                            S_pv * model[region].parameters.template get<InfectedFactorPartiallyImmune>()[i] *
                            model[region].populations[{i, InfectionState::InfectedPartiallyImmune}] * denom_I;
                        model[region].populations[{i, InfectionState::InfectedFullyImmune}] =
                            S_v * model[region].parameters.template get<InfectedFactorFullyImmune>()[i] *
                            model[region].populations[{i, InfectionState::InfectedFullyImmune}] * denom_I;

                        model[region].populations[{i, InfectionState::Hospitalized}] =
                            S * model[region].populations[{i, InfectionState::Hospitalized}] * denom_HU;
                        model[region].populations[{i, InfectionState::HospitalizedPartiallyImmune}] =
                            S_pv * model[region].parameters.template get<HospitalizedFactorPartiallyImmune>()[i] *
                            model[region].populations[{i, InfectionState::HospitalizedPartiallyImmune}] * denom_HU;
                        model[region].populations[{i, InfectionState::HospitalizedFullyImmune}] =
                            S_v * model[region].parameters.template get<HospitalizedFactorFullyImmune>()[i] *
                            model[region].populations[{i, InfectionState::HospitalizedFullyImmune}] * denom_HU;

                        model[region].populations[{i, InfectionState::ICU}] =
                            S * model[region].populations[{i, InfectionState::ICU}] * denom_HU;
                        model[region].populations[{i, InfectionState::ICUPartiallyImmune}] =
                            S_pv * model[region].parameters.template get<HospitalizedFactorPartiallyImmune>()[i] *
                            model[region].populations[{i, InfectionState::ICU}] * denom_HU;
                        model[region].populations[{i, InfectionState::ICUFullyImmune}] =
                            S_v * model[region].parameters.template get<HospitalizedFactorFullyImmune>()[i] *
                            model[region].populations[{i, InfectionState::ICU}] * denom_HU;

                        model[region].populations[{i, InfectionState::Recovered}] =
                            model[region].parameters.template get<DailyFullVaccination>()[{i, SimulationDay(0)}] +
                            model[region].populations[{i, InfectionState::Recovered}] -
                            (model[region].populations[{i, InfectionState::Infected}] +
                             model[region].populations[{i, InfectionState::InfectedPartiallyImmune}] +
                             model[region].populations[{i, InfectionState::InfectedFullyImmune}] +
                             model[region].populations[{i, InfectionState::Hospitalized}] +
                             model[region].populations[{i, InfectionState::HospitalizedPartiallyImmune}] +
                             model[region].populations[{i, InfectionState::HospitalizedFullyImmune}] +
                             model[region].populations[{i, InfectionState::ICU}] +
                             model[region].populations[{i, InfectionState::ICUPartiallyImmune}] +
                             model[region].populations[{i, InfectionState::ICUFullyImmune}] +
                             model[region].populations[{i, InfectionState::Dead}]);

                        model[region].populations[{i, InfectionState::Recovered}] =
                            std::min(S + S_pv + S_v,
                                     std::max(0.0, double(model[region].populations[{i, InfectionState::Recovered}])));

                        model[region].populations[{i, InfectionState::SusceptiblePartiallyImmune}] = std::max(
                            0.0, S_pv - model[region].populations[{i, InfectionState::ExposedPartiallyImmune}] -
                                     model[region].populations[{i, InfectionState::CarrierPartiallyImmune}] -
                                     model[region].populations[{i, InfectionState::InfectedPartiallyImmune}] -
                                     model[region].populations[{i, InfectionState::HospitalizedPartiallyImmune}] -
                                     model[region].populations[{i, InfectionState::ICUPartiallyImmune}]);

                        model[region].populations.template set_difference_from_group_total<AgeGroup>(
                            {i, InfectionState::Susceptible}, num_population[region][size_t(i)]);
                    }

                    for (auto i = AgeGroup(0); i < AgeGroup(6); i++) {
                        for (auto j = Index<InfectionState>(0); j < InfectionState::Count; ++j) {
                            if (model[region].populations[{i, j}] < 0) {
                                log_warning("Compartment at age group {}, infection state {}, is negative: {}",
                                            size_t(i), size_t(j),
                                            model[region].populations[{i, j}] / num_population[region][size_t(i)]);
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

        IOResult<void> set_vaccination_data(std::vector<SecirModel>& model, const std::string& path, Date date,
                                            const std::vector<int>& vregion, int num_days);
    } // namespace details

    /**
    * @brief Exports the time series of extrapolated real data according to
    *   the extrapolation / approximation method used to initialize the model from
    *   real world data.
        (This is the vector-valued functionality of set_rki_data())
    * @param model vector of objects in which the data is set
    * @param data_dir Path to RKI files
    * @param results_dir Path to result files
    * @param date Start date of the time series to be exported.
    * @param region vector of keys of the region of interest
    * @param scaling_factor_inf factors by which to scale the confirmed cases of rki data
    * @param scaling_factor_icu factors by which to scale the intensive care data
    * @param num_days Number of days for which the time series is exported.
    */
    template <class Model>
    IOResult<void> export_input_data_county_timeseries(std::vector<Model>& model, const std::string& data_dir,
                                                       const std::string& results_dir, std::vector<int> const& region,
                                                       Date date, const std::vector<double>& scaling_factor_inf,
                                                       double scaling_factor_icu, int num_days)
    {
        auto num_age_groups = (size_t)model[0].parameters.get_num_groups();
        assert(scaling_factor_inf.size() == num_age_groups);
        assert(num_age_groups == StringRkiAgeGroup::age_group_names.size());
        assert(model.size() == region.size());

        BOOST_OUTCOME_TRY(rki_data, read_rki_data(path_join(data_dir, "all_county_age_ma7_rki.json")));
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
        std::vector<std::vector<int>> t_hosp_to_icu_uv{model.size()}; // R7
        std::vector<std::vector<int>> t_icu_to_dead_uv{model.size()}; // R10
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
                t_car_to_rec_uv[county].push_back(static_cast<int>(
                    t_car_to_inf_uv[county][group] +
                    0.5 * model[county].parameters.template get<InfectiousTimeMild>()[(AgeGroup)group]));
                t_exp_to_car_uv[county].push_back(
                    static_cast<int>(2 * model[county].parameters.template get<SerialInterval>()[(AgeGroup)group] -
                                     model[county].parameters.template get<IncubationTime>()[(AgeGroup)group]));
                t_inf_to_rec_uv[county].push_back(
                    static_cast<int>(model[county].parameters.template get<InfectiousTimeMild>()[(AgeGroup)group]));
                t_inf_to_hosp_uv[county].push_back(
                    static_cast<int>(model[county].parameters.template get<HomeToHospitalizedTime>()[(AgeGroup)group]));
                t_hosp_to_rec_uv[county].push_back(
                    static_cast<int>(model[county].parameters.template get<HospitalizedToHomeTime>()[(AgeGroup)group]));
                t_hosp_to_icu_uv[county].push_back(
                    static_cast<int>(model[county].parameters.template get<HospitalizedToICUTime>()[(AgeGroup)group]));
                t_icu_to_dead_uv[county].push_back(
                    static_cast<int>(model[county].parameters.template get<ICUToDeathTime>()[(AgeGroup)group]));
                t_icu_to_rec_uv[county].push_back(
                    static_cast<int>(model[county].parameters.template get<ICUToHomeTime>()[(AgeGroup)group]));

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
        std::vector<std::vector<int>> t_hosp_to_icu_pv{model.size()}; // R7
        std::vector<std::vector<int>> t_icu_to_dead_pv{model.size()}; // R10
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
                    0.5 * model[county].parameters.template get<InfectiousTimeMild>()[(AgeGroup)group] * reduc_t));
                t_exp_to_car_pv[county].push_back(
                    static_cast<int>(2 * model[county].parameters.template get<SerialInterval>()[(AgeGroup)group] -
                                     model[county].parameters.template get<IncubationTime>()[(AgeGroup)group]));
                t_inf_to_rec_pv[county].push_back(static_cast<int>(
                    model[county].parameters.template get<InfectiousTimeMild>()[(AgeGroup)group] * reduc_t));
                t_inf_to_hosp_pv[county].push_back(
                    static_cast<int>(model[county].parameters.template get<HomeToHospitalizedTime>()[(AgeGroup)group]));
                t_hosp_to_rec_pv[county].push_back(
                    static_cast<int>(model[county].parameters.template get<HospitalizedToHomeTime>()[(AgeGroup)group]));
                t_hosp_to_icu_pv[county].push_back(
                    static_cast<int>(model[county].parameters.template get<HospitalizedToICUTime>()[(AgeGroup)group]));
                t_icu_to_dead_pv[county].push_back(
                    static_cast<int>(model[county].parameters.template get<ICUToDeathTime>()[(AgeGroup)group]));
                t_icu_to_rec_pv[county].push_back(
                    static_cast<int>(model[county].parameters.template get<ICUToHomeTime>()[(AgeGroup)group]));

                double reduc_vacc_exp  = model[county].parameters.template get<ExposedFactorPartiallyImmune>()[(AgeGroup)group];
                double reduc_vacc_inf  = model[county].parameters.template get<InfectedFactorPartiallyImmune>()[(AgeGroup)group];
                double reduc_vacc_hosp = model[county].parameters.template get<HospitalizedFactorPartiallyImmune>()[(AgeGroup)group];
                double reduc_vacc_icu  = model[county].parameters.template get<HospitalizedFactorPartiallyImmune>()[(AgeGroup)group];
                double reduc_vacc_dead = model[county].parameters.template get<HospitalizedFactorPartiallyImmune>()[(AgeGroup)group];
                mu_C_R_pv[county].push_back(
                    (1 -
                     reduc_vacc_inf / reduc_vacc_exp *
                         (1 - model[county].parameters.template get<AsymptoticCasesPerInfectious>()[(AgeGroup)group])));
                mu_I_H_pv[county].push_back(
                    reduc_vacc_hosp / reduc_vacc_inf *
                    model[county].parameters.template get<HospitalizedCasesPerInfectious>()[(AgeGroup)group]);
                // transfer from H to U, D unchanged.
                mu_H_U_pv[county].push_back(
                    reduc_vacc_icu / reduc_vacc_hosp *
                    model[county].parameters.template get<ICUCasesPerHospitalized>()[(AgeGroup)group]);
                mu_U_D_pv[county].push_back(reduc_vacc_dead / reduc_vacc_icu *
                                            model[county].parameters.template get<DeathsPerICU>()[(AgeGroup)group]);

                /* begin: NOT in set_rki_data() */
                sum_mu_I_U_pv[county] +=
                    reduc_vacc_icu / reduc_vacc_hosp *
                    model[county].parameters.template get<ICUCasesPerHospitalized>()[AgeGroup(group)] *
                    reduc_vacc_hosp / reduc_vacc_inf *
                    model[county].parameters.template get<HospitalizedCasesPerInfectious>()[AgeGroup(group)];
                mu_I_U_pv[county].push_back(
                    reduc_vacc_icu / reduc_vacc_hosp *
                    model[county].parameters.template get<ICUCasesPerHospitalized>()[AgeGroup(group)] *
                    reduc_vacc_hosp / reduc_vacc_inf *
                    model[county].parameters.template get<HospitalizedCasesPerInfectious>()[AgeGroup(group)]);
                /* end: NOT in set_rki_data() */
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
        std::vector<std::vector<int>> t_hosp_to_icu_fv{model.size()}; // R7
        std::vector<std::vector<int>> t_icu_to_dead_fv{model.size()}; // R10
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
                    0.5 * model[county].parameters.template get<InfectiousTimeMild>()[(AgeGroup)group] * reduc_t));
                t_exp_to_car_fv[county].push_back(
                    static_cast<int>(2 * model[county].parameters.template get<SerialInterval>()[(AgeGroup)group] -
                                     model[county].parameters.template get<IncubationTime>()[(AgeGroup)group]));
                t_inf_to_rec_fv[county].push_back(static_cast<int>(
                    model[county].parameters.template get<InfectiousTimeMild>()[(AgeGroup)group] * reduc_t));
                t_inf_to_hosp_fv[county].push_back(
                    static_cast<int>(model[county].parameters.template get<HomeToHospitalizedTime>()[(AgeGroup)group]));
                t_hosp_to_rec_fv[county].push_back(
                    static_cast<int>(model[county].parameters.template get<HospitalizedToHomeTime>()[(AgeGroup)group]));
                t_hosp_to_icu_fv[county].push_back(
                    static_cast<int>(model[county].parameters.template get<HospitalizedToICUTime>()[(AgeGroup)group]));
                t_icu_to_dead_fv[county].push_back(
                    static_cast<int>(model[county].parameters.template get<ICUToDeathTime>()[(AgeGroup)group]));
                t_icu_to_rec_fv[county].push_back(
                    static_cast<int>(model[county].parameters.template get<ICUToHomeTime>()[(AgeGroup)group]));

                double reduc_immune_exp  = model[county].parameters.template get<ExposedFactorFullyImmune>()[(AgeGroup)group];
                double reduc_immune_inf  = model[county].parameters.template get<InfectedFactorFullyImmune>()[(AgeGroup)group];
                double reduc_immune_hosp = model[county].parameters.template get<HospitalizedFactorFullyImmune>()[(AgeGroup)group];
                double reduc_immune_icu  = model[county].parameters.template get<HospitalizedFactorFullyImmune>()[(AgeGroup)group];
                double reduc_immune_dead = model[county].parameters.template get<HospitalizedFactorFullyImmune>()[(AgeGroup)group];
                mu_C_R_fv[county].push_back(
                    (1 -
                     reduc_immune_inf / reduc_immune_exp *
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

                /* begin: NOT in set_rki_data() */
                sum_mu_I_U_fv[county] +=
                    reduc_immune_icu / reduc_immune_hosp *
                    model[county].parameters.template get<ICUCasesPerHospitalized>()[AgeGroup(group)] *
                    reduc_immune_hosp / reduc_immune_inf *
                    model[county].parameters.template get<HospitalizedCasesPerInfectious>()[AgeGroup(group)];
                mu_I_U_fv[county].push_back(
                    reduc_immune_icu / reduc_immune_hosp *
                    model[county].parameters.template get<ICUCasesPerHospitalized>()[AgeGroup(group)] *
                    reduc_immune_hosp / reduc_immune_inf *
                    model[county].parameters.template get<HospitalizedCasesPerInfectious>()[AgeGroup(group)]);
                /* end: NOT in set_rki_data() */
            }
        }
        // TODO: Restart here!
        /* begin: similar functionality in set_rki_data(), here only for vector of TimeSeries */
        /* object and with additions for read_divi and read_population... */
        std::vector<TimeSeries<double>> extrapolated_rki(
            model.size(), TimeSeries<double>::zero(num_days, (size_t)InfectionState::Count * num_age_groups));

        for (size_t day = 0; day < static_cast<size_t>(num_days); day++) {

            // unvaccinated
            std::vector<std::vector<double>> num_exp_uv(model.size(), std::vector<double>(num_age_groups, 0.0));
            std::vector<std::vector<double>> num_car_uv(model.size(), std::vector<double>(num_age_groups, 0.0));
            std::vector<std::vector<double>> num_inf_uv(model.size(), std::vector<double>(num_age_groups, 0.0));
            std::vector<std::vector<double>> num_rec_uv(model.size(), std::vector<double>(num_age_groups, 0.0));
            std::vector<std::vector<double>> num_hosp_uv(model.size(), std::vector<double>(num_age_groups, 0.0));
            std::vector<std::vector<double>> num_death_uv(model.size(), std::vector<double>(num_age_groups, 0.0));
            std::vector<std::vector<double>> dummy_icu(model.size(), std::vector<double>(num_age_groups, 0.0));
            BOOST_OUTCOME_TRY(details::read_rki_data(
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
            BOOST_OUTCOME_TRY(details::read_rki_data(
                rki_data, region, date, num_exp_pv, num_car_pv, num_inf_pv, num_hosp_pv, dummy_icu, dummy_death,
                dummy_rec, t_car_to_rec_pv, t_car_to_inf_pv, t_exp_to_car_pv, t_inf_to_rec_pv, t_inf_to_hosp_pv,
                t_hosp_to_rec_pv, t_hosp_to_icu_pv, t_icu_to_dead_pv, t_icu_to_rec_pv, mu_C_R_pv, mu_I_H_pv, mu_H_U_pv,
                mu_U_D_pv, scaling_factor_inf));

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
            BOOST_OUTCOME_TRY(details::read_rki_data(
                rki_data, region, date, num_exp_fv, num_car_fv, num_inf_fv, num_hosp_fv, dummy_icu, dummy_death,
                dummy_rec, t_car_to_rec_fv, t_car_to_inf_fv, t_exp_to_car_fv, t_inf_to_rec_fv, t_inf_to_hosp_fv,
                t_hosp_to_rec_fv, t_hosp_to_icu_fv, t_icu_to_dead_fv, t_icu_to_rec_fv, mu_C_R_fv, mu_I_H_fv, mu_H_U_fv,
                mu_U_D_fv, scaling_factor_inf));

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
            BOOST_OUTCOME_TRY(details::read_rki_data_confirmed_to_recovered(rki_data, region, date, num_rec, 14.));

            for (size_t county = 0; county < region.size(); county++) {
                if (std::accumulate(num_population[county].begin(), num_population[county].end(), 0.0) > 0) {
                    for (size_t age = 0; age < num_age_groups; age++) {
                        auto age_group_offset = age * (size_t)InfectionState::Count;
                        double S_v =
                            std::min(model[county].parameters.template get<DailyFullVaccination>()[{AgeGroup(age), SimulationDay(day)}] +
                                         num_rec[county][age],
                                     num_population[county][age]);
                        double S_pv = std::max(
                            model[county].parameters.template get<DailyFirstVaccination>()[{AgeGroup(age), SimulationDay(day)}] -
                                model[county].parameters.template get<DailyFullVaccination>()[{AgeGroup(age), SimulationDay(day)}],
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
                            1 / (S + S_pv * model[county].parameters.template get<ExposedFactorPartiallyImmune>()[AgeGroup(age)] +
                                 S_v * model[county].parameters.template get<ExposedFactorFullyImmune>()[AgeGroup(age)]);
                        double denom_C = 1 / (S + S_pv + S_v);
                        double denom_I =
                            1 / (S + S_pv * model[county].parameters.template get<InfectedFactorPartiallyImmune>()[AgeGroup(age)] +
                                 S_v * model[county].parameters.template get<InfectedFactorFullyImmune>()[AgeGroup(age)]);
                        double denom_HU =
                            1 / (S + S_pv * model[county].parameters.template get<HospitalizedFactorPartiallyImmune>()[AgeGroup(age)] +
                                 S_v * model[county].parameters.template get<HospitalizedFactorFullyImmune>()[AgeGroup(age)]);

                        extrapolated_rki[county][day]((size_t)InfectionState::Exposed + age_group_offset) =
                            S * denom_E * num_exp_uv[county][age];
                        extrapolated_rki[county][day]((size_t)InfectionState::ExposedPartiallyImmune +
                                                      age_group_offset) =
                            S_pv * model[county].parameters.template get<ExposedFactorPartiallyImmune>()[AgeGroup(age)] * denom_E *
                            num_exp_pv[county][age];
                        extrapolated_rki[county][day]((size_t)InfectionState::ExposedFullyImmune + age_group_offset) =
                            S_v * model[county].parameters.template get<ExposedFactorFullyImmune>()[AgeGroup(age)] * denom_E *
                            num_exp_fv[county][age];

                        extrapolated_rki[county][day]((size_t)InfectionState::Carrier + age_group_offset) =
                            S * denom_C * num_car_uv[county][age];
                        extrapolated_rki[county][day]((size_t)InfectionState::CarrierPartiallyImmune +
                                                      age_group_offset) = S_pv * denom_C * num_car_pv[county][age];
                        extrapolated_rki[county][day]((size_t)InfectionState::CarrierFullyImmune + age_group_offset) =
                            S_v * denom_C * num_car_fv[county][age];

                        extrapolated_rki[county][day]((size_t)InfectionState::Infected + age_group_offset) =
                            S * denom_I * num_inf_uv[county][age];
                        extrapolated_rki[county][day]((size_t)InfectionState::InfectedPartiallyImmune +
                                                      age_group_offset) =
                            S_pv * model[county].parameters.template get<InfectedFactorPartiallyImmune>()[AgeGroup(age)] * denom_I *
                            num_inf_pv[county][age];
                        extrapolated_rki[county][day]((size_t)InfectionState::InfectedFullyImmune + age_group_offset) =
                            S_v * model[county].parameters.template get<InfectedFactorFullyImmune>()[AgeGroup(age)] * denom_I *
                            num_inf_fv[county][age];

                        extrapolated_rki[county][day]((size_t)InfectionState::Hospitalized + age_group_offset) =
                            S * denom_HU * num_hosp_uv[county][age];
                        extrapolated_rki[county][day]((size_t)InfectionState::HospitalizedPartiallyImmune +
                                                      age_group_offset) =
                            S_pv * model[county].parameters.template get<HospitalizedFactorPartiallyImmune>()[AgeGroup(age)] * denom_HU *
                            num_hosp_pv[county][age];
                        extrapolated_rki[county][day]((size_t)InfectionState::HospitalizedFullyImmune +
                                                      age_group_offset) =
                            S_v * model[county].parameters.template get<HospitalizedFactorFullyImmune>()[AgeGroup(age)] * denom_HU *
                            num_hosp_fv[county][age];

                        extrapolated_rki[county][day]((size_t)InfectionState::ICU + age_group_offset) =
                            S * denom_HU * num_icu[county][age];
                        extrapolated_rki[county][day]((size_t)InfectionState::ICUPartiallyImmune + age_group_offset) =
                            S_pv * model[county].parameters.template get<HospitalizedFactorPartiallyImmune>()[AgeGroup(age)] * denom_HU *
                            num_icu[county][age];
                        extrapolated_rki[county][day]((size_t)InfectionState::ICUFullyImmune + age_group_offset) =
                            S_v * model[county].parameters.template get<HospitalizedFactorFullyImmune>()[AgeGroup(age)] * denom_HU *
                            num_icu[county][age];

                        // in set_rki_data initilization, deaths are now set to 0. In order to visualize
                        // the extrapolated real number of deaths, they have to be set here. In the comparison of data
                        // it has to be paid attention to the fact, the the simulation starts with deaths=0
                        // while this method starts with deaths=number of reported deaths so far...
                        extrapolated_rki[county][day]((size_t)InfectionState::Dead + age_group_offset) =
                            num_death_uv[county][age];

                        extrapolated_rki[county][day]((size_t)InfectionState::Recovered + age_group_offset) =
                            model[county].parameters.template get<DailyFullVaccination>()[{AgeGroup(age), SimulationDay(day)}] +
                            num_rec_uv[county][age] -
                            (extrapolated_rki[county][day]((size_t)InfectionState::Infected + age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::InfectedPartiallyImmune +
                                                           age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::InfectedFullyImmune +
                                                           age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::Hospitalized + age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::HospitalizedPartiallyImmune +
                                                           age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::HospitalizedFullyImmune +
                                                           age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::ICU + age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::ICUPartiallyImmune +
                                                           age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::ICUFullyImmune + age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::Dead + age_group_offset));

                        extrapolated_rki[county][day]((size_t)InfectionState::Recovered + age_group_offset) = std::min(
                            S + S_pv + S_v, std::max(0.0, double(extrapolated_rki[county][day](
                                                              (size_t)InfectionState::Recovered + age_group_offset))));

                        extrapolated_rki[county][day]((size_t)InfectionState::SusceptiblePartiallyImmune +
                                                      age_group_offset) =
                            std::max(0.0,
                                     S_pv -
                                         extrapolated_rki[county][day]((size_t)InfectionState::ExposedPartiallyImmune +
                                                                       age_group_offset) -
                                         extrapolated_rki[county][day]((size_t)InfectionState::CarrierPartiallyImmune +
                                                                       age_group_offset) -
                                         extrapolated_rki[county][day]((size_t)InfectionState::InfectedPartiallyImmune +
                                                                       age_group_offset) -
                                         extrapolated_rki[county][day](
                                             (size_t)InfectionState::HospitalizedPartiallyImmune + age_group_offset) -
                                         extrapolated_rki[county][day]((size_t)InfectionState::ICUPartiallyImmune +
                                                                       age_group_offset));

                        extrapolated_rki[county][day]((size_t)InfectionState::Susceptible + age_group_offset) =
                            num_population[county][age] -
                            (extrapolated_rki[county][day]((size_t)InfectionState::SusceptiblePartiallyImmune +
                                                           age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::Recovered + age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::Exposed + age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::ExposedPartiallyImmune +
                                                           age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::ExposedFullyImmune +
                                                           age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::Carrier + age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::CarrierPartiallyImmune +
                                                           age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::CarrierFullyImmune +
                                                           age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::Infected + age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::InfectedPartiallyImmune +
                                                           age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::InfectedFullyImmune +
                                                           age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::Hospitalized + age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::HospitalizedPartiallyImmune +
                                                           age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::HospitalizedFullyImmune +
                                                           age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::ICU + age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::ICUPartiallyImmune +
                                                           age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::ICUFullyImmune + age_group_offset) +
                             extrapolated_rki[county][day]((size_t)InfectionState::Dead + age_group_offset));
                    }
                }
                else {
                    log_warning("No population data available for region " + std::to_string(county) +
                                ". Population data has not been set.");
                }
            }
            log_info("extrapolated real data for date: {}-{}-{}", date.day, date.month, date.year);
            date = offset_date_by_days(date, 1);
        }
        /* end: similar functionality in set_rki_data(), here only for vector of TimeSeries */
        auto num_groups = (int)(size_t)model[0].parameters.get_num_groups();
        BOOST_OUTCOME_TRY(save_result(extrapolated_rki, region, num_groups, path_join(results_dir, "Results_rki.h5")));

        auto extrapolated_rki_data_sum = sum_nodes(std::vector<std::vector<TimeSeries<double>>>{extrapolated_rki});
        BOOST_OUTCOME_TRY(save_result({extrapolated_rki_data_sum[0][0]}, {0}, num_groups,
                                      path_join(results_dir, "Results_rki_sum.h5")));

        return success();
    }

    template <class Model>
    IOResult<void> read_input_data_county(std::vector<Model>& model, Date date, const std::vector<int>& county,
                                          const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                          const std::string& dir, int num_days)
    {
        BOOST_OUTCOME_TRY(details::set_vaccination_data(model, path_join(dir, "all_county_ageinf_vacc_ma7.json"), date,
                                                        county, num_days));

        // TODO: set_divi_data and set_rki_data to be merged! Possibly also with set_population_data.
        // set_divi_data and a potential set_divi_data_vaccmodel only need a different ModelType (InfectionState vs InfectionState)
        if (date > Date(2020, 4, 23)) {
            BOOST_OUTCOME_TRY(details::set_divi_data(model, path_join(dir, "county_divi_ma7.json"), county, date,
                                                     scaling_factor_icu));
        }
        else {
            log_warning("No DIVI data available for this date");
        }

        BOOST_OUTCOME_TRY(details::set_rki_data(model, path_join(dir, "all_county_age_ma7_rki.json"), county, date,
                                                scaling_factor_inf));
        BOOST_OUTCOME_TRY(details::set_population_data(model, path_join(dir, "county_current_population.json"),
                                                       path_join(dir, "all_county_age_ma7_rki.json"), county, date));

        // Use only if extrapolated real data is needed for comparison. EXPENSIVE !
        // Run time equals run time of the previous functions times the num_days !
        // (This only represents the vectorization of the previous function over all simulation days...)
        // BOOST_OUTCOME_TRY(export_input_data_county_timeseries(model, dir, dir, county, date, scaling_factor_inf,
        //                                                       scaling_factor_icu, num_days));

        return success();
    }

} // namespace secirv
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // SECIRV_PARAMETERS_IO_H
