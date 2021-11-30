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
#include <epidemiology_io/secir_parameters_io.h>
#include <epidemiology_io/secir_result_io.h>
#include <epidemiology_io/io.h>
#include <epidemiology/utils/memory.h>
#include <epidemiology/utils/uncertain_value.h>
#include <epidemiology/utils/stl_util.h>
#include <epidemiology/utils/graph.h>
#include <epidemiology/migration/migration.h>
#include <epidemiology/secir/damping.h>
#include <epidemiology/model/populations.h>
#include <epidemiology/secir/uncertain_matrix.h>
#include <epidemiology/secir/secir.h>
#include <epidemiology/utils/compiler_diagnostics.h>
#include <epidemiology/utils/date.h>

#include <json/json.h>

#include <boost/filesystem.hpp>

#include <numeric>
#include <vector>
#include <iostream>
#include <string>
#include <random>
#include <fstream>

namespace epi
{

namespace details
{
    IOResult<Date> get_max_date(const Json::Value& root)
    {
        Date max_date = Date(0, 1, 1);
        for (unsigned int i = 0; i < root.size(); i++) {
            auto js_date = root[i]["Date"];
            if (!js_date.isString()) {
                log_error("Date element must be a string.");
                return failure(StatusCode::InvalidType, "Date element must be a string.");
            }
            auto date_temp = parse_date(js_date.asString());
            if (date_temp > max_date) {
                max_date = date_temp;
            }
        }

        return success(max_date);
    }

    void interpolate_ages(const std::vector<double>& age_ranges, std::vector<std::vector<double>>& interpolation,
                          std::vector<bool>& carry_over)
    {
        std::vector<double> param_ranges = {5., 10., 20., 25., 20., 20.};

        //counter for parameter age groups
        size_t counter = 0;

        //residual of param age groups
        double res = 0.0;
        for (size_t i = 0; i < age_ranges.size(); i++) {

            // if current param age group didn't fit into previous rki age group, transfer residual to current age group
            if (res < 0) {
                interpolation[i].push_back(std::min(-res / age_ranges[i], 1.0));
            }

            if (counter < param_ranges.size() - 1) {
                res += age_ranges[i];
                if (std::abs(res) < age_ranges[i]) {
                    counter++;
                }
                // iterate over param age groups while there is still room in the current rki age group
                while (res > 0) {
                    res -= param_ranges[counter];
                    interpolation[i].push_back((param_ranges[counter] + std::min(res, 0.0)) / age_ranges[i]);
                    if (res >= 0) {
                        counter++;
                    }
                }
                if (res < 0) {
                    carry_over.push_back(true);
                }
                else if (res == 0) {
                    carry_over.push_back(false);
                }
            }
            // if last param age group is reached
            else {
                interpolation[i].push_back((age_ranges[i] + res) / age_ranges[i]);
                if (res < 0 || counter == 0) {
                    carry_over.push_back(true);
                }
                else if (res == 0) {
                    carry_over.push_back(false);
                }
                res = 0;
            }
        }
        // last entries for "unknown" age group
        interpolation.push_back({1.0});
        carry_over.push_back(true);
    }

    IOResult<void> read_rki_data(
        std::string const& path, const std::string& id_name, std::vector<int> const& vregion, Date date,
        std::vector<std::vector<double>>& vnum_exp, std::vector<std::vector<double>>& vnum_car,
        std::vector<std::vector<double>>& vnum_inf, std::vector<std::vector<double>>& vnum_hosp,
        std::vector<std::vector<double>>& vnum_icu, std::vector<std::vector<double>>& vnum_death,
        std::vector<std::vector<double>>& vnum_rec, const std::vector<std::vector<int>>& vt_car_to_rec,
        const std::vector<std::vector<int>>& vt_car_to_inf, const std::vector<std::vector<int>>& vt_exp_to_car,
        const std::vector<std::vector<int>>& vt_inf_to_rec, const std::vector<std::vector<int>>& vt_inf_to_hosp,
        const std::vector<std::vector<int>>& vt_hosp_to_rec, const std::vector<std::vector<int>>& vt_hosp_to_icu,
        const std::vector<std::vector<int>>& vt_icu_to_dead, const std::vector<std::vector<double>>& vmu_C_R,
        const std::vector<std::vector<double>>& vmu_I_H, const std::vector<std::vector<double>>& vmu_H_U,
        const std::vector<double>& scaling_factor_inf)
    {
        if (!boost::filesystem::exists(path)) {
            log_error("RKI data file not found: {}.", path);
            return failure(StatusCode::FileNotFound, path);
        }

        Json::Reader reader;
        Json::Value root;

        std::ifstream rki(path);
        if (!reader.parse(rki, root)) {
            log_error(reader.getFormattedErrorMessages());
            return failure(StatusCode::UnknownError, path + ", " + reader.getFormattedErrorMessages());
        }

        BOOST_OUTCOME_TRY(max_date, get_max_date(root));
        if (max_date == Date(0, 1, 1)) {
            log_error("RKI data file is empty.");
            return failure(StatusCode::InvalidFileFormat, path + ", file is empty.");
        }
        if (max_date < date) {
            log_error("Specified date does not exist in RKI data");
            return failure(StatusCode::OutOfRange, path + ", specified date does not exist in RKI data.");
        }
        auto days_surplus = get_offset_in_days(max_date, date) - 6;

        if (days_surplus > 0) {
            days_surplus = 0;
        }

        std::vector<std::string> age_names = {"A00-A04", "A05-A14", "A15-A34", "A35-A59", "A60-A79", "A80+", "unknown"};
        std::vector<double> age_ranges     = {5., 10., 20., 25., 20., 20.};

        for (unsigned int i = 0; i < root.size(); i++) {
            auto it = std::find_if(vregion.begin(), vregion.end(), [&root, i, &id_name](auto r) {
                return r == 0 || root[i][id_name].asInt() == r;
            });
            if (it != vregion.end()) {
                auto region_idx = size_t(it - vregion.begin());

                auto& t_exp_to_car  = vt_exp_to_car[region_idx];
                auto& t_car_to_rec  = vt_car_to_rec[region_idx];
                auto& t_car_to_inf  = vt_car_to_inf[region_idx];
                auto& t_inf_to_rec  = vt_inf_to_rec[region_idx];
                auto& t_inf_to_hosp = vt_inf_to_hosp[region_idx];
                auto& t_hosp_to_rec = vt_hosp_to_rec[region_idx];
                auto& t_hosp_to_icu = vt_hosp_to_icu[region_idx];
                auto& t_icu_to_dead = vt_icu_to_dead[region_idx];

                auto& num_car   = vnum_car[region_idx];
                auto& num_inf   = vnum_inf[region_idx];
                auto& num_rec   = vnum_rec[region_idx];
                auto& num_exp   = vnum_exp[region_idx];
                auto& num_hosp  = vnum_hosp[region_idx];
                auto& num_death = vnum_death[region_idx];
                auto& num_icu   = vnum_icu[region_idx];

                auto& mu_C_R = vmu_C_R[region_idx];
                auto& mu_I_H = vmu_I_H[region_idx];
                auto& mu_H_U = vmu_H_U[region_idx];

                auto date_df = parse_date(root[i]["Date"].asString());

                auto it_age = std::find(age_names.begin(), age_names.end() - 1, root[i]["Age_RKI"].asString());
                if (it_age != age_names.end() - 1) {
                    auto age = size_t(it_age - age_names.begin());

                    bool read_icu = false; //params.populations.get({age, SecirCompartments::U}) == 0;

                    if (date_df == offset_date_by_days(date, 0)) {
                        num_inf[age] += (1 - mu_C_R[age]) * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                        num_rec[age] += root[i]["Confirmed"].asDouble();
                    }
                    if (date_df == offset_date_by_days(date, days_surplus)) {
                        num_car[age] +=
                            (2 * mu_C_R[age] - 1) * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // -R9
                    if (date_df == offset_date_by_days(date, -t_car_to_rec[age] + days_surplus)) {
                        num_car[age] -= mu_C_R[age] * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // +R2
                    if (date_df == offset_date_by_days(date, +t_exp_to_car[age] + days_surplus)) {
                        num_exp[age] += mu_C_R[age] * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // +R3
                    if (date_df == offset_date_by_days(date, +t_car_to_inf[age] + days_surplus)) {
                        num_car[age] += (1 - mu_C_R[age]) * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                        num_exp[age] -= (1 - mu_C_R[age]) * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // R2 - R9
                    if (date_df == offset_date_by_days(date, t_exp_to_car[age] - t_car_to_rec[age] + days_surplus)) {
                        num_exp[age] -= mu_C_R[age] * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // R2 + R3
                    if (date_df == offset_date_by_days(date, t_exp_to_car[age] + t_car_to_inf[age] + days_surplus)) {
                        num_exp[age] += (1 - mu_C_R[age]) * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // -R4
                    if (date_df == offset_date_by_days(date, -t_inf_to_rec[age])) {
                        num_inf[age] -= (1 - mu_C_R[age]) * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // -R6
                    if (date_df == offset_date_by_days(date, -t_inf_to_hosp[age])) {
                        num_inf[age] -= mu_I_H[age] * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                        num_hosp[age] += mu_I_H[age] * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // -R6 - R7
                    if (date_df == offset_date_by_days(date, -t_inf_to_hosp[age] - t_hosp_to_icu[age])) {
                        num_inf[age] +=
                            mu_I_H[age] * mu_H_U[age] * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                        num_hosp[age] -=
                            mu_I_H[age] * mu_H_U[age] * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                        if (read_icu) {
                            num_icu[age] +=
                                mu_H_U[age] * mu_I_H[age] * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                        }
                    }
                    // -R6 - R5
                    if (date_df == offset_date_by_days(date, -t_inf_to_hosp[age] - t_hosp_to_rec[age])) {
                        num_inf[age] +=
                            mu_I_H[age] * (1 - mu_H_U[age]) * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                        num_hosp[age] -=
                            mu_I_H[age] * (1 - mu_H_U[age]) * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // -R10 - R6 - R7
                    if (date_df ==
                        offset_date_by_days(date, -t_icu_to_dead[age] - t_inf_to_hosp[age] - t_hosp_to_icu[age])) {
                        num_death[age] += root[i]["Deaths"].asDouble();
                    }
                    if (read_icu) {
                        // -R6 - R7 - R7
                        if (date_df == offset_date_by_days(date, -t_inf_to_hosp[age] - 2 * t_hosp_to_icu[age])) {
                            num_icu[age] -= mu_I_H[age] * mu_H_U[age] * mu_H_U[age] * scaling_factor_inf[age] *
                                            root[i]["Confirmed"].asDouble();
                        }
                        // -R6 - R5 - R7
                        if (date_df == offset_date_by_days(date, -t_inf_to_hosp[age] - t_hosp_to_icu[age])) {
                            num_icu[age] -= mu_I_H[age] * mu_H_U[age] * (1 - mu_H_U[age]) * scaling_factor_inf[age] *
                                            root[i]["Confirmed"].asDouble();
                        }
                    }
                }
            }
        }

        for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
            auto region = vregion[region_idx];

            auto& num_car   = vnum_car[region_idx];
            auto& num_inf   = vnum_inf[region_idx];
            auto& num_rec   = vnum_rec[region_idx];
            auto& num_exp   = vnum_exp[region_idx];
            auto& num_hosp  = vnum_hosp[region_idx];
            auto& num_death = vnum_death[region_idx];
            auto& num_icu   = vnum_icu[region_idx];

            size_t num_groups = age_ranges.size();
            for (size_t i = 0; i < num_groups; i++) {
                num_rec[i] -=
                    (num_inf[i] + num_car[i] + num_exp[i] + num_hosp[i]) / scaling_factor_inf[i] + num_death[i];
                auto try_fix_constraints = [region, &age_names, i](double& value, double error, auto str) {
                    if (value < error) {
                        //this should probably return a failure
                        //but the algorithm is not robust enough to avoid large negative values and there are tests that rely on it
                        log_error("{:s} for age group {:s} is {:.4f} for region {:d}, exceeds expected negative value.",
                                  str, age_names[i], value, region);
                        value = 0.0;
                    }
                    else if (value < 0) {
                        log_info("{:s} for age group {:s} is {:.4f} for region {:d}, automatically corrected", str,
                                 age_names[i], value, region);
                        value = 0.0;
                    }
                };

                try_fix_constraints(num_inf[i], -5, "Infected");
                try_fix_constraints(num_car[i], -5, "Carrier");
                try_fix_constraints(num_exp[i], -5, "Exposed");
                try_fix_constraints(num_hosp[i], -5, "Hospitalized");
                try_fix_constraints(num_death[i], -5, "Dead");
                try_fix_constraints(num_icu[i], -5, "ICU");
                try_fix_constraints(num_rec[i], -20, "Recovered");
            }
        }

        return success();
    }

    IOResult<void> read_divi_data(const std::string& path, const std::string& id_name, const std::vector<int>& vregion,
                                  Date date, std::vector<double>& vnum_icu)
    {
        if (!boost::filesystem::exists(path)) {
            log_error("DIVI data file not found: {}.", path);
            return failure(StatusCode::FileNotFound, path);
        }

        Json::Reader reader;
        Json::Value root;

        std::ifstream divi(path);
        if (!reader.parse(divi, root)) {
            log_error(reader.getFormattedErrorMessages());
            return failure(StatusCode::UnknownError, path + ", " + reader.getFormattedErrorMessages());
        }

        BOOST_OUTCOME_TRY(max_date, get_max_date(root));
        if (max_date == Date(0, 1, 1)) {
            log_error("DIVI data file is empty.");
            return failure(StatusCode::InvalidFileFormat, path + ", file is empty.");
        }
        if (max_date < date) {
            log_error("Specified date does not exist in DIVI data.");
            return failure(StatusCode::OutOfRange, path + ", specified date does not exist in DIVI data.");
        }

        for (unsigned int i = 0; i < root.size(); i++) {
            auto it      = std::find_if(vregion.begin(), vregion.end(), [&root, i, &id_name](auto r) {
                return r == 0 || r == root[i][id_name].asInt();
            });
            auto date_df = parse_date(root[i]["Date"].asString());
            if (it != vregion.end() && date_df == date) {
                auto region_idx = size_t(it - vregion.begin());
                auto& num_icu   = vnum_icu[region_idx];
                num_icu         = root[i]["ICU"].asDouble();
            }
        }

        return success();
    }

    IOResult<std::vector<std::vector<double>>> read_population_data(const std::string& path, const std::string& id_name,
                                                                    const std::vector<int>& vregion)
    {
        if (!boost::filesystem::exists(path)) {
            log_error("Population data file not found: {}.", path);
            return failure(StatusCode::FileNotFound, path);
        }

        Json::Reader reader;
        Json::Value root;

        std::ifstream census(path);
        if (!reader.parse(census, root)) {
            log_error(reader.getFormattedErrorMessages());
            return failure(StatusCode::UnknownError, path + ", " + reader.getFormattedErrorMessages());
        }

        std::vector<std::string> age_names = {"<3 years",    "3-5 years",   "6-14 years",  "15-17 years",
                                              "18-24 years", "25-29 years", "30-39 years", "40-49 years",
                                              "50-64 years", "65-74 years", ">74 years"};
        std::vector<double> age_ranges     = {3., 3., 9., 3., 7., 5., 10., 10., 15., 10., 25.};

        std::vector<std::vector<double>> interpolation(age_names.size());
        std::vector<bool> carry_over;

        interpolate_ages(age_ranges, interpolation, carry_over);

        std::vector<std::vector<double>> vnum_population(vregion.size(), std::vector<double>(age_names.size(), 0.0));

        for (unsigned int i = 0; i < root.size(); i++) {
            auto it = std::find_if(vregion.begin(), vregion.end(), [&root, i, &id_name](auto r) {
                return r == 0 || (int)root[i][id_name].asDouble() / 1000 == r || root[i][id_name] == r;
            });
            if (it != vregion.end()) {
                auto region_idx      = size_t(it - vregion.begin());
                auto& num_population = vnum_population[region_idx];
                for (size_t age = 0; age < age_names.size(); age++) {
                    num_population[age] += root[i][age_names[age]].asDouble();
                }
            }
        }

        std::vector<std::vector<double>> interpol_population(vregion.size(),
                                                             std::vector<double>(age_names.size(), 0.0));
        for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
            auto& num_population = vnum_population[region_idx];

            int counter = 0;
            for (size_t i = 0; i < interpolation.size() - 1; i++) {
                for (size_t j = 0; j < interpolation[i].size(); j++) {
                    interpol_population[region_idx][counter] += interpolation[i][j] * num_population[i];
                    if (j < interpolation[i].size() - 1 || !carry_over[i]) {
                        counter++;
                    }
                }
            }
        }

        return success(interpol_population);
    }

    void set_vaccine_data(std::vector<SecirModelV>& model, const std::string& path, const std::string& id_name,
                          const std::vector<int>& vregion)
    {
        if (!boost::filesystem::exists(path)) {
            log_error("Vaccine data file missing: {}.", path);
            return;
        }
        Json::Reader reader;
        Json::Value root;

        std::ifstream vaccine(path);
        reader.parse(vaccine, root);

        std::vector<std::vector<double>> vnum_vaccine(model.size(), std::vector<double>(6, 0.0));
        for (unsigned int i = 0; i < root.size(); i++) {
            auto it = std::find_if(vregion.begin(), vregion.end(), [&root, i, &id_name](auto r) {
                return r == 0 || r == root[i][id_name].asInt();
            });
            if (it != vregion.end()) {
                auto region_idx                 = size_t(it - vregion.begin());
                auto& num_first_vac             = vnum_vaccine[region_idx][0];
                auto& num_full_vac              = vnum_vaccine[region_idx][1];
                auto& num_vac_ratio_young_first = vnum_vaccine[region_idx][2];
                auto& num_vac_ratio_old_first   = vnum_vaccine[region_idx][3];
                auto& num_vac_ratio_young_full  = vnum_vaccine[region_idx][4];
                auto& num_vac_ratio_old_full    = vnum_vaccine[region_idx][5];
                num_first_vac                   = root[i]["First_Shot"].asDouble();
                num_full_vac                    = root[i]["Full_Vaccination"].asDouble();
                num_vac_ratio_young_first       = root[i]["Ratio_Young_First"].asDouble();
                num_vac_ratio_old_first         = root[i]["Ratio_Old_First"].asDouble();
                num_vac_ratio_young_full        = root[i]["Ratio_Young_Full"].asDouble();
                num_vac_ratio_old_full          = root[i]["Ratio_Old_Full"].asDouble();
            }
        }
        for (size_t region_idx = 0; region_idx < model.size(); ++region_idx) {
            auto& params      = model[region_idx];
            auto& num_vaccine = vnum_vaccine[region_idx];

            auto num_groups = params.parameters.get_num_groups();

            double r_young_first = num_vaccine[2] / 100;
            double r_old_first   = num_vaccine[3] / 100;
            double r_young_full  = num_vaccine[4] / 100;
            double r_old_full    = num_vaccine[5] / 100;
            std::vector<double> num_vaccine_first((size_t)num_groups, 0);
            std::vector<double> num_vaccine_full((size_t)num_groups, 0);
            double sum_population = 0;
            for (auto i = AgeGroup(0); i < num_groups; i++) {
                sum_population += params.populations.get_group_total(i);
            }
            for (auto i = AgeGroup(0); i < num_groups; i++) {
                if ((size_t)i < 4 && (size_t)i > 1) {
                    num_vaccine_first[(size_t)i] = r_young_first * params.populations.get_group_total(i);
                    num_vaccine_full[(size_t)i]  = r_young_full * params.populations.get_group_total(i);
                }
                else if ((size_t)i > 3) {
                    num_vaccine_first[(size_t)i] = r_old_first * params.populations.get_group_total(i);
                    num_vaccine_full[(size_t)i]  = r_old_full * params.populations.get_group_total(i);
                }
                if (num_vaccine_first[(size_t)i] < 0) {
                    num_vaccine_first[(size_t)i] = 0;
                }
                if (num_vaccine_full[(size_t)i] < 0) {
                    num_vaccine_full[(size_t)i] = 0;
                }
            }

            for (auto i = AgeGroup(0); i < num_groups; i++) {
                if (params.populations[{i, InfectionStateV::Susceptible}] < num_vaccine_first[(size_t)i]) {
                    num_vaccine_first[(size_t)i] = 0.9 * params.populations[{i, InfectionStateV::Susceptible}];
                }
                params.populations[{i, InfectionStateV::Susceptible}] =
                    params.populations[{i, InfectionStateV::Susceptible}] - num_vaccine_first[(size_t)i];
                params.populations[{i, InfectionStateV::SusceptibleV1}] =
                    params.populations[{i, InfectionStateV::SusceptibleV1}] + num_vaccine_first[(size_t)i] -
                    num_vaccine_full[(size_t)i];
                params.populations[{i, InfectionStateV::Recovered}] =
                    params.populations[{i, epi::InfectionStateV::Recovered}] + num_vaccine_full[(size_t)i];
            }
        }
    }

    void get_new_vaccine_growth(std::vector<SecirModelV>& model, const std::string& path, Date date,
                                const std::string& id_name, const std::vector<int>& vregion, int num_days)
    {
        if (!boost::filesystem::exists(path)) {
            log_error("Vaccine data file missing: {}.", path);
            return;
        }

        auto num_groups = model[0].parameters.get_num_groups();

        Json::Reader reader;
        Json::Value root;

        auto days_until_effective1 = (int)(double)model[0].parameters.get<DaysUntilEffective>()[AgeGroup(0)];
        auto days_until_effective2 =
            (int)(double)model[0].parameters.get<DaysUntilEffectiveFull>()[AgeGroup(0)];
        auto vaccination_gap =  
            (int)(double)model[0].parameters.get<VaccinationGap>()[AgeGroup(0)];

        std::vector<std::string> age_names = {"0-4", "5-14", "15-34", "35-59", "60-79", "80-99"};

        std::ifstream vaccine(path);
        reader.parse(vaccine, root);

        for (size_t i = 0; i < model.size(); ++i) {
            for (auto g = AgeGroup(0); g < num_groups; ++g) {

                model[i].parameters.template get<DailyFirstVaccination>()[g] = {};
                model[i].parameters.template get<DailyFullVaccination>()[g]  = {};
                for (int d = 0; d < num_days + 1; ++d) {
                    model[i].parameters.template get<DailyFirstVaccination>()[g].push_back(0.0);
                    model[i].parameters.template get<DailyFullVaccination>()[g].push_back(0.0);
                }
            }
        }

        Date max_date = Date(2020, 1, 1);
        for (unsigned int i = 0; i < root.size(); i++) {
            auto date_df = parse_date(root[i]["Date"].asString());
            if (date_df > max_date){
                max_date = date_df;
            }
        }
        auto max_full_date = offset_date_by_days(max_date,  days_until_effective1 - days_until_effective2 - vaccination_gap);

        for (unsigned int i = 0; i < root.size(); i++) {
            auto it = std::find_if(vregion.begin(), vregion.end(), [&root, i, &id_name](auto r) {
                return r == 0 || root[i][id_name].asInt() == r;
            });
            auto date_df = parse_date(root[i]["Date"].asString());
            if (it != vregion.end()) {
                auto region_idx = size_t(it - vregion.begin());


                auto it_age = std::find(age_names.begin(), age_names.end() - 1, root[i]["Age_RKI"].asString());
                if (it_age != age_names.end() - 1) {
                    auto age = size_t(it_age - age_names.begin());

                    for (int d = 0; d < num_days + 1; ++d) {

                        auto offset_first_date = offset_date_by_days(date, d - days_until_effective1 + vaccination_gap);
                        if (date_df == offset_first_date && max_date >= offset_first_date) {
                            model[region_idx].parameters.template get<DailyFirstVaccination>()[AgeGroup(age)][d] =
                                root[i]["Vacc_completed"].asDouble();
                        }
                        else if(max_date < offset_first_date && date_df == offset_date_by_days(max_date, -1)){
                            model[region_idx].parameters.template get<DailyFirstVaccination>()[AgeGroup(age)][d] -=
                                root[i]["Vacc_completed"].asDouble();
                        }
                        else if(max_date < offset_first_date && date_df == max_date){
                            model[region_idx].parameters.template get<DailyFirstVaccination>()[AgeGroup(age)][d] +=
                                2*root[i]["Vacc_completed"].asDouble();
                        }
                        
                        auto offset_full_date = offset_date_by_days(date, d - days_until_effective2);
                        if (date_df == offset_full_date && max_full_date >= offset_full_date) {
                            model[region_idx].parameters.template get<DailyFullVaccination>()[AgeGroup(age)][d] =
                                root[i]["Vacc_completed"].asDouble();
                        }
                        else if (max_full_date < offset_full_date && date_df == offset_date_by_days(max_full_date, -1)) {
                            model[region_idx].parameters.template get<DailyFullVaccination>()[AgeGroup(age)][d] -=
                                root[i]["Vacc_completed"].asDouble();
                        }
                        else if (max_full_date < offset_full_date && date_df == max_full_date) {
                            model[region_idx].parameters.template get<DailyFullVaccination>()[AgeGroup(age)][d] +=
                                2*root[i]["Vacc_completed"].asDouble();
                        }
                    }
                }
            }
        }
    }

    void get_vaccine_growth(std::vector<SecirModelV>& model, const std::string& path, int window)
    {

        if (!boost::filesystem::exists(path)) {
            log_error("Vaccine data file missing: {}.", path);
            return;
        }

        auto num_groups = model[0].parameters.get_num_groups();

        Json::Reader reader;
        Json::Value root;

        size_t vaccination_gap      = (size_t)(int)(double)model[0].parameters.get<VaccinationGap>()[AgeGroup(0)];
        size_t days_until_effective = (size_t)(int)(double)model[0].parameters.get<DaysUntilEffective>()[AgeGroup(0)];

        std::ifstream vaccine(path);
        reader.parse(vaccine, root);
        std::vector<double> vaccine_first(root.size(), 0);
        std::vector<double> vaccine_full(root.size(), 0);
        std::vector<double> daily_first_vaccination;
        std::vector<double> daily_full_vaccination;
        for (unsigned int i = 0; i < root.size(); i++) {
            vaccine_first[i] = root[i]["Erstimpfung"].asDouble();
            vaccine_full[i]  = root[i]["Zweitimpfung"].asDouble();
        }

        double vaccine_first_mean;
        double vaccine_full_mean;
        for (size_t j = 0; j < vaccination_gap + days_until_effective; ++j) {
            vaccine_first_mean = 0;
            vaccine_full_mean  = 0;
            for (size_t i = root.size() - window + j - vaccination_gap - days_until_effective;
                 i < root.size() + j - vaccination_gap - days_until_effective; ++i) {
                vaccine_first_mean += vaccine_first[i] / window;
                vaccine_full_mean += vaccine_full[i] / window;
            }
            daily_first_vaccination.push_back(vaccine_first_mean);
            daily_full_vaccination.push_back(vaccine_full_mean);
        }
        double sum_population = 0;
        for (auto& param : model) {
            for (auto i = AgeGroup(2); i < num_groups; ++i) {
                sum_population += param.populations.get_group_total(i);
            }
        }

        for (auto& param : model) {

            for (auto i = AgeGroup(2); i < num_groups; ++i) {
                param.parameters.get<DailyFirstVaccination>()[i] = {};
                param.parameters.get<DailyFullVaccination>()[i]  = {};
                for (size_t j = 0; j < daily_first_vaccination.size(); ++j) {
                    param.parameters.get<DailyFirstVaccination>()[i].push_back(
                        daily_first_vaccination[j] * (param.populations.get_group_total(i) / sum_population));
                    param.parameters.get<DailyFullVaccination>()[i].push_back(
                        daily_full_vaccination[j] * (param.populations.get_group_total(i) / sum_population));
                }

                double vaccine_growth_first =
                    daily_first_vaccination.back() * (param.populations.get_group_total(i) / sum_population);
                double vaccine_growth_full =
                    daily_full_vaccination.back() * (param.populations.get_group_total(i) / sum_population);
                param.parameters.get<VaccineGrowthFirst>()[i] = vaccine_growth_first;
                param.parameters.get<VaccineGrowthFull>()[i]  = vaccine_growth_full;
            }
        }
    }

} // namespace details

IOResult<std::vector<int>> get_county_ids(const std::string& path)
{
    Json::Reader reader;
    Json::Value root;

    std::vector<int> id;

    auto filename = path_join(path, "county_current_population.json");
    std::ifstream census(filename);
    if (!reader.parse(census, root)) {
        log_error(reader.getFormattedErrorMessages());
        return failure(StatusCode::UnknownError, filename + ", " + reader.getFormattedErrorMessages());
    }

    for (unsigned int i = 0; i < root.size(); i++) {
        auto val = root[i]["ID_County"];
        if (!val.isInt()) {
            log_error("ID_County field must be an integer.");
            return failure(StatusCode::InvalidFileFormat, filename + ", ID_County field must be an integer.");
        }
        id.push_back(val.asInt());
    }

    return success(id);
}

} // namespace epi
