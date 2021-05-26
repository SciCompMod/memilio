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

#include <tixi.h>

#include <json/json.h>
#include <json/value.h>

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
    Date get_max_date(Json::Value root)
    {
        Date max_date = Date(0, 1, 1);
        for (unsigned int i = 0; i < root.size(); i++) {
            auto date_temp = parse_date(root[i]["Date"].asString());
            if (date_temp > max_date) {
                max_date = date_temp;
            }
        }

        return max_date;
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

    void read_rki_data(std::string const& path, const std::string& id_name, std::vector<int> const& vregion, Date date,
                       std::vector<std::vector<double>>& vnum_exp, std::vector<std::vector<double>>& vnum_car,
                       std::vector<std::vector<double>>& vnum_inf, std::vector<std::vector<double>>& vnum_hosp,
                       std::vector<std::vector<double>>& vnum_icu, std::vector<std::vector<double>>& vnum_death,
                       std::vector<std::vector<double>>& vnum_rec, const std::vector<std::vector<int>>& vt_car_to_rec,
                       const std::vector<std::vector<int>>& vt_car_to_inf,
                       const std::vector<std::vector<int>>& vt_exp_to_car,
                       const std::vector<std::vector<int>>& vt_inf_to_rec,
                       const std::vector<std::vector<int>>& vt_inf_to_hosp,
                       const std::vector<std::vector<int>>& vt_hosp_to_rec,
                       const std::vector<std::vector<int>>& vt_hosp_to_icu,
                       const std::vector<std::vector<int>>& vt_icu_to_dead,
                       const std::vector<std::vector<double>>& vmu_C_R, const std::vector<std::vector<double>>& vmu_I_H,
                       const std::vector<std::vector<double>>& vmu_H_U, const std::vector<double>& scaling_factor_inf)
    {
        if (!boost::filesystem::exists(path)) {
            log_error("RKI data file missing: {}.", path);
            return;
        }

        Json::Reader reader;
        Json::Value root;

        std::ifstream rki(path);
        if (!reader.parse(rki, root)) {
            log_error(reader.getFormattedErrorMessages());
            return;
        }

        Date max_date = get_max_date(root);
        assert(max_date != Date(0, 1, 1) && "File is epmty");
        if (max_date == Date(0, 1, 1)) {
            log_error("File is empty");
            return;
        }
        assert(max_date >= date && "Specified date does not exist in divi data");
        if (max_date < date) {
            log_error("Specified date does not exist in divi data");
            return;
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
                auto try_fix_constraints = [region, &age_names, i](double& value, double error, auto str) {
                    if (value < error) {
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
    }


    void set_rki_data(std::vector<SecirModel>& model, const std::string& path, const std::string& id_name,
                      std::vector<int> const& region, Date date, const std::vector<double>& scaling_factor_inf)
    {

        std::vector<double> age_ranges = {5., 10., 20., 25., 20., 20.};
        assert(scaling_factor_inf.size() == age_ranges.size());

        std::vector<std::vector<int>> t_car_to_rec{model.size()}; // R9
        std::vector<std::vector<int>> t_car_to_inf{model.size()}; // R3
        std::vector<std::vector<int>> t_exp_to_car{model.size()}; // R2
        std::vector<std::vector<int>> t_inf_to_rec{model.size()}; // R4
        std::vector<std::vector<int>> t_inf_to_hosp{model.size()}; // R6
        std::vector<std::vector<int>> t_hosp_to_rec{model.size()}; // R5
        std::vector<std::vector<int>> t_hosp_to_icu{model.size()}; // R7
        std::vector<std::vector<int>> t_icu_to_dead{model.size()}; // R10

        std::vector<std::vector<double>> mu_C_R{model.size()};
        std::vector<std::vector<double>> mu_I_H{model.size()};
        std::vector<std::vector<double>> mu_H_U{model.size()};

        for (size_t county = 0; county < model.size(); county++) {
            for (size_t group = 0; group < age_ranges.size(); group++) {

                t_car_to_inf[county].push_back(
                    static_cast<int>(2 * (model[county].parameters.get<epi::IncubationTime>()[(epi::AgeGroup)group] -
                                          model[county].parameters.get<epi::SerialInterval>()[(epi::AgeGroup)group])));
                t_car_to_rec[county].push_back(static_cast<int>(
                    t_car_to_inf[county][group] + 0.5 * model[county].parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)group]));
                t_exp_to_car[county].push_back(
                    static_cast<int>(2 * model[county].parameters.get<epi::SerialInterval>()[(epi::AgeGroup)group] -
                                     model[county].parameters.get<epi::IncubationTime>()[(epi::AgeGroup)group]));
                t_inf_to_rec[county].push_back(
                    static_cast<int>(model[county].parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)group]));
                t_inf_to_hosp[county].push_back(
                    static_cast<int>(model[county].parameters.get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)group]));
                t_hosp_to_rec[county].push_back(
                    static_cast<int>(model[county].parameters.get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)group]));
                t_hosp_to_icu[county].push_back(
                    static_cast<int>(model[county].parameters.get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)group]));
                t_icu_to_dead[county].push_back(
                    static_cast<int>(model[county].parameters.get<epi::ICUToDeathTime>()[(epi::AgeGroup)group]));

                mu_C_R[county].push_back(model[county].parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)group]);
                mu_I_H[county].push_back(
                    model[county].parameters.get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)group]);
                mu_H_U[county].push_back(model[county].parameters.get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)group]);
            }
        }
        std::vector<std::vector<double>> num_inf(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_death(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_rec(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_exp(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_car(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_hosp(model.size(), std::vector<double>(age_ranges.size(), 0.0));
        std::vector<std::vector<double>> num_icu(model.size(), std::vector<double>(age_ranges.size(), 0.0));

        read_rki_data(path, id_name, region, date, num_exp, num_car, num_inf, num_hosp, num_icu, num_death, num_rec,
                      t_car_to_rec, t_car_to_inf, t_exp_to_car, t_inf_to_rec, t_inf_to_hosp, t_hosp_to_rec,
                      t_hosp_to_icu, t_icu_to_dead, mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf);

        for (size_t county = 0; county < model.size(); county++) {
            if (std::accumulate(num_inf[county].begin(), num_inf[county].end(), 0.0) > 0) {
                size_t num_groups = (size_t)model[county].parameters.get_num_groups();
                for (size_t i = 0; i < num_groups; i++) {
                    model[county].populations[{AgeGroup(i), InfectionState::Exposed}] =
                        num_exp[county][i];
                    model[county].populations[{AgeGroup(i), InfectionState::Carrier}] =
                        num_car[county][i];
                    model[county].populations[{AgeGroup(i), InfectionState::Infected}] =
                        num_inf[county][i];
                    model[county].populations[{AgeGroup(i), InfectionState::Hospitalized}] =
                        num_hosp[county][i];
                    model[county].populations[{AgeGroup(i), InfectionState::Dead}] =
                        num_death[county][i];
                    model[county].populations[{AgeGroup(i), InfectionState::Recovered}] =
                        num_rec[county][i];
                }
            }
            else {
                log_warning("No infections reported on date " + std::to_string(date.year) + "-" +
                            std::to_string(date.month) + "-" + std::to_string(date.day) + " for region " +
                            std::to_string(region[county]) + ". Population data has not been set.");
            }
        }
    }

    void read_divi_data(const std::string& path, const std::string& id_name, const std::vector<int>& vregion, Date date,
                        std::vector<double>& vnum_icu)
    {
        if (!boost::filesystem::exists(path)) {
            log_error("DIVI data file missing: {}.", path);
            return;
        }

        Json::Reader reader;
        Json::Value root;

        std::ifstream divi(path);
        reader.parse(divi, root);

        Date max_date = get_max_date(root);
        assert(max_date != Date(0, 1, 1) && "File is epmty");
        if (max_date == Date(0, 1, 1)) {
            log_error("File is empty");
            return;
        }
        assert(max_date >= date && "Specified date does not exist in divi data");
        if (max_date < date) {
            log_error("Specified date does not exist in divi data");
            return;
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
    }

    std::vector<std::vector<double>> read_population_data(const std::string& path, const std::string& id_name,
                                                          const std::vector<int>& vregion)
    {
        if (!boost::filesystem::exists(path)) {
            log_error("Population data file missing: {}.", path);
            std::vector<std::vector<double>> num_population;
            return num_population;
        }

        Json::Reader reader;
        Json::Value root;

        std::ifstream census(path);
        reader.parse(census, root);

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

        return interpol_population;
    }

    void set_population_data(std::vector<SecirModel>& model, const std::string& path, const std::string& id_name,
                             const std::vector<int> vregion)
    {
        std::vector<std::vector<double>> num_population = read_population_data(path, id_name, vregion);

        for (size_t region = 0; region < vregion.size(); region++) {
            if (std::accumulate(num_population[region].begin(), num_population[region].end(), 0.0) > 0) {
                auto num_groups = model[region].parameters.get_num_groups();
                for (auto i = AgeGroup(0); i < num_groups; i++) {
                    model[region].populations.set_difference_from_group_total<epi::AgeGroup>(
                    {i, InfectionState::Susceptible}, num_population[region][size_t(i)]);
                }
            }
            else {
                log_warning("No population data available for region " + std::to_string(region) +
                            ". Population data has not been set.");
            }
        }
    }


    void set_divi_data(std::vector<SecirModel>& model, const std::string& path, const std::string& id_name,
                       const std::vector<int> vregion, Date date, double scaling_factor_icu)
    {
        std::vector<double> sum_mu_I_U(vregion.size(), 0);
        std::vector<std::vector<double>> mu_I_U{model.size()};
        for (size_t region = 0; region < vregion.size(); region++) {
            auto num_groups = model[region].parameters.get_num_groups();
            for (auto i = epi::AgeGroup(0); i < num_groups; i++) {
                sum_mu_I_U[region] += model[region].parameters.get<epi::ICUCasesPerHospitalized>()[i] *
                                      model[region].parameters.get<epi::HospitalizedCasesPerInfectious>()[i];
                mu_I_U[region].push_back(model[region].parameters.get<epi::ICUCasesPerHospitalized>()[i] *
                                         model[region].parameters.get<epi::HospitalizedCasesPerInfectious>()[i]);
            }
        }
        std::vector<double> num_icu(model.size(), 0.0);
        read_divi_data(path, id_name, vregion, date, num_icu);

        for (size_t region = 0; region < vregion.size(); region++) {
            auto num_groups = model[region].parameters.get_num_groups();
            for (auto i = epi::AgeGroup(0); i < num_groups; i++) {
                model[region].populations[{i, epi::InfectionState::ICU}] =
                    scaling_factor_icu * num_icu[region] * mu_I_U[region][(size_t)i] / sum_mu_I_U[region];
            }
        }

    }

} // namespace details

void write_element(const TixiDocumentHandle& handle, const std::string& path, const std::string& element_name,
                   const UncertainValue& element, int io_mode, int num_runs)
{
    auto element_path = path_join(path, element_name);

    if (io_mode == 0) {
        tixiAddDoubleElement(handle, path.c_str(), element_name.c_str(), (double)element, "%g");
    }
    else if (io_mode == 1 || io_mode == 2 || io_mode == 3) {
        assert(element.get_distribution().get() && ("No Distribution detected for " + element_name +
                                                    ". Either define a distribution or choose a different io_mode.")
                                                       .c_str());
        auto distribution = element.get_distribution().get();
        write_distribution(handle, path, element_name, *distribution);
        if (io_mode == 2) {
            tixiAddDoubleElement(handle, element_path.c_str(), "Value", (double)element, "%g");
        }
    }
    else {
        // TODO error handling
        epi::log_error("Wrong input for io_mode.");
    }

    if (io_mode == 3) {
        std::vector<double> predef_sample;
        for (int run = 0; run < num_runs; run++) {
            predef_sample.push_back((double)element);
        }
        write_predef_sample(handle, element_path, predef_sample);
    }
}

void write_distribution(const TixiDocumentHandle& handle, const std::string& path, const std::string& element_name,
                        const ParameterDistribution& distribution)
{

    struct WriteDistVisitor : public ConstParameterDistributionVisitor {
        WriteDistVisitor(const std::string& xml_path, TixiDocumentHandle tixi_handle)
            : handle(tixi_handle)
            , element_path(xml_path)
        {
        }

        void visit(const ParameterDistributionNormal& normal_distribution) override
        {
            tixiAddTextElement(handle, element_path.c_str(), "Distribution", "Normal");
            tixiAddDoubleElement(handle, element_path.c_str(), "Mean", normal_distribution.get_mean(), "%g");
            tixiAddDoubleElement(handle, element_path.c_str(), "Deviation", normal_distribution.get_standard_dev(),
                                 "%g");
            tixiAddDoubleElement(handle, element_path.c_str(), "Min", normal_distribution.get_lower_bound(), "%g");
            tixiAddDoubleElement(handle, element_path.c_str(), "Max", normal_distribution.get_upper_bound(), "%g");
        }

        void visit(const ParameterDistributionUniform& uniform_distribution) override
        {
            tixiAddTextElement(handle, element_path.c_str(), "Distribution", "Uniform");
            tixiAddDoubleElement(handle, element_path.c_str(), "Min", uniform_distribution.get_lower_bound(), "%g");
            tixiAddDoubleElement(handle, element_path.c_str(), "Max", uniform_distribution.get_upper_bound(), "%g");
        }

        TixiDocumentHandle handle;
        std::string element_path;
    };

    tixiCreateElement(handle, path.c_str(), element_name.c_str());
    auto element_path = path_join(path, element_name);

    WriteDistVisitor visitor(element_path, handle);
    distribution.accept(visitor);

    tixiAddFloatVector(handle, element_path.c_str(), "PredefinedSamples", distribution.get_predefined_samples().data(),
                       static_cast<int>(distribution.get_predefined_samples().size()), "%g");
}

std::unique_ptr<UncertainValue> read_element(TixiDocumentHandle handle, const std::string& path, int io_mode)
{
    std::unique_ptr<UncertainValue> value;
    ReturnCode status;
    unused(status);

    if (io_mode == 0) {
        double read_buffer;
        status = tixiGetDoubleElement(handle, path.c_str(), &read_buffer);
        assert(status == SUCCESS && ("Failed to read value at " + path).c_str());
        value = std::make_unique<UncertainValue>(read_buffer);
    }
    else if (io_mode == 1 || io_mode == 2 || io_mode == 3) {
        std::unique_ptr<ParameterDistribution> distribution = read_distribution(handle, path);

        if (io_mode == 2) {
            double read_buffer;
            status = tixiGetDoubleElement(handle, path_join(path, "Value").c_str(), &read_buffer);
            assert(status == SUCCESS && ("Failed to read value at " + path).c_str());
            value = std::make_unique<UncertainValue>(read_buffer);
        }
        value->set_distribution(*distribution.get());
    }
    else {
        // TODO error handling
        epi::log_error("Wrong input for io_mode.");
    }
    return value;
}

std::unique_ptr<ParameterDistribution> read_distribution(TixiDocumentHandle handle, const std::string& path)
{
    ReturnCode status;
    unused(status);
    std::unique_ptr<ParameterDistribution> distribution;

    char* distri_str;
    status = tixiGetTextElement(handle, path_join(path, "Distribution").c_str(), &distri_str);
    assert(status == SUCCESS && ("Failed to read distribution type at " + path).c_str());
    if (strcmp("Normal", distri_str) == 0) {
        double mean;
        double dev;
        double min;
        double max;
        status = tixiGetDoubleElement(handle, path_join(path, "Mean").c_str(), &mean);
        assert(status == SUCCESS && ("Failed to read mean at " + path).c_str());

        status = tixiGetDoubleElement(handle, path_join(path, "Deviation").c_str(), &dev);
        assert(status == SUCCESS && ("Failed to read deviation at " + path).c_str());

        status = tixiGetDoubleElement(handle, path_join(path, "Min").c_str(), &min);
        assert(status == SUCCESS && ("Failed to read min value at " + path).c_str());

        status = tixiGetDoubleElement(handle, path_join(path, "Max").c_str(), &max);
        assert(status == SUCCESS && ("Failed to read max value at " + path).c_str());

        distribution = std::make_unique<ParameterDistributionNormal>(min, max, mean, dev);
    }
    else if (strcmp("Uniform", distri_str) == 0) {
        double min;
        double max;
        status = tixiGetDoubleElement(handle, path_join(path, "Min").c_str(), &min);
        assert(status == SUCCESS && ("Failed to read min value at " + path).c_str());

        status = tixiGetDoubleElement(handle, path_join(path, "Max").c_str(), &max);
        assert(status == SUCCESS && ("Failed to read max value at " + path).c_str());

        distribution = std::make_unique<ParameterDistributionUniform>(min, max);
    }
    else {
        // TODO: true error handling
        epi::log_error("Unknown distribution.");
        assert(false && "Unknown distribution.");
    }

    auto predef_path = path_join(path, "PredefinedSamples");
    int n_predef;
    tixiGetVectorSize(handle, predef_path.c_str(), &n_predef);

    double* predef = nullptr;
    tixiGetFloatVector(handle, predef_path.c_str(), &predef, n_predef);

    for (int i = 0; i < n_predef; i++) {
        distribution->add_predefined_sample(predef[i]);
    }

    return distribution;
}

void write_predef_sample(TixiDocumentHandle handle, const std::string& path, const std::vector<double>& samples)
{
    tixiRemoveElement(handle, path_join(path, "PredefinedSamples").c_str());
    tixiAddFloatVector(handle, path.c_str(), "PredefinedSamples", samples.data(), static_cast<int>(samples.size()),
                       "%g");
}

SecirModel read_parameter_space(TixiDocumentHandle handle, const std::string& path, int io_mode)
{
    ReturnCode status;
    unused(status);

    int num_groups;
    status = tixiGetIntegerElement(handle, path_join(path, "NumberOfGroups").c_str(), &num_groups);
    assert(status == SUCCESS && ("Failed to read num_groups at " + path).c_str());

    status = tixiGetIntegerElement(handle, path_join(path, "NumberOfGroups").c_str(), &num_groups);

    SecirModel model(num_groups);
    double read_buffer;
    status = tixiGetDoubleElement(handle, path_join(path, "StartDay").c_str(), &read_buffer);
    assert(status == SUCCESS && ("Failed to read StartDay at " + path).c_str());

    model.parameters.set<epi::StartDay>(read_buffer);
    model.parameters.set<epi::Seasonality>(* read_element(handle, path_join(path, "Seasonality"), io_mode));
    model.parameters.set<epi::ICUCapacity>(* read_element(handle, path_join(path, "ICUCapacity"), io_mode));
    model.parameters.set<epi::TestAndTraceCapacity>(* read_element(handle, path_join(path, "TestAndTraceCapacity"), io_mode));
    model.parameters.set<epi::ContactPatterns>(read_contact(handle, path_join(path, "ContactFreq"), io_mode));

    for (auto i = AgeGroup(0); i < AgeGroup(num_groups); i++) {
        auto group_name = "Group" + std::to_string((size_t)i + 1);
        auto group_path = path_join(path, group_name);

        // populations
        auto population_path = path_join(group_path, "Population");

        status = tixiGetDoubleElement(handle, path_join(population_path, "Dead").c_str(), &read_buffer);
        assert(status == SUCCESS && ("Failed to read number of deaths at " + path).c_str());

        model.populations[{i, epi::InfectionState::Dead}] = read_buffer;

        model.populations[{i, epi::InfectionState::Exposed}] =
            *read_element(handle, path_join(population_path, "Exposed"), io_mode);
        model.populations[{i, epi::InfectionState::Carrier}] =
            *read_element(handle, path_join(population_path, "Carrier"), io_mode);
        model.populations[{i, epi::InfectionState::Infected}] =
            *read_element(handle, path_join(population_path, "Infectious"), io_mode);
        model.populations[{i, epi::InfectionState::Hospitalized}] =
            *read_element(handle, path_join(population_path, "Hospitalized"), io_mode);
        model.populations[{i, epi::InfectionState::ICU}] =
            *read_element(handle, path_join(population_path, "ICU"), io_mode);
        model.populations[{i, epi::InfectionState::Recovered}] =
            *read_element(handle, path_join(population_path, "Recovered"), io_mode);

        status = tixiGetDoubleElement(handle, path_join(population_path, "Total").c_str(), &read_buffer);
        assert(status == SUCCESS && ("Failed to read total population at " + path).c_str());

        model.populations.set_difference_from_group_total<AgeGroup>({i, epi::InfectionState::Susceptible}, read_buffer);

        // times
        auto times_path = path_join(group_path, "StageTimes");

        model.parameters.get<IncubationTime>()[i] = * read_element(handle, path_join(times_path, "Incubation"), io_mode);
        model.parameters.get<InfectiousTimeMild>()[i] = * read_element(handle, path_join(times_path, "InfectiousMild"), io_mode);
        model.parameters.get<SerialInterval>()[i] = * read_element(handle, path_join(times_path, "SerialInterval"), io_mode);
        model.parameters.get<HospitalizedToHomeTime>()[i] = * read_element(handle, path_join(times_path, "HospitalizedToRecovered"), io_mode);
        model.parameters.get<HomeToHospitalizedTime>()[i] = * read_element(handle, path_join(times_path, "InfectiousToHospitalized"), io_mode);
        model.parameters.get<InfectiousTimeAsymptomatic>()[i] = * read_element(handle, path_join(times_path, "InfectiousAsympt"), io_mode);
        model.parameters.get<HospitalizedToICUTime>()[i] = * read_element(handle, path_join(times_path, "HospitalizedToICU"), io_mode);
        model.parameters.get<ICUToHomeTime>()[i] = * read_element(handle, path_join(times_path, "ICUToRecovered"), io_mode);
        model.parameters.get<ICUToDeathTime>()[i] = * read_element(handle, path_join(times_path, "ICUToDead"), io_mode);

        // probabilities
        auto probabilities_path = path_join(group_path, "Probabilities");

        model.parameters.get<InfectionProbabilityFromContact>()[i] = * read_element(handle, path_join(probabilities_path, "InfectedFromContact"), io_mode);
        model.parameters.get<RelativeCarrierInfectability>()[i] = * read_element(handle, path_join(probabilities_path, "Carrierinfectability"), io_mode);
        model.parameters.get<AsymptoticCasesPerInfectious>()[i] = * read_element(handle, path_join(probabilities_path, "AsympPerInfectious"), io_mode);
        model.parameters.get<RiskOfInfectionFromSympomatic>()[i] = * read_element(handle, path_join(probabilities_path, "RiskFromSymptomatic"), io_mode);
        model.parameters.get<MaxRiskOfInfectionFromSympomatic>()[i] = * read_element(handle, path_join(probabilities_path, "TestAndTraceMaxRiskFromSymptomatic"), io_mode);
        model.parameters.get<DeathsPerHospitalized>()[i] = * read_element(handle, path_join(probabilities_path, "DeadPerICU"), io_mode);
        model.parameters.get<HospitalizedCasesPerInfectious>()[i] = * read_element(handle, path_join(probabilities_path, "HospitalizedPerInfectious"), io_mode);
        model.parameters.get<ICUCasesPerHospitalized>()[i] = * read_element(handle, path_join(probabilities_path, "ICUPerHospitalized"), io_mode);
    }

    return model;

}

void write_parameter_space(TixiDocumentHandle handle, const std::string& path, SecirModel const& model, int num_runs,
                           int io_mode)
{
    auto num_groups = model.parameters.get_num_groups();
    tixiAddIntegerElement(handle, path.c_str(), "NumberOfGroups", (int)(size_t)num_groups, "%d");

    tixiAddDoubleElement(handle, path.c_str(), "StartDay", model.parameters.get<epi::StartDay>(), "%g");
    write_element(handle, path, "Seasonality", model.parameters.get<epi::Seasonality>(), io_mode, num_runs);
    write_element(handle, path, "ICUCapacity", model.parameters.get<epi::ICUCapacity>(), io_mode, num_runs);
    write_element(handle, path, "TestAndTraceCapacity", model.parameters.get<epi::TestAndTraceCapacity>(), io_mode,
                  num_runs);

    for (auto i = AgeGroup(0); i < AgeGroup(num_groups); i++) {
        auto group_name = "Group" + std::to_string((size_t)i + 1);
        auto group_path = path_join(path, group_name);

        tixiCreateElement(handle, path.c_str(), group_name.c_str());

        // populations
        auto population_path = path_join(group_path, "Population");
        tixiCreateElement(handle, group_path.c_str(), "Population");

        tixiAddDoubleElement(handle, population_path.c_str(), "Total",
                             model.populations.get_group_total(i), "%g");
        tixiAddDoubleElement(handle, population_path.c_str(), "Dead",
                             model.populations[{i, InfectionState::Dead}], "%g");
        write_element(handle, population_path, "Exposed",
                      model.populations[{i, InfectionState::Exposed}], io_mode, num_runs);
        write_element(handle, population_path, "Carrier",
                      model.populations[{i, InfectionState::Carrier}], io_mode, num_runs);
        write_element(handle, population_path, "Infectious",
                      model.populations[{i, InfectionState::Infected}], io_mode, num_runs);
        write_element(handle, population_path, "Hospitalized",
                      model.populations[{i, InfectionState::Hospitalized}], io_mode, num_runs);
        write_element(handle, population_path, "ICU",
                      model.populations[{i, InfectionState::ICU}], io_mode, num_runs);
        write_element(handle, population_path, "Recovered",
                      model.populations[{i, InfectionState::Recovered}], io_mode, num_runs);

        // times
        auto times_path = path_join(group_path, "StageTimes");
        tixiCreateElement(handle, group_path.c_str(), "StageTimes");

        write_element(handle, times_path, "Incubation", model.parameters.get<IncubationTime>()[i], io_mode, num_runs);
        write_element(handle, times_path, "InfectiousMild", model.parameters.get<InfectiousTimeMild>()[i], io_mode,
                      num_runs);
        write_element(handle, times_path, "SerialInterval", model.parameters.get<SerialInterval>()[i], io_mode,
                      num_runs);
        write_element(handle, times_path, "HospitalizedToRecovered",
                      model.parameters.get<HospitalizedToHomeTime>()[i], io_mode, num_runs);
        write_element(handle, times_path, "InfectiousToHospitalized",
                      model.parameters.get<HomeToHospitalizedTime>()[i], io_mode, num_runs);
        write_element(handle, times_path, "InfectiousAsympt", model.parameters.get<InfectiousTimeAsymptomatic>()[i], io_mode,
                      num_runs);
        write_element(handle, times_path, "HospitalizedToICU", model.parameters.get<HospitalizedToICUTime>()[i],
                      io_mode, num_runs);
        write_element(handle, times_path, "ICUToRecovered", model.parameters.get<ICUToHomeTime>()[i], io_mode,
                      num_runs);
        write_element(handle, times_path, "ICUToDead", model.parameters.get<ICUToDeathTime>()[i], io_mode, num_runs);

        // probabilities
        auto probabilities_path = path_join(group_path, "Probabilities");
        tixiCreateElement(handle, group_path.c_str(), "Probabilities");

        write_element(handle, probabilities_path, "InfectedFromContact",
                      model.parameters.get<InfectionProbabilityFromContact>()[i], io_mode, num_runs);
        write_element(handle, probabilities_path, "Carrierinfectability",
                      model.parameters.get<RelativeCarrierInfectability>()[i], io_mode, num_runs);
        write_element(handle, probabilities_path, "AsympPerInfectious",
                      model.parameters.get<AsymptoticCasesPerInfectious>()[i], io_mode, num_runs);
        write_element(handle, probabilities_path, "RiskFromSymptomatic",
                      model.parameters.get<RiskOfInfectionFromSympomatic>()[i], io_mode, num_runs);
        write_element(handle, probabilities_path, "TestAndTraceMaxRiskFromSymptomatic",
                      model.parameters.get<MaxRiskOfInfectionFromSympomatic>()[i], io_mode,
                      num_runs);
        write_element(handle, probabilities_path, "DeadPerICU", model.parameters.get<DeathsPerHospitalized>()[i],
                      io_mode, num_runs);
        write_element(handle, probabilities_path, "HospitalizedPerInfectious",
                      model.parameters.get<HospitalizedCasesPerInfectious>()[i], io_mode, num_runs);
        write_element(handle, probabilities_path, "ICUPerHospitalized",
                      model.parameters.get<ICUCasesPerHospitalized>()[i], io_mode, num_runs);
    }

    write_contact(handle, path, model.parameters.get<epi::ContactPatterns>(), io_mode);
}

ParameterStudy<SecirModel> read_parameter_study(TixiDocumentHandle handle, const std::string& path)
{
    ReturnCode status;

    int io_mode;
    int num_runs;
    double t0;
    double tmax;

    status = tixiGetIntegerElement(handle, path_join(path, "IOMode").c_str(), &io_mode);
    assert(status == SUCCESS && ("Failed to read io_mode at " + path).c_str());

    status = tixiGetIntegerElement(handle, path_join(path, "Runs").c_str(), &num_runs);
    assert(status == SUCCESS && ("Failed to read num_runs at " + path).c_str());

    status = tixiGetDoubleElement(handle, path_join(path, "T0").c_str(), &t0);
    assert(status == SUCCESS && ("Failed to read t0 at " + path).c_str());

    status = tixiGetDoubleElement(handle, path_join(path, "TMax").c_str(), &tmax);
    assert(status == SUCCESS && ("Failed to read tmax at " + path).c_str());

    unused(status);

    SecirModel model = read_parameter_space(handle, path, io_mode);
    return ParameterStudy<SecirModel>(model, t0, tmax, num_runs);
}

void write_parameter_study(TixiDocumentHandle handle, const std::string& path,
                           const ParameterStudy<SecirModel>& parameter_study, int io_mode)
{
    tixiAddIntegerElement(handle, path.c_str(), "IOMode", io_mode, "%d");
    tixiAddIntegerElement(handle, path.c_str(), "Runs", parameter_study.get_num_runs(), "%d");
    tixiAddDoubleElement(handle, path.c_str(), "T0", parameter_study.get_t0(), "%g");
    tixiAddDoubleElement(handle, path.c_str(), "TMax", parameter_study.get_tmax(), "%g");

    write_parameter_space(handle, path, parameter_study.get_model(), parameter_study.get_num_runs(), io_mode);
}

void write_single_run_params(const int run,
                             epi::Graph<epi::ModelNode<epi::Simulation<SecirModel>>, epi::MigrationEdge> graph, double t0,
                             double tmax)
{
    assert(graph.nodes().size() > 0 && "Graph Nodes are empty");

    std::string abs_path;
    bool created = create_directory("results", abs_path);

    if (created) {
        log_info("Results are stored in {:s}/results.", epi::get_current_dir_name());
    }
    else if (run == 0) {
        log_info(
            "Directory '{:s}' already exists. Results are stored in {:s}/results. Files from previous runs will be "
            "overwritten",
            epi::get_current_dir_name());
    }
    std::vector<TimeSeries<double>> all_results;
    std::vector<int> ids;

    ids.reserve(graph.nodes().size());
    all_results.reserve(graph.nodes().size());
    std::transform(graph.nodes().begin(), graph.nodes().end(), std::back_inserter(all_results), [](auto& node) {
        return node.property.get_result();
    });
    std::transform(graph.nodes().begin(), graph.nodes().end(), std::back_inserter(ids), [](auto& node) {
        return node.id;
    });

    int node_id = 0;
    for (auto& node : graph.nodes()) {
        int num_runs     = 1;
        std::string path = "/Parameters";
        TixiDocumentHandle handle;
        tixiCreateDocument("Parameters", &handle);
        ParameterStudy<SecirModel> study(node.property.get_simulation().get_model(), t0, tmax, num_runs);

        write_parameter_study(handle, path, study);

        tixiSaveDocument(handle, path_join(abs_path, ("Parameters_run" + std::to_string(run) + "_node" +
                                                      std::to_string(node_id) + ".xml"))
                                     .c_str());
        tixiCloseDocument(handle);
        node_id++;
    }

    save_result(all_results, ids,
                path_join(abs_path, ("Results_run" + std::to_string(run) + std::to_string(node_id) + ".h5")));
}

void write_contact_frequency_matrix_collection(TixiDocumentHandle handle, const std::string& path,
                                               const ContactMatrixGroup& cfmc)
{
    tixiCreateElement(handle, path.c_str(), "ContactMatrixGroup");
    auto collection_path = path_join(path, "ContactMatrixGroup");
    for (size_t i = 0; i < cfmc.get_num_matrices(); ++i) {
        auto& cfm     = cfmc[i];
        auto cfm_name = "ContactMatrix" + std::to_string(i + 1);
        tixiCreateElement(handle, collection_path.c_str(), cfm_name.c_str());
        auto cfm_path = path_join(collection_path, cfm_name);
        write_matrix(handle, cfm_path, "Baseline", cfm.get_baseline());
        write_matrix(handle, cfm_path, "Minimum", cfm.get_minimum());
        tixiCreateElement(handle, cfm_path.c_str(), "Dampings");
        auto dampings_path = path_join(cfm_path, "Dampings");
        for (size_t j = 0; j < cfm.get_dampings().size(); ++j) {
            auto& damping     = cfm.get_dampings()[j];
            auto damping_name = "Damping" + std::to_string(j + 1);
            tixiCreateElement(handle, dampings_path.c_str(), damping_name.c_str());
            auto damping_path = path_join(dampings_path, damping_name);
            tixiAddDoubleAttribute(handle, damping_path.c_str(), "Time", double(damping.get_time()), "%.18g");
            tixiAddIntegerAttribute(handle, damping_path.c_str(), "Type", int(damping.get_type()), "%d");
            tixiAddIntegerAttribute(handle, damping_path.c_str(), "Level", int(damping.get_level()), "%d");
            write_matrix(handle, damping_path, "Values", damping.get_coeffs());
        }
    }
}

ContactMatrixGroup read_contact_frequency_matrix_collection(TixiDocumentHandle handle, const std::string& path)
{
    auto status = SUCCESS;
    unused(status);

    auto collection_path = path_join(path, "ContactMatrixGroup");
    int num_matrices;
    status = tixiGetNumberOfChilds(handle, collection_path.c_str(), &num_matrices);
    assert(status == SUCCESS && "Failed to read ContactMatrixGroup.");
    ContactMatrixGroup cfmc{1, size_t(num_matrices)};
    for (size_t i = 0; i < cfmc.get_num_matrices(); ++i) {
        auto cfm_path      = path_join(collection_path, "ContactMatrix" + std::to_string(i + 1));
        cfmc[i]            = ContactMatrix(read_matrix<>(handle, path_join(cfm_path, "Baseline")),
                                read_matrix<>(handle, path_join(cfm_path, "Minimum")));
        auto dampings_path = path_join(cfm_path, "Dampings");
        int num_dampings;
        status = tixiGetNumberOfChilds(handle, dampings_path.c_str(), &num_dampings);
        assert(status == SUCCESS && "Failed to read Dampings from ContactMatrix.");
        for (int j = 0; j < num_dampings; ++j) {
            auto damping_path = path_join(dampings_path, "Damping" + std::to_string(j + 1));
            double t;
            status = tixiGetDoubleAttribute(handle, damping_path.c_str(), "Time", &t);
            assert(status == SUCCESS && "Failed to read Damping Time.");
            int type;
            status = tixiGetIntegerAttribute(handle, damping_path.c_str(), "Type", &type);
            assert(status == SUCCESS && "Failed to read Damping Type.");
            int level;
            status = tixiGetIntegerAttribute(handle, damping_path.c_str(), "Level", &level);
            assert(status == SUCCESS && "Failed to read Damping Level.");
            cfmc[i].add_damping(read_matrix<>(handle, path_join(damping_path, "Values")), DampingLevel(level),
                                DampingType(type), SimulationTime(t));
        }
    }
    return cfmc;
}

void write_contact(TixiDocumentHandle handle, const std::string& path, const UncertainContactMatrix& contact_pattern,
                   int io_mode)
{
    tixiCreateElement(handle, path.c_str(), "ContactFreq");
    auto contact_path = path_join(path, "ContactFreq");

    write_damping_matrix_expression_collection(handle, contact_path, contact_pattern.get_cont_freq_mat());

    if (io_mode == 1 || io_mode == 2 || io_mode == 3) {
        write_distribution(handle, contact_path, "NumDampings", *contact_pattern.get_distribution_damp_nb().get());
        write_distribution(handle, contact_path, "DampingDay", *contact_pattern.get_distribution_damp_days().get());
        write_distribution(handle, contact_path, "DampingDiagBase",
                           *contact_pattern.get_distribution_damp_diag_base().get());
        write_distribution(handle, contact_path, "DampingDiagRel",
                           *contact_pattern.get_distribution_damp_diag_rel().get());
        write_distribution(handle, contact_path, "DampingOffdiagRel",
                           *contact_pattern.get_distribution_damp_offdiag_rel().get());
    }
}

UncertainContactMatrix read_contact(TixiDocumentHandle handle, const std::string& path, int io_mode)
{
    UncertainContactMatrix contact_patterns{read_damping_matrix_expression_collection(handle, path)};

    if (io_mode == 1 || io_mode == 2 || io_mode == 3) {
        contact_patterns.set_distribution_damp_nb(*read_distribution(handle, path_join(path, "NumDampings")));
        contact_patterns.set_distribution_damp_days(*read_distribution(handle, path_join(path, "DampingDay")));
        contact_patterns.set_distribution_damp_diag_base(
            *read_distribution(handle, path_join(path, "DampingDiagBase")));
        contact_patterns.set_distribution_damp_diag_rel(*read_distribution(handle, path_join(path, "DampingDiagRel")));
        contact_patterns.set_distribution_damp_offdiag_rel(
            *read_distribution(handle, path_join(path, "DampingOffdiagRel")));
    }
    return contact_patterns;
}

std::vector<int> get_county_ids(const std::string& path)
{
    Json::Reader reader;
    Json::Value root;

    std::vector<int> id;

    std::ifstream census(path_join(path, "county_current_population.json"));
    if (!reader.parse(census, root)) {
        log_error(reader.getFormattedErrorMessages());
        return id;
    }

    for (unsigned int i = 0; i < root.size(); i++) {
        id.push_back(root[i]["ID_County"].asInt());
    }

    return id;
}

} // namespace epi
