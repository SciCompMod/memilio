/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

#include "ode_secirvvs/parameters_io.h"
#include "ode_secirvvs/parameters.h"
#include "memilio/epidemiology/regions.h"
#include "memilio/io/io.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "memilio/io/epi_data.h"
#include "memilio/utils/memory.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/utils/stl_util.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/mobility.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/date.h"
#include "memilio/io/mobility_io.h"

#include <boost/filesystem.hpp>

#include <numeric>
#include <vector>
#include <iostream>
#include <string>
#include <random>
#include <fstream>

namespace mio
{
namespace osecirvvs
{
namespace details
{
//gets the county or state id of the entry if available, 0 (for whole country) otherwise
//used for comparisons of entry to integer region id
template <class EpiDataEntry>
int get_region_id(const EpiDataEntry& rki_entry)
{
    return rki_entry.county_id ? rki_entry.county_id->get() : (rki_entry.state_id ? rki_entry.state_id->get() : 0);
}

IOResult<void> read_confirmed_cases_data(
    std::string const& path, std::vector<int> const& vregion, Date date, std::vector<std::vector<double>>& vnum_Exposed,
    std::vector<std::vector<double>>& vnum_InfectedNoSymptoms, std::vector<std::vector<double>>& vnum_InfectedSymptoms,
    std::vector<std::vector<double>>& vnum_InfectedSevere, std::vector<std::vector<double>>& vnum_icu,
    std::vector<std::vector<double>>& vnum_death, std::vector<std::vector<double>>& vnum_rec,
    std::vector<std::vector<double>>& vnum_timm_i, const std::vector<std::vector<int>>& vt_Exposed,
    const std::vector<std::vector<int>>& vt_InfectedNoSymptoms,
    const std::vector<std::vector<int>>& vt_InfectedSymptoms, const std::vector<std::vector<int>>& vt_InfectedSevere,
    const std::vector<std::vector<int>>& vt_InfectedCritical, const std::vector<std::vector<int>>& vt_imm_interval_i,
    const std::vector<std::vector<double>>& vmu_C_R, const std::vector<std::vector<double>>& vmu_I_H,
    const std::vector<std::vector<double>>& vmu_H_U, const std::vector<double>& scaling_factor_inf)
{
    BOOST_OUTCOME_TRY(rki_data, mio::read_confirmed_cases_data(path));
    return read_confirmed_cases_data(
        rki_data, vregion, date, vnum_Exposed, vnum_InfectedNoSymptoms, vnum_InfectedSymptoms, vnum_InfectedSevere,
        vnum_icu, vnum_death, vnum_rec, vnum_timm_i, vt_Exposed, vt_InfectedNoSymptoms, vt_InfectedSymptoms,
        vt_InfectedSevere, vt_InfectedCritical, vt_imm_interval_i, vmu_C_R, vmu_I_H, vmu_H_U, scaling_factor_inf);
}

IOResult<void> read_confirmed_cases_data(
    const std::vector<ConfirmedCasesDataEntry>& rki_data, std::vector<int> const& vregion, Date date,
    std::vector<std::vector<double>>& vnum_Exposed, std::vector<std::vector<double>>& vnum_InfectedNoSymptoms,
    std::vector<std::vector<double>>& vnum_InfectedSymptoms, std::vector<std::vector<double>>& vnum_InfectedSevere,
    std::vector<std::vector<double>>& vnum_icu, std::vector<std::vector<double>>& vnum_death,
    std::vector<std::vector<double>>& vnum_rec, std::vector<std::vector<double>>& vnum_timm_i,
    const std::vector<std::vector<int>>& vt_Exposed, const std::vector<std::vector<int>>& vt_InfectedNoSymptoms,
    const std::vector<std::vector<int>>& vt_InfectedSymptoms, const std::vector<std::vector<int>>& vt_InfectedSevere,
    const std::vector<std::vector<int>>& vt_InfectedCritical, const std::vector<std::vector<int>>& vt_imm_interval_i,
    const std::vector<std::vector<double>>& vmu_C_R, const std::vector<std::vector<double>>& vmu_I_H,
    const std::vector<std::vector<double>>& vmu_H_U, const std::vector<double>& scaling_factor_inf)
{
    auto max_date_entry = std::max_element(rki_data.begin(), rki_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == rki_data.end()) {
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

    auto vnum_InfectedSymptoms_timm_i = vnum_InfectedSymptoms;
    auto vnum_Recovered_timm_i        = vnum_rec;
    auto vnum_Exposed_timm_i          = vnum_Exposed;
    auto vnum_InfectedSevere_timm_i   = vnum_InfectedSevere;
    auto vnum_Death_timm_i            = vnum_death;
    auto vnum_InfectedCritical_timm_i = vnum_icu;

    for (auto&& entry : rki_data) {
        auto it = std::find_if(vregion.begin(), vregion.end(), [&entry](auto r) {
            return r == 0 || get_region_id(entry) == r;
        });
        if (it != vregion.end()) {
            auto region_idx = size_t(it - vregion.begin());

            auto& t_Exposed            = vt_Exposed[region_idx];
            auto& t_InfectedNoSymptoms = vt_InfectedNoSymptoms[region_idx];
            auto& t_InfectedSymptoms   = vt_InfectedSymptoms[region_idx];
            auto& t_InfectedSevere     = vt_InfectedSevere[region_idx];
            auto& t_InfectedCritical   = vt_InfectedCritical[region_idx];
            auto& t_imm_interval_i     = vt_imm_interval_i[region_idx];

            auto& num_InfectedNoSymptoms = vnum_InfectedNoSymptoms[region_idx];
            auto& num_InfectedSymptoms   = vnum_InfectedSymptoms[region_idx];
            auto& num_rec                = vnum_rec[region_idx];
            auto& num_Exposed            = vnum_Exposed[region_idx];
            auto& num_InfectedSevere     = vnum_InfectedSevere[region_idx];
            auto& num_death              = vnum_death[region_idx];
            auto& num_icu                = vnum_icu[region_idx];

            auto& num_InfectedSymptoms_timm_i = vnum_InfectedSymptoms_timm_i[region_idx];
            auto& num_Recovered_timm_i        = vnum_Recovered_timm_i[region_idx];
            auto& num_InfectedSevere_timm_i   = vnum_InfectedSevere_timm_i[region_idx];
            auto& num_Death_timm_i            = vnum_Death_timm_i[region_idx];
            auto& num_InfectedCritical_timm_i = vnum_InfectedCritical_timm_i[region_idx];

            auto& mu_C_R = vmu_C_R[region_idx];
            auto& mu_I_H = vmu_I_H[region_idx];
            auto& mu_H_U = vmu_H_U[region_idx];

            bool read_icu = false; // params.populations.get({age, SecirCompartments::U}) == 0;

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
                if (read_icu) {
                    num_icu[age] += mu_I_H[age] * mu_H_U[age] * scaling_factor_inf[age] * entry.num_confirmed;
                }
            }
            if (entry.date ==
                offset_date_by_days(date, -t_InfectedSymptoms[age] - t_InfectedSevere[age] - t_InfectedCritical[age])) {
                num_death[age] += entry.num_deaths;
                if (read_icu) {
                    num_icu[age] -= mu_I_H[age] * mu_H_U[age] * scaling_factor_inf[age] * entry.num_confirmed;
                }
            }
            // calc all values again for t0 - t_imm_interval_i
            if (entry.date == offset_date_by_days(date, 0 - t_imm_interval_i[age])) {
                num_InfectedSymptoms_timm_i[age] += scaling_factor_inf[age] * entry.num_confirmed;
                num_Recovered_timm_i[age] += entry.num_confirmed;
            }
            if (entry.date == offset_date_by_days(date, -t_InfectedSymptoms[age] - t_imm_interval_i[age])) {
                num_InfectedSymptoms_timm_i[age] -= scaling_factor_inf[age] * entry.num_confirmed;
                num_InfectedSevere_timm_i[age] += mu_I_H[age] * scaling_factor_inf[age] * entry.num_confirmed;
            }
            if (entry.date ==
                offset_date_by_days(date, -t_InfectedSymptoms[age] - t_InfectedSevere[age] - t_imm_interval_i[age])) {
                num_InfectedSevere_timm_i[age] -= mu_I_H[age] * scaling_factor_inf[age] * entry.num_confirmed;
                if (read_icu) {
                    num_InfectedCritical_timm_i[age] +=
                        mu_I_H[age] * mu_H_U[age] * scaling_factor_inf[age] * entry.num_confirmed;
                }
            }
            if (entry.date == offset_date_by_days(date, -t_InfectedSymptoms[age] - t_InfectedSevere[age] -
                                                            t_InfectedCritical[age] - t_imm_interval_i[age])) {
                num_Death_timm_i[age] += entry.num_deaths;
                if (read_icu) {
                    num_InfectedCritical_timm_i[age] -=
                        mu_I_H[age] * mu_H_U[age] * scaling_factor_inf[age] * entry.num_confirmed;
                }
            }
        }
    }

    for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
        auto region = vregion[region_idx];

        auto& num_InfectedNoSymptoms = vnum_InfectedNoSymptoms[region_idx];
        auto& num_InfectedSymptoms   = vnum_InfectedSymptoms[region_idx];
        auto& num_rec                = vnum_rec[region_idx];
        auto& num_Exposed            = vnum_Exposed[region_idx];
        auto& num_InfectedSevere     = vnum_InfectedSevere[region_idx];
        auto& num_death              = vnum_death[region_idx];
        auto& num_icu                = vnum_icu[region_idx];

        auto& num_InfectedSymptoms_timm_i = vnum_InfectedSymptoms_timm_i[region_idx];
        auto& num_Recovered_timm_i        = vnum_Recovered_timm_i[region_idx];
        auto& num_InfectedSevere_timm_i   = vnum_InfectedSevere_timm_i[region_idx];
        auto& num_Death_timm_i            = vnum_Death_timm_i[region_idx];
        auto& num_InfectedCritical_timm_i = vnum_InfectedCritical_timm_i[region_idx];

        auto& num_timm_i = vnum_timm_i[region_idx];

        size_t num_groups = ConfirmedCasesDataEntry::age_group_names.size();
        for (size_t i = 0; i < num_groups; i++) {
            // subtract infected confirmed cases which are not yet recovered
            // and remove dark figure scaling factor
            num_rec[i] -= num_InfectedSymptoms[i] / scaling_factor_inf[i];
            num_rec[i] -= num_InfectedSevere[i] / scaling_factor_inf[i];
            num_rec[i] -=
                num_icu[i] /
                scaling_factor_inf[i]; // TODO: this has to be adapted for scaling_factor_inf != 1 or != ***_icu
            num_rec[i] -= num_death[i] / scaling_factor_inf[i];

            num_Recovered_timm_i[i] -= num_InfectedSymptoms_timm_i[i] / scaling_factor_inf[i];
            num_Recovered_timm_i[i] -= num_InfectedSevere_timm_i[i] / scaling_factor_inf[i];
            num_Recovered_timm_i[i] -=
                num_InfectedCritical_timm_i[i] /
                scaling_factor_inf[i]; // TODO: this has to be adapted for scaling_factor_inf != 1 or != ***_icu
            num_Recovered_timm_i[i] -= num_Death_timm_i[i] / scaling_factor_inf[i];

            num_timm_i[i] = num_rec[i] - num_Recovered_timm_i[i];

            // Since we are calculating the compartments over the Confirmed tests, if the incidence decreases, it may result in a negative value here.
            // if this happens, set zero.
            if (num_timm_i[i] < 0) {
                num_timm_i[i] = 0;
            }

            auto try_fix_constraints = [region, i](double& value, double error, auto str) {
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
    }

    return success();
}

IOResult<void> read_confirmed_cases_data_fix_recovered(std::string const& path, std::vector<int> const& vregion,
                                                       Date date, std::vector<std::vector<double>>& vnum_rec,
                                                       double delay)
{
    BOOST_OUTCOME_TRY(rki_data, mio::read_confirmed_cases_data(path));
    return read_confirmed_cases_data_fix_recovered(rki_data, vregion, date, vnum_rec, delay);
}

IOResult<void> read_confirmed_cases_data_fix_recovered(const std::vector<ConfirmedCasesDataEntry>& rki_data,
                                                       std::vector<int> const& vregion, Date date,
                                                       std::vector<std::vector<double>>& vnum_rec, double delay)
{
    auto max_date_entry = std::max_element(rki_data.begin(), rki_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == rki_data.end()) {
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

    for (auto&& rki_entry : rki_data) {
        auto it = std::find_if(vregion.begin(), vregion.end(), [&rki_entry](auto r) {
            return r == 0 || get_region_id(rki_entry) == r;
        });
        if (it != vregion.end()) {
            auto region_idx = size_t(it - vregion.begin());
            if (rki_entry.date == offset_date_by_days(date, int(-delay))) {
                vnum_rec[region_idx][size_t(rki_entry.age_group)] = rki_entry.num_confirmed;
            }
        }
    }

    for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
        auto region   = vregion[region_idx];
        auto& num_rec = vnum_rec[region_idx];

        size_t num_groups = ConfirmedCasesDataEntry::age_group_names.size();
        for (size_t i = 0; i < num_groups; i++) {
            auto try_fix_constraints = [region, i](double& value, double error, auto str) {
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
    }

    return success();
}

IOResult<void> read_divi_data(const std::string& path, const std::vector<int>& vregion, Date date,
                              std::vector<double>& vnum_icu)
{
    BOOST_OUTCOME_TRY(divi_data, mio::read_divi_data(path));
    return read_divi_data(divi_data, vregion, date, vnum_icu);
}

IOResult<void> read_divi_data(const std::vector<DiviEntry>& divi_data, const std::vector<int>& vregion, Date date,
                              std::vector<double>& vnum_icu)
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

IOResult<std::vector<std::vector<double>>> read_immunity_population(const std::string& path, const int& num_age_groups)
{
    std::vector<std::vector<double>> ans(num_age_groups, std::vector<double>(3, 0.0));
    std::fstream immunity_file;
    immunity_file.open(path, std::ios::in);
    if (immunity_file.fail()) { // checks to see if file opended
        return mio::failure(mio::StatusCode::InvalidValue, "Failed to open immunity_population.txt.");
    }
    if (immunity_file.is_open()) {
        // std::vector<double> tp;
        std::string tp;
        int linenumber = 0;
        while (linenumber < num_age_groups) {
            getline(immunity_file, tp);
            // delete /r at the end by delete last entry.
            tp.erase(tp.size() - 1);
            auto line          = split(tp, ' ');
            ans[linenumber][0] = std::stod(line[0]);
            ans[linenumber][1] = std::stod(line[1]);
            ans[linenumber][2] = std::stod(line[2]);
            linenumber++;
        }
        immunity_file.close(); //close the file object.
    }
    return ans;
}

IOResult<std::vector<std::vector<double>>> read_population_data(const std::string& path,
                                                                const std::vector<int>& vregion)
{
    BOOST_OUTCOME_TRY(population_data, mio::read_population_data(path));
    return read_population_data(population_data, vregion);
}

IOResult<std::vector<std::vector<double>>> read_population_data(const std::vector<PopulationDataEntry>& population_data,
                                                                const std::vector<int>& vregion)
{
    std::vector<std::vector<double>> vnum_population(
        vregion.size(), std::vector<double>(ConfirmedCasesDataEntry::age_group_names.size(), 0.0));

    for (auto&& county_entry : population_data) {
        //accumulate population of states or country from population of counties
        if (!county_entry.county_id) {
            return failure(StatusCode::InvalidFileFormat, "File with county population expected.");
        }
        //find region that this county belongs to
        //all counties belong to the country (id = 0)
        auto it = std::find_if(vregion.begin(), vregion.end(), [&county_entry](auto r) {
            return r == 0 ||
                   (county_entry.county_id &&
                    regions::de::StateId(r) == regions::de::get_state_id(*county_entry.county_id)) ||
                   (county_entry.county_id && regions::de::CountyId(r) == *county_entry.county_id);
        });
        if (it != vregion.end()) {
            auto region_idx      = size_t(it - vregion.begin());
            auto& num_population = vnum_population[region_idx];
            for (size_t age = 0; age < num_population.size(); age++) {
                num_population[age] += county_entry.population[AgeGroup(age)];
            }
        }
    }

    return success(vnum_population);
}

IOResult<void> set_vaccination_data(std::vector<Model>& model, const std::string& path, Date date,
                                    const std::vector<int>& vregion, int num_days)
{
    BOOST_OUTCOME_TRY(vacc_data, read_vaccination_data(path));

    auto num_groups = model[0].parameters.get_num_groups();

    auto days_until_effective1 = (int)(double)model[0].parameters.get<DaysUntilEffectivePartialImmunity>()[AgeGroup(0)];
    auto days_until_effective2 =
        (int)(double)model[0].parameters.get<DaysUntilEffectiveImprovedImmunity>()[AgeGroup(0)];
    auto vaccination_distance = (int)(double)model[0].parameters.get<VaccinationGap>()[AgeGroup(0)];

    // calculating the current rate of vaccinations based on the last (num_days_vacc_rate) days of simulation.
    // scaling by the number of population is done later, when population is set
    const int num_days_vacc_rate = 30;

    // iterate over regions (e.g., counties)
    for (size_t i = 0; i < model.size(); ++i) {
        // iterate over age groups in region
        for (auto g = AgeGroup(0); g < num_groups; ++g) {

            model[i].parameters.template get<DailyPartialVaccination>().resize(SimulationDay(num_days + 1));
            model[i].parameters.template get<DailyFullVaccination>().resize(SimulationDay(num_days + 1));
            for (auto d = SimulationDay(0); d < SimulationDay(num_days + 1); ++d) {
                model[i].parameters.template get<DailyPartialVaccination>()[{g, d}] = 0.0;
                model[i].parameters.template get<DailyFullVaccination>()[{g, d}]    = 0.0;
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
    auto max_full_date =
        offset_date_by_days(max_date, days_until_effective1 - days_until_effective2 - vaccination_distance);

    for (auto&& vacc_data_entry : vacc_data) {
        auto it      = std::find_if(vregion.begin(), vregion.end(), [&vacc_data_entry](auto&& r) {
            return r == 0 || (vacc_data_entry.county_id && vacc_data_entry.county_id == regions::de::CountyId(r)) ||
                   (vacc_data_entry.state_id && vacc_data_entry.state_id == regions::de::StateId(r));
        });
        auto date_df = vacc_data_entry.date;
        if (it != vregion.end()) {
            auto region_idx = size_t(it - vregion.begin());
            auto age        = vacc_data_entry.age_group;

            if (date_df == offset_date_by_days(date, -days_until_effective1)) {
                model[region_idx].parameters.template get<VaccinationTemporaryImm1>()[{age}] =
                    model[region_idx].parameters.template get<VaccinationTemporaryImm1>()[{age}] +
                    vacc_data_entry.num_first_vaccinations_completed +
                    vacc_data_entry.num_second_vaccinations_completed;
            }

            if (date_df ==
                offset_date_by_days(date, -days_until_effective1 -
                                              model[region_idx].parameters.template get<ImmunityInterval2>()[{age}])) {
                model[region_idx].parameters.template get<VaccinationTemporaryImm1>()[{age}] =
                    model[region_idx].parameters.template get<VaccinationTemporaryImm1>()[{age}] -
                    vacc_data_entry.num_first_vaccinations_completed -
                    vacc_data_entry.num_second_vaccinations_completed;
            }

            if (date_df == offset_date_by_days(date, -days_until_effective2)) {
                model[region_idx].parameters.template get<VaccinationTemporaryImm2>()[{age}] +=
                    vacc_data_entry.num_vaccinations_refreshed + vacc_data_entry.num_vaccinations_refreshed_2;
            }

            if (date_df ==
                offset_date_by_days(date, -days_until_effective2 -
                                              model[region_idx].parameters.template get<ImmunityInterval2>()[{age}])) {
                model[region_idx].parameters.template get<VaccinationTemporaryImm2>()[{age}] -=
                    vacc_data_entry.num_vaccinations_refreshed - vacc_data_entry.num_vaccinations_refreshed_2;
            }

            if (date_df == offset_date_by_days(date, -num_days_vacc_rate)) {
                model[region_idx].parameters.template get<RateOfDailyBoosterVaccinations>()[{age}] -=
                    vacc_data_entry.num_vaccinations_refreshed_2 / num_days_vacc_rate;
                model[region_idx].parameters.template get<RateOfDailyImprovedVaccinations>()[{age}] -=
                    vacc_data_entry.num_vaccinations_refreshed / num_days_vacc_rate;
                model[region_idx].parameters.template get<RateOfDailyPartialVaccinations>()[{age}] -=
                    (vacc_data_entry.num_first_vaccinations_completed +
                     vacc_data_entry.num_second_vaccinations_completed) /
                    num_days_vacc_rate;
            }

            if (date_df == date) {
                model[region_idx].parameters.template get<RateOfDailyBoosterVaccinations>()[{age}] +=
                    vacc_data_entry.num_vaccinations_refreshed_2 / num_days_vacc_rate;
                model[region_idx].parameters.template get<RateOfDailyImprovedVaccinations>()[{age}] +=
                    vacc_data_entry.num_vaccinations_refreshed / num_days_vacc_rate;
                model[region_idx].parameters.template get<RateOfDailyPartialVaccinations>()[{age}] +=
                    (vacc_data_entry.num_first_vaccinations_completed +
                     vacc_data_entry.num_second_vaccinations_completed) /
                    num_days_vacc_rate;
            }

            for (size_t d = 0; d < (size_t)num_days + 1; ++d) {
                // int days_plus;
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
                auto offset_first_date = offset_date_by_days(date, (int)d - days_until_effective1);
                if (max_full_date >= offset_first_date) {
                    // Option 1: considered offset_first_date is available in input data frame
                    if (date_df == offset_first_date) {
                        model[region_idx].parameters.template get<DailyPartialVaccination>()[{age, SimulationDay(d)}] =
                            vacc_data_entry.num_first_vaccinations_completed +
                            vacc_data_entry.num_second_vaccinations_completed;
                    }
                }
                // else { // offset_first_date > max_date
                //     // Option 2: considered offset_first_date is NOT available in input data frame
                //     // Here, a constant number of first and second doses is assumed, i.e.,
                //     // the the number of vaccinationes at day d (N days after max_date) will be:
                //     // total number of vaccinations up to day max_date + N * number of vaccinations ON max_date
                //     // (where the latter is computed as the difference between the total number at max_date and max_date-1)
                //     days_plus = get_offset_in_days(offset_first_date, max_date);
                //     if (date_df == offset_date_by_days(max_date, -1)) {
                //         model[region_idx].parameters.template get<DailyFirstVaccination>()[{age, SimulationDay(d)}] -=
                //             (days_plus - 1) * vacc_data_entry.num_vaccinations_completed;
                //     }
                //     else if (date_df == max_date) {
                //         model[region_idx].parameters.template get<DailyFirstVaccination>()[{age, SimulationDay(d)}] +=
                //             days_plus * vacc_data_entry.num_vaccinations_completed;
                //     }
                // }

                // a person whose second dose is reported at start_date + simulation_day - days_until_effective2
                // has the full protection of the second dose at day X = start_date + simulation_day
                // Storing its value in get<DailyFullVaccination>() will eventually (in the simulation)
                // transfer the difference (between get<DailyFullVaccination>() at d and d-1) of
                // N susceptible, partially vaccinated individuals to 'SusceptibleImprovedImmunity' state at day d; see secir_vaccinated.h
                auto offset_full_date = offset_date_by_days(date, (int)d - days_until_effective2);
                if (max_full_date >= offset_full_date) {
                    // Option 1: considered offset_full_date is available in input data frame
                    if (date_df == offset_full_date) {
                        model[region_idx].parameters.template get<DailyFullVaccination>()[{age, SimulationDay(d)}] =
                            vacc_data_entry.num_vaccinations_refreshed + vacc_data_entry.num_vaccinations_refreshed_2;
                    }
                }
                // else { // offset_full_date > max_full_date
                //     // Option 2: considered offset_full_date is NOT available in input data frame
                //     days_plus = get_offset_in_days(offset_full_date, max_date);
                //     if (date_df == offset_date_by_days(max_full_date, -1)) {
                //         model[region_idx].parameters.template get<DailyFullVaccination>()[{age, SimulationDay(d)}] -=
                //             (days_plus - 1) * vacc_data_entry.num_vaccinations_completed;
                //     }
                //     else if (date_df == max_full_date) {
                //         model[region_idx].parameters.template get<DailyFullVaccination>()[{age, SimulationDay(d)}] +=
                //             days_plus * vacc_data_entry.num_vaccinations_completed;
                //     }
                // }
            }
        }
    }
    return success();
}

IOResult<void> set_vaccination_rate(std::vector<Model>& model, int num_days, int vacc_campaign_szenario)
{
    auto num_groups = model[0].parameters.get_num_groups();

    // the average number of vaccinations in the last 30 days was already set_vaccination_data and stored in RateOfDaily...Vaccinations.
    // Since we called the function before the population was set, we still need to divide by the total numbers of people in the age group.
    for (size_t i = 0; i < model.size(); ++i) {
        for (auto g = AgeGroup(0); g < num_groups; ++g) {
            auto group_total = model[i].populations.get_group_total(g);
            model[i].parameters.template get<RateOfDailyBoosterVaccinations>()[{g}] /= group_total;
            model[i].parameters.template get<RateOfDailyImprovedVaccinations>()[{g}] /= group_total;
            model[i].parameters.template get<RateOfDailyPartialVaccinations>()[{g}] /= group_total;

            if (model[i].parameters.template get<RateOfDailyBoosterVaccinations>()[{g}] > 0.001 ||
                model[i].parameters.template get<RateOfDailyImprovedVaccinations>()[{g}] > 0.001 ||
                model[i].parameters.template get<RateOfDailyPartialVaccinations>()[{g}] > 0.001) {
                log_info("Rate of daily vaccinations unusually high");
            }
        }
    }

    // iterate over regions (e.g., counties)
    for (size_t i = 0; i < model.size(); ++i) {
        // iterate over age groups in region
        for (auto g = AgeGroup(0); g < num_groups; ++g) {
            model[i].parameters.template get<DailyPartialVaccination>().resize(SimulationDay(num_days + 1));
            model[i].parameters.template get<DailyFullVaccination>().resize(SimulationDay(num_days + 1));
            model[i].parameters.template get<DailyBoosterVaccination>().resize(SimulationDay(num_days + 1));

            // Vaccine campaing start
            auto start_simulation = Date(2022, 6, 1);
            auto campaign_start   = Date(2022, 6, 3);

            auto days_till_start = get_offset_in_days(campaign_start, start_simulation);
            // to calculate the number of daily vaccinations when the campaign begins, we need to know the number of possible vaccinations.
            auto num_partial_vaccinations =
                model[i].populations[{g, InfectionState::SusceptibleNaive}] *
                pow((1 - model[i].parameters.template get<RateOfDailyPartialVaccinations>()[{g}]), days_till_start);
            auto num_possible_booster =
                model[i].populations[{g, InfectionState::SusceptiblePartialImmunity}] *
                    pow((1 - model[i].parameters.template get<RateOfDailyImprovedVaccinations>()[{g}]),
                        days_till_start) +
                model[i].populations[{g, InfectionState::SusceptibleImprovedImmunity}] *
                    pow((1 - model[i].parameters.template get<RateOfDailyBoosterVaccinations>()[{g}]), days_till_start);
            auto num_possible_vacc = num_partial_vaccinations + num_possible_booster;
            // number of vaccinations per day, assuming we vaccinate 5% of the possible population per week
            auto vacc_per_day         = num_possible_vacc * 0.05 / 7;
            auto booster_vacc_per_day = num_possible_booster * 0.05 / 7;
            double campaign_duration  = 0;
            // get the number of days, we need to either vaccinate half of the boostered persons or all already boostered persons.
            if (vacc_campaign_szenario == 2) {
                campaign_duration = (0.5 * num_possible_booster) / booster_vacc_per_day;
            }
            else if (vacc_campaign_szenario == 3) {
                campaign_duration = num_possible_booster / booster_vacc_per_day;
            }
            auto campaign_end = offset_date_by_days(campaign_start, campaign_duration);

            auto vacc_daily_partially =
                model[i].populations[{g, InfectionState::SusceptibleNaive}] / num_possible_vacc * vacc_per_day;
            auto vacc_daily_first_booster = model[i].populations[{g, InfectionState::SusceptiblePartialImmunity}] /
                                            num_possible_booster * booster_vacc_per_day;
            auto vacc_daily_more_booster = model[i].populations[{g, InfectionState::SusceptibleImprovedImmunity}] /
                                           num_possible_booster * booster_vacc_per_day;

            auto day_campaign_idx = 0;
            for (size_t d = 0; d < (size_t)num_days + 1; ++d) {
                auto date = offset_date_by_days(start_simulation, d);
                if (vacc_campaign_szenario > 1 && date >= campaign_start && date <= campaign_end) {

                    // get number of people in compartment left
                    auto number_S_left = model[i].populations[{g, InfectionState::SusceptibleNaive}] -
                                         day_campaign_idx * vacc_daily_partially;
                    auto number_Sp_left = model[i].populations[{g, InfectionState::SusceptiblePartialImmunity}] -
                                          day_campaign_idx * vacc_daily_first_booster;
                    auto number_R_left = model[i].populations[{g, InfectionState::SusceptibleImprovedImmunity}] -
                                         day_campaign_idx * vacc_daily_more_booster;

                    model[i].parameters.template get<DailyPartialVaccination>()[{g, SimulationDay{d}}] =
                        vacc_daily_partially / number_S_left;
                    model[i].parameters.template get<DailyFullVaccination>()[{g, SimulationDay{d}}] =
                        vacc_daily_first_booster / number_Sp_left;
                    model[i].parameters.template get<DailyBoosterVaccination>()[{g, SimulationDay{d}}] =
                        vacc_daily_more_booster / number_R_left;
                    day_campaign_idx++;
                }
                else {
                    model[i].parameters.template get<DailyPartialVaccination>()[{g, SimulationDay{d}}] =
                        model[i].parameters.template get<RateOfDailyPartialVaccinations>()[{g}];
                    model[i].parameters.template get<DailyFullVaccination>()[{g, SimulationDay{d}}] =
                        model[i].parameters.template get<RateOfDailyImprovedVaccinations>()[{g}];
                    model[i].parameters.template get<DailyBoosterVaccination>()[{g, SimulationDay{d}}] =
                        model[i].parameters.template get<RateOfDailyBoosterVaccinations>()[{g}];
                }
            }
        }
    }
    return success();
}
} // namespace details
} // namespace osecirvvs
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP
