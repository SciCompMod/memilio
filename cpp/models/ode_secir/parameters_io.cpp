/* 
* Copyright (C) 2020-2024 MEmilio
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

#include "memilio/utils/compiler_diagnostics.h"

//see below for line that causes this warning
GCC_CLANG_DIAGNOSTIC(push)
GCC_CLANG_DIAGNOSTIC(ignored "-Wmaybe-uninitialized")

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "ode_secir/parameters_io.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"
#include "memilio/utils/stl_util.h"
#include "memilio/utils/date.h"

namespace mio
{

namespace osecir
{

namespace details
{
//district, county or state id of a data entry if available, 0 (for whole country) otherwise
//used to compare data entries to integer ids in STL algorithms
template <class EpiDataEntry>
int get_region_id(const EpiDataEntry& entry)
{
    return entry.county_id
               ? entry.county_id->get()
               : (entry.state_id ? entry.state_id->get() : (entry.district_id ? entry.district_id->get() : 0));
}
//overload for integers, so the comparison of data entry to integers is symmetric (required by e.g. equal_range)
int get_region_id(int id)
{
    return id;
}

IOResult<void> read_confirmed_cases_data(
    std::vector<ConfirmedCasesDataEntry>& rki_data, std::vector<int> const& vregion, Date date,
    std::vector<std::vector<double>>& vnum_Exposed, std::vector<std::vector<double>>& vnum_InfectedNoSymptoms,
    std::vector<std::vector<double>>& vnum_InfectedSymptoms, std::vector<std::vector<double>>& vnum_InfectedSevere,
    std::vector<std::vector<double>>& vnum_icu, std::vector<std::vector<double>>& vnum_death,
    std::vector<std::vector<double>>& vnum_rec, const std::vector<std::vector<int>>& vt_Exposed,
    const std::vector<std::vector<int>>& vt_InfectedNoSymptoms,
    const std::vector<std::vector<int>>& vt_InfectedSymptoms, const std::vector<std::vector<int>>& vt_InfectedSevere,
    const std::vector<std::vector<int>>& vt_InfectedCritical, const std::vector<std::vector<double>>& vmu_C_R,
    const std::vector<std::vector<double>>& vmu_I_H, const std::vector<std::vector<double>>& vmu_H_U,
    const std::vector<double>& scaling_factor_inf)
{
    auto max_date_entry = std::max_element(rki_data.begin(), rki_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == rki_data.end()) {
        log_error("RKI data file is empty.");
        return failure(StatusCode::InvalidFileFormat, "RKI file is empty.");
    }
    auto max_date = max_date_entry->date;
    if (max_date < date) {
        log_error("Specified date does not exist in RKI data");
        return failure(StatusCode::OutOfRange, "Specified date does not exist in RKI data.");
    }
    auto days_surplus = std::min(get_offset_in_days(max_date, date) - 6, 0);

    //this statement causes maybe-uninitialized warning for some versions of gcc.
    //the error is reported in an included header, so the warning is disabled for the whole file
    std::sort(rki_data.begin(), rki_data.end(), [](auto&& a, auto&& b) {
        return std::make_tuple(get_region_id(a), a.date) < std::make_tuple(get_region_id(b), b.date);
    });

    for (auto region_idx = size_t(0); region_idx < vregion.size(); ++region_idx) {
        auto region_entry_range_it =
            std::equal_range(rki_data.begin(), rki_data.end(), vregion[region_idx], [](auto&& a, auto&& b) {
                return get_region_id(a) < get_region_id(b);
            });
        auto region_entry_range = make_range(region_entry_range_it);
        if (region_entry_range.begin() == region_entry_range.end()) {
            log_error("No entries found for region {}", vregion[region_idx]);
            return failure(StatusCode::InvalidFileFormat,
                           "No entries found for region " + std::to_string(vregion[region_idx]));
        }
        for (auto&& region_entry : region_entry_range) {

            auto& t_Exposed            = vt_Exposed[region_idx];
            auto& t_InfectedNoSymptoms = vt_InfectedNoSymptoms[region_idx];
            auto& t_InfectedSymptoms   = vt_InfectedSymptoms[region_idx];
            auto& t_InfectedSevere     = vt_InfectedSevere[region_idx];
            auto& t_InfectedCritical   = vt_InfectedCritical[region_idx];

            auto& num_InfectedNoSymptoms = vnum_InfectedNoSymptoms[region_idx];
            auto& num_InfectedSymptoms   = vnum_InfectedSymptoms[region_idx];
            auto& num_rec                = vnum_rec[region_idx];
            auto& num_Exposed            = vnum_Exposed[region_idx];
            auto& num_InfectedSevere     = vnum_InfectedSevere[region_idx];
            auto& num_death              = vnum_death[region_idx];
            auto& num_icu                = vnum_icu[region_idx];

            auto& mu_C_R = vmu_C_R[region_idx];
            auto& mu_I_H = vmu_I_H[region_idx];
            auto& mu_H_U = vmu_H_U[region_idx];

            auto date_df = region_entry.date;
            auto age     = size_t(region_entry.age_group);

            bool read_icu = false; //params.populations.get({age, SecirCompartments::U}) == 0;

            if (date_df == offset_date_by_days(date, 0)) {
                num_InfectedSymptoms[age] += scaling_factor_inf[age] * region_entry.num_confirmed;
                num_rec[age] += region_entry.num_confirmed;
            }
            if (date_df == offset_date_by_days(date, days_surplus)) {
                num_InfectedNoSymptoms[age] -=
                    1 / (1 - mu_C_R[age]) * scaling_factor_inf[age] * region_entry.num_confirmed;
            }
            if (date_df == offset_date_by_days(date, t_InfectedNoSymptoms[age] + days_surplus)) {
                num_InfectedNoSymptoms[age] +=
                    1 / (1 - mu_C_R[age]) * scaling_factor_inf[age] * region_entry.num_confirmed;
                num_Exposed[age] -= 1 / (1 - mu_C_R[age]) * scaling_factor_inf[age] * region_entry.num_confirmed;
            }
            if (date_df == offset_date_by_days(date, t_Exposed[age] + t_InfectedNoSymptoms[age] + days_surplus)) {
                num_Exposed[age] += 1 / (1 - mu_C_R[age]) * scaling_factor_inf[age] * region_entry.num_confirmed;
            }
            if (date_df == offset_date_by_days(date, -t_InfectedSymptoms[age])) {
                num_InfectedSymptoms[age] -= scaling_factor_inf[age] * region_entry.num_confirmed;
                num_InfectedSevere[age] += mu_I_H[age] * scaling_factor_inf[age] * region_entry.num_confirmed;
            }
            if (date_df == offset_date_by_days(date, -t_InfectedSymptoms[age] - t_InfectedSevere[age])) {
                num_InfectedSevere[age] -= mu_I_H[age] * scaling_factor_inf[age] * region_entry.num_confirmed;
                if (read_icu) {
                    num_icu[age] += mu_I_H[age] * mu_H_U[age] * scaling_factor_inf[age] * region_entry.num_confirmed;
                }
            }
            if (date_df ==
                offset_date_by_days(date, -t_InfectedSymptoms[age] - t_InfectedSevere[age] - t_InfectedCritical[age])) {
                num_death[age] += region_entry.num_deaths;
                if (read_icu) {
                    num_icu[age] -= mu_I_H[age] * mu_H_U[age] * scaling_factor_inf[age] * region_entry.num_confirmed;
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

        for (size_t i = 0; i < ConfirmedCasesDataEntry::age_group_names.size(); i++) {
            auto try_fix_constraints = [region, i](double& value, double error, auto str) {
                if (value < error) {
                    //this should probably return a failure
                    //but the algorithm is not robust enough to avoid large negative values and there are tests that rely on it
                    log_error("{:s} for age group {:s} is {:.4f} for region {:d}, exceeds expected negative value.",
                              str, ConfirmedCasesDataEntry::age_group_names[i], value, region);
                    value = 0.0;
                }
                else if (value < 0) {
                    log_info("{:s} for age group {:s} is {:.4f} for region {:d}, automatically corrected", str,
                             ConfirmedCasesDataEntry::age_group_names[i], value, region);
                    value = 0.0;
                }
            };

            try_fix_constraints(num_InfectedSymptoms[i], -5, "InfectedSymptoms");
            try_fix_constraints(num_InfectedNoSymptoms[i], -5, "InfectedNoSymptoms");
            try_fix_constraints(num_Exposed[i], -5, "Exposed");
            try_fix_constraints(num_InfectedSevere[i], -5, "InfectedSevere");
            try_fix_constraints(num_death[i], -5, "Dead");
            try_fix_constraints(num_icu[i], -5, "InfectedCritical");
            try_fix_constraints(num_rec[i], -20, "Recovered");
        }
    }

    return success();
}

IOResult<void> read_divi_data(const std::string& path, const std::vector<int>& vregion, Date date,
                              std::vector<double>& vnum_icu)
{
    BOOST_OUTCOME_TRY(auto&& divi_data, mio::read_divi_data(path));

    auto max_date_entry = std::max_element(divi_data.begin(), divi_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == divi_data.end()) {
        log_error("DIVI data file is empty.");
        return failure(StatusCode::InvalidFileFormat, path + ", file is empty.");
    }
    auto max_date = max_date_entry->date;
    if (max_date < date) {
        log_error("Specified date does not exist in DIVI data.");
        return failure(StatusCode::OutOfRange, path + ", specified date does not exist in DIVI data.");
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

IOResult<std::vector<std::vector<double>>>
read_population_data(const std::string& path, const std::vector<int>& vregion, bool accumulate_age_groups)
{
    BOOST_OUTCOME_TRY(auto&& population_data, mio::read_population_data(path, !accumulate_age_groups));
    //if we set up the model for one age group, the population data should be read in with the
    //age groups given in the population data json file and are accumulated later
    //otherwise the populations are directly saved for the correct model age groups
    size_t age_group_size = accumulate_age_groups ? PopulationDataEntry::age_group_names.size()
                                                  : ConfirmedCasesDataEntry::age_group_names.size();
    std::vector<std::vector<double>> vnum_population(vregion.size(), std::vector<double>(age_group_size, 0.0));

    for (auto&& entry : population_data) {
        auto it = std::find_if(vregion.begin(), vregion.end(), [&entry](auto r) {
            return r == 0 || (entry.county_id && regions::StateId(r) == regions::get_state_id(int(*entry.county_id))) ||
                   (entry.county_id && regions::CountyId(r) == *entry.county_id) ||
                   (entry.district_id && regions::DistrictId(r) == *entry.district_id);
        });
        if (it != vregion.end()) {
            auto region_idx      = size_t(it - vregion.begin());
            auto& num_population = vnum_population[region_idx];
            for (size_t age = 0; age < num_population.size(); age++) {
                num_population[age] += entry.population[AgeGroup(age)];
            }
        }
    }
    if (accumulate_age_groups) {
        std::vector<std::vector<double>> vnum_pop_acc(vregion.size(), std::vector<double>(1, 0));
        for (size_t region = 0; region < vregion.size(); ++region) {
            vnum_pop_acc[region][0] =
                std::accumulate(vnum_population[region].begin(), vnum_population[region].end(), 0.0);
        }
        return success(vnum_pop_acc);
    }
    else {
        return success(vnum_population);
    }
}

} // namespace details
} // namespace osecir
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

GCC_CLANG_DIAGNOSTIC(pop)
