/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Henrik Zunker
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
#ifndef MEMILIO_IO_PARAMETER_H
#define MEMILIO_IO_PARAMETER_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "memilio/io/epi_data.h"
#include "memilio/io/result_io.h"
#include "json/value.h"
#include <string>
#include <vector>

namespace mio
{
/**
 * @brief Gets the region ID (county, state, or district) of an EpiDataEntry.
 * 
 * If none are available, it defaults to 0 which is representing the whole country.
 * 
 * @tparam EpiDataEntry The type of the data entry.
 * @param data_entry The (RKI) data entry to extract the region ID from.
 * @return The region ID as integer, or 0 if no specific region information is available.
 */
template <class EpiDataEntry>
int get_region_id(const EpiDataEntry& data_entry)
{
    return data_entry.county_id ? data_entry.county_id->get()
                                : (data_entry.state_id ? data_entry.state_id->get()
                                                       : (data_entry.district_id ? data_entry.district_id->get() : 0));
}

/**
 * @brief Extracts the number of individuals in critical condition (ICU) for each region 
 * on a specified date from the provided DIVI data.
 *
 * @tparam FP Floating point type (default: double).
 *
 * @param[in] divi_data Vector of DIVI data entries containing date, region, and ICU information.
 * @param[in] vregion Vector of region IDs for which the data is computed.
 * @param[in] date Date for which the ICU data is computed.
 * @param[out] vnum_icu Output vector containing the number of ICU cases for each region.
 *
 * @return An IOResult indicating success or failure.
 */
template <typename FP = double>
IOResult<void> compute_divi_data(const std::vector<DiviEntry>& divi_data, const std::vector<int>& vregion, Date date,
                                 std::vector<FP>& vnum_icu)
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

/**
 * @brief Reads DIVI data from a file and computes the ICU data for specified regions and date.
 *
 * @tparam FP Floating point type (default: double).
 *
 * @param[in] path Path to the file containing DIVI data.
 * @param[in] vregion Vector of region IDs for which the data is computed.
 * @param[in] date Date for which the ICU data is computed.
 * @param[out] vnum_icu Output vector containing the number of ICU cases for each region.
 *
 * @return An IOResult indicating success or failure.
 */
template <typename FP = double>
IOResult<void> read_divi_data(const std::string& path, const std::vector<int>& vregion, Date date,
                              std::vector<FP>& vnum_icu)
{
    BOOST_OUTCOME_TRY(auto&& divi_data, mio::read_divi_data(path));
    return compute_divi_data(divi_data, vregion, date, vnum_icu);
}

/**
 * @brief Reads population data from a vector of population data entries.
 * 
 * @param[in] population_data Vector of population data entries.
 * @param[in] vregion Vector of keys representing the regions of interest.
 * @return An IOResult containing a vector of vectors, where each inner vector represents the population
 *         distribution across age groups for a specific region, or an error if the function fails.
 */
IOResult<std::vector<std::vector<double>>> read_population_data(const std::vector<PopulationDataEntry>&,
                                                                const std::vector<int>&);

/**
 * @brief Reads population data from census data.
 * 
 * @param[in] path Path to the population data file.
 * @param[in] vregion Vector of keys representing the regions of interest.
 * @return An IOResult containing a vector of vectors, where each inner vector represents the population
 *         distribution across age groups for a specific region, or an error if the function fails.
 */
IOResult<std::vector<std::vector<double>>> read_population_data(const std::string&, const std::vector<int>&);

} // namespace mio

#endif //MEMILIO_HAS_JSONCPP

#endif //MEMILIO_IO_PARAMETER_H
