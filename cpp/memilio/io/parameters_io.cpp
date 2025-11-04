/* 
* Copyright (C) 2020-2025 MEmilio
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

#include "memilio/config.h"
#include "memilio/io/parameters_io.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "memilio/io/epi_data.h"
#include "memilio/io/result_io.h"
#include "json/value.h"
#include <string>
#include <vector>

namespace mio
{
IOResult<std::vector<ScalarType>> compute_divi_data(const std::vector<DiviEntry>& divi_data, const std::vector<int>& vregion, Date date)
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

    std::vector<ScalarType> vnum_icu(vregion.size(), 0.0);


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

    return success(vnum_icu);
}

IOResult<std::vector<ScalarType>> read_divi_data(const std::string& path, const std::vector<int>& vregion, Date date)
{
    BOOST_OUTCOME_TRY(auto&& divi_data, mio::read_divi_data(path));
    return compute_divi_data(divi_data, vregion, date);
}

IOResult<std::vector<std::vector<ScalarType>>>
read_population_data(const std::vector<PopulationDataEntry>& population_data, const std::vector<int>& vregion)
{
    std::vector<std::vector<ScalarType>> vnum_population(
        vregion.size(), std::vector<ScalarType>(ConfirmedCasesDataEntry::age_group_names.size(), 0.0));

    for (auto&& county_entry : population_data) {
        //accumulate population of states or country from population of counties
        if (!county_entry.county_id && !county_entry.district_id) {
            return failure(StatusCode::InvalidFileFormat, "File with county population expected.");
        }
        //find region that this county belongs to
        //all counties belong to the country (id = 0)
        auto it = std::find_if(vregion.begin(), vregion.end(), [&county_entry](auto r) {
            return r == 0 ||
                   (county_entry.county_id &&
                    regions::de::StateId(r) == regions::de::get_state_id(int(*county_entry.county_id))) ||
                   (county_entry.county_id && regions::de::CountyId(r) == *county_entry.county_id) ||
                   (county_entry.district_id && regions::de::DistrictId(r) == *county_entry.district_id);
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

IOResult<std::vector<std::vector<ScalarType>>> read_population_data(const std::string& path,
                                                                    const std::vector<int>& vregion)
{
    BOOST_OUTCOME_TRY(auto&& population_data, mio::read_population_data(path));
    return read_population_data(population_data, vregion);
}
} // namespace mio
#endif //MEMILIO_HAS_JSONCPP
