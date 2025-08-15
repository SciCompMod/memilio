/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele
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

#include "memilio/io/epi_data.h"

#ifdef MEMILIO_HAS_JSONCPP

namespace mio
{

std::vector<const char*> ConfirmedCasesDataEntry::age_group_names  = {"Population"};
std::vector<const char*> PopulationDataEntry::age_group_names      = {"Population"};
// std::vector<const char*> PopulationDataEntrySpain::age_group_names = {"Population"};

std::vector<const char*> VaccinationDataEntry::age_group_names = {"0-4", "5-14", "15-34", "35-59", "60-79", "80-99"};

IOResult<std::vector<int>> get_node_ids(const std::string& path, bool is_node_for_county, bool rki_age_groups)
{
    BOOST_OUTCOME_TRY(auto&& population_data, read_population_data(path, rki_age_groups));
    std::vector<int> id;
    id.reserve(population_data.size());
    for (auto&& entry : population_data) {
        if (is_node_for_county) {
            if (entry.county_id) {
                id.push_back(entry.county_id->get());
            }
            else {
                return failure(StatusCode::InvalidValue, "Population data file is missing county ids.");
            }
        }
        else {
            if (entry.state_id) {
                id.push_back(entry.state_id->get());
            }
            else if (entry.district_id) {
                id.push_back(entry.district_id->get());
            }
            else {
                return failure(StatusCode::InvalidValue, "Population data file is missing district and state ids.");
            }
        }
    }

    //remove duplicate node ids
    id.erase(std::unique(id.begin(), id.end()), id.end());
    return success(id);
}

IOResult<std::vector<int>> get_country_id(const std::string& /*path*/, bool /*is_node_for_county*/,
                                          bool /*rki_age_groups*/)
{
    std::vector<int> id = {0};
    return success(id);
}
} // namespace mio

#endif //MEMILIO_HAS_JSONCPP
