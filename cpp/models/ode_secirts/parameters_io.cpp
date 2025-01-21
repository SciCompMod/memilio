/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker, Wadim Koslow, Daniel Abele, Martin J. Kühn
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

#include "ode_secirts/parameters_io.h"
#include "memilio/geography/regions.h"
#include "memilio/io/io.h"
#include "ode_secirts/parameters.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "memilio/io/epi_data.h"
#include "memilio/utils/memory.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/utils/stl_util.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/date.h"

#include <boost/filesystem.hpp>

#include <numeric>
#include <vector>
#include <iostream>
#include <string>
#include <random>
#include <fstream>

namespace mio
{
namespace osecirts
{
namespace details
{

IOResult<std::vector<std::vector<double>>> read_population_data(const std::string& path,
                                                                const std::vector<int>& vregion)
{
    BOOST_OUTCOME_TRY(auto&& population_data, mio::read_population_data(path));
    return read_population_data(population_data, vregion);
}

IOResult<std::vector<std::vector<double>>> read_population_data(const std::vector<PopulationDataEntry>& population_data,
                                                                const std::vector<int>& vregion)
{
    std::vector<std::vector<double>> vnum_population(
        vregion.size(), std::vector<double>(ConfirmedCasesDataEntry::age_group_names.size(), 0.0));

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
                    regions::StateId(r) == regions::get_state_id(int(*county_entry.county_id))) ||
                   (county_entry.county_id && regions::CountyId(r) == *county_entry.county_id) ||
                   (county_entry.district_id && regions::DistrictId(r) == *county_entry.district_id);
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

} // namespace details
} // namespace osecirts
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP
