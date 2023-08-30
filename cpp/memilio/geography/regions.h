/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#ifndef MIO_EPI_REGIONS_H
#define MIO_EPI_REGIONS_H

#include "memilio/utils/date.h"
#include "memilio/utils/stl_util.h"
#include "memilio/utils/type_safe.h"

#include "boost/filesystem.hpp"

//is used to provide some paths as function arguments
namespace fs = boost::filesystem;

namespace mio
{
/**
 * Contains utilities that depend on geographical regions.
 */
namespace regions
{

/**
         * Id of a state.
         * For Germany the Ids are:
         * 1 = Schleswig-Holstein
         * 2 = Hamburg
         * 3 = Niedersachsen
         * 4 = Bremen
         * 5 = Nordrhein-Westfalen
         * 6 = Hessen
         * 7 = Rheinland-Pfalz
         * 8 = Baden-Württemberg
         * 9 = Bayern
         * 10 = Saarland
         * 11 = Berlin
         * 12 = Brandenburg
         * 13 = Mecklenburg-Vorpommern
         * 14 = Sachsen
         * 15 = Sachsen-Anhalt
         * 16 = Thüringen
         */
DECL_TYPESAFE(int, StateId);

/**
         * Id of a county.
         * Format ssxxx where ss is the id of the state that the county is in (first s may be 0) and xxx are other digits.
         * Ids are generally not consecutive, even within one state.
         */
DECL_TYPESAFE(int, CountyId);

DECL_TYPESAFE(int, DistrictId);

/**
         * get the id of the state that the specified county is in. 
         * @param[in, out] county a county id.
         */
StateId get_state_id(int county);

/**
         * get the holidays in a german state.
         * @param[in] state id of the state.
         * @return range of pairs of start and end dates of holiday periods, sorted by start date.
         */
Range<std::pair<std::vector<std::pair<Date, Date>>::const_iterator, std::vector<std::pair<Date, Date>>::const_iterator>>
get_holidays(StateId state);

/**
         * get the holidays in a german state in a given time period.
         * The returned periods may not be completely included in the queried period,
         * they may only partially overlap.
         * @param[in] state id of the state.
         * @param[in] start_date start of the queried period.
         * @param[in] end_date end of the queried period.
         * @return range of pairs of start and end dates of holiday periods, sorted by start date.
         */
Range<std::pair<std::vector<std::pair<Date, Date>>::const_iterator, std::vector<std::pair<Date, Date>>::const_iterator>>
get_holidays(StateId state, Date start_date, Date end_date);

} // namespace regions
} // namespace mio

#endif //MIO_EPI_REGIONS_H
