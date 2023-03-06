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
#include "memilio/geography/regions.h"
#include "memilio/geography/holiday_data.ipp"

#include <tuple>

namespace mio
{

namespace regions
{

StateId get_state_id(int county)
{
    return StateId(county / 1000);
}

Range<std::pair<std::vector<std::pair<Date, Date>>::const_iterator, std::vector<std::pair<Date, Date>>::const_iterator>>
get_holidays(StateId state)
{
    static const std::vector<std::pair<mio::Date, mio::Date>> def;
    assert(int(state) >= 1 && int(state) <= 16 && "invalid state_id");

    auto data = &def;
    switch (int(state)) {
    case 1:
        data = &holidays::de::holidays_01_sh;
        break;
    case 2:
        data = &holidays::de::holidays_02_hh;
        break;
    case 3:
        data = &holidays::de::holidays_03_ni;
        break;
    case 4:
        data = &holidays::de::holidays_04_hb;
        break;
    case 5:
        data = &holidays::de::holidays_05_nw;
        break;
    case 6:
        data = &holidays::de::holidays_06_he;
        break;
    case 7:
        data = &holidays::de::holidays_07_rp;
        break;
    case 8:
        data = &holidays::de::holidays_08_bw;
        break;
    case 9:
        data = &holidays::de::holidays_09_by;
        break;
    case 10:
        data = &holidays::de::holidays_10_sl;
        break;
    case 11:
        data = &holidays::de::holidays_11_be;
        break;
    case 12:
        data = &holidays::de::holidays_12_bb;
        break;
    case 13:
        data = &holidays::de::holidays_13_mv;
        break;
    case 14:
        data = &holidays::de::holidays_14_sn;
        break;
    case 15:
        data = &holidays::de::holidays_15_st;
        break;
    case 16:
        data = &holidays::de::holidays_16_th;
        break;
    default:
        break;
    }
    return make_range(data->cbegin(), data->cend());
}

Range<std::pair<std::vector<std::pair<Date, Date>>::const_iterator, std::vector<std::pair<Date, Date>>::const_iterator>>
get_holidays(StateId state, Date start_date, Date end_date)
{
    auto all = get_holidays(state);

    //all holiday periods that overlap with the specified period
    auto holidays_in_range =
        std::equal_range(all.begin(), all.end(), std::make_pair(start_date, end_date), [](auto& p1, auto& p2) {
            return std::tie(p1.second.year, p1.second.month, p1.second.day) <
                   std::tie(p2.first.year, p2.first.month, p2.first.day);
        });

    return {holidays_in_range};
}

} // namespace regions

} // namespace mio
