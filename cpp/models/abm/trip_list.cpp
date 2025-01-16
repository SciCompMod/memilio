/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Elisabeth Kluth, Daniel Abele
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
#include "abm/trip_list.h"
#include "abm/random_events.h"
#include "memilio/utils/stl_util.h"

namespace mio
{
namespace abm
{

TripList::TripList()
    : m_trips_weekday({})
    , m_trips_weekend({})
    , m_current_index(0)
{
}

const Trip& TripList::get_next_trip(bool weekend) const
{
    return weekend ? m_trips_weekend[m_current_index] : m_trips_weekday[m_current_index];
}

TimePoint TripList::get_next_trip_time(bool weekend) const
{
    return weekend ? m_trips_weekend[m_current_index].time : m_trips_weekday[m_current_index].time;
}

void TripList::use_weekday_trips_on_weekend()
{
    m_trips_weekend = m_trips_weekday;
}

void TripList::add_trip(Trip trip, bool weekend)
{
    //Trips are sorted by time.
    //Also include the person id in the comparison so different persons can make trips at the same time.
    //The same person can only make one trip at the same time.
    if (!weekend) {
        insert_sorted_replace(m_trips_weekday, trip, [](auto& trip1, auto& trip2) {
            return std::tie(trip1.time, trip1.person_id) < std::tie(trip2.time, trip2.person_id);
        });
    }
    else {
        insert_sorted_replace(m_trips_weekend, trip, [](auto& trip1, auto& trip2) {
            return std::tie(trip1.time, trip1.person_id) < std::tie(trip2.time, trip2.person_id);
        });
    }
}

} // namespace abm
} // namespace mio
