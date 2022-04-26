/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#include "abm/location.h"
#include "abm/random_events.h"

#include <numeric>

namespace mio
{

TripList::TripList()
    : m_trips({})
    , m_current_index(0)
{
}

const Trip& TripList::get_next_trip() const
{
    return m_trips[m_current_index];
}

TimePoint TripList::get_next_trip_time() const
{
    return m_trips[m_current_index].time;
}

void TripList::add_trip(Trip trip)
{
    //trips are sorted by time
    //also include the person id in the comparison so different persons can make trips at the same time
    //same person can only make one trip at the same time
    insert_sorted_replace(m_trips, trip, [](auto& trip1, auto& trip2) {
        return std::tie(trip1.time, trip1.person_id) < std::tie(trip2.time, trip2.person_id);
    });
}

} // namespace mio
