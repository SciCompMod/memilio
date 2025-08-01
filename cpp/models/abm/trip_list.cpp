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

const Trip& TripList::get_next_trip() const
{
    return m_trips[m_current_index];
}

TimePoint TripList::get_next_trip_time() const
{
    return m_trips[m_current_index].trip_time;
}

void TripList::add_trips(std::vector<Trip> trip)
{
    //Trips are sorted by time.
    //Also include the person id in the comparison so different persons can make trips at the same time.
    //The same person can only make one trip at the same time.

    std::sort(trip.begin(), trip.end(), [](auto& trip1, auto& trip2) {
        return std::tie(trip1.trip_time, trip1.person_id) < std::tie(trip2.trip_time, trip2.person_id);
    });
    // Avoid storage duplication by using in-place merge
    const size_t original_size = m_trips.size();

    // Append the sorted new trips to the end
    m_trips.insert(m_trips.end(), std::make_move_iterator(trip.begin()), std::make_move_iterator(trip.end()));

    // Use in-place merge to merge the two sorted ranges
    std::inplace_merge(
        m_trips.begin(), m_trips.begin() + original_size, m_trips.end(), [](const auto& trip1, const auto& trip2) {
            return std::tie(trip1.trip_time, trip1.person_id) < std::tie(trip2.trip_time, trip2.person_id);
        });
}

} // namespace abm
} // namespace mio
