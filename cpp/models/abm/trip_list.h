/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Elisabeth Kluth, Khoa Nguyen
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
#ifndef EPI_ABM_TRIP_LIST_H
#define EPI_ABM_TRIP_LIST_H

#include "abm/parameters.h"
#include "abm/location.h"
#include "abm/infection.h"
#include "abm/location_type.h"

#include "memilio/math/eigen.h"
#include <array>
#include <random>

namespace mio
{
namespace abm
{

/**
 * @brief A trip describes a migration from one Location to another Location.
 */
struct Trip {
    uint32_t person_id; /**< Person that makes the trip and corresponds to the index into the structure m_persons from
    World, where all Person%s are saved.*/
    TimePoint time; ///< Time at which a Person changes the Location.
    LocationId migration_destination; ///< Location where the Person migrates to.
    LocationId migration_origin; ///< Location where the Person starts the Trip.
    std::vector<uint32_t> cells; /**< If migration_destination consists of different Cell%s, this gives the index of the
    Cell%s the Person migrates to.*/
    TransportMode
        trip_mode; ///< Mode of transportation. 1:Bike, 2:Car (Driver), 3:Car (Co-Driver)), 4:Public Transport, 5:Walking, 6:Other/Unknown
    ActivityType
        activity_type; ///< Type of activity. 1:Workplace, 2:Education, 3:Shopping, 4:Leisure, 5:Private Matters, 6:Other Activity, 7:Home, 8:Unknown Activity

    /**
     * @brief Construct a new Trip.
     * @param[in] id ID of the Person that makes the Trip.
     * @param[in] time_new Time at which a Person changes the Location this currently cant be set for s specific day just a timepoint in a day.
     * @param[in] destination Location where the Person migrates to.
     * @param[in] origin Location where the person starts the Trip.
     * @param[in] input_cells The index of the Cell%s the Person migrates to.
     */
    Trip(uint32_t id, TimePoint time_new, LocationId destination, LocationId origin, TransportMode mode_of_transport,
         ActivityType type_of_activity, const std::vector<uint32_t>& input_cells = {})
        : person_id(id)
        , time(mio::abm::TimePoint(time_new.time_since_midnight().seconds()))
        , migration_destination(destination)
        , migration_origin(origin)
        , cells(input_cells)
        , trip_mode(mode_of_transport)
        , activity_type(type_of_activity)
    {
    }

    Trip(uint32_t id, TimePoint time_new, LocationId destination, const std::vector<uint32_t>& input_cells = {})
        : Trip(id, time_new, destination, destination, mio::abm::TransportMode::Unknown,
               mio::abm::ActivityType::UnknownActivity, input_cells)
    {
    }

    Trip(uint32_t id, TimePoint time_new, LocationId destination, LocationId origin,
         const std::vector<uint32_t>& input_cells = {})
        : Trip(id, time_new, destination, origin, mio::abm::TransportMode::Unknown,
               mio::abm::ActivityType::UnknownActivity, input_cells)
    {
    }

    /**
     * @brief Compare two Trip%s.
     */
    bool operator==(const Trip& other) const
    {
        return (person_id == other.person_id) && (time == other.time) &&
               (migration_destination == other.migration_destination) && (migration_origin == other.migration_origin);
    }

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Trip");
        obj.add_element("person_id", person_id);
        obj.add_element("time", time.seconds());
        obj.add_element("destination_index", migration_destination.index);
        obj.add_element("destination_type", migration_destination.type);
        obj.add_element("origin_index", migration_origin.index);
        obj.add_element("origin_type", migration_origin.type);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Trip> deserialize(IOContext& io)
    {
        auto obj               = io.expect_object("Trip");
        auto person_id         = obj.expect_element("person_id", Tag<uint32_t>{});
        auto time              = obj.expect_element("time", Tag<int>{});
        auto destination_index = obj.expect_element("destination_index", Tag<uint32_t>{});
        auto destination_type  = obj.expect_element("destination_type", Tag<uint32_t>{});
        auto origin_index      = obj.expect_element("origin_index", Tag<uint32_t>{});
        auto origin_type       = obj.expect_element("origin_type", Tag<uint32_t>{});
        return apply(
            io,
            [](auto&& person_id_, auto&& time_, auto&& destination_index_, auto&& destination_type_,
               auto&& origin_index_, auto&& origin_type_) {
                return Trip(person_id_, TimePoint(time_),
                            LocationId{destination_index_, LocationType(destination_type_)},
                            LocationId{origin_index_, LocationType(origin_type_)});
            },
            person_id, time, destination_index, destination_type, origin_index, origin_type);
    }
};

/**
 * @brief A list of Trip%s a Person follows.
 */
class TripList
{
public:
    /**
     * @brief Construct empty TripList.
     */
    TripList();

    /**
     * @brief Get the next Trip.
     * @param weekend Whether the Trip%s during the week or on the weekend are used.
     */
    const Trip& get_next_trip(bool weekend) const;

    /**
     * @brief Get the time at which the next Trip will happen.
     * @param weekend Whether the Trip%s during the week or on the weekend are used.
     */
    TimePoint get_next_trip_time(bool weekend) const;

    /**
     * @brief Add a Trip to migration data.
     * @param[in] trip The Trip to be added.
     * @param[in] weekend If the Trip is made on a weekend day.     
     */
    void add_trip(Trip trip, bool weekend = false);

    /**
     * @brief Use the same TripList for weekend and weekday.
     */
    void use_weekday_trips_on_weekend();

    /**
     * @brief Increment the current index to select the next Trip.
     */
    void increase_index()
    {
        m_current_index++;
    }

    /**
     * @brief Reset the current index to 0.
     */
    void reset_index()
    {
        m_current_index = 0;
    }

    /**
     * @brief Get the length of the TripList.
     * @param weekend Whether the Trip%s during the week or on the weekend are used.
     */
    size_t num_trips(bool weekend = false) const
    {
        return weekend ? m_trips_weekend.size() : m_trips_weekday.size();
    }

    /**
     * @brief Get the current index.
     */
    uint32_t get_current_index() const
    {
        return m_current_index;
    }

private:
    std::vector<Trip> m_trips_weekday; ///< The list of Trip%s a Person makes on a weekday.
    std::vector<Trip> m_trips_weekend; ///< The list of Trip%s a Person makes on a weekend day.
    uint32_t m_current_index; ///< The index of the Trip a Person makes next.
};

} // namespace abm
} // namespace mio

#endif
