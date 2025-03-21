/* 
* Copyright (C) 2020-2025 MEmilio
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
#ifndef MIO_ABM_TRIP_LIST_H
#define MIO_ABM_TRIP_LIST_H

#include "abm/location_id.h"
#include "abm/mobility_data.h"
#include "abm/person.h"
#include "abm/person_id.h"
#include "abm/time.h"
#include "abm/location_type.h"
#include "memilio/io/io.h"
#include "memilio/io/default_serialize.h"
#include <cstdint>
#include <vector>

namespace mio
{
namespace abm
{

/**
 * @brief A trip describes a change of Location from one Location to another Location.
 */
struct Trip {
    //TODO: Origin is currently not used for the trips. Should we delete it then?
    PersonId person_id; /**< Person that makes the trip and corresponds to the index into the structure m_persons from
    Model, where all Person%s are saved.*/
    TimePoint time; ///< Daytime at which a Person changes the Location.
    LocationId destination; ///< Location where the Person changes to.
    int destination_model_id; ///< Model id of destination Location.
    LocationId origin; ///< Location where the Person starts the Trip.
    int origin_model_id; ///< Model id of origin Location.
    std::vector<uint32_t> cells; /**< If destination consists of different Cell%s, this gives the index of the
    Cell%s the Person changes to.*/
    TransportMode
        trip_mode; ///< Mode of transportation. 1:Bike, 2:Car (Driver), 3:Car (Co-Driver)), 4:Public Transport, 5:Walking, 6:Other/Unknown
    LocationType destination_type; ///< Type of destination Location.

    /**
     * @brief Construct a new Trip.
     * @param[in] id ID of the Person that makes the Trip.
     * @param[in] time_new Time at which a Person changes the Location this currently cant be set for s specific day just a timepoint in a day.
     * @param[in] destination Location where the Person changes to.
     * @param[in] destination_model_id Model the Person changes to.
     * @param[in] origin Location where the person starts the Trip.
     * @param[in] origin_model_id Model the Person starts the Trip.
     * @param[in] input_cells The index of the Cell%s the Person changes to.
     */
    Trip(PersonId id, TimePoint time_new, LocationId dest, int dest_model_id, LocationId orig, int orig_model_id,
         TransportMode mode_of_transport, LocationType type_of_activity, const std::vector<uint32_t>& input_cells = {})
        : person_id(id)
        , time(mio::abm::TimePoint(time_new.time_since_midnight().seconds()))
        , destination(dest)
        , destination_model_id(dest_model_id)
        , origin(orig)
        , origin_model_id(orig_model_id)
        , cells(input_cells)
        , trip_mode(mode_of_transport)
        , destination_type(type_of_activity)
    {
    }

    Trip(PersonId id, TimePoint time_new, LocationId dest, LocationId orig, TransportMode mode_of_transport,
         LocationType type_of_activity, const std::vector<uint32_t>& input_cells = {})
        : person_id(id)
        , time(mio::abm::TimePoint(time_new.time_since_midnight().seconds()))
        , destination(dest)
        , destination_model_id(0)
        , origin(orig)
        , origin_model_id(0)
        , cells(input_cells)
        , trip_mode(mode_of_transport)
        , destination_type(type_of_activity)
    {
    }

    Trip(PersonId id, TimePoint time_new, LocationId dest, LocationId orig, LocationType type_of_activity,
         const std::vector<uint32_t>& input_cells = {})
        : Trip(id, time_new, dest, orig, mio::abm::TransportMode::Unknown, type_of_activity, input_cells)
    {
    }

    Trip(PersonId id, TimePoint time_new, LocationId dest, LocationType type_of_activity,
         const std::vector<uint32_t>& input_cells = {})
        : Trip(id, time_new, dest, dest, mio::abm::TransportMode::Unknown, type_of_activity, input_cells)
    {
    }

    /**
     * @brief Compare two Trip%s.
     */
    bool operator==(const Trip& other) const
    {
        return (person_id == other.person_id) && (time == other.time) && (destination == other.destination) &&
               (origin == other.origin);
    }

    auto default_serialize()
    {
        return Members("Trip")
            .add("person_id", person_id)
            .add("time", time)
            .add("destination", destination)
            .add("origin", origin);
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
     * @brief Add a Trip to the trip list.
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

    /// This method is used by the default serialization feature.
    auto default_serialize()
    {
        return Members("TestingScheme")
            .add("trips_weekday", m_trips_weekday)
            .add("trips_weekend", m_trips_weekend)
            .add("index", m_current_index);
    }

private:
    std::vector<Trip> m_trips_weekday; ///< The list of Trip%s a Person makes on a weekday.
    std::vector<Trip> m_trips_weekend; ///< The list of Trip%s a Person makes on a weekend day.
    uint32_t m_current_index; ///< The index of the Trip a Person makes next.
};

} // namespace abm

/// @brief Creates an instance of abm::Trip for default serialization.
template <>
struct DefaultFactory<abm::Trip> {
    static abm::Trip create()
    {
        return abm::Trip{abm::PersonId{}, abm::TimePoint{}, abm::LocationId{}, abm::LocationType{}};
    }
};

} // namespace mio

#endif
