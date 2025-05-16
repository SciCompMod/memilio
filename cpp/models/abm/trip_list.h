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
    PersonId person_id; /**< Person that makes the trip and corresponds to the index into the structure m_persons from
    Model, where all Person%s are saved.*/
    TimePoint time; ///< Daytime at which a Person changes the Location.
    LocationId destination; ///< Location where the Person changes to.
    int destination_model_id; ///< Model id of destination Location.
    TransportMode trip_mode; ///< Mode of transportation. See TransportMode for all possible modes of transportation.
    std::vector<uint32_t> cells; /**< If destination consists of different Cell%s, this gives the index of the
    Cell%s the Person changes to.*/

    /**
     * @brief Construct a new Trip.
     * @param[in] id ID of the Person that makes the Trip.
     * @param[in] time Time at which a Person changes the Location this currently cant be set for s specific day just a timepoint in a day.
     * @param[in] destination Location where the Person changes to.
     * @param[in] destination_model_id Model the Person changes to.
     * @param[in] origin Location where the person starts the Trip.
     * @param[in] origin_model_id Model the Person starts the Trip.
     * @param[in] input_cells The index of the Cell%s the Person changes to.
     */
    Trip(PersonId id, const TimePoint time, const LocationId dest, const int dest_model_id = 0,
         const TransportMode mode_of_transport    = mio::abm::TransportMode::Unknown,
         const std::vector<uint32_t>& input_cells = {})
        : person_id(id)
        , time(mio::abm::TimePoint(time.time_since_midnight().seconds()))
        , destination(dest)
        , destination_model_id(dest_model_id)
        , trip_mode(mode_of_transport)
        , cells(input_cells)

    {
    }

    /**
     * @brief Compare two Trip%s.
     */
    bool operator==(const Trip& other) const
    {
        return (person_id == other.person_id) && (time == other.time) && (destination == other.destination) &&
               (destination_model_id == other.destination_model_id) && (trip_mode == other.trip_mode);
    }

    auto default_serialize()
    {
        return Members("Trip")
            .add("person_id", person_id)
            .add("time", time)
            .add("destination", destination)
            .add("model_id", destination_model_id)
            .add("trip_mode", trip_mode);
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
    TripList() = default;

    /**
     * @brief Get the next Trip.
     * @param weekend Whether the Trip%s during the week or on the weekend are used.
     */
    const Trip& get_next_trip() const;

    /**
     * @brief Get the time at which the next Trip will happen.
     * @param weekend Whether the Trip%s during the week or on the weekend are used.
     */
    TimePoint get_next_trip_time() const;

    /**
     * @brief Adds Trips to the trip list.
     * @param[in] trips The Trips to be added.
     * @param[in] weekend If the Trip is made on a weekend day.     
     */
    void add_trips(std::vector<Trip> trip);

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
    size_t num_trips() const
    {
        return m_trips.size();
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
        return Members("TestingScheme").add("trips", m_trips).add("index", m_current_index);
    }

private:
    std::vector<Trip> m_trips; ///< The list of Trip%s a Person makes on a weekday.
    uint32_t m_current_index = 0; ///< The index of the Trip a Person makes next.
};

} // namespace abm

/// @brief Creates an instance of abm::Trip for default serialization.
template <>
struct DefaultFactory<abm::Trip> {
    static abm::Trip create()
    {
        return abm::Trip{abm::PersonId{}, abm::TimePoint{}, abm::LocationId{}};
    }
};

} // namespace mio

#endif
