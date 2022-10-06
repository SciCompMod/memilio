/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Elisabeth Kluth
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
#include "abm/state.h"
#include "abm/location_type.h"

#include "memilio/math/eigen.h"
#include <array>
#include <random>

namespace mio
{
namespace abm
{

/**
 * A trip describes a migration from one location to another location.
 */
struct Trip {
    /** person that makes the trip and corresponds to the index into the structure m_persons from world, where all people are saved*/
    uint32_t person_id;
    /** time at which a person changes the location*/
    TimePoint time;
    /**location where the person migrates to */
    LocationId migration_destination;
    /**location where te person starts the trip*/
    LocationId migration_origin;
    /**If migration_destination consists of different cells, this gives the index of the cells the person migrates to.*/
    std::vector<uint32_t> cells;

    Trip(uint32_t id, TimePoint time_new, LocationId destination, LocationId origin,
         const std::vector<uint32_t>& input_cells = {})
    {
        person_id             = id;
        time                  = time_new;
        migration_destination = destination;
        migration_origin      = origin;
        cells                 = input_cells;
    }
};

class TripList
{
public:
    /**
     * Construct empty trip list.
     */
    TripList();

    /*
     * Get the next trip.
     */
    const Trip& get_next_trip() const;

    /*
     * Get the time at which the next trip will happen.
     */
    TimePoint get_next_trip_time() const;

    /**
     * Add a trip to migration data.
     */
    void add_trip(Trip trip);

    /*
     * Increment the current index to select the next trip.
     */
    void increase_index()
    {
        m_current_index++;
    }

    /* 
     * Get the length of the TripList.
     */
    size_t num_trips() const
    {
        return m_trips.size();
    }

    /* 
     * Get the current index.
     */
    uint32_t get_current_index() const
    {
        return m_current_index;
    }

private:
    std::vector<Trip> m_trips;
    uint32_t m_current_index;
};

} // namespace abm
} // namespace mio

#endif
