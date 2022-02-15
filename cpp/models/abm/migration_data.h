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
#ifndef EPI_ABM_MIGRATION_DATA_H
#define EPI_ABM_MIGRATION_DATA_H

#include "abm/parameters.h"
#include "abm/location.h"
#include "abm/state.h"
#include "abm/location_type.h"

#include "memilio/math/eigen.h"
#include <array>
#include <random>

namespace mio
{

/**
 * A trip describes a migration from one location to another location.
 * @param person_index describes the person that makes the trip and corresponds to the index into the structure m_persons from world, where all people are saved.
 * @param migration_time is the time at which a person changes the location.
 * @param migration_target is the location where the person migrates to.
 * @param migration_start is the location where te person starts the trip.
 * @param cells If migration_target consists of different cells, this gives the index of the cells the person migrates to.
 */
struct Trip {
    uint32_t person_id;
    TimePoint migration_time;
    LocationId migration_target;
    LocationId migration_start;
    std::vector<uint32_t> cells;

    Trip(uint32_t id, TimePoint time, LocationId target, LocationId start,
         const std::vector<uint32_t>& input_cells = {})
    {
        person_id        = id;
        migration_time   = time;
        migration_target = target;
        migration_start  = start;
        cells            = input_cells;
    }
};

class MigrationData
{
public:
    /**
     * construct Migration data
     */
    MigrationData();

    /*
     * returns the next trip
     */
    Trip& get_next_trip();

    /*
     * returns the time at which the next trip will happen
     */
    TimePoint get_next_trip_time();

    /**
     * sort the trips by migration time
     */
    void sort_trips();

    /**
     * add a trip to migration data
     */
    void add_trip(Trip trip);

    /*
     * increases m_current_index by one
     */
    void increase_index()
    {
        m_current_index++;
    }

    /* 
     * returns vector of trips
     */
    std::vector<Trip>& get_trips()
    {
        return m_trips;
    }

    /* 
     * returns currents index
     */
    uint32_t get_current_index()
    {
        return m_current_index;
    }

private:
    std::vector<Trip> m_trips;
    uint32_t m_current_index;
};
} // namespace mio

#endif
