/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Elisabeth Kluth, Khoa Nguyen
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
#ifndef EPI_ABM_LOCATION_H
#define EPI_ABM_LOCATION_H

#include "abm/parameters.h"
#include "abm/state.h"
#include "abm/location_type.h"

#include "memilio/math/eigen.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/time_series.h"
#include <array>
#include <random>

namespace mio
{
namespace abm
{
class Person;

/**
 * LocationCapacity describes the size of a location. 
 * It consists of a volume and a capacity in persons which is an upper bound for the number
 * of people that can be at the location at the same time.
 */
struct LocationCapacity {
    LocationCapacity()
        : volume(0)
        , persons(std::numeric_limits<int>::max())
    {
    }
    int volume;
    int persons;
};

/**
 * LocationId identifies a Location uniquely. It consists of the LocationType of the Location and an Index.
 * The index corresponds to the index into the structure m_locations from world, where all Locations are saved.
 */
struct LocationId {
    uint32_t index;
    LocationType type;

    bool operator==(const LocationId& rhs) const
    {
        return (index == rhs.index && type == rhs.type);
    }

    bool operator!=(const LocationId& rhs) const
    {
        return !(index == rhs.index && type == rhs.type);
    }
};

/**
 * The location can be split up into several cells. This allows a finer division of the people in public transport.
 */
struct Cell {
    uint32_t num_people;
    uint32_t num_carriers;
    uint32_t num_infected;
    CustomIndexArray<double, AgeGroup, VaccinationState> cached_exposure_rate;

    Cell()
        : num_people(0)
        , num_carriers(0)
        , num_infected(0)
        , cached_exposure_rate({{AgeGroup::Count, VaccinationState::Count}, 0.})
    {
    }

    Cell(uint32_t num_p, uint32_t num_c, uint32_t num_i,
         CustomIndexArray<double, AgeGroup, VaccinationState> cached_exposure_rate_new)
        : num_people(num_p)
        , num_carriers(num_c)
        , num_infected(num_i)
        , cached_exposure_rate(cached_exposure_rate_new)
    {
    }

}; // namespace mio

/**
 * all locations in the simulated world where persons gather.
 */
class Location
{
public:
    /**
     * construct a Location of a certain type.
     * @param type the type of the location
     * @param index the index of the location
     * @param num_cells the number of cells in which the location is divided
     */
    Location(LocationType type, uint32_t index, uint32_t num_cells = 0);

    /**
     * get the type of this location.
     */
    LocationType get_type() const
    {
        return m_type;
    }

    /**
     *get the index of this location.
     */
    unsigned get_index() const
    {
        return m_index;
    }

    /** 
     * a person interacts with the population at this location, may change infection state.
     * @param person the person that interacts with the population
     * @param dt length of the current simulation time step
     * @param global_params global infection parameters
     * @return new infection state of the person
     */
    InfectionState interact(const Person& person, TimeSpan dt, const GlobalInfectionParameters& global_params) const;

    /** 
     * add a person to the population at this location.
     * @param person the person arriving
    */
    void add_person(const Person& person);

    /** 
     * remove a person from the population of this location.
     * @param person the person leaving
     */
    void remove_person(const Person& person);

    /** 
     * notification that one person in this location changed infection state.
     * @param person the person that changed infection state
     * @param old_state the previous infection state of the person
     */
    void changed_state(const Person& person, InfectionState old_infection_state);

    /** 
     * prepare the location for the next simulation step.
     * @param dt the duration of the simulation step
     * @param global_params global infection parameters
     */
    void begin_step(TimeSpan dt, const GlobalInfectionParameters& global_params);

    /** 
     * number of persons at this location in one infection state.
     * @return number of persons at this location that are in the specified infection state
     */
    unsigned get_subpopulation(InfectionState s) const;

    /**
     * @return parameters of the infection that are specific to this location
     */
    LocalInfectionParameters& get_infection_parameters()
    {
        return m_parameters;
    }

    const LocalInfectionParameters& get_infection_parameters() const
    {
        return m_parameters;
    }

    const std::vector<Cell>& get_cells() const
    {
        return m_cells;
    }

    /**
     * get the number of persons at the location
     * @return number of persons
     */
    int get_population()
    {
        return m_num_persons;
    }

    /**
     * get the exposure rate of the location
     */
    CustomIndexArray<double, AgeGroup, VaccinationState> get_cached_exposure_rate()
    {
        return m_cached_exposure_rate;
    }

    /**
    * Set the capacity of the location in person and volume
    * @param persons maximum number of people that can visit the location at the same time
    * @param volume volume of the location in m^3
    */
    void set_capacity(int persons, int volume)
    {
        m_capacity.persons = persons;
        m_capacity.volume  = volume;
    }

    /**
    * @return the capacity of the location in person and volume
    */
    LocationCapacity get_capacity()
    {
        return m_capacity;
    }

    /**
    * computes a relative transmission risk factor for the location
    * @return the relative risk factor for the location
    */
    double compute_relative_transmission_risk();

    /**
    * Set the capacity adapted transmission risk flag
    * @param consider_capacity if true considers the capacity of the location for the computation of relative 
    * transmission risk
    */
    void set_capacity_adapted_transmission_risk_flag(bool consider_capacity)
    {
        m_capacity_adapted_transmission_risk = consider_capacity;
    }

    /**
     * Store the current number of individuals in the each infection state at a specific time step
     * @param t the time step to save
    */
    void store_current_population(const TimePoint& t);

    /**
     * Return the time series object of the current number of individuals in the each infection state
     * @return the time series object of the current number of individuals in the each infection state
    */
    const TimeSeries<unsigned>& get_time_series_subpopulations() const
    {
        return m_subpopulations_time_series;
    }

private:
    void change_subpopulation(InfectionState s, int delta);
private:
    LocationType m_type;
    uint32_t m_index;
    int m_num_persons = 0;
    LocationCapacity m_capacity;
    bool m_capacity_adapted_transmission_risk;
    TimeSeries<unsigned> m_subpopulations_time_series;
    LocalInfectionParameters m_parameters;
    CustomIndexArray<double, AgeGroup, VaccinationState> m_cached_exposure_rate;
    std::vector<Cell> m_cells;
};

} // namespace abm
} // namespace mio

#endif
