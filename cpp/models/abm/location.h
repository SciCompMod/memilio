/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Elisabeth Kluth, David Kerkmann
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

#include "abm/person.h"
#include "abm/parameters.h"
#include "abm/location_type.h"
#include "abm/infection_state.h"
#include "abm/vaccine.h"

#include "memilio/math/eigen.h"
#include "memilio/utils/custom_index_array.h"
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
    std::vector<Person> m_persons;
    CustomIndexArray<double, VirusVariant, AgeGroup> m_cached_exposure_rate;

    Cell(std::vector<Person> persons)
        : m_persons(std::move(persons))
        , m_cached_exposure_rate({{VirusVariant::Count, AgeGroup::Count}, 0.})
    {
    }

    Cell()
        : m_persons{}
        , m_cached_exposure_rate({{VirusVariant::Count, AgeGroup::Count}, 0.})
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
     * @return new infection of the person
     */
    VirusVariant interact(const Person& person, const TimePoint& t, const TimeSpan& dt) const;

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
     * prepare the location for the next simulation step.
     * @param dt the duration of the simulation step
     * @param global_params global infection parameters
     */
    void begin_step(TimePoint t, TimeSpan dt);

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
        return m_persons.size();
    }

    /**
     * get the exposure rate of the location
     */
    CustomIndexArray<double, VirusVariant, AgeGroup> get_cached_exposure_rate()
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

    int get_subpopulation(TimePoint t, InfectionState state = InfectionState::Infected) const;

    Eigen::Ref<const Eigen::VectorXi> get_subpopulations(TimePoint t) const;

    int get_number_infected_total(TimePoint t) const;

private:
    LocationType m_type;
    uint32_t m_index;
    LocationCapacity m_capacity;
    bool m_capacity_adapted_transmission_risk;
    std::vector<Person> m_persons;
    LocalInfectionParameters m_parameters;
    CustomIndexArray<double, VirusVariant, AgeGroup> m_cached_exposure_rate;
    std::vector<Cell> m_cells{};
};

} // namespace abm
} // namespace mio

#endif
