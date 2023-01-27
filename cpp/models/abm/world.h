/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: Daniel Abele, Majid Abedi, Elisabeth Kluth, David Kerkmann, Sascha Korf, Martin J. Kuehn
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
#ifndef EPI_ABM_WORLD_H
#define EPI_ABM_WORLD_H

#include "abm/parameters.h"
#include "abm/location.h"
#include "abm/person.h"
#include "abm/lockdown_rules.h"
#include "abm/trip_list.h"
#include "abm/testing_strategy.h"
#include "memilio/utils/pointer_dereferencing_iterator.h"
#include "memilio/utils/stl_util.h"

#include <vector>
#include <memory>

namespace mio
{
namespace abm
{

/**
 * The world of the simulation.
 * consists of Locations and Persons (Agents)
 */
class World
{
public:
    using LocationIterator      = PointerDereferencingIterator<std::vector<std::shared_ptr<Location>>::iterator>;
    using ConstLocationIterator = PointerDereferencingIterator<std::vector<std::shared_ptr<Location>>::const_iterator>;
    using PersonIterator        = PointerDereferencingIterator<std::vector<std::shared_ptr<Person>>::iterator>;
    using ConstPersonIterator   = PointerDereferencingIterator<std::vector<std::shared_ptr<Person>>::const_iterator>;

    /**
     * create a World.
     * @param params parameters of the infection that are the same everywhere in the world.
     */
    World(const GlobalInfectionParameters& params = {})
        : m_locations((uint32_t)LocationType::Count)
        , m_infection_parameters(params)
        , m_migration_parameters()
        , m_trip_list()
    {
        use_migration_rules(true);
    }

    //type is move-only for stable references of persons/locations
    World(World&& other)            = default;
    World& operator=(World&& other) = default;
    World(const World&)             = delete;
    World& operator=(const World&)  = delete;

    /** 
     * Prepare the World for the next simulation step.
     * @param[in] t Current time.
     * @param[in] dt Length of the time step.
     */
    void begin_step(TimePoint t, TimeSpan dt);

    /** 
     * Evolve the world one time step.
     * @param[in] t Current time.
     * @param[in] dt Length of the time step.
     */
    void evolve(TimePoint t, TimeSpan dt);

    /** 
     * Add a Location to the world.
     * @param[in] type Type of Location to add.
     * @param[in] num_cells [Default: 1] Number of Cell%s that the Location is divided into.
     * @return Index and type of the newly created Location.
     */
    LocationId add_location(LocationType type, uint32_t num_cells = 1);

    /** 
     * @brief Add a Person to the world.
     * @param[in] id Index and type of the initial Location of the Person.
     * @param[in] age AgeGroup of the person.
     * @return Reference to the newly created Person.
     */
    Person& add_person(const LocationId id, AgeGroup age);

    /**
     * Sets the current infection state of the person.
     * Warning: Use only during setup, may distort the simulation results!
     * @param[in,out] person Person to set infection state.
     * @param[in] inf_state Infection state to set.
     * @param[in] t [Default: 0] TimePoint of initialization of the infection.
     */
    void set_infection_state(Person& person, InfectionState inf_state, TimePoint t = TimePoint(0));

    /**
     * get a range of all locations in the world.
     * @return a range of all locations.
     */
    Range<std::pair<std::vector<std::vector<std::shared_ptr<Location>>>::const_iterator,
                    std::vector<std::vector<std::shared_ptr<Location>>>::const_iterator>>
    get_locations() const;

    /**
     * get a range of all persons in the world.
     * @return a range of all persons.
     */
    Range<std::pair<ConstPersonIterator, ConstPersonIterator>> get_persons() const;

    /**
     * get an individualized location
     * @param id LocationId of the location
     * @return reference to the location
     */
    const std::shared_ptr<Location>& get_individualized_location(LocationId id) const;

    std::shared_ptr<Location>& get_individualized_location(LocationId id);

    /**
     * get the current location of a person
     * @return reference to the current location of the person
     */
    const std::shared_ptr<Location>& get_location(const Person& person) const;

    std::shared_ptr<Location>& get_location(Person& person);

    /**
     * find an assigned location of a person
     * @param type the location type that specifies the assigned location
     * @return pointer to the assigned location
     */
    std::shared_ptr<Location> find_location(LocationType type, const Person& person);

    /** 
     * number of persons in one infection state at all locations of a type.
     * @param type specified location type
     * @return number of persons that are in the specified infection state
     */
    int get_subpopulation_combined(TimePoint t, InfectionState s, LocationType type) const;

    /** 
     * get migration parameters
     */
    MigrationParameters& get_migration_parameters();

    const MigrationParameters& get_migration_parameters() const;

    /** 
     * get global infection parameters
     */
    GlobalInfectionParameters& get_global_infection_parameters();

    const GlobalInfectionParameters& get_global_infection_parameters() const;

    /**
     * get migration data
     */
    TripList& get_trip_list();

    const TripList& get_trip_list() const;

    /** 
     * decide if migration rules (like go to school/work) are used or not
     * the migration rules regarding hospitalization/ICU/quarantine are always used
     */
    void use_migration_rules(bool param);
    bool use_migration_rules() const;

    /** 
     * get testing strategy
     */
    TestingStrategy& get_testing_strategy();

    const TestingStrategy& get_testing_strategy() const;

private:
    void interaction(TimePoint t, TimeSpan dt);
    void migration(TimePoint t, TimeSpan dt);

    std::vector<std::shared_ptr<Person>> m_persons;
    std::vector<std::vector<std::shared_ptr<Location>>> m_locations;
    TestingStrategy m_testing_strategy;
    GlobalInfectionParameters m_infection_parameters;
    MigrationParameters m_migration_parameters;
    TripList m_trip_list;
    bool m_use_migration_rules;
    std::vector<std::pair<LocationType (*)(const Person&, TimePoint, TimeSpan, const MigrationParameters&),
                          std::vector<LocationType>>>
        m_migration_rules;
};

} // namespace abm
} // namespace mio

#endif
