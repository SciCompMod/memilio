/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: Daniel Abele, Majid Abedi, Elisabeth Kluth
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
#include "abm/testing_scheme.h"
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
    using LocationIterator      = PointerDereferencingIterator<std::vector<std::unique_ptr<Location>>::iterator>;
    using ConstLocationIterator = PointerDereferencingIterator<std::vector<std::unique_ptr<Location>>::const_iterator>;
    using PersonIterator        = PointerDereferencingIterator<std::vector<std::unique_ptr<Person>>::iterator>;
    using ConstPersonIterator   = PointerDereferencingIterator<std::vector<std::unique_ptr<Person>>::const_iterator>;

    /**
     * create a World.
     * @param params parameters of the infection that are the same everywhere in the world.
     */
    World(const GlobalInfectionParameters& params = {})
        : m_locations((uint32_t)LocationType::Count)
        , m_infection_parameters(params)
        , m_migration_parameters()
        , m_testing_parameters()
        , m_trip_list()
        , m_use_migration_rules(true)
    {
    }

    //type is move-only for stable references of persons/locations
    World(World&& other)            = default;
    World& operator=(World&& other) = default;
    World(const World&)             = delete;
    World& operator=(const World&)  = delete;

    /** 
     * prepare the world for the next simulation step.
     * @param dt length of the time step 
     */
    void begin_step(TimePoint t, TimeSpan dt);

    /** 
     * evolve the world one time step.
     * @param dt length of the time step
     */
    void evolve(TimePoint t, TimeSpan dt);

    /** 
     * add a location to the world.
     * @param type type of location to add
     * @param num_cells number of cells that the location is divided into
     * @return index and type of the newly created location
     */
    LocationId add_location(LocationType type, uint32_t num_cells = 0);

    /** add a person to the world 
     * @param id index and type of the initial location of the person
     * @param state initial infection state of the person
     * @return reference to the newly created person
     */
    Person& add_person(LocationId id, InfectionState infection_state, AgeGroup age = AgeGroup::Age15to34);

    /**
     * Sets the current infection state of the person.
     * Use only during setup, may distort the simulation results
     * @param person
     * @param inf_state
     */
    void set_infection_state(Person& person, InfectionState inf_state);

    /**
     * get a range of all locations in the world.
     * @return a range of all locations.
     */
    Range<std::pair<std::vector<std::vector<Location>>::const_iterator,
                    std::vector<std::vector<Location>>::const_iterator>>
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
    const Location& get_individualized_location(LocationId id) const;

    Location& get_individualized_location(LocationId id);

    /**
     * get the current location of a person
     * @return reference to the current location of the person
     */
    const Location& get_location(const Person& person) const;

    Location& get_location(const Person& person);

    /**
     * find an assigned location of a person
     * @param type the location type that specifies the assigned location
     * @return pointer to the assigned location
     */
    Location* find_location(LocationType type, const Person& person);

    /** 
     * number of persons in one infection state at all locations of a type.
     * @param type specified location type
     * @return number of persons that are in the specified infection state
     */
    int get_subpopulation_combined(InfectionState s, LocationType type) const;

    /** 
     *get migration parameters
     */
    MigrationParameters& get_migration_parameters();

    const MigrationParameters& get_migration_parameters() const;

    /** 
     *get global infection parameters
     */
    GlobalInfectionParameters& get_global_infection_parameters();

    const GlobalInfectionParameters& get_global_infection_parameters() const;

    /** 
     *get global testing parameters
     */
    GlobalTestingParameters& get_global_testing_parameters();

    const GlobalTestingParameters& get_global_testing_parameters() const;

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

    void add_testing_scheme(const TestingScheme& testing_scheme);
    const std::vector<TestingScheme>& get_testing_schemes() const;
    void set_testing_schemes(const std::vector<TestingScheme> testing_schemes);
    
private:
    void interaction(TimePoint t, TimeSpan dt);
    void migration(TimePoint t, TimeSpan dt);
    
    std::vector<std::unique_ptr<Person>> m_persons;
    std::vector<std::vector<Location>> m_locations;
    std::vector<TestingScheme> m_testing_schemes;
    GlobalInfectionParameters m_infection_parameters;
    MigrationParameters m_migration_parameters;
    GlobalTestingParameters m_testing_parameters;
    TripList m_trip_list;
    bool m_use_migration_rules;
};

} // namespace abm
} // namespace mio

#endif
