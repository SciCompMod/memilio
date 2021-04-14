#ifndef EPI_ABM_WORLD_H
#define EPI_ABM_WORLD_H

#include "epidemiology/abm/age.h"
#include "epidemiology/abm/parameters.h"
#include "epidemiology/abm/location.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/utils/pointer_dereferencing_iterator.h"
#include "epidemiology/utils/stl_util.h"

#include <vector>
#include <memory>

namespace epi
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
    World(const GlobalInfectionParameters& params = {AbmAgeGroup::Count})
        : m_infection_parameters(params)
    {
    }

    //type is move-only for stable references of persons/locations
    World(World&& other) = default;
    World& operator=(World&& other) = default;
    World(const World&)             = delete;
    World& operator=(const World&) = delete;

    /** 
     * prepare the world for the next simulation step.
     * @param dt length of the time step 
     */
    void begin_step(double dt);

    /** 
     * evolve the world one time step.
     * @param dt length of the time step
     */
    void evolve(double dt);

    /** 
     * add a location to the world.
     * @param type type of location to add
     * @return reference to the newly created location
     */
    Location& add_location(LocationType type);

    /** add a person to the world 
     * @param location initial location of the person
     * @param state initial infection state of the person
     * @return reference to the newly created person
     */
    Person& add_person(Location& location, InfectionState state, Index<AbmAgeGroup> age = AbmAgeGroup::Age15to34);

    /**
     * get a range of all locations in the world.
     * @return a range of all locations.
     */
    Range<std::pair<ConstLocationIterator, ConstLocationIterator>> get_locations() const;

    /**
     * get a range of all persons in the world.
     * @return a range of all persons.
     */
    Range<std::pair<ConstPersonIterator, ConstPersonIterator>> get_persons() const;

private:
    void interaction(double dt);
    void migration(double dt);

    std::vector<std::unique_ptr<Person>> m_persons;
    std::vector<std::unique_ptr<Location>> m_locations;
    GlobalInfectionParameters m_infection_parameters;
};

} // namespace epi

#endif
