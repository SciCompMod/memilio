#ifndef EPI_ABM_PERSON_H
#define EPI_ABM_PERSON_H

#include "epidemiology/abm/state.h"
#include <functional>

namespace epi
{

class Location;
class GlobalInfectionParameters;

/**
 * Agents in the simulated world that can carry and spread the infection.
 */
class Person
{
public:
    /**
     * create a Person.
     * @param location the initial location of the person
     * @param state the initial infection state of the person
     */
    Person(Location& location, InfectionState state);

    /** 
     * interact with the population at its current location.
     * The person might change infection state 
     * @param dt length of the current simulation time step
     * @param global_infection_parameters infection parameters that are the same in all locations
     */
    void interact(double dt, const GlobalInfectionParameters& global_infection_parameters);

    /** 
     * migrate to a different location.
     * @param location the new location of the person.
     * */
    void migrate_to(Location& location);

    /**
     * get the current infection state of the person.
     * @returns the current infection state of the person
     */
    InfectionState get_infection_state() const
    {
        return m_state;
    }

    /**
     * get the current location of the person.
     * @returns the current location of the person
     */
    const Location& get_location() const
    {
        return m_location;
    }

private:
    InfectionState m_state;
    float m_time_until_carrier;
    std::reference_wrapper<Location> m_location;
    //age, times, ...
};

} // namespace epi

#endif