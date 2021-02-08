#ifndef EPI_ABM_PERSON_H
#define EPI_ABM_PERSON_H

#include "epidemiology/abm/state.h"
#include "epidemiology/abm/age.h"
#include "epidemiology/abm/time.h"

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
    Person(Location& location, InfectionState state, AbmAgeGroup age);

    /** 
     * interact with the population at its current location.
     * The person might change infection state 
     * @param dt length of the current simulation time step
     * @param global_infection_parameters infection parameters that are the same in all locations
     */
    void interact(TimeSpan dt, const GlobalInfectionParameters& global_infection_parameters);

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
     * get the age of this person.
     * @return age.
     */
    AbmAgeGroup get_age() const
    {
        return m_age;
    }

    /**
     * get the current location of the person.
     * @returns the current location of the person
     */
    const Location& get_location() const
    {
        return m_location;
    }

    /**
     * how long has the person been at its current location.
     * @return time span.
     */
    TimeSpan get_time_at_location() const
    {
        return m_time_at_location;
    }

private:
    std::reference_wrapper<Location> m_location;
    InfectionState m_state;
    TimeSpan m_time_until_carrier;
    AbmAgeGroup m_age;
    TimeSpan m_time_at_location;
    //age, times, ...
};

} // namespace epi

#endif
