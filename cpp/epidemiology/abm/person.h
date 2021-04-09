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
     * Create a Person.
     * @param location the initial location of the person
     * @param state the initial infection state of the person
     */
    Person(Location& location, InfectionState state, AbmAgeGroup age);

    /** 
     * Time passes and the person interacts with the population at its current location.
     * The person might change infection state.
     * @param dt length of the current simulation time step
     * @param global_infection_parameters infection parameters that are the same in all locations
     */
    void interact(TimeSpan dt, const GlobalInfectionParameters& global_infection_parameters);

    /** 
     * Migrate to a different location.
     * @param location the new location of the person.
     * */
    void migrate_to(Location& location);

    /**
     * Get the current infection state of the person.
     * @returns the current infection state of the person
     */
    InfectionState get_infection_state() const
    {
        return m_state;
    }

    /**
     * Get the age group of this person.
     * @return age.
     */
    AbmAgeGroup get_age() const
    {
        return m_age;
    }

    /**
     * Get the current location of the person.
     * @returns the current location of the person
     */
    const Location& get_location() const
    {
        return m_location;
    }

    /**
     * Get the time the person has been at its current location.
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
};

} // namespace epi

#endif
