#ifndef EPI_ABM_PERSON_H
#define EPI_ABM_PERSON_H

#include "epidemiology/abm/state.h"
#include "epidemiology/abm/age.h"
#include "epidemiology/abm/time.h"
#include "epidemiology/abm/parameters.h"
#include "epidemiology/abm/world.h"
#include "epidemiology/abm/time.h"

#include <functional>

namespace epi
{

class Location;



/**
 * Agents in the simulated world that can carry and spread the infection.
 */
class Person
{
public:
    /**
     * create a Person.
     * @param id index and type of the initial location of the person
     * @param state the initial infection state of the person
     * @param age the age group of the person
     * @param global_params the global infection parameters
     */
    Person(LocationId id, InfectionState state, AbmAgeGroup age, const GlobalInfectionParameters& global_params);
    
    /**
     * create a Person.
     * @param location the initial location of the person
     * @param state the initial infection state of the person
     * @param age the age group of the person
     * @param global_params the global infection parameters
     */
    Person(Location& location, InfectionState state, AbmAgeGroup age, const GlobalInfectionParameters& global_params);

    /** 
     * Time passes and the person interacts with the population at its current location.
     * The person might change infection state.
     * @param dt length of the current simulation time step
     * @param global_infection_parameters infection parameters that are the same in all locations
     */
    void interact(TimeSpan dt, const GlobalInfectionParameters& global_infection_parameters, Location& loc);

    /** 
     * migrate to a different location.
     * @param loc_new the new location of the person.
     * */
    void migrate_to(Location& loc_old, Location& loc_new);

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
    epi::Index<AbmAgeGroup> get_age() const
    {
        return m_age;
    }

    /**
     * get index and type of the current location of the person.
     * @returns index and type of the current location of the person
     */
    LocationId get_location_id() const
    {
        return m_location_id;
    } 

    /**
     * Get the time the person has been at its current location.
     * @return time span.
     */
    TimeSpan get_time_at_location() const
    {
        return m_time_at_location;
    }
    /**
     * set an assigned location of the person. The assigned location is saved by its index.
     * Assume that a person has at most one assigned location per location type.
     * @param location the new assigned location
     */
    void set_assigned_location (Location& location);

    /**
     * set an assigned location of the person. The assigned location is saved by an index.
     * Assume that a person has at most one assigned location of a certain location type
     * @param location_type location type of the new assigned location
     * @param index index of the new assigned location
     */
    void set_assigned_location (LocationId id);

    /**
     * returns the index of a assigned location of the person.
     * Assume that a person has at most one assigned location of a certrain location type.
     * @param type location type of the assigned location
     */
    uint32_t get_assigned_location_index (LocationType type) const;

    /**
     *returns the assigned locations of the person.
     */
    const std::vector<uint32_t>& get_assigned_locations () const{
        return m_assigned_locations;
    }

    /**
     * Every person has a random number.
     * Depending on this number and the time, the person works from home in case of a lockdown.
     * @return if the person works from home
     */
    bool goes_to_work(TimePoint t, const AbmMigrationParameters& params) const;

    /**
     * Every person has a random number that determines if they go to school in case of a lockdown.
     * @return if the person goes to school
     */
    bool goes_to_school(TimePoint t, const AbmMigrationParameters& params) const;

    /**
     * Answers the question if a person is currently in quarantine.
     * @return if the person is in quarantine
     */
    bool is_in_quarantine () const {
        return m_quarantine;
    }

private:
    LocationId m_location_id;
    std::vector<uint32_t> m_assigned_locations;
    InfectionState m_state;
    TimeSpan m_time_until_carrier;
    bool m_quarantine;
    epi::Index<AbmAgeGroup> m_age;
    TimeSpan m_time_at_location;
    double m_random_workgroup;
    double m_random_schoolgroup;
};

} // namespace epi

#endif
