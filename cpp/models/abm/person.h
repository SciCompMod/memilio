/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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
#ifndef EPI_ABM_PERSON_H
#define EPI_ABM_PERSON_H

#include "abm/state.h"
#include "abm/age.h"
#include "abm/time.h"
#include "abm/parameters.h"
#include "abm/world.h"
#include "abm/time.h"

#include <functional>

namespace mio
{

class Location;

/**
 * LocationId identifies a Location uniquely. It consists of the LocationType of the Location and an Index.
 * The index corresponds to the index into the structure m_locations from world, where all Locations are saved.
 */
struct InfectionProperties {
    InfectionProperties(InfectionState infection_state, bool infection_detected = false)
        : state(infection_state)
        , detected(infection_detected)
    {
    }
    InfectionState state;
    bool detected;
};

/**
 * Agents in the simulated world that can carry and spread the infection.
 */
class Person
{
public:
    
    /**
     * create a Person.
     * @param id index and type of the initial location of the person
     * @param infection_properties the initial infection state of the person and if infection is detected
     * @param vaccination_state the initial infection state of the person
     * @param age the age group of the person
     * @param global_params the global infection parameters
     */
    Person(LocationId id, InfectionProperties infection_properties, VaccinationState vaccination_state, AbmAgeGroup age, const GlobalInfectionParameters& global_params);
    /**
     * create a Person.
     * @param id index and type of the initial location of the person
     * @param infection_properties the initial infection state of the person and if infection is detected
     * @param age the age group of the person
     * @param global_params the global infection parameters
     */
    Person(LocationId id, InfectionProperties infection_properties, AbmAgeGroup age, const GlobalInfectionParameters& global_params);
    
    /**
     * create a Person.
     * @param location the initial location of the person
     * @param infection_properties the initial infection state of the person and if infection is detected
     * @param age the age group of the person
     * @param global_params the global infection parameters
     */
    Person(Location& location, InfectionProperties infection_properties, AbmAgeGroup age, const GlobalInfectionParameters& global_params);
    
    /**
     * create a Person.
     * @param location the initial location of the person
     * @param infection_properties the initial infection state of the person and if infection is detected
     * @param age the age group of the person
     * @param global_params the global infection parameters
     */
    Person(Location& location, InfectionProperties infection_properties, VaccinationState vaccination_state, AbmAgeGroup age, const GlobalInfectionParameters& global_params);

    /** 
     * Time passes and the person interacts with the population at its current location.
     * The person might change infection state.
     * @param dt length of the current simulation time step
     * @param global_infection_parameters infection parameters that are the same in all locations
     */
    void interact(TimeSpan dt, const GlobalInfectionParameters& global_infection_parameters, Location& loc,
                  const GlobalTestingParameters& global_testing_params);

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
        return m_infection_state;
    }
    
    /**
     * Get the current vaccination state of the person.
     * @returns the current vaccination state of the person
     */
    VaccinationState get_vaccination_state() const
    {
        return m_vaccination_state;
    }
    
    /**
     * Sets the current infection state of the person.
     */
    void set_infection_state(InfectionState inf_state);
    
    /**
     * Get the age group of this person.
     * @return age.
     */
    mio::Index<AbmAgeGroup> get_age() const
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
     * Get the time since the person has been testes.
     * @return time span.
     */
    TimeSpan get_time_since_negative_test() const
    {
        return m_time_since_negative_test;
    }
    /**
     * set an assigned location of the person. The assigned location is saved by its index.
     * Assume that a person has at most one assigned location per location type.
     * @param location the new assigned location
     */
    void set_assigned_location(Location& location);

    /**
     * set an assigned location of the person. The assigned location is saved by an index.
     * Assume that a person has at most one assigned location of a certain location type
     * @param location_type location type of the new assigned location
     * @param index index of the new assigned location
     */
    void set_assigned_location(LocationId id);

    /**
     * returns the index of a assigned location of the person.
     * Assume that a person has at most one assigned location of a certrain location type.
     * @param type location type of the assigned location
     */
    uint32_t get_assigned_location_index(LocationType type) const;

    /**
     *returns the assigned locations of the person.
     */
    const std::vector<uint32_t>& get_assigned_locations() const
    {
        return m_assigned_locations;
    }

    /**
     * Every person has a random number.
     * Depending on this number and the time, the person works from home in case of a lockdown.
     * @return if the person works from home
     */
    bool goes_to_work(TimePoint t, const AbmMigrationParameters& params) const;

    /**
     * Every person has a random number to determine what time to go to work.
     * Depending on this number person decides what time has to go to work;
     * @return the hour of going to work
     */
    int get_go_to_work_hour(const AbmMigrationParameters& params) const;

    /**
     * Every person has a random number that determines if they go to school in case of a lockdown.
     * @return if the person goes to school
     */
    bool goes_to_school(TimePoint t, const AbmMigrationParameters& params) const;

    /**
     * Answers the question if a person is currently in quarantine.
     * @return if the person is in quarantine
     */
    bool is_in_quarantine() const
    {
        return m_quarantine;
    }

    /**
     * Simulates a Corona test and returns the test result of the person.
     * If the test is positive, the person has to quarantine.
     * If the test is negative, quarantine ends.
     * @param params sensitivity and specificity of the test method
     * @return true if the test result of the person is positive
     */
    bool get_tested(const TestParameters& params);

private:
    LocationId m_location_id;
    std::vector<uint32_t> m_assigned_locations;
    InfectionState m_infection_state;
    VaccinationState m_vaccination_state;
    TimeSpan m_time_until_carrier;
    bool m_quarantine;
    mio::Index<AbmAgeGroup> m_age;
    TimeSpan m_time_at_location;
    double m_random_workgroup;
    double m_random_schoolgroup;
    double m_random_goto_work_hour;
    TimeSpan m_time_since_negative_test;
};

} // namespace mio

#endif
