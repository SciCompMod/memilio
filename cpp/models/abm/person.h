/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Elisabeth Kluth
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

#include "abm/location_type.h"
#include "abm/mask_type.h"
#include "abm/infection_state.h"
#include "abm/age.h"
#include "abm/time.h"
#include "abm/parameters.h"
#include "abm/location.h"
#include "abm/mask.h"

#include <functional>
#include <vector>

namespace mio
{
namespace abm
{

class Location;

/**
 * Infection properties describe the infection state of a person and if a infection is detected
 */
struct InfectionProperties {
    InfectionProperties(InfectionState infection_state, bool infection_detected = false)
        : state(infection_state)
        , detected(infection_detected)
    {
    }
    InfectionState state;
    bool detected;

    bool operator==(const InfectionProperties& other) const
    {
        return std::tie(state, detected) == std::tie(other.state, detected);
    }
    bool operator!=(const InfectionProperties& other) const
    {
        return std::tie(state, detected) != std::tie(other.state, detected);
    }
};

static constexpr uint32_t INVALID_PERSON_ID = std::numeric_limits<uint32_t>::max();

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
     * @param person_id index of the person
     */
    Person(LocationId id, InfectionProperties infection_properties, AgeGroup age,
           const GlobalInfectionParameters& global_params,
           VaccinationState vaccination_state = VaccinationState::Unvaccinated, uint32_t person_id = INVALID_PERSON_ID);

    /**
     * create a Person.
     * @param location the initial location of the person
     * @param infection_properties the initial infection state of the person and if infection is detected
     * @param age the age group of the person
     * @param global_params the global infection parameters
     * @param person_id index of the person
     */
    Person(Location& location, InfectionProperties infection_properties, AgeGroup age,
           const GlobalInfectionParameters& global_params,
           VaccinationState vaccination_state = VaccinationState::Unvaccinated, uint32_t person_id = INVALID_PERSON_ID);

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
     * @param cells_new the new cells of the person.
     * */
    void migrate_to(Location& loc_old, Location& loc_new, const std::vector<uint32_t>& cells_new = {});

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
    AgeGroup get_age() const
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
    bool goes_to_work(TimePoint t, const MigrationParameters& params) const;

    /**
     * Every person has a random number to determine what time to go to work.
     * Depending on this number person decides what time has to go to work;
     * @return the time of going to work
     */
    TimeSpan get_go_to_work_time(const MigrationParameters& params) const;

    /**
     * Every person has a random number that determines if they go to school in case of a lockdown.
     * @return if the person goes to school
     */
    bool goes_to_school(TimePoint t, const MigrationParameters& params) const;

    /**
     * Every person has a random number to determine what time to go to school.
     * Depending on this number person decides what time has to go to school;
     * @return the time of going to school
     */
    TimeSpan get_go_to_school_time(const MigrationParameters& params) const;

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

    /**
     * get the person id of the person
     * the person id should correspondet to the index in m_persons in world
     */
    uint32_t get_person_id();

    /**
     * get index of cells of the person
     */
    std::vector<uint32_t>& get_cells();

    const std::vector<uint32_t>& get_cells() const;

    /**
     * @brief Get the mask of the person.
     * @return Current mask of the person.
     */
    Mask& get_mask()
    {
        return m_mask;
    }

    const Mask& get_mask() const
    {
        return m_mask;
    }

    /**
     * @brief Get the protection of the mask. A value of 1 represents full protection and a value of 0 means no protection.
     * @return The protection factor of the mask.
     */
    double get_protective_factor(const GlobalInfectionParameters& params) const;

    /**
     * @brief For every LocationType a person has a compliance value between -1 and 1.
     * -1 means that the Person never complies to any mask duty at the given LocationType.
     * 1 means that the Person always wears a Mask a the LocationType even if it is not required.
     * @param preferences The vector of mask compliance values for all LocationTypes.
     */
    void set_mask_preferences(std::vector<double> preferences)
    {
        m_mask_compliance = preferences;
    }

    /**
     * @brief Get the mask compliance of the person for the current location.
     * @param location the current location of the person
     * @return The probability that the person does not comply to any mask duty/wears a
     * mask even if it is not required.
     */
    double get_mask_compliance(LocationType location) const
    {
        return m_mask_compliance[static_cast<int>(location)];
    }

    /**
     * @brief Checks whether the person wears a mask at the target location.
     * @param target The target location.
     */
    bool apply_mask_intervention(const Location& target);

    /**
     * @brief Decide if a person is currently wearing a mask.
     * @param wear_mask if true the protection of the mask is considered when
     * computing the exposure rate
     */
    void set_wear_mask(bool wear_mask)
    {
        m_wears_mask = wear_mask;
    }

    /**
     * @return true if the person is currently wearing a mask
     */
    bool get_wear_mask() const
    {
        return m_wears_mask;
    }

private:
    LocationId m_location_id;
    std::vector<uint32_t> m_assigned_locations;
    InfectionState m_infection_state;
    VaccinationState m_vaccination_state;
    TimeSpan m_time_until_infected_no_symptoms;
    bool m_quarantine;
    AgeGroup m_age;
    TimeSpan m_time_at_location;
    double m_random_workgroup;
    double m_random_schoolgroup;
    double m_random_goto_work_hour;
    double m_random_goto_school_hour;
    TimeSpan m_time_since_negative_test;
    Mask m_mask;
    bool m_wears_mask;
    std::vector<double> m_mask_compliance;
    uint32_t m_person_id;
    std::vector<uint32_t> m_cells;
};

} // namespace abm
} // namespace mio

#endif
