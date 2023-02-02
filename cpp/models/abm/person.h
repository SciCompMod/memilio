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
#ifndef EPI_ABM_PERSON_H
#define EPI_ABM_PERSON_H

#include "abm/age.h"
#include "abm/location_type.h"
#include "abm/parameters.h"
#include "abm/time.h"
#include "abm/infection.h"
#include "abm/vaccine.h"
#include "abm/mask_type.h"
#include "abm/mask.h"

#include <functional>

namespace mio
{
namespace abm
{

struct LocationId;
class Location;

static constexpr uint32_t INVALID_PERSON_ID = std::numeric_limits<uint32_t>::max();

/**
 * Agents in the simulated world that can carry and spread the infection.
 */
class Person : public std::enable_shared_from_this<Person>
{
public:
    /**
     * @brief Create a Person.
     * @param location Initial location of the person.
     * @param age The age group of the person.
     * @param person_id Index of the person.
     */
    explicit Person(Location& location, AgeGroup age, uint32_t person_id = INVALID_PERSON_ID);

    /**
    * compare two persons
    */
    bool operator==(const Person& other) const
    {
        return (m_person_id == other.m_person_id);
    }

    /** 
     * @brief Time passes and the person interacts with the population at its current location.
     * The person might become infected.
     * @param[in] t Current time.
     * @param[in] dt Length of the current simulation time step.
     * @param[in,out] global_infection_parameters Infection parameters that are the same in all locations.
     */
    void interact(TimePoint t, TimeSpan dt, GlobalInfectionParameters& params);

    /** 
     * @brief Migrate to a different location.
     * @param[in] loc_new The new location of the person.
     * @param[in] cells_new The new cells of the person.
     * */
    void migrate_to(Location& loc_old, Location& loc_new, const std::vector<uint32_t>& cells_new = {});

    /**
     * @brief Get the latest Infection of the Person.
     * @returns The latest Infection of the Person.
     */
    Infection& get_infection()
    {
        return m_infections.back();
    }

    const Infection& get_infection() const
    {
        return m_infections.back();
    }

    /** 
     * @returns All vaccinations.
    */
    std::vector<Vaccination>& get_vaccinations()
    {
        return m_vaccinations;
    }

    const std::vector<Vaccination>& get_vaccinations() const
    {
        return m_vaccinations;
    }

    /**
     * @brief Returns if the person is infected at the time point.
     * @param[in] t Time point of querry. Usually the current time of the simulation.
     * @returns True if the person is infected at the time point.
    */
    bool is_infected(TimePoint t) const;

    /**
     * @param[in] t Time point of querry. Usually the current time of the simulation.
     * @returns The infection state of the latest infection at time t.
    */
    InfectionState get_infection_state(TimePoint t) const;

    /**
     * @brief Adds a new infection to the list of infections.
     * @param[in] inf The new infection.
    */
    void add_new_infection(Infection&& inf);

    /**
     * Get the age group of this person.
     * @return Age.
     */
    AgeGroup get_age() const
    {
        return m_age;
    }

    /**
     * @brief Get the current Location of the Person.
     * @returns Current Location of the Person.
     */
    Location& get_location();

    const Location& get_location() const;

    /**
     * @brief Get the time the person has been at its current location.
     * @return Time span.
     */
    TimeSpan get_time_at_location() const
    {
        return m_time_at_location;
    }

    /**
     * @brief Get the time since the person has been testes.
     * @return Time span.
     */
    TimeSpan get_time_since_negative_test() const
    {
        return m_time_since_negative_test;
    }
    /**
     * @brief Set an assigned location of the person. The assigned location is saved by its index.
     * Assume that a person has at most one assigned location per location type.
     * @param[in] location The new assigned location.
     */
    void set_assigned_location(Location& location);

    /**
     * Set an assigned location of the person. The assigned location is saved by an index.
     * Assume that a person has at most one assigned location of a certain location type.
     * @param[in] location_type Location type of the new assigned location.
     * @param[in] index Index of the new assigned location.
     */
    void set_assigned_location(LocationId id);

    /**
     * @brief Returns the index of a assigned location of the person.
     * Assume that a person has at most one assigned location of a certrain location type.
     * @param[in] type Location type of the assigned location.
     */
    uint32_t get_assigned_location_index(LocationType type) const;

    /**
     * @brief Returns the assigned locations of the person.
     */
    const std::vector<uint32_t>& get_assigned_locations() const
    {
        return m_assigned_locations;
    }

    /**
     * Every person has a random number.
     * Depending on this number and the time, the person works from home in case of a lockdown.
     * @return If the person works from home.
     */
    bool goes_to_work(TimePoint t, const MigrationParameters& params) const;

    /**
     * Every person has a random number to determine what time to go to work.
     * Depending on this number person decides what time has to go to work;
     * @return The time of going to work.
     */
    TimeSpan get_go_to_work_time(const MigrationParameters& params) const;

    /**
     * Every person has a random number that determines if they go to school in case of a lockdown.
     * @return If the person goes to school.
     */
    bool goes_to_school(TimePoint t, const MigrationParameters& params) const;

    /**
     * Every person has a random number to determine what time to go to school.
     * Depending on this number person decides what time has to go to school;
     * @return The time of going to school.
     */
    TimeSpan get_go_to_school_time(const MigrationParameters& params) const;

    /**
     * @brief Answers the question if a person is currently in quarantine.
     * @return If the person is in quarantine.
     */
    bool is_in_quarantine() const
    {
        return m_quarantine;
    }

    /**
     * @brief Simulates a Corona test and returns the test result of the person.
     * If the test is positive, the person has to quarantine.
     * If the test is negative, quarantine ends.
     * @param[in] t TimePoint of the test.
     * @param[in] params Sensitivity and specificity of the test method.
     * @return True if the test result of the person is positive.
     */
    bool get_tested(TimePoint t, const TestParameters& params);

    /**
     * @brief Get the person id of the person.
     * The person id should correspond to the index in m_persons in world.
     */
    uint32_t get_person_id();

    /**
     * @brief Get index of cells of the person.
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
    ScalarType get_mask_protective_factor(const GlobalInfectionParameters& params) const;

    /**
     * @brief For every LocationType a person has a compliance value between -1 and 1.
     * -1 means that the Person never complies to any mask duty at the given LocationType.
     * 1 means that the Person always wears a Mask a the LocationType even if it is not required.
     * @param preferences The vector of mask compliance values for all LocationTypes.
     */
    void set_mask_preferences(std::vector<ScalarType> preferences)
    {
        m_mask_compliance = preferences;
    }

    /**
     * @brief Get the mask compliance of the person for the current location.
     * @param[in] location The current location of the person.
     * @return The probability that the person does not comply to any mask duty/wears a
     * mask even if it is not required.
     */
    ScalarType get_mask_compliance(LocationType location) const
    {
        return m_mask_compliance[static_cast<int>(location)];
    }

    /**
     * @brief Checks whether the person wears a mask at the target location.
     * @param[in] target The target location.
     */
    bool apply_mask_intervention(const Location& target);

    /**
     * @brief Decide if a person is currently wearing a mask.
     * @param[in] wear_mask If true, the protection of the mask is considered when
     * computing the exposure rate.
     */
    void set_wear_mask(bool wear_mask)
    {
        m_wears_mask = wear_mask;
    }

    /**
     * @return True if the person is currently wearing a mask.
     */
    bool get_wear_mask() const
    {
        return m_wears_mask;
    }

private:
    struct ImmunityLevel {
        /**
         * @brief Get the multiplicative factor on how likely an infection is due to the immune system.
         * @param[in] v VirusVariant to take into consideration.
         * @param[in] t TimePoint of check.
         * @returns Protection factor of the immune system to the given VirusVariant at the given TimePoint.
        */
        ScalarType get_protection_factor(VirusVariant /*v*/, TimePoint /*t*/) const
        {
            return 1.; // put implementation in .cpp
        }

        /**
         * @brief Get the multiplicative factor on how severe a new infection is due to the immune system.
         * @param[in] v VirusVariant to take into consideration.
         * @param[in] t TimePoint of check.
         * @returns Severity factor of a new infection with the given VirusVariant at the given TimePoint.
        */
        ScalarType get_severity_factor(VirusVariant /*v*/, TimePoint /*t*/) const
        {
            return 1.; // put implementation in .cpp
        }
    };

public:
    /**
     * @brief Get the multiplicative factor on how likely an infection is due to the immune system.
     * @param[in] v VirusVariant to take into consideration.
     * @param[in] t TimePoint of check.
     * @returns Protection factor of the immune system to the given VirusVariant at the given TimePoint.
    */
    ScalarType get_protection_factor(VirusVariant v, TimePoint t) const
    {
        return m_immunity_level.get_protection_factor(v, t);
    }
    //ScalarType get_severity_factor = ImmunityLevel::get_severity_factor;

private:
    std::shared_ptr<Location> m_location;
    std::vector<uint32_t> m_assigned_locations;
    std::vector<Vaccination> m_vaccinations;
    std::vector<Infection> m_infections;
    ImmunityLevel m_immunity_level; // do we need this? Or can we call p.ImmunityLevel.get_...() ?
    bool m_quarantine = false;
    AgeGroup m_age;
    TimeSpan m_time_at_location;
    ScalarType m_random_workgroup;
    ScalarType m_random_schoolgroup;
    ScalarType m_random_goto_work_hour;
    ScalarType m_random_goto_school_hour;
    TimeSpan m_time_since_negative_test;
    Mask m_mask;
    bool m_wears_mask = false;
    std::vector<ScalarType> m_mask_compliance;
    uint32_t m_person_id;
    std::vector<uint32_t> m_cells;
};

} // namespace abm
} // namespace mio

#endif
