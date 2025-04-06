/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Elisabeth Kluth, David Kerkmann, Khoa Nguyen
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
#include "abm/infection_state.h"
#include "abm/parameters.h"
#include "abm/time.h"
#include "abm/vaccine.h"
#include "abm/mask_type.h"
#include "abm/mask.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/memory.h"
#include "abm/movement_data.h"
#include <functional>

namespace mio
{
namespace abm
{

struct LocationId;
class Location;
class Infection;

static constexpr uint32_t INVALID_PERSON_ID = std::numeric_limits<uint32_t>::max();

/**
 * @brief Agents in the simulated World that can carry and spread the Infection.
 */
class Person
{
public:
    /**
    * Random number generator of individual persons.
    * Increments the random number generator counter of the person when used.
    * Does not store its own key or counter.
    * Instead the key needs to be provided from the outside, so that the RNG
    * for all persons share the same key.
    * The counter is taken from the person.
    * Person::RandomNumberGenerator is cheap to construct and transparent
    * for the compiler to optimize, so we don't store the RNG persistently, only the 
    * counter, so we don't need to store the key in each person. This increases
    * consistency (if the key is changed after the person is created) and 
    * reduces the memory required per person.
    * @see mio::RandomNumberGeneratorBase
    */
    class RandomNumberGenerator : public RandomNumberGeneratorBase<RandomNumberGenerator>
    {
    public:
        /**
        * Creates a RandomNumberGenerator for a person.
        * @param key Key to be used by the generator.
        * @param id Id of the Person.
        * @param counter Reference to the Person's RNG Counter. 
        */
        RandomNumberGenerator(Key<uint64_t> key, uint32_t id, Counter<uint32_t>& counter)
            : m_key(key)
            , m_person_id(id)
            , m_counter(counter)
        {
        }

        /**
        * Creates a RandomNumberGenerator for a person.
        * Uses the same key as another RandomNumberGenerator.
        * @param rng RandomNumberGenerator who's key will be used.
        * @param person Reference to the Person who's counter will be used. 
        */
        RandomNumberGenerator(const mio::RandomNumberGenerator& rng, Person& person)
            : RandomNumberGenerator(rng.get_key(), person.get_person_id(), person.get_rng_counter())
        {
        }

        /**
        * @return Get the key.
        */
        Key<uint64_t> get_key() const
        {
            return m_key;
        }

        /**
        * @return Get the current counter.
        */
        Counter<uint64_t> get_counter() const
        {
            return rng_totalsequence_counter<uint64_t>(m_person_id, m_counter);
        }

        /**
        * Increment the counter.
        */
        void increment_counter()
        {
            ++m_counter;
        }

    private:
        Key<uint64_t> m_key; ///< Global RNG Key
        uint32_t m_person_id; ///< Id of the Person
        Counter<uint32_t>& m_counter; ///< Reference to the Person's rng counter
    };

    /**
     * @brief Create a Person.
     * @param[in, out] rng RandomNumberGenerator.
     * @param[in, out] location Initial Location of the Person.
     * @param[in] age The AgeGroup of the Person.
     * @param[in] person_id Index of the Person.
     */
    explicit Person(mio::RandomNumberGenerator& rng, Location& location, AgeGroup age,
                    uint32_t person_id = INVALID_PERSON_ID);

    /**
     * @brief Create a copy of this #Person object with a new Location.
     * @param[in, out] location The new #Location of the Person.
     */
    Person copy_person(Location& location);

    /**
     * @brief Compare two Person%s.
     */
    bool operator==(const Person& other) const
    {
        return (m_person_id == other.m_person_id);
    }

    /** 
     * @brief Time passes and the Person interacts with the population at its current Location.
     * The Person might become infected.
     * @param[in] t Current time.
     * @param[in] dt Length of the current Simulation TimeStep.
     * @param[in, out] params Infection parameters that are the same in all Location%s.
     */
    void interact(RandomNumberGenerator& rng, TimePoint t, TimeSpan dt, const Parameters& params);

    /** 
     * @brief Migrate to a different Location.
     * @param[in, out] loc_new The new Location of the Person.
     * @param[in] cells_new The Cell%s that the Person visits at the new Location.
     * */
    void migrate_to(Location& loc_new, const std::vector<uint32_t>& cells_new = {0})
    {
        migrate_to(loc_new, TransportMode::Unknown, cells_new);
    }

    /** 
     * @brief Migrate to a different Location.
     * @param[in] loc_new The new Location of the Person.
     * @param[in] transport_mode The TransportMode the Person used to get to the new Location.
     * @param[in] cells_new The Cell%s that the Person visits at the new Location.
     * */
    void migrate_to(Location& loc_new, mio::abm::TransportMode transport_mode,
                    const std::vector<uint32_t>& cells = {0});

    /**
     * @brief Get the latest #Infection of the Person.
     * @return The latest #Infection of the Person.
     */
    Infection& get_infection();
    const Infection& get_infection() const;

    /** 
     * @brief Get all Vaccination%s of the Person.
     * @return A vector with all Vaccination%s.
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
     * @brief Returns if the Person is infected at the TimePoint.
     * @param[in] t TimePoint of querry. Usually the current time of the Simulation.
     * @return True if the Person is infected at the TimePoint.
     */
    bool is_infected(TimePoint t) const;

    /**
     * @brief Get the InfectionState of the Person at a specific TimePoint.
     * @param[in] t TimePoint of querry. Usually the current time of the Simulation.
     * @return The InfectionState of the latest Infection at time t.
     */
    InfectionState get_infection_state(TimePoint t) const;

    /**
     * @brief Adds a new Infection to the list of Infection%s.
     * @param[in] inf The new Infection.
     */
    void add_new_infection(Infection&& inf);

    /**
     * @brief Adds a new Vaccination to the list of Vaccination%s that is sorted in vaccination time.
     * @param[in] vacc The new Vaccination.
     */
    void add_new_vaccination(Vaccination&& vacc);

    /**
     * @brief Get the AgeGroup of this Person.
     * @return AgeGroup of the Person.
     */
    AgeGroup get_age() const
    {
        return m_age;
    }

    /**
     * @brief Get the current Location of the Person.
     * @return Current Location of the Person.
     */
    Location& get_location();

    const Location& get_location() const;

    /**
     * @brief Get the time the Person has been at its current Location.
     * @return TimeSpan the Person has been at the Location.
     */
    TimeSpan get_time_at_location() const
    {
        return m_time_at_location;
    }

    /**
     * @brief Get the TimePoint of the last negative test.
     * @return TimePoint since the last test.
     */
    TimePoint get_time_of_last_test() const
    {
        return m_time_of_last_test;
    }
    /**
     * @brief Set an assigned Location of the Person. 
     * The assigned Location is saved by its index of its LocationId. Assume that a Person has at most one assigned
     * Location per #LocationType.
     * @param[in] location The new assigned Location.
     */
    void set_assigned_location(Location& location);

    /**
     * @brief Set an assigned Location of the Person. 
     * The assigned Location is saved by the index of its LocationId. Assume that a Person has at most one assigned
     * Location of a certain #LocationType.
     * @param[in] id The LocationId of the Location.
     */
    void set_assigned_location(LocationId id);

    /**
     * @brief Returns the index of an assigned Location of the Person.
     * Assume that a Person has at most one assigned Location of a certain #LocationType.
     * @param[in] type #LocationType of the assigned Location.
     * @return The index in the LocationId of the assigned Location.
     */
    uint32_t get_assigned_location_index(LocationType type) const;

    /**
     * @brief Get the assigned Location%s of the Person.
     * @return A vector with the indices of the assigned Location%s of the Person.
     */
    const std::vector<uint32_t>& get_assigned_locations() const
    {
        return m_assigned_locations;
    }

    /**
     * @brief Draw if the Person goes to work or is in home office during lockdown at a specific TimePoint.
     * Every Person has a random number. Depending on this number and the time, the Person works from home in case of a
     * lockdown.
     * @param[in] t The TimePoint of interest. Usually the current time of the Simulation.
     * @param[in] params Parameters that describe the migration between Location%s.
     * @return True the Person works from home.
     */
    bool goes_to_work(TimePoint t, const Parameters& params) const;

    /**
     * @brief Draw at what time the Person goes to work.
     * Every Person has a random number to determine what time to go to work.
     * Depending on this number Person decides what time has to go to work.
     * @param[in] params Parameters that describe the migration between Location%s.
     * @return The time of going to work.
     */
    TimeSpan get_go_to_work_time(const Parameters& params) const;

    /**
     * @brief Draw if the Person goes to school or stays at home during lockdown.
     * Every Person has a random number that determines if they go to school in case of a lockdown.
     * @param[in] t The TimePoint of interest. Usually the current time of the Simulation.
     * @param[in] params Parameters that describe the migration between Location%s.
     * @return True if the Person goes to school.
     */
    bool goes_to_school(TimePoint t, const Parameters& params) const;

    /**
     * @brief Draw at what time the Person goes to work.
     * Every Person has a random number to determine what time to go to school.
     * Depending on this number Person decides what time has to go to school.
     * @param[in] params Parameters that describe the migration between Location%s.
     * @return The time of going to school.
     */
    TimeSpan get_go_to_school_time(const Parameters& params) const;

    /**
     * @brief Answers the question if a Person is currently in quarantine.
     * If a Person is in quarantine this Person cannot migrate to Location%s other than Home or the Hospital.
     * @param[in] t The TimePoint of interest. Usually the current time of the Simulation.
     * @param[in] params Parameter that includes the length of a quarantine.
     * @return True if the Person is in quarantine.
     */
    bool is_in_quarantine(TimePoint t, const Parameters& params) const
    {
        return t < m_quarantine_start + params.get<mio::abm::QuarantineDuration>();
    }

    /**
     * @brief Removes the quarantine status of the Person.
     */
    void remove_quarantine();

    /**
     * @brief Simulates a viral test and returns the test result of the Person.
     * If the test is positive, the Person has to quarantine.
     * If the test is negative, quarantine ends.
     * @param[inout] rng RandomNumberGenerator of the Person.
     * @param[in] t TimePoint of the test.
     * @param[in] params Sensitivity and specificity of the test method.
     * @return True if the test result of the Person is positive.
     */
    bool get_tested(RandomNumberGenerator& rng, TimePoint t, const TestParameters& params);

    /**
     * @brief Get the PersonID of the Person.
     * The PersonID should correspond to the index in m_persons in world.
     * @return The PersonID.
     */
    uint32_t get_person_id() const;

    /**
     * @brief Get index of Cell%s of the Person.
     * @return A vector of all Cell indices the Person visits at the current Location.
     */
    std::vector<uint32_t>& get_cells();

    const std::vector<uint32_t>& get_cells() const;

    /**
     * @brief Get the current Mask of the Person.
     * @return Reference to the Mask object of the Person.
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
     * @brief Get the protection of the Mask.
     * A value of 1 represents full protection and a value of 0 means no protection. This depends on the MaskType of the
     * Mask the Person is wearing.
     * @param[in] params The parameters of the Infection that are the same everywhere within the World.
     * @return The reduction factor of getting an Infection when wearing the Mask.
     */
    ScalarType get_mask_protective_factor(const Parameters& params) const;

    /**
     * @brief For every #LocationType a Person has a compliance value between -1 and 1.
     * -1 means that the Person never complies to any Mask duty at the given #LocationType.
     * 1 means that the Person always wears a Mask a the #LocationType even if it is not required.
     * @param[in] preferences The vector of Mask compliance values for all #LocationType%s.
     */
    void set_mask_preferences(std::vector<ScalarType> preferences)
    {
        m_mask_compliance = preferences;
    }

    /**
     * @brief Get the Mask compliance of the Person for the current Location.
     * @param[in] location The current Location of the Person.
     * @return The probability that the Person does not comply to any Mask duty/wears a Mask even if it is not required.
     */
    ScalarType get_mask_compliance(LocationType location) const
    {
        return m_mask_compliance[static_cast<int>(location)];
    }

    /**
     * @brief Checks whether the Person wears a Mask at the target Location.
     * @param[inout] rng RandomNumberGenerator of the Person.
     * @param[in] target The target Location.
     * @return Whether a Person wears a Mask at the Location.
     */
    bool apply_mask_intervention(RandomNumberGenerator& rng, const Location& target);

    /**
     * @brief Decide if a Person is currently wearing a Mask.
     * @param[in] wear_mask If true, the protection of the Mask is considered when computing the exposure rate.
     */
    void set_wear_mask(bool wear_mask)
    {
        m_wears_mask = wear_mask;
    }

    /**
     * @brief Get the information if the Person is currently wearing a Mask.
     * @return True if the Person is currently wearing a Mask.
     */
    bool get_wear_mask() const
    {
        return m_wears_mask;
    }

    /**
     * @brief Get the multiplicative factor on how likely an #Infection is due to the immune system.
     * @param[in] t TimePoint of check.
     * @param[in] virus VirusVariant to check
     * @param[in] params Parameters in the model.
     * @returns Protection factor for general #Infection of the immune system to the given VirusVariant at the given TimePoint.
     */
    ScalarType get_protection_factor(TimePoint t, VirusVariant virus, const Parameters& params) const;

    /**
     * @brief Get the transport mode the Person used to get to its current Location.
     * @return TransportMode the Person used to get to its current Location.
    */
    mio::abm::TransportMode get_last_transport_mode() const
    {
        return m_last_transport_mode;
    }

    /**
     * @brief Sett if a Person should be logged
     */
    void set_should_be_logged(bool should_be_logged)
    {
        m_should_be_logged = should_be_logged;
    }

    /**
     * @brief Get if a Person should be logged
     */
    bool get_should_be_logged() const
    {
        return m_should_be_logged;
    }

    /**
    * Get this persons RandomNumberGenerator counter.
    * @see mio::abm::Person::RandomNumberGenerator.
    */
    Counter<uint32_t>& get_rng_counter()
    {
        return m_rng_counter;
    }

    /**
     * @brief Get the latest #Infection or #Vaccination and its initial TimePoint of the Person. 
    */
    std::pair<ExposureType, TimePoint> get_latest_protection(TimePoint t) const;

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Person");
        obj.add_element("Location", *m_location);
        obj.add_element("age", m_age);
        obj.add_element("id", m_person_id);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Person> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Person");
        auto loc = obj.expect_element("Location", mio::Tag<Location>{});
        auto age = obj.expect_element("age", Tag<uint32_t>{});
        auto id  = obj.expect_element("id", Tag<uint32_t>{});
        return apply(
            io,
            [](auto&& loc_, auto&& age_, auto&& id_) {
                return Person{mio::RandomNumberGenerator(), loc_, AgeGroup(age_), id_};
            },
            loc, age, id);
    }

private:
    observer_ptr<Location> m_location; ///< Current Location of the Person.
    std::vector<uint32_t> m_assigned_locations; /**! Vector with the indices of the assigned Locations so that the 
    Person always visits the same Home or School etc. */
    std::vector<Vaccination> m_vaccinations; ///< Vector with all Vaccination%s the Person has received.
    std::vector<Infection> m_infections; ///< Vector with all Infection%s the Person had.
    TimePoint m_quarantine_start; ///< TimePoint when the Person started quarantine.
    AgeGroup m_age; ///< AgeGroup the Person belongs to.
    TimeSpan m_time_at_location; ///< Time the Person has spent at its current Location so far.
    double m_random_workgroup; ///< Value to determine if the Person goes to work or works from home during lockdown.
    double m_random_schoolgroup; ///< Value to determine if the Person goes to school or stays at home during lockdown.
    double m_random_goto_work_hour; ///< Value to determine at what time the Person goes to work.
    double m_random_goto_school_hour; ///< Value to determine at what time the Person goes to school.
    TimePoint m_time_of_last_test; ///< TimePoint of the last negative test.
    Mask m_mask; ///< The Mask of the Person.
    bool m_wears_mask = false; ///< Whether the Person currently wears a Mask.
    std::vector<ScalarType> m_mask_compliance; ///< Vector of Mask compliance values for all #LocationType%s.
    uint32_t m_person_id; ///< Id of the Person.
    std::vector<uint32_t> m_cells; ///< Vector with all Cell%s the Person visits at its current Location.
    mio::abm::TransportMode m_last_transport_mode; ///< TransportMode the Person used to get to its current Location.
    Counter<uint32_t> m_rng_counter{0}; ///< counter for RandomNumberGenerator.
    bool m_should_be_logged = true; ///< Flag to determine if the Person should be logged.
    bool m_tested_positive = false;
};

} // namespace abm
} // namespace mio

#endif
