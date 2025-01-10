/* 
* Copyright (C) 2020-2025 MEmilio
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
#ifndef MIO_ABM_PERSON_H
#define MIO_ABM_PERSON_H

#include "abm/infection.h"
#include "abm/infection_state.h"
#include "abm/location_id.h"
#include "abm/location.h"
#include "abm/location_type.h"
#include "abm/parameters.h"
#include "abm/person_id.h"
#include "abm/personal_rng.h"
#include "memilio/io/default_serialize.h"
#include "abm/time.h"
#include "abm/test_type.h"
#include "abm/protection_event.h"
#include "abm/intervention_type.h"
#include "abm/mask.h"
#include "abm/mobility_data.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/random_number_generator.h"

namespace mio
{
namespace abm
{

static constexpr uint32_t INVALID_PERSON_ID = std::numeric_limits<uint32_t>::max();

/**
 * @brief Agents in the simulated Model that can carry and spread the Infection.
 */
class Person
{
public:
    /**
     * @brief Create a Person.
     * @param[in, out] rng RandomNumberGenerator.
     * @param[in, out] location Initial Location of the Person.
     * @param[in] age The AgeGroup of the Person.
     * @param[in] person_id Index of the Person.
     */
    explicit Person(mio::RandomNumberGenerator& rng, LocationType location_type, LocationId location_id, AgeGroup age,
                    PersonId person_id = PersonId::invalid_id());

    explicit Person(const Person& other, PersonId id);

    /**
     * @brief Compare two Person%s.
     */
    bool operator==(const Person& other) const
    {
        return (m_person_id == other.m_person_id);
    }

    /**
     * @brief Get the latest #Infection of the Person.
     * @return The latest #Infection of the Person.
     */
    Infection& get_infection();
    const Infection& get_infection() const;

    /**
     * @brief Get all vaccinations of the Person.
     * @return A vector with all vaccinations.
     * @{
     */
    std::vector<ProtectionEvent>& get_vaccinations()
    {
        return m_vaccinations;
    }

    const std::vector<ProtectionEvent>& get_vaccinations() const
    {
        return m_vaccinations;
    }
    /** @} */

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
    void add_new_infection(Infection&& inf, TimePoint current_time = TimePoint(0));

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
    LocationId get_location() const;

    LocationType get_location_type() const
    {
        return m_location_type;
    }

    /**
     * @brief Change the location of the person.
     * @param[in] id The new location.
     */
    void set_location(LocationType type, LocationId id);

    /**
     * @brief Get the time the Person has been at its current Location.
     * @return TimeSpan the Person has been at the Location.
     */
    TimeSpan get_time_at_location() const
    {
        return m_time_at_location;
    }

    /**
     * @brief Add to the time the Person has been at its current Location.
     * @param[in] dt TimeSpan the Person has spent at the Location.
     */
    void add_time_at_location(const TimeSpan dt)
    {
        m_time_at_location += dt;
    }

    /**
     * @brief Set an assigned Location of the Person.
     *
     * Important: Setting incorrect values will cause issues during simulation. It is preferable to use
     *            Model::assign_location with a valid LocationId, obtained e.g. through Model::add_location.
     *
     * The assigned Location is saved by the index of its LocationId. Assume that a Person has at most one assigned
     * Location of a certain #LocationType.
     * @param[in] type The LocationType of the Location.
     * @param[in] id The LocationId of the Location.
     */
    void set_assigned_location(LocationType type, LocationId id);

    /**
     * @brief Returns the index of an assigned Location of the Person.
     * Assume that a Person has at most one assigned Location of a certain #LocationType.
     * @param[in] type #LocationType of the assigned Location.
     * @return The index in the LocationId of the assigned Location.
     */
    LocationId get_assigned_location(LocationType type) const;

    /**
     * @brief Get the assigned Location%s of the Person.
     * @return A vector with the indices of the assigned Location%s of the Person.
     */
    const std::vector<LocationId>& get_assigned_locations() const
    {
        return m_assigned_locations;
    }

    /**
     * @brief Draw if the Person goes to work or is in home office during lockdown at a specific TimePoint.
     * Every Person has a random number. Depending on this number and the time, the Person works from home in case of a
     * lockdown.
     * @param[in] t The TimePoint of interest. Usually the current time of the Simulation.
     * @param[in] params Parameters that describe the mobility between Location%s.
     * @return True the Person works from home.
     */
    bool goes_to_work(TimePoint t, const Parameters& params) const;

    /**
     * @brief Draw at what time the Person goes to work.
     * Every Person has a random number to determine what time to go to work.
     * Depending on this number Person decides what time has to go to work.
     * @param[in] params Parameters that describe the mobility between Location%s.
     * @return The time of going to work.
     */
    TimeSpan get_go_to_work_time(const Parameters& params) const;

    /**
     * @brief Draw at what time the Person returns from work.
     * Every Person has a random number to determine what time to return from work.
     * Depending on this number Person decides what time it has to return from work.
     * @param[in] params Parameters that describe the migration between Location%s.
     * @return The time of returning from work.
     */
    TimeSpan get_return_from_work_time(const Parameters& params) const;

    /**
     * @brief Draw if the Person goes to school or stays at home during lockdown.
     * Every Person has a random number that determines if they go to school in case of a lockdown.
     * @param[in] t The TimePoint of interest. Usually the current time of the Simulation.
     * @param[in] params Parameters that describe the mobility between Location%s.
     * @return True if the Person goes to school.
     */
    bool goes_to_school(TimePoint t, const Parameters& params) const;

    /**
     * @brief Draw at what time the Person goes to school.
     * Every Person has a random number to determine what time to go to school.
     * Depending on this number Person decides what time has to go to school.
     * @param[in] params Parameters that describe the mobility between Location%s.
     * @return The time of going to school.
     */
    TimeSpan get_go_to_school_time(const Parameters& params) const;

    /**
     * @brief Draw at what time the Person return from school.
     * Every Person has a random number to determine what time to return from school.
     * Depending on this number Person decides what time it has to return from school.
     * @param[in] params Parameters that describe the migration between Location%s.
     * @return The time of returning from school.
     */
    TimeSpan get_return_from_school_time(const Parameters& params) const;

    /**
     * @brief Answers the question if a Person is currently in quarantine.
     * If a Person is in quarantine this Person cannot change to Location%s other than Home or the Hospital.
     * @param[in] t The TimePoint of interest. Usually the current time of the Simulation.
     * @param[in] params Parameter that includes the length of a quarantine.
     * @return True if the Person is in quarantine.
     */
    bool is_in_quarantine(TimePoint t, const Parameters& params) const
    {
        return t < m_home_isolation_start + params.get<mio::abm::QuarantineDuration>();
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
    bool get_tested(PersonalRandomNumberGenerator& rng, TimePoint t, const TestParameters& params);

    /**
     * @brief Get the PersonId of the Person.
     * The PersonId should correspond to the index in m_persons in the Model.
     * @return The PersonId.
     */
    PersonId get_id() const;

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
     * @param[in] params The parameters of the Infection that are the same everywhere within the Model.
     * @return The reduction factor of getting an Infection when wearing the Mask.
     */
    ScalarType get_mask_protective_factor(const Parameters& params) const;

    /**
     * @brief For every #InterventionType a Person has a compliance value between 0 and 1.
     * 0 means that the Person never complies to the Intervention.
     * 1 means that the Person always complies to the Intervention.
     * @param[in] intervention_type The #InterventionType.
     * @param[in] value The compliance value.
     */
    void set_compliance(InterventionType intervention_type, ScalarType value)
    {
        m_compliance[static_cast<uint32_t>(intervention_type)] = value;
    }

    /**
     * @brief Get the compliance of the Person for an Intervention.
     * @param[in] intervention_type The #InterventionType.
     * @return The probability that the Person complies to an Intervention.
     */
    ScalarType get_compliance(InterventionType intervention_type) const
    {
        return m_compliance[static_cast<uint32_t>(intervention_type)];
    }

    /**
     * @brief Checks whether the Person complies an Intervention.
     * @param[inout] rng PersonalRandomNumberGenerator of the Person.
     * @param[in] intervention The #InterventionType.
     * @return Checks whether the Person complies an Intervention.
     */
    bool is_compliant(PersonalRandomNumberGenerator& rng, InterventionType intervention) const;

    /**
     * @brief Change the mask to new type.
     * @param[in] type The required #MaskType.
     * @param[in] t The TimePoint of mask change.
     */
    void set_mask(MaskType type, TimePoint t);

    /**
     * @brief Get the multiplicative factor on how likely an #Infection is due to the immune system.
     * @param[in] t TimePoint of check.
     * @param[in] virus VirusVariant to check
     * @param[in] params Parameters in the model.
     * @returns Protection factor for general #Infection of the immune system to the given VirusVariant at the given TimePoint.
     */
    ScalarType get_protection_factor(TimePoint t, VirusVariant virus, const Parameters& params) const;

    /**
     * @brief Add a new vaccination
     * @param[in] v ProtectionType (i. e. vaccine) the person takes.
     * @param[in] t TimePoint of the vaccination.
     */
    void add_new_vaccination(ProtectionType v, TimePoint t)
    {
        m_vaccinations.push_back(ProtectionEvent(v, t));
    }

    /**
     * @brief Get the transport mode the Person used to get to its current Location.
     * @return TransportMode the Person used to get to its current Location.
     */
    mio::abm::TransportMode get_last_transport_mode() const
    {
        return m_last_transport_mode;
    }

    /**
     * @brief Set the transport mode the Person used to get to its current Location.
     * @param[in] mode TransportMode the Person used to get to its current Location.
     */
    void set_last_transport_mode(const mio::abm::TransportMode mode)
    {
        m_last_transport_mode = mode;
    }

    /**
     * @brief Get this persons RandomNumberGenerator counter.
     * @see mio::abm::PersonalRandomNumberGenerator.
     */
    Counter<uint32_t>& get_rng_counter()
    {
        return m_rng_counter;
    }

    /**
     * @brief Get the latest #ProtectionType and its initial TimePoint of the Person.
     */
    ProtectionEvent get_latest_protection() const;

    /// This method is used by the default serialization feature.
    auto default_serialize()
    {
        return Members("Person")
            .add("location", m_location)
            .add("location_type", m_location_type)
            .add("assigned_locations", m_assigned_locations)
            .add("vaccinations", m_vaccinations)
            .add("infections", m_infections)
            .add("home_isolation_start", m_home_isolation_start)
            .add("age_group", m_age)
            .add("time_at_location", m_time_at_location)
            .add("rnd_workgroup", m_random_workgroup)
            .add("rnd_schoolgroup", m_random_schoolgroup)
            .add("rnd_go_to_work_hour", m_random_goto_work_hour)
            .add("rnd_go_to_school_hour", m_random_goto_school_hour)
            .add("mask", m_mask)
            .add("compliance", m_compliance)
            .add("id", m_person_id)
            .add("cells", m_cells)
            .add("last_transport_mode", m_last_transport_mode)
            .add("rng_counter", m_rng_counter)
            .add("test_results", m_test_results);
    }

    /**
     * @brief Add TestResult to the Person
     * @param[in] t The TimePoint of the test.
     * @param[in] type The TestType of the test.
     * @param[in] result The result of the test.
    */
    void add_test_result(TimePoint t, TestType type, bool result);

    /**
     * @brief Get the most recent TestResult performed from the Person based on the TestType.
     * If time_of_testing == TimePoint(std::numeric_limits<int>::min()), there is no previous TestResult.
     * @param[in] type The TestType of the test.
     * @return The latest TestResult of the given Type.
    */
    TestResult get_test_result(TestType type) const;

    TimeSpan get_time_since_transmission() const
    {
        return m_time_since_transmission;
    };

private:
    void change_time_since_transmission(const TimeSpan dt, TimePoint t);

    LocationId m_location; ///< Current Location of the Person.
    LocationType m_location_type; ///< Type of the current Location.
    std::vector<LocationId> m_assigned_locations; /**! Vector with the indices of the assigned Locations so that the
    Person always visits the same Home or School etc. */
    std::vector<ProtectionEvent> m_vaccinations; ///< Vector with all vaccinations the Person has received.
    std::vector<Infection> m_infections; ///< Vector with all Infection%s the Person had.
    TimeSpan m_time_since_transmission;
    TimePoint m_home_isolation_start; ///< TimePoint when the Person started isolation at home.
    AgeGroup m_age; ///< AgeGroup the Person belongs to.
    TimeSpan m_time_at_location; ///< Time the Person has spent at its current Location so far.
    double m_random_workgroup; ///< Value to determine if the Person goes to work or works from home during lockdown.
    double m_random_schoolgroup; ///< Value to determine if the Person goes to school or stays at home during lockdown.
    double m_random_goto_work_hour; ///< Value to determine at what time the Person goes to work.
    double m_random_return_from_work_hour; ///< Value to determine at what time the Person returns from work.
    double m_random_goto_school_hour; ///< Value to determine at what time the Person goes to school.
    double m_random_return_from_school_hour; ///< Value to determine at what time the Person returns from school.
    Mask m_mask; ///< The Mask of the Person.
    std::vector<ScalarType>
        m_compliance; ///< Vector of compliance values for all #InterventionType%s. Values from 0 to 1.
    PersonId m_person_id; ///< Id of the Person.
    std::vector<uint32_t> m_cells; ///< Vector with all Cell%s the Person visits at its current Location.
    mio::abm::TransportMode m_last_transport_mode; ///< TransportMode the Person used to get to its current Location.
    Counter<uint32_t> m_rng_counter{0}; ///< counter for RandomNumberGenerator.
    CustomIndexArray<TestResult, TestType> m_test_results; ///< CustomIndexArray for TestResults.
};

} // namespace abm

/// @brief Creates an instance of abm::Person for default serialization.
template <>
struct DefaultFactory<abm::Person> {
    static abm::Person create()
    {
        return abm::Person(thread_local_rng(), abm::LocationType::Count, abm::LocationId(), AgeGroup(0),
                           abm::PersonId());
    }
};

} // namespace mio

#endif
