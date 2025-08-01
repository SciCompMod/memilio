/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Majid Abedi, Elisabeth Kluth, David Kerkmann, Sascha Korf, Martin J. Kuehn, Khoa Nguyen
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
#ifndef MIO_ABM_MODEL_H
#define MIO_ABM_MODEL_H

#include "abm/model_functions.h"
#include "abm/location_type.h"
#include "abm/mobility_data.h"
#include "abm/parameters.h"
#include "abm/location.h"
#include "abm/person.h"
#include "abm/person_id.h"
#include "abm/time.h"
#include "abm/trip_list.h"
#include "abm/random_events.h"
#include "abm/testing_strategy.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/stl_util.h"

#include <bitset>
#include <cstdint>
#include <vector>

namespace mio
{
namespace abm
{

/**
 * @brief The Model of the Simulation.
 * It consists of Location%s and Person%s (Agents).
 */
class Model
{
public:
    using LocationIterator        = std::vector<Location>::iterator;
    using ConstLocationIterator   = std::vector<Location>::const_iterator;
    using PersonIterator          = std::vector<Person>::iterator;
    using ConstPersonIterator     = std::vector<Person>::const_iterator;
    using ActivenessIterator      = std::vector<bool>::iterator;
    using ConstActivenessIterator = std::vector<bool>::const_iterator;
    using MobilityRuleType        = LocationType (*)(PersonalRandomNumberGenerator&, const Person&, TimePoint, TimeSpan,
                                              const Parameters&);

    /**
     * @brief Create a Model.
     * @param[in] num_agegroups The number of AgeGroup%s in the simulated Model. Must be less than MAX_NUM_AGE_GROUPS.
     */
    Model(size_t num_agegroups, int id = 0)
        : parameters(num_agegroups)
        , m_id(id)
        , m_trip_list()
        , m_use_mobility_rules(true)
        , m_cemetery_id(add_location(LocationType::Cemetery))
        , m_person_ids_equal_index(true)
    {
        assert(num_agegroups < MAX_NUM_AGE_GROUPS && "MAX_NUM_AGE_GROUPS exceeded.");
    }

    /**
     * @brief Create a Model.
     * @param[in] params Initial simulation parameters.
     */
    Model(const Parameters& params, int id = 0)
        : parameters(params.get_num_groups())
        , m_id(id)
        , m_trip_list()
        , m_use_mobility_rules(true)
        , m_cemetery_id(add_location(LocationType::Cemetery))
        , m_person_ids_equal_index(true)
    {
        parameters = params;
    }

    Model(const Model& other, int id = 0)
        : parameters(other.parameters)
        , m_local_population_cache()
        , m_air_exposure_rates_cache()
        , m_contact_exposure_rates_cache()
        , m_is_local_population_cache_valid(false)
        , m_are_exposure_caches_valid(false)
        , m_exposure_caches_need_rebuild(true)
        , m_id(id)
        , m_persons(other.m_persons)
        , m_locations(other.m_locations)
        , m_activeness_statuses(other.m_activeness_statuses)
        , m_has_locations(other.m_has_locations)
        , m_testing_strategy(other.m_testing_strategy)
        , m_trip_list(other.m_trip_list)
        , m_use_mobility_rules(other.m_use_mobility_rules)
        , m_mobility_rules(other.m_mobility_rules)
        , m_cemetery_id(other.m_cemetery_id)
        , m_rng(other.m_rng)
        , m_person_ids_equal_index(true)
    {
    }
    Model& operator=(const Model&) = default;
    Model(Model&&)                 = default;
    Model& operator=(Model&&)      = default;

    /**
     * serialize this.
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Model");
        obj.add_element("parameters", parameters);
        // skip caches, they are rebuild by the deserialized model
        obj.add_list("persons", get_persons().begin(), get_persons().end());
        obj.add_list("locations", get_locations().begin(), get_locations().end());
        obj.add_element("location_types", m_has_locations.to_ulong());
        obj.add_element("testing_strategy", m_testing_strategy);
        obj.add_element("trip_list", m_trip_list);
        obj.add_element("use_mobility_rules", m_use_mobility_rules);
        obj.add_element("cemetery_id", m_cemetery_id);
        obj.add_element("rng", m_rng);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Model> deserialize(IOContext& io)
    {
        auto obj                = io.expect_object("Model");
        auto params             = obj.expect_element("parameters", Tag<Parameters>{});
        auto persons            = obj.expect_list("persons", Tag<Person>{});
        auto locations          = obj.expect_list("locations", Tag<Location>{});
        auto location_types     = obj.expect_element("location_types", Tag<unsigned long>{});
        auto trip_list          = obj.expect_element("trip_list", Tag<TripList>{});
        auto use_mobility_rules = obj.expect_element("use_mobility_rules", Tag<bool>{});
        auto cemetery_id        = obj.expect_element("cemetery_id", Tag<LocationId>{});
        auto rng                = obj.expect_element("rng", Tag<RandomNumberGenerator>{});
        return apply(
            io,
            [](auto&& params_, auto&& persons_, auto&& locations_, auto&& location_types_, auto&& trip_list_,
               auto&& use_mobility_rules_, auto&& cemetery_id_, auto&& rng_) {
                Model model{params_};
                model.m_persons.assign(persons_.cbegin(), persons_.cend());
                model.m_locations.assign(locations_.cbegin(), locations_.cend());
                model.m_has_locations      = location_types_;
                model.m_trip_list          = trip_list_;
                model.m_use_mobility_rules = use_mobility_rules_;
                model.m_cemetery_id        = cemetery_id_;
                model.m_rng                = rng_;
                return model;
            },
            params, persons, locations, location_types, trip_list, use_mobility_rules, cemetery_id, rng);
    }

    /** 
     * @brief Prepare the Model for the next Simulation step.
     * @param[in] t Current time.
     * @param[in] dt Length of the time step.
     */
    void begin_step(TimePoint t, TimeSpan dt);

    /** 
     * @brief Evolve the Model one time step.
     * @param[in] t Current time.
     * @param[in] dt Length of the time step.
     */
    void evolve(TimePoint t, TimeSpan dt);

    /** 
     * @brief Add a Location to the Model.
     * @param[in] type Type of Location to add.
     * @param[in] num_cells [Default: 1] Number of Cell%s that the Location is divided into.
     * @return ID of the newly created Location.
     */
    LocationId add_location(LocationType type, uint32_t num_cells = 1);

    /** 
     * @brief Add a Person to the Model.
     * @param[in] id The LocationID of the initial Location of the Person.
     * @param[in] age AgeGroup of the person.
     * @return Id of the newly created Person.
     */
    PersonId add_person(const LocationId id, AgeGroup age);

    /**
     * @brief Adds a copy of a given Person to the Model.
     * @param[in] person The Person to copy from. 
     * @return Id of the newly created Person.
     */
    PersonId add_person(Person&& person);

    /**
     * @brief Get a range of all Location%s in the Model.
     * @return A range of all Location%s.
     * @{
     */
    Range<std::pair<ConstLocationIterator, ConstLocationIterator>> get_locations() const;
    Range<std::pair<LocationIterator, LocationIterator>> get_locations();
    /** @} */

    /**
     * @brief Get a range of all Person%s in the Model.
     * @return A range of all Person%s.
     * @{
     */
    Range<std::pair<ConstPersonIterator, ConstPersonIterator>> get_persons() const;
    Range<std::pair<PersonIterator, PersonIterator>> get_persons();
    /** @} */

    /**
     * @brief Get a range of all Person%s activeness statuses in the Model.
     * @return A range of all Person%s activeness statuses.
     * @{
     */
    Range<std::pair<ConstActivenessIterator, ConstActivenessIterator>> get_activeness_statuses() const;
    Range<std::pair<ActivenessIterator, ActivenessIterator>> get_activeness_statuses();
    /** @} */

    /**
     * @brief Find an assigned Location of a Person.
     * @param[in] type The #LocationType that specifies the assigned Location.
     * @param[in] person Id of the Person.
     * @return ID of the Location of LocationType type assigend to person.
     */
    LocationId find_location(LocationType type, const PersonId person) const;

    /**
     * @brief Assign a Location to a Person.
     * A Person can have at most one assigned Location of a certain LocationType.
     * Assigning another Location of an already assigned LocationType will replace the prior assignment.  
     * @param[in] person The Id of the person this location will be assigned to.
     * @param[in] location The LocationId of the Location.
     */
    void assign_location(PersonId person, LocationId location)
    {
        assign_location(get_person(person), location);
    }

    /**
     * @brief Get the number of Persons in one #InfectionState at all Location%s.
     * @param[in] t Specified #TimePoint.
     * @param[in] s Specified #InfectionState.
     */
    size_t get_subpopulation_combined(TimePoint t, InfectionState s) const;

    /**
     * @brief Get the number of Persons in one #InfectionState at all Location%s of a type.
     * @param[in] t Specified #TimePoint.
     * @param[in] s Specified #InfectionState.
     * @param[in] type Specified #LocationType.
     */
    size_t get_subpopulation_combined_per_location_type(TimePoint t, InfectionState s, LocationType type) const;

    /**
     * @brief Get the mobility data.
     * @return Reference to the list of Trip%s that the Person%s make.
     */
    TripList& get_trip_list();

    const TripList& get_trip_list() const;

    /**
     * @brief Decide if mobility rules (like go to school/work) are used or not;
     * The mobility rules regarding hospitalization/ICU/quarantine are always used.
     * @param[in] param If true uses the mobility rules for changing location to school/work etc., else only the rules
     * regarding hospitalization/ICU/quarantine.
     */
    void use_mobility_rules(bool param);
    bool use_mobility_rules() const;

    /**
    * @brief Check if at least one Location with a specified LocationType exists.
    * @return True if there is at least one Location of LocationType `type`. False otherwise.
    */
    bool has_location(LocationType type) const
    {
        return m_has_locations[size_t(type)];
    }

    /**
    * @brief Check if at least one Location of every specified LocationType exists.
    * @tparam C A type of container of LocationType.
    * @param location_types A container of LocationType%s.
    * @return True if there is at least one Location of every LocationType in `location_types`. False otherwise.
    */
    template <class C = std::initializer_list<LocationType>>
    bool has_locations(const C& location_types) const
    {
        return std::all_of(location_types.begin(), location_types.end(), [&](auto loc) {
            return has_location(loc);
        });
    }

    /**
     * @brief Get the TestingStrategy.
     * @return Reference to the list of TestingScheme%s that are checked for testing.
     */
    TestingStrategy& get_testing_strategy();

    const TestingStrategy& get_testing_strategy() const;

    /** 
     * @brief The simulation parameters of the Model.
     */
    Parameters parameters;

    /**
    * Get the RandomNumberGenerator used by this Model for random events.
    * Persons use their own generators with the same key as the global one. 
    * @return The random number generator.
    */
    RandomNumberGenerator& get_rng()
    {
        return m_rng;
    }

    /**
     * Get the model id. Is only relevant for graph abm or hybrid model.
     * @return The model id
     */
    int get_id() const
    {
        return m_id;
    }

    /**
     * @brief Get a reference to a Person from this Model.
     * @param[in] person_id A Person's PersonId.
     * @return A reference to the Person.
     */
    const Person& get_person(PersonId person_id) const
    {
        return get_person_impl(*this, person_id);
    }

    Person& get_person(PersonId person_id)
    {
        return get_person_impl(*this, person_id);
    }

    /**
     * @brief Get the number of Person%s of a particular #InfectionState for all Cell%s.
     * @param[in] location A LocationId from the Model.
     * @param[in] t TimePoint of querry.
     * @param[in] state #InfectionState of interest.
     * @return Amount of Person%s of the #InfectionState in all Cell%s of the Location.
     */
    size_t get_subpopulation(LocationId location, TimePoint t, InfectionState state) const
    {
        return std::count_if(m_persons.begin(), m_persons.end(), [&](auto&& p) {
            return p.get_location_model_id() == m_id && p.get_location() == location &&
                   p.get_infection_state(t) == state;
        });
    }

    /**
     * @brief Get the total number of Person%s at the Location.
     * @param[in] location A LocationId from the Model.
     * @return Number of Person%s in the location.
     */
    size_t get_number_persons(LocationId location) const
    {
        if (!m_is_local_population_cache_valid) {
            build_compute_local_population_cache();
        }
        return m_local_population_cache[location.get()];
    }

    // Change the Location of a Person. this requires that Location is part of this Model.
    /**
     * @brief Let a Person change to another Location.
     * @param[in] person Id of a person from this Model.
     * @param[in] destination LocationId of the Location in this Model, which the Person should change to.
     * @param[in] mode The transport mode the person uses to change the Location.
     * @param[in] cells The cells within the destination the person should be in.
     */
    inline void change_location(PersonId person, LocationId destination, TransportMode mode = TransportMode::Unknown,
                                const std::vector<uint32_t>& cells = {0})
    {
        change_location(get_person(person), destination, mode, cells);
    }

    /**
     * @brief Let a person interact with the population at its current location.
     * @param[in] person Id of a person from this Model.
     * @param[in] t Time step of the simulation.
     * @param[in] dt Step size of the simulation.
     */
    inline void interact(PersonId person, TimePoint t, TimeSpan dt)
    {
        interact(get_person(person), t, dt);
    }

    /**
     * @brief Get a reference to a location in this Model.
     * @param[in] id LocationId of the Location.
     * @return Reference to the Location.
     * @{
     */
    const Location& get_location(LocationId id) const
    {
        assert(id != LocationId::invalid_id() && "Given LocationId must be valid.");
        assert(id < LocationId((uint32_t)m_locations.size()) && "Given LocationId is not in this Model.");
        return m_locations[id.get()];
    }

    /**
     * @brief Assign a Location to a Person.
     * A Person can have at most one assigned Location of a certain LocationType.
     * Assigning another Location of an already assigned LocationType will replace the prior assignment.  
     * @param[in] person reference to the Person the location will be assigned to.
     * @param[in] location The LocationId of the Location.
     */
    void assign_location(Person& person, LocationId location)
    {
        person.set_assigned_location(get_location(location).get_type(), location, m_id);
    }

    Location& get_location(LocationId id)
    {
        assert(id != LocationId::invalid_id() && "Given LocationId must be valid.");
        assert(id < LocationId((uint32_t)m_locations.size()) && "Given LocationId is not in this Model.");
        return m_locations[id.get()];
    }
    /** @} */

    /**
     * @brief Get a reference to the location of a person.
     * @param[in] id Id of a person.
     * @return Reference to the Location.
     * @{
     */
    inline Location& get_location(PersonId id)
    {
        return get_location(get_person(id).get_location());
    }

    inline const Location& get_location(PersonId id) const
    {
        return get_location(get_person(id).get_location());
    }
    /** @} */

    /**
     * @brief Get index of person in m_persons.
     * @param[in] person_id A person's unique PersonId. 
     * First 32 bit are the Person's individual id and second 32 bit the Persons's home model id. 
     * @return Index of Person in m_persons vector.
     * @{
     */
    uint32_t get_person_index(PersonId person_id) const
    {
        mio::log_debug("get_person_index is used leading to a search in m_persons.");
        auto it = std::find_if(m_persons.begin(), m_persons.end(), [person_id](auto& person) {
            return person.get_id() == person_id;
        });
        if (it == m_persons.end()) {
            log_error("Given PersonId is not in this Model.");
            return std::numeric_limits<uint32_t>::max();
        }
        else {
            return static_cast<uint32_t>(std::distance(m_persons.begin(), it));
        }
    }

protected:
    /**
     * @brief Person%s interact at their Location and may become infected.
     * @param[in] t The current TimePoint.
     * @param[in] dt The length of the time step of the Simulation.
     */
    void interaction(TimePoint t, TimeSpan dt);
    /**
     * @brief Person%s change location in the Model according to rules.
     * @param[in] t The current TimePoint.
     * @param[in] dt The length of the time step of the Simulation.
     */
    void perform_mobility(TimePoint t, TimeSpan dt);

    /// @brief Shape the cache and store how many Person%s are at any Location. Use from single thread!
    void build_compute_local_population_cache() const;

    /// @brief Shape the air and contact exposure cache according to the current Location%s.
    void build_exposure_caches();

    /**
     * @brief Store all air/contact exposures for the current simulation step.
     * @param[in] t Current TimePoint of the simulation.
     * @param[in] dt The duration of the simulation step.
     */
    void compute_exposure_caches(TimePoint t, TimeSpan dt);

    // Change the Location of a Person. this requires that Location is part of this Model.
    /**
     * @brief Let a Person change to another Location.
     * @param[in] person Reference to Person.
     * @param[in] destination LocationId of the Location in this Model, which the Person should change to.
     * @param[in] mode The transport mode the person uses to change the Location.
     * @param[in] cells The cells within the destination the person should be in.
     */
    inline void change_location(Person& person, LocationId destination, TransportMode mode = TransportMode::Unknown,
                                const std::vector<uint32_t>& cells = {0})
    {
        LocationId origin               = get_location(person).get_id();
        const bool has_changed_location = mio::abm::change_location(person, get_location(destination), mode, cells);
        // if the person has changed location, invalidate exposure caches but keep population caches valid
        if (has_changed_location) {
            m_are_exposure_caches_valid = false;
            if (m_is_local_population_cache_valid) {
                --m_local_population_cache[origin.get()];
                ++m_local_population_cache[destination.get()];
            }
        }
    }

    /**
     * @brief Get a reference to the location of a person.
     * @param[in] person Reference to a Person.
     * @return Reference to the Location.
     * @{
     */
    inline Location& get_location(Person& person)
    {
        return get_location(person.get_location());
    }

    inline const Location& get_location(Person& person) const
    {
        return get_location(person.get_location());
    }

    /**
     * @brief Find an assigned Location of a Person.
     * @param[in] type The #LocationType that specifies the assigned Location.
     * @param[in] person Reference to Person.
     * @return ID of the Location of LocationType type assigend to person.
     */
    LocationId find_location(LocationType type, const Person& person) const
    {
        auto location_id = person.get_assigned_location(type);
        assert(location_id != LocationId::invalid_id() && "The person has no assigned location of that type.");
        return location_id;
    }

    /**
     * @brief Let a person interact with the population at its current location.
     * @param[in] person Reference to Person.
     * @param[in] t Time step of the simulation.
     * @param[in] dt Step size of the simulation.
     */
    inline void interact(Person& person, TimePoint t, TimeSpan dt)
    {
        if (!m_are_exposure_caches_valid) {
            // checking caches is only needed for external calls
            // during simulation (i.e. in evolve()), the caches are computed in begin_step
            compute_exposure_caches(t, dt);
            m_are_exposure_caches_valid = true;
        }
        auto personal_rng = PersonalRandomNumberGenerator(person);
        mio::abm::interact(personal_rng, person, get_location(person.get_location()),
                           m_air_exposure_rates_cache[person.get_location().get()],
                           m_contact_exposure_rates_cache[person.get_location().get()], t, dt, parameters);
    }

    /**
     * @brief Implementation of Model::get_person. 
     * This function needs to use a template to deduce whether the model and returned person should be const.
     * @param[in] m A reference to `*this`, so we can access m_persons.
     * @param[in] person_id A Person's PersonId.
     * @return A reference to the Person with matching ID.
     */
    template <class M>
    static std::conditional_t<std::is_const_v<M>, const Person&, Person&> get_person_impl(M& m, PersonId person_id)
    {
        if (m.m_person_ids_equal_index) {
            assert(static_cast<uint32_t>(person_id.get()) < m.m_persons.size() &&
                   "Given PersonId is not in this Model.");
            return m.m_persons[static_cast<uint32_t>(person_id.get())];
        }
        else {
            mio::log_warning("get_person is accessed by PersonId which does not align with the index of the person due "
                             "to former removal of persons. Therefore m_persons is searched.");
            auto it = std::find_if(m.m_persons.begin(), m.m_persons.end(), [person_id](auto& person) {
                return person.get_id() == person_id;
            });
            assert(it != m.m_persons.end() && "Given PersonId is not in this Model.");
            return *it;
        };
    }

    mutable Eigen::Matrix<std::atomic_int_fast32_t, Eigen::Dynamic, 1>
        m_local_population_cache; ///< Current number of Persons in a given location.
    Eigen::Matrix<AirExposureRates, Eigen::Dynamic, 1>
        m_air_exposure_rates_cache; ///< Cache for local exposure through droplets in #transmissions/day.
    Eigen::Matrix<ContactExposureRates, Eigen::Dynamic, 1>
        m_contact_exposure_rates_cache; ///< Cache for local exposure through contacts in #transmissions/day.
    bool m_is_local_population_cache_valid = false;
    bool m_are_exposure_caches_valid       = false;
    bool m_exposure_caches_need_rebuild    = true;

    int m_id; ///< Model id. Is only used for abm graph model or hybrid model.
    std::vector<Person> m_persons; ///< Vector of every Person.
    std::vector<Location> m_locations; ///< Vector of every Location.
    std::vector<bool>
        m_activeness_statuses; ///< Vector with activeness status for every person. Is only used for abm graph model or hybrid model.
    std::bitset<size_t(LocationType::Count)>
        m_has_locations; ///< Flags for each LocationType, set if a Location of that type exists.
    TestingStrategy m_testing_strategy; ///< List of TestingScheme%s that are checked for testing.
    TripList m_trip_list; ///< List of all Trip%s the Person%s do.
    bool m_use_mobility_rules; ///< Whether mobility rules are considered.
    std::vector<MobilityRuleType> m_mobility_rules; ///< Rules that govern the mobility between Location%s.
    LocationId m_cemetery_id; // Central cemetery for all dead persons.
    RandomNumberGenerator m_rng; ///< Global random number generator
    bool m_person_ids_equal_index;
};

} // namespace abm
} // namespace mio

#endif
