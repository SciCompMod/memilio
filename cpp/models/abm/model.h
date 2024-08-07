/* 
* Copyright (C) 2020-2024 MEmilio
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
    using LocationIterator      = std::vector<Location>::iterator;
    using ConstLocationIterator = std::vector<Location>::const_iterator;
    using PersonIterator        = std::vector<Person>::iterator;
    using ConstPersonIterator   = std::vector<Person>::const_iterator;

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
        obj.add_element("num_agegroups", parameters.get_num_groups());
        std::vector<Trip> trips;
        TripList trip_list = get_trip_list();
        for (size_t i = 0; i < trip_list.num_trips(false); i++) {
            trips.push_back(trip_list.get_next_trip(false));
            trip_list.increase_index();
        }
        trip_list.reset_index();
        for (size_t i = 0; i < trip_list.num_trips(true); i++) {
            trips.push_back(trip_list.get_next_trip(true));
            trip_list.increase_index();
        }
        obj.add_list("trips", trips.begin(), trips.end());
        obj.add_list("locations", get_locations().begin(), get_locations().end());
        obj.add_list("persons", get_persons().begin(), get_persons().end());
        obj.add_element("use_mobility_rules", m_use_mobility_rules);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Model> deserialize(IOContext& io)
    {
        auto obj                = io.expect_object("Model");
        auto size               = obj.expect_element("num_agegroups", Tag<size_t>{});
        auto locations          = obj.expect_list("locations", Tag<Location>{});
        auto trip_list          = obj.expect_list("trips", Tag<Trip>{});
        auto persons            = obj.expect_list("persons", Tag<Person>{});
        auto use_mobility_rules = obj.expect_element("use_mobility_rules", Tag<bool>{});
        return apply(
            io,
            [](auto&& size_, auto&& locations_, auto&& trip_list_, auto&& persons_, auto&& use_mobility_rule_) {
                return Model{size_, locations_, trip_list_, persons_, use_mobility_rule_};
            },
            size, locations, trip_list, persons, use_mobility_rules);
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
     * @return ID of the newly created Person.
     */
    PersonId add_person(const LocationId id, AgeGroup age);

    /**
     * @brief Adds a copy of a given Person to the Model.
     * @param[in] person The Person to copy from. 
     * @return ID of the newly created Person.
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
     * @brief Find an assigned Location of a Person.
     * @param[in] type The #LocationType that specifies the assigned Location.
     * @param[in] person_index Index of the Person.
     * @return ID of the Location of LocationType type assigend to person.
     */
    LocationId find_location(LocationType type, const uint32_t person_index) const;

    /**
     * @brief Assign a Location to a Person.
     * A Person can have at most one assigned Location of a certain LocationType.
     * Assigning another Location of an already assigned LocationType will replace the prior assignment.  
     * @param[in] person_index The Index of the person this location will be assigned to.
     * @param[in] location The LocationId of the Location.
     */
    void assign_location(uint32_t person_index, LocationId location)
    {
        get_person(person_index).set_assigned_location(get_location(location).get_type(), location, m_id);
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
     * Get activeness status of all persons in the model.
     * @return Activeness vector
     */
    std::vector<bool>& get_activeness_statuses()
    {
        return m_activeness_statuses;
    }

    /**
     * @brief Add a TestingScheme to the set of schemes that are checked for testing at all Locations that have
     * the LocationType.
     * @param[in] loc_type LocationId key for TestingScheme to be added.
     * @param[in] scheme TestingScheme to be added.
     */
    void add_testing_scheme(const LocationType& loc_type, const TestingScheme& scheme);

    /**
     * @brief Remove a TestingScheme from the set of schemes that are checked for testing at all Locations that have
     * the LocationType.
     * @param[in] loc_type LocationId key for TestingScheme to be added.
     * @param[in] scheme TestingScheme to be added.
     */
    void remove_testing_scheme(const LocationType& loc_type, const TestingScheme& scheme);

    /**
     * @brief Get a reference to a Person from this Model.
     * @param[in] index A Person's index in m_persons.
     * @return A reference to the Person.
     * @{
     */
    Person& get_person(uint32_t index)
    {
        assert(index < m_persons.size() && "Given PersonId is not in this Model.");
        return m_persons[index];
    }

    const Person& get_person(uint32_t index) const
    {
        assert(index < m_persons.size() && "Given PersonId is not in this Model.");
        return m_persons[index];
    }

    /**
     * @brief Get a reference to a Person from this Model.
     * @param[in] id A Person's id.
     * @return A reference to the Person.
     */
    Person& get_person(PersonId id)
    {
        auto it = std::find_if(m_persons.begin(), m_persons.end(), [id](auto& person) {
            return person.get_id() == id;
        });
        if (it == m_persons.end()) {
            log_error("Given PersonId is not in this Model.");
        }
        return *it;
    }

    const Person& get_person(PersonId id) const
    {
        auto it = std::find_if(m_persons.begin(), m_persons.end(), [id](auto& person) {
            return person.get_id() == id;
        });
        if (it == m_persons.end()) {
            log_error("Given PersonId is not in this Model.");
        }
        return *it;
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
     * @param[in] person_index Index of a Person in m_persons vector of this Model.
     * @param[in] destination LocationId of the Location in this Model, which the Person should change to.
     * @param[in] mode The transport mode the person uses to change the Location.
     * @param[in] cells The cells within the destination the person should be in.
     */
    inline void change_location(uint32_t person_index, LocationId destination,
                                TransportMode mode = TransportMode::Unknown, const std::vector<uint32_t>& cells = {0})
    {
        LocationId origin = get_location(person_index).get_id();
        const bool has_changed_location =
            mio::abm::change_location(get_person(person_index), get_location(destination), mode, cells);
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
     * @brief Let a person interact with the population at its current location.
     * @param[in] person_index Index of a person in m_persons vector of this Model.
     * @param[in] t Time step of the simulation.
     * @param[in] dt Step size of the simulation.
     */
    inline void interact(uint32_t person_index, TimePoint t, TimeSpan dt)
    {
        if (!m_are_exposure_caches_valid) {
            // checking caches is only needed for external calls
            // during simulation (i.e. in evolve()), the caches are computed in begin_step
            compute_exposure_caches(t, dt);
            m_are_exposure_caches_valid = true;
        }
        auto personal_rng = PersonalRandomNumberGenerator(m_rng, get_person(person_index));
        mio::abm::interact(personal_rng, get_person(person_index), get_location(person_index),
                           m_air_exposure_rates_cache[get_location(person_index).get_id().get()],
                           m_contact_exposure_rates_cache[get_location(person_index).get_id().get()], t, dt,
                           parameters);
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

    Location& get_location(LocationId id)
    {
        assert(id != LocationId::invalid_id() && "Given LocationId must be valid.");
        assert(id < LocationId((uint32_t)m_locations.size()) && "Given LocationId is not in this Model.");
        return m_locations[id.get()];
    }
    /** @} */

    /**
     * @brief Get a reference to the location of a person.
     * @param[in] index Index of a person in m_persons.
     * @return Reference to the Location.
     * @{
     */
    inline Location& get_location(uint32_t index)
    {
        return get_location(get_person(index).get_location());
    }

    inline const Location& get_location(uint32_t index) const
    {
        return get_location(get_person(index).get_location());
    }
    /** @} */

    /**
     * @brief Flip activeness status of a person in the model.
     * @param[in] person_id PersonId of Person whose activeness status is fipped.
     */
    void change_activeness(PersonId person_id)
    {
        m_activeness_statuses[person_id.get()] = !m_activeness_statuses[person_id.get()];
    }

    void set_activeness(PersonId person_id)
    {
        m_activeness_statuses[person_id.get()] = false;
    }

    // /**
    //  * @brief Copy the persons from another Model to this Model.
    //  * If the persons are at a location in this model they are activated, otherwise they are deactivated.
    //  * If necessary the person ids are changed such that they correspond to the index in this model's m_persons vector.
    //  * @param[in] other The Model the Person%s are copied from.
    //  */
    // void copy_persons_from_other_model(const Model& other)
    // {
    //     for (auto& p : other.get_persons()) {
    //         if (p.get_id() != static_cast<uint32_t>(m_persons.size())) {
    //             mio::log_debug("In model.copy_persons_from_other_model: PersonId does not correspond to index in "
    //                            "m_persons vector. Person is copied with adapted Id");
    //         }
    //         PersonId new_id{static_cast<uint32_t>(m_persons.size())};
    //         m_persons.emplace_back(p, new_id);
    //         if (p.get_location_model_id() == m_id) {
    //             m_activeness_statuses.push_back(true);
    //         }
    //         else {
    //             m_activeness_statuses.push_back(false);
    //         }
    //     }
    // }

    /**
    * @brief Set the Person%s of the Model.
    * @param[in] persons The Person%s of the Model.
    */
    void set_persons(std::vector<Person>& persons)
    {
        m_is_local_population_cache_valid = false;
        m_are_exposure_caches_valid       = false;
        //first clear old person vector and corresponding activeness vector
        m_persons.clear();
        m_activeness_statuses.clear();
        for (auto& p : persons) {
            if (p.get_id() != static_cast<uint32_t>(m_persons.size())) {
                mio::log_debug("In model.copy_persons_from_other_model: PersonId does not correspond to index in "
                               "m_persons vector. Person is copied with adapted Id");
            }
            PersonId new_id{static_cast<uint32_t>(m_persons.size())};
            m_persons.emplace_back(p, new_id);
            if (p.get_location_model_id() == m_id) {
                m_activeness_statuses.push_back(true);
            }
            else {
                m_activeness_statuses.push_back(false);
            }
        }
    }

    /**
     * @brief Invalidate local population and exposure rate cache.
     */
    void invalidate_cache()
    {
        m_are_exposure_caches_valid       = false;
        m_is_local_population_cache_valid = false;
    }

    /**
     * @brief Get index of person in m_persons.
     * @param[in] id A person's PersonId. 
     * First 32 bit are the Person's individual id and second 32 bit the Persons's home model id. 
     * @return Index of Person in m_persons vector.
     * @{
     */
    uint32_t get_person_index(PersonId id) const
    {
        auto it = std::find_if(m_persons.begin(), m_persons.end(), [id](auto& person) {
            return person.get_id() == id;
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
    std::vector<bool> m_activeness_statuses; ///< Vector with activeness status for every person
    std::bitset<size_t(LocationType::Count)>
        m_has_locations; ///< Flags for each LocationType, set if a Location of that type exists.
    TestingStrategy m_testing_strategy; ///< List of TestingScheme%s that are checked for testing.
    TripList m_trip_list; ///< List of all Trip%s the Person%s do.
    bool m_use_mobility_rules; ///< Whether mobility rules are considered.
    std::vector<std::pair<LocationType (*)(PersonalRandomNumberGenerator&, const Person&, TimePoint, TimeSpan,
                                           const Parameters&),
                          std::vector<LocationType>>>
        m_mobility_rules; ///< Rules that govern the mobility between Location%s.
    LocationId m_cemetery_id; // Central cemetery for all dead persons.
    RandomNumberGenerator m_rng; ///< Global random number generator
};

} // namespace abm
} // namespace mio

#endif
