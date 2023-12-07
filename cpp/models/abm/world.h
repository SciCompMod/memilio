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
#ifndef EPI_ABM_WORLD_H
#define EPI_ABM_WORLD_H

#include "abm/location_type.h"
#include "abm/movement_data.h"
#include "abm/parameters.h"
#include "abm/location.h"
#include "abm/person.h"
#include "abm/lockdown_rules.h"
#include "abm/trip_list.h"
#include "abm/random_events.h"
#include "abm/testing_strategy.h"
#include "memilio/utils/pointer_dereferencing_iterator.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/stl_util.h"

#include <bitset>
#include <cassert>
#include <initializer_list>
#include <vector>
#include <memory>

namespace mio
{
namespace abm
{

/**
 * @brief The World of the Simulation.
 * It consists of Location%s and Person%s (Agents).
 */
class World
{
public:
    using LocationIterator      = PointerDereferencingIterator<std::vector<std::unique_ptr<Location>>::iterator>;
    using ConstLocationIterator = PointerDereferencingIterator<std::vector<std::unique_ptr<Location>>::const_iterator>;
    using PersonIterator        = PointerDereferencingIterator<std::vector<std::unique_ptr<Person>>::iterator>;
    using ConstPersonIterator   = PointerDereferencingIterator<std::vector<std::unique_ptr<Person>>::const_iterator>;

    /**
     * @brief Create a World.
     * @param[in] num_agegroups The number of AgeGroup%s in the simulated World.
     */
    World(size_t num_agegroups)
        : parameters(num_agegroups)
        , m_trip_list()
        , m_use_migration_rules(true)
        , m_cemetery_id(add_location(LocationType::Cemetery))
    {
    }

    /**
     * @brief Create a World.
     * @param[in] params Initial simulation parameters.
     */
    World(const Parameters& params)
        : parameters(params.get_num_groups())
        , m_trip_list()
        , m_use_migration_rules(true)
        , m_cemetery_id(add_location(LocationType::Cemetery))
    {
        parameters = params;
    }

    /**
     * @brief Create a copied World.
     * @param[in] other The World that needs to be copied. 
     */
    World(const World& other)
        : parameters(other.parameters)
        , m_persons()
        , m_locations()
        , m_trip_list(other.m_trip_list)
        , m_cemetery_id(add_location(LocationType::Cemetery))
    {
        for (auto& origin_loc : other.get_locations()) {
            if (origin_loc.get_type() != LocationType::Cemetery) {
                // Copy a location
                m_locations.emplace_back(
                    std::make_unique<Location>(origin_loc.copy_location_without_persons(parameters.get_num_groups())));
            }
            for (auto& person : other.get_persons()) {
                // If a person is in this location, copy this person and add it to this location.
                if (person.get_location() == origin_loc.get_id()) {
                    m_persons.push_back(
                        std::make_unique<Person>(person.copy_person(get_individualized_location(origin_loc.get_id()))));
                }
            }
        }
        use_migration_rules(other.m_use_migration_rules);
    }

    //type is move-only for stable references of persons/locations
    World(World&& other)            = default;
    World& operator=(World&& other) = default;
    World& operator=(const World&)  = delete;

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("World");
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
        obj.add_element("use_migration_rules", m_use_migration_rules);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<World> deserialize(IOContext& io)
    {
        auto obj                 = io.expect_object("World");
        auto size                = obj.expect_element("num_agegroups", Tag<size_t>{});
        auto locations           = obj.expect_list("locations", Tag<Location>{});
        auto trip_list           = obj.expect_list("trips", Tag<Trip>{});
        auto persons             = obj.expect_list("persons", Tag<Person>{});
        auto use_migration_rules = obj.expect_element("use_migration_rules", Tag<bool>{});
        return apply(
            io,
            [](auto&& size_, auto&& locations_, auto&& trip_list_, auto&& persons_, auto&& use_migration_rule_) {
                return World{size_, locations_, trip_list_, persons_, use_migration_rule_};
            },
            size, locations, trip_list, persons, use_migration_rules);
    }

    /** 
     * @brief Prepare the World for the next Simulation step.
     * @param[in] t Current time.
     * @param[in] dt Length of the time step.
     */
    void begin_step(TimePoint t, TimeSpan dt);

    /** 
     * @brief Evolve the world one time step.
     * @param[in] t Current time.
     * @param[in] dt Length of the time step.
     */
    void evolve(TimePoint t, TimeSpan dt);

    /** 
     * @brief Add a Location to the World.
     * @param[in] type Type of Location to add.
     * @param[in] num_cells [Default: 1] Number of Cell%s that the Location is divided into.
     * @return Index and type of the newly created Location.
     */
    LocationId add_location(LocationType type, uint32_t num_cells = 1);

    /** 
     * @brief Add a Person to the World.
     * @param[in] id Index and type of the initial Location of the Person.
     * @param[in] age AgeGroup of the person.
     * @return Reference to the newly created Person.
     */
    Person& add_person(const LocationId id, AgeGroup age);

    /**
     * @brief Get a range of all Location%s in the World.
     * @return A range of all Location%s.
     */
    Range<std::pair<ConstLocationIterator, ConstLocationIterator>> get_locations() const;

    /**
     * @brief Get a range of all Person%s in the World.
     * @return A range of all Person%s.
     */
    Range<std::pair<ConstPersonIterator, ConstPersonIterator>> get_persons() const;

    /**
     * @brief Get an individualized Location.
     * @param[in] id LocationId of the Location.
     * @return Reference to the Location.
     */
    const Location& get_individualized_location(LocationId id) const;

    Location& get_individualized_location(LocationId id);

    /**
     * @brief Find an assigned Location of a Person.
     * @param[in] type The #LocationType that specifies the assigned Location.
     * @param[in] person The Person.
     * @return Reference to the assigned Location.
     */
    const Location& find_location(LocationType type, const Person& person) const;

    Location& find_location(LocationType type, const Person& person);

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
     * @brief Get the migration data.
     * @return Reference to the list of Trip%s that the Person%s make.
     */
    TripList& get_trip_list();

    const TripList& get_trip_list() const;

    /** 
     * @brief Decide if migration rules (like go to school/work) are used or not;
     * The migration rules regarding hospitalization/ICU/quarantine are always used.
     * @param[in] param If true uses the migration rules for migration to school/work etc., else only the rules 
     * regarding hospitalization/ICU/quarantine.
     */
    void use_migration_rules(bool param);
    bool use_migration_rules() const;

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
     * @brief The simulation parameters of the world.
     */
    Parameters parameters;

    /**
    * Get the RandomNumberGenerator used by this world for random events.
    * Persons use their own generators with the same key as the global one. 
    * @return The random number generator.
    */
    RandomNumberGenerator& get_rng()
    {
        return m_rng;
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

    // std::unordered_map<LocationId, std::unordered_set<PersonID>> m_cache_local_population;

    // move a person to another location. this requires that location is part of this world.
    void migrate(Person& person, Location& destination, TransportMode mode = TransportMode::Unknown,
                 const std::vector<uint32_t>& cells = {0})
    {
        assert(get_location(destination.get_id()) == destination &&
               "Destination is outside of World."); // always true, but may trigger asserts in get_location
        if (person.get_location() == destination.get_id()) {
            return;
        }
        get_location(person).remove_person(person);
        person.set_location(destination);
        person.get_cells() = cells;
        destination.add_person(person, cells);
        person.set_last_transport_mode(mode);
    }

    // let a person interact with its current location
    void interact(Person& person, TimePoint t, TimeSpan dt)
    {
        auto personal_rng = Person::RandomNumberGenerator(m_rng, person);
        interact(person, t, dt, personal_rng, parameters);
    }

    // let a person interact with its current location
    void interact(Person& person, TimePoint t, TimeSpan dt, Person::RandomNumberGenerator& personal_rng,
                  const Parameters& global_parameters)
    {
        interact(person, get_location(person), t, dt, personal_rng, global_parameters);
    }

    // let a person interact with a location for and at some time
    static void interact(Person& person, Location& location, TimePoint t, TimeSpan dt,
                         Person::RandomNumberGenerator& personal_rng, const Parameters& global_parameters)
    {
        if (person.get_infection_state(t) == InfectionState::Susceptible) {
            auto& local_parameters = location.get_infection_parameters();
            // TODO: we need to define what a cell is used for, as the loop may lead to incorrect results for multiple cells
            auto age_receiver          = person.get_age();
            ScalarType mask_protection = person.get_mask_protective_factor(global_parameters);
            assert(person.get_cells().size() &&
                   "Person is in multiple cells. Interact logic is incorrect at the moment.");
            for (auto cell_index :
                 person.get_cells()) { // TODO: the logic here is incorrect in case a person is in multiple cells
                std::pair<VirusVariant, ScalarType> local_indiv_trans_prob[static_cast<uint32_t>(VirusVariant::Count)];
                for (uint32_t v = 0; v != static_cast<uint32_t>(VirusVariant::Count); ++v) {
                    VirusVariant virus = static_cast<VirusVariant>(v);
                    ScalarType local_indiv_trans_prob_v =
                        (std::min(local_parameters.get<MaximumContacts>(),
                                  location.transmission_contacts_per_day(cell_index, virus, age_receiver,
                                                                         global_parameters.get_num_groups())) +
                         location.transmission_air_per_day(cell_index, virus, global_parameters)) *
                        (1 - mask_protection) * dt.days() *
                        (1 - person.get_protection_factor(t, virus, global_parameters));

                    local_indiv_trans_prob[v] = std::make_pair(virus, local_indiv_trans_prob_v);
                }
                VirusVariant virus =
                    random_transition(personal_rng, VirusVariant::Count, dt,
                                      local_indiv_trans_prob); // use VirusVariant::Count for no virus submission
                if (virus != VirusVariant::Count) {
                    person.add_new_infection(Infection(personal_rng, virus, age_receiver, global_parameters, t + dt / 2,
                                                       mio::abm::InfectionState::Exposed,
                                                       person.get_latest_protection(),
                                                       false)); // Starting time in first approximation
                }
            }
        }
        person.add_time_at_location(dt);
    }

    // get location by id
    Location& get_location(LocationId id)
    {
        assert(id.index != INVALID_LOCATION_INDEX);
        for (auto&& location : m_locations) {
            if (location->get_id() == id) {
                return *location;
            }
        }
        assert(false && "Location id not found");
    }

    // get current location of the Person
    Location& get_location(const Person& p)
    {
        return get_location(p.get_location());
    }

private:
    /**
     * @brief Person%s interact at their Location and may become infected.
     * @param[in] t The current TimePoint.
     * @param[in] dt The length of the time step of the Simulation.
     */
    void interaction(TimePoint t, TimeSpan dt);
    /**
     * @brief Person%s move in the World according to rules.
     * @param[in] t The current TimePoint.
     * @param[in] dt The length of the time step of the Simulation.
     */
    void migration(TimePoint t, TimeSpan dt);

    std::vector<std::unique_ptr<Person>> m_persons; ///< Vector with pointers to every Person.
    std::vector<std::unique_ptr<Location>> m_locations; ///< Vector with pointers to every Location.
    std::bitset<size_t(LocationType::Count)>
        m_has_locations; ///< Flags for each LocationType, set if a Location of that type exists.
    TestingStrategy m_testing_strategy; ///< List of TestingScheme%s that are checked for testing.
    TripList m_trip_list; ///< List of all Trip%s the Person%s do.
    bool m_use_migration_rules; ///< Whether migration rules are considered.
    std::vector<std::pair<LocationType (*)(Person::RandomNumberGenerator&, const Person&, TimePoint, TimeSpan,
                                           const Parameters&),
                          std::vector<LocationType>>>
        m_migration_rules; ///< Rules that govern the migration between Location%s.
    LocationId m_cemetery_id; // Central cemetery for all dead persons.
    RandomNumberGenerator m_rng; ///< Global random number generator
};

} // namespace abm
} // namespace mio

#endif
