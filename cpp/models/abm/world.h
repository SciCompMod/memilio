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

#include "abm/config.h" // IWYU pragma: keep
#include "abm/location_type.h"
#include "abm/parameters.h"
#include "abm/location.h"
#include "abm/person.h"
#include "abm/trip_list.h"
#include "abm/testing_strategy.h"
#include "memilio/utils/pointer_dereferencing_iterator.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/stl_util.h"
#include "abm/migration_rules.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/mioomp.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/stl_util.h"

#include <bitset>
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
template <typename FP = double>
class World
{
public:
    using LocationIterator =
        PointerDereferencingIterator<typename std::vector<std::unique_ptr<Location<FP>>>::iterator>;
    using ConstLocationIterator =
        PointerDereferencingIterator<typename std::vector<std::unique_ptr<Location<FP>>>::const_iterator>;
    using PersonIterator = PointerDereferencingIterator<typename std::vector<std::unique_ptr<Person<FP>>>::iterator>;
    using ConstPersonIterator =
        PointerDereferencingIterator<typename std::vector<std::unique_ptr<Person<FP>>>::const_iterator>;

    /**
     * @brief Create a World.
     * @param[in] num_agegroups The number of AgeGroup%s in the simulated World. Must be less than MAX_NUM_AGE_GROUPS.
     */
    World(size_t num_agegroups)
        : parameters(num_agegroups)
        , m_trip_list()
        , m_use_migration_rules(true)
        , m_cemetery_id(add_location(LocationType::Cemetery))
    {
        assert(num_agegroups < MAX_NUM_AGE_GROUPS && "MAX_NUM_AGE_GROUPS exceeded.");
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
                    std::make_unique<Location<FP>>(origin_loc.copy_location_without_persons(parameters.get_num_groups())));
            }
            for (auto& person : other.get_persons()) {
                // If a person is in this location, copy this person and add it to this location.
                if (person.get_location() == origin_loc) {
                    LocationId origin_id = {origin_loc.get_index(), origin_loc.get_type()};
                    m_persons.push_back(
                        std::make_unique<Person<FP>>(person.copy_person(get_individualized_location(origin_id))));
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
        auto locations           = obj.expect_list("locations", Tag<Location<FP>>{});
        auto trip_list           = obj.expect_list("trips", Tag<Trip>{});
        auto persons             = obj.expect_list("persons", Tag<Person<FP>>{});
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
    void begin_step(TimePoint t, TimeSpan dt)
    {
        m_testing_strategy.update_activity_status(t);
    PRAGMA_OMP(parallel for)
    for (auto i = size_t(0); i < m_locations.size(); ++i) {
        auto&& location = m_locations[i];
        location->cache_exposure_rates(t, dt, parameters.get_num_groups());
    }
    }

    /** 
     * @brief Evolve the world one time step.
     * @param[in] t Current time.
     * @param[in] dt Length of the time step.
     */
    void evolve(TimePoint t, TimeSpan dt)
    {
        begin_step(t, dt);
        log_info("ABM World interaction.");
        interaction(t, dt);
        log_info("ABM World migration.");
        migration(t, dt);
    }

    /** 
     * @brief Add a Location to the World.
     * @param[in] type Type of Location to add.
     * @param[in] num_cells [Default: 1] Number of Cell%s that the Location is divided into.
     * @return Index and type of the newly created Location.
     */
    LocationId add_location(LocationType type, uint32_t num_cells = 1)
    {
        LocationId id = {static_cast<uint32_t>(m_locations.size()), type};
        m_locations.emplace_back(std::make_unique<Location<FP>>(id, parameters.get_num_groups(), num_cells));
        m_has_locations[size_t(type)] = true;
        return id;
    }

    /** 
     * @brief Add a Person to the World.
     * @param[in] id Index and type of the initial Location of the Person.
     * @param[in] age AgeGroup of the person.
     * @return Reference to the newly created Person.
     */
    Person<FP>& add_person(const LocationId id, AgeGroup age)
    {
        assert(age.get() < parameters.get_num_groups());
        uint32_t person_id = static_cast<uint32_t>(m_persons.size());
        m_persons.push_back(std::make_unique<Person<FP>>(m_rng, get_individualized_location(id), age, person_id));
        auto& person = *m_persons.back();
        person.set_assigned_location(m_cemetery_id);
        get_individualized_location(id).add_person(person);
        return person;
    }

    /**
     * @brief Get a range of all Location%s in the World.
     * @return A range of all Location%s.
     */
    Range<std::pair<ConstLocationIterator, ConstLocationIterator>> get_locations() const
    {
        return std::make_pair(ConstLocationIterator(m_locations.begin()), ConstLocationIterator(m_locations.end()));
    }

    /**
     * @brief Get a range of all Person%s in the World.
     * @return A range of all Person%s.
     */
    Range<std::pair<ConstPersonIterator, ConstPersonIterator>> get_persons() const
    {
        return std::make_pair(ConstPersonIterator(m_persons.begin()), ConstPersonIterator(m_persons.end()));
    }

    /**
     * @brief Get an individualized Location.
     * @param[in] id LocationId of the Location.
     * @return Reference to the Location.
     */
    const Location<FP>& get_individualized_location(LocationId id) const
    {
        return *m_locations[id.index];
    }

    Location<FP>& get_individualized_location(LocationId id)
    {
        return *m_locations[id.index];
    }

    /**
     * @brief Find an assigned Location of a Person.
     * @param[in] type The #LocationType that specifies the assigned Location.
     * @param[in] person The Person.
     * @return Reference to the assigned Location.
     */
    const Location<FP>& find_location(LocationType type, const Person<FP>& person) const
    {
        auto index = person.get_assigned_location_index(type);
        assert(index != INVALID_LOCATION_INDEX && "unexpected error.");
        return get_individualized_location({index, type});
    }

    Location<FP>& find_location(LocationType type, const Person<FP>& person)
    {
        auto index = person.get_assigned_location_index(type);
        assert(index != INVALID_LOCATION_INDEX && "unexpected error.");
        return get_individualized_location({index, type});
    }

    /** 
     * @brief Get the number of Persons in one #InfectionState at all Location%s.
     * @param[in] t Specified #TimePoint.
     * @param[in] s Specified #InfectionState.
     */
    size_t get_subpopulation_combined(TimePoint t, InfectionState s) const
    {
        return std::accumulate(m_locations.begin(), m_locations.end(), (size_t)0,
                               [t, s](size_t running_sum, const std::unique_ptr<Location<FP>>& loc) {
                                   return running_sum + loc->get_subpopulation(t, s);
                               });
    }

    /** 
     * @brief Get the number of Persons in one #InfectionState at all Location%s of a type.
     * @param[in] t Specified #TimePoint.
     * @param[in] s Specified #InfectionState.
     * @param[in] type Specified #LocationType.
     */
    size_t get_subpopulation_combined_per_location_type(TimePoint t, InfectionState s, LocationType type) const
    {
        return std::accumulate(m_locations.begin(), m_locations.end(), (size_t)0,
                               [t, s, type](size_t running_sum, const std::unique_ptr<Location<FP>>& loc) {
                                   return loc->get_type() == type ? running_sum + loc->get_subpopulation(t, s)
                                                                  : running_sum;
                               });
    }

    /**
     * @brief Get the migration data.
     * @return Reference to the list of Trip%s that the Person%s make.
     */
    TripList& get_trip_list()
    {
        return m_trip_list;
    }

    const TripList& get_trip_list() const
    {
        return m_trip_list;
    }

    /** 
     * @brief Decide if migration rules (like go to school/work) are used or not;
     * The migration rules regarding hospitalization/ICU/quarantine are always used.
     * @param[in] param If true uses the migration rules for migration to school/work etc., else only the rules 
     * regarding hospitalization/ICU/quarantine.
     */
    void use_migration_rules(bool param)
    {
        m_use_migration_rules = param;
    }
    bool use_migration_rules() const
    {
        return m_use_migration_rules;
    }

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
    TestingStrategy<FP>& get_testing_strategy()
    {
        return m_testing_strategy;
    }

    const TestingStrategy<FP>& get_testing_strategy() const
    {
        return m_testing_strategy;
    }

    /** 
     * @brief The simulation parameters of the world.
     */
    Parameters<FP> parameters;

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
    void add_testing_scheme(const LocationType& loc_type, const TestingScheme<FP>& scheme);

    /**
     * @brief Remove a TestingScheme from the set of schemes that are checked for testing at all Locations that have 
     * the LocationType.
     * @param[in] loc_type LocationId key for TestingScheme to be added.
     * @param[in] scheme TestingScheme to be added.
     */
    void remove_testing_scheme(const LocationType& loc_type, const TestingScheme<FP>& scheme);

private:
    /**
     * @brief Person%s interact at their Location and may become infected.
     * @param[in] t The current TimePoint.
     * @param[in] dt The length of the time step of the Simulation.
     */
    void interaction(TimePoint t, TimeSpan dt)
    {
    PRAGMA_OMP(parallel for)
    for (auto i = size_t(0); i < m_persons.size(); ++i) {
        auto&& person     = m_persons[i];
        auto personal_rng = typename Person<FP>::RandomNumberGenerator(m_rng, *person);
        person->interact(personal_rng, t, dt, parameters);
    }
    }
    /**
     * @brief Person%s move in the World according to rules.
     * @param[in] t The current TimePoint.
     * @param[in] dt The length of the time step of the Simulation.
     */
    void migration(TimePoint t, TimeSpan dt)
{
    PRAGMA_OMP(parallel for)
    for (auto i = size_t(0); i < m_persons.size(); ++i) {
        auto&& person     = m_persons[i];
        auto personal_rng = typename Person<FP>::RandomNumberGenerator(m_rng, *person);

        auto try_migration_rule = [&](auto rule) -> bool {
            //run migration rule and check if migration can actually happen
            auto target_type       = rule(personal_rng, *person, t, dt, parameters);
            auto& target_location  = find_location(target_type, *person);
            auto& current_location = person->get_location();
            if (m_testing_strategy.run_strategy(personal_rng, *person, target_location, t)) {
                if (target_location != current_location &&
                    target_location.get_number_persons() < target_location.get_capacity().persons) {
                    bool wears_mask = person->apply_mask_intervention(personal_rng, target_location);
                    if (wears_mask) {
                        person->migrate_to(target_location);
                    }
                    return true;
                }
            }
            return false;
        };

        //run migration rules one after the other if the corresponding location type exists
        //shortcutting of bool operators ensures the rules stop after the first rule is applied
        if (m_use_migration_rules) {
            (has_locations({LocationType::Cemetery}) && try_migration_rule(&get_buried<FP>)) ||
                (has_locations({LocationType::Home}) && try_migration_rule(&return_home_when_recovered<FP>)) ||
                (has_locations({LocationType::Hospital}) && try_migration_rule(&go_to_hospital<FP>)) ||
                (has_locations({LocationType::ICU}) && try_migration_rule(&go_to_icu<FP>)) ||
                (has_locations({LocationType::School, LocationType::Home}) && try_migration_rule(&go_to_school<FP>)) ||
                (has_locations({LocationType::Work, LocationType::Home}) && try_migration_rule(&go_to_work<FP>)) ||
                (has_locations({LocationType::BasicsShop, LocationType::Home}) && try_migration_rule(&go_to_shop<FP>)) ||
                (has_locations({LocationType::SocialEvent, LocationType::Home}) && try_migration_rule(&go_to_event<FP>)) ||
                (has_locations({LocationType::Home}) && try_migration_rule(&go_to_quarantine<FP>));
        }
        else {
            //no daily routine migration, just infection related
            (has_locations({LocationType::Cemetery}) && try_migration_rule(&get_buried<FP>)) ||
                (has_locations({LocationType::Home}) && try_migration_rule(&return_home_when_recovered<FP>)) ||
                (has_locations({LocationType::Hospital}) && try_migration_rule(&go_to_hospital<FP>)) ||
                (has_locations({LocationType::ICU}) && try_migration_rule(&go_to_icu<FP>)) ||
                (has_locations({LocationType::Home}) && try_migration_rule(&go_to_quarantine<FP>));
        }
    }

    // check if a person makes a trip
    bool weekend     = t.is_weekend();
    size_t num_trips = m_trip_list.num_trips(weekend);

    if (num_trips != 0) {
        while (m_trip_list.get_current_index() < num_trips &&
               m_trip_list.get_next_trip_time(weekend).seconds() < (t + dt).time_since_midnight().seconds()) {
            auto& trip        = m_trip_list.get_next_trip(weekend);
            auto& person      = m_persons[trip.person_id];
            auto personal_rng = typename Person<FP>::RandomNumberGenerator(m_rng, *person);
            if (!person->is_in_quarantine(t, parameters) && person->get_infection_state(t) != InfectionState::Dead) {
                auto& target_location = get_individualized_location(trip.migration_destination);
                if (m_testing_strategy.run_strategy(personal_rng, *person, target_location, t)) {
                    person->apply_mask_intervention(personal_rng, target_location);
                    person->migrate_to(target_location, trip.trip_mode);
                }
            }
            m_trip_list.increase_index();
        }
    }
    if (((t).days() < std::floor((t + dt).days()))) {
        m_trip_list.reset_index();
    }
}

    std::vector<std::unique_ptr<Person<FP>>> m_persons; ///< Vector with pointers to every Person.
    std::vector<std::unique_ptr<Location<FP>>> m_locations; ///< Vector with pointers to every Location.
    std::bitset<size_t(LocationType::Count)>
        m_has_locations; ///< Flags for each LocationType, set if a Location of that type exists.
    TestingStrategy<FP> m_testing_strategy; ///< List of TestingScheme%s that are checked for testing.
    TripList m_trip_list; ///< List of all Trip%s the Person%s do.
    bool m_use_migration_rules; ///< Whether migration rules are considered.
    std::vector<std::pair<LocationType (*)(typename Person<FP>::RandomNumberGenerator&, const Person<FP>&, TimePoint,
                                           TimeSpan, const Parameters<FP>&),
                          std::vector<LocationType>>>
        m_migration_rules; ///< Rules that govern the migration between Location%s.
    LocationId m_cemetery_id; // Central cemetery for all dead persons.
    RandomNumberGenerator m_rng; ///< Global random number generator
};

} // namespace abm
} // namespace mio

#endif
