/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

#include "abm/parameters.h"
#include "abm/location.h"
#include "abm/person.h"
#include "abm/lockdown_rules.h"
#include "abm/trip_list.h"
#include "abm/testing_strategy.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/pointer_dereferencing_iterator.h"
#include "memilio/utils/stl_util.h"
#include "memilio/epidemiology/populations.h"

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
        , m_cemetery_id(add_location(LocationType::Cemetery))
    {
        use_migration_rules(true);
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
                m_locations.emplace_back(std::make_unique<Location>(origin_loc));
            }
            for (auto& person : other.get_persons()) {
                // If a person is in this location, copy this person and add it to this location.
                if (person.get_location() == origin_loc) {
                    LocationId origin_id = {origin_loc.get_index(), origin_loc.get_type()};
                    m_persons.push_back(std::make_unique<Person>(person, get_individualized_location(origin_id)));
                    auto& copied_person = *m_persons.back();
                    get_individualized_location(origin_id).add_person(copied_person);
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
        for (size_t i = 0; i < trip_list.num_trips(); i++) {
            trips.push_back(trip_list.get_next_trip());
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
     * @brief Follow up on the World after the Simulation step.
     * @param[in] t Current time.
     * @param[in] dt Length of the time step.
     */
    void end_step(TimePoint t, TimeSpan dt);

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
    Location& find_location(LocationType type, const Person& person);

    /** 
     * @brief Get the number of Persons in one #InfectionState at all Location%s of a type.
     * @param[in] t Specified #TimePoint.
     * @param[in] s Specified #InfectionState.
     * @param[in] type Specified #LocationType.
     */
    size_t get_subpopulation_combined(TimePoint t, InfectionState s, LocationType type) const;

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
     * @brief Get the TestingStrategy.
     * @return Reference to the list of TestingScheme%s that are checked for testing.
     */
    TestingStrategy& get_testing_strategy();

    const TestingStrategy& get_testing_strategy() const;

    /** 
     * @brief The simulation parameters of the world.
     */
    Parameters parameters;

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
    TestingStrategy m_testing_strategy; ///< List of TestingScheme%s that are checked for testing.
    TripList m_trip_list; ///< List of all Trip%s the Person%s do.
    bool m_use_migration_rules; ///< Whether migration rules are considered.
    std::vector<std::pair<LocationType (*)(const Person&, TimePoint, TimeSpan, const Parameters&),
                          std::vector<LocationType>>>
        m_migration_rules; ///< Rules that govern the migration between Location%s.
    LocationId m_cemetery_id; // Central cemetery for all dead persons.
};

} // namespace abm
} // namespace mio

#endif