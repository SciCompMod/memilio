/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Majid Abedi, Elisabeth Kluth, Carlotta Gerstein, Martin J. Kuehn , David Kerkmann, Khoa Nguyen
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
#include "abm/world.h"
#include "abm/mask_type.h"
#include "abm/person.h"
#include "abm/location.h"
#include "abm/migration_rules.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/stl_util.h"
#include "abm/infection.h"
#include "abm/vaccine.h"

namespace mio
{
namespace abm
{

LocationId World::add_location(LocationType type, uint32_t num_cells)
{
    LocationId id = {static_cast<uint32_t>(m_locations.size()), type};
    m_locations.emplace_back(std::make_unique<Location>(id, num_cells));
    return id;
}

Person& World::add_person(const LocationId id, AgeGroup age)
{
    uint32_t person_id = static_cast<uint32_t>(m_persons.size());
    m_persons.push_back(std::make_unique<Person>(get_individualized_location(id), age, person_id));
    auto& person = *m_persons.back();
    person.set_assigned_location(m_cemetery_id);
    get_individualized_location(id).add_person(person);
    return person;
}

void World::evolve(TimePoint t, TimeSpan dt)
{
    begin_step(t, dt);
    interaction(t, dt);
    m_testing_strategy.update_activity_status(t);
    migration(t, dt);
    end_step(t, dt);
}

void World::interaction(TimePoint t, TimeSpan dt)
{
    for (auto&& person : m_persons) {
        person->interact(t, dt, m_infection_parameters);
    }
}

void World::migration(TimePoint t, TimeSpan dt)
{
    std::vector<std::pair<LocationType (*)(const Person&, TimePoint, TimeSpan, const MigrationParameters&),
                          std::vector<LocationType>>>
        m_enhanced_migration_rules;
    for (auto rule : m_migration_rules) {
        //check if transition rule can be applied
        bool nonempty         = false;
        const auto& loc_types = rule.second;
        for (auto loc_type : loc_types) {
            nonempty = std::find_if(m_locations.begin(), m_locations.end(),
                                    [loc_type](const std::unique_ptr<Location>& location) {
                                        return location->get_type() == loc_type;
                                    }) != m_locations.end();
        }

        if (nonempty) {
            m_enhanced_migration_rules.push_back(rule);
        }
    }
    for (auto& person : m_persons) {
        for (auto rule : m_enhanced_migration_rules) {
            //check if transition rule can be applied
            auto target_type      = rule.first(*person, t, dt, m_migration_parameters);
            auto& target_location = find_location(target_type, *person);
            auto current_location = person->get_location();
            if (m_testing_strategy.run_strategy(*person, target_location, t)) {
                if (target_location != current_location &&
                    target_location.get_number_persons() < target_location.get_capacity().persons) {
                    bool wears_mask = person->apply_mask_intervention(target_location);
                    if (wears_mask) {
                        person->migrate_to(target_location);
                    }
                    break;
                }
            }
        }
    }
    // check if a person makes a trip
    size_t num_trips = m_trip_list.num_trips();
    if (num_trips != 0) {
        while (m_trip_list.get_current_index() < num_trips && m_trip_list.get_next_trip_time() < t + dt) {
            auto& trip            = m_trip_list.get_next_trip();
            auto& person          = m_persons[trip.person_id];
            auto current_location = person->get_location();
            if (!person->is_in_quarantine() && person->get_infection_state(t) != InfectionState::Dead &&
                current_location == get_individualized_location(trip.migration_origin)) {
                auto& target_location = get_individualized_location(trip.migration_destination);
                if (m_testing_strategy.run_strategy(*person, target_location, t)) {
                    person->apply_mask_intervention(target_location);
                    person->migrate_to(target_location);
                }
            }
            m_trip_list.increase_index();
        }
    }
}

void World::begin_step(TimePoint t, TimeSpan dt)
{
    for (auto& location : m_locations) {
        location->cache_exposure_rates(t, dt);
    }
}

void World::end_step(TimePoint t, TimeSpan dt)
{
    for (auto& location : m_locations) {
        location->store_subpopulations(t + dt);
    }
}

auto World::get_locations() const -> Range<std::pair<ConstLocationIterator, ConstLocationIterator>>
{
    return std::make_pair(ConstLocationIterator(m_locations.begin()), ConstLocationIterator(m_locations.end()));
}

auto World::get_persons() const -> Range<std::pair<ConstPersonIterator, ConstPersonIterator>>
{
    return std::make_pair(ConstPersonIterator(m_persons.begin()), ConstPersonIterator(m_persons.end()));
}

const Location& World::get_individualized_location(LocationId id) const
{
    return *m_locations[id.index];
}

Location& World::get_individualized_location(LocationId id)
{
    return *m_locations[id.index];
}

Location& World::find_location(LocationType type, const Person& person)
{
    auto index = person.get_assigned_location_index(type);
    assert(index != INVALID_LOCATION_INDEX && "unexpected error.");
    return get_individualized_location({index, type});
}

size_t World::get_subpopulation_combined(TimePoint t, InfectionState s, LocationType type) const
{
    return std::accumulate(m_locations.begin(), m_locations.end(), (size_t)0,
                           [t, s, type](size_t running_sum, const std::unique_ptr<Location>& loc) {
                               return loc->get_type() == type ? running_sum + loc->get_subpopulation(t, s)
                                                              : running_sum;
                           });
}

MigrationParameters& World::get_migration_parameters()
{
    return m_migration_parameters;
}

const MigrationParameters& World::get_migration_parameters() const
{
    return m_migration_parameters;
}

GlobalInfectionParameters& World::get_global_infection_parameters()
{
    return m_infection_parameters;
}

const GlobalInfectionParameters& World::get_global_infection_parameters() const
{
    return m_infection_parameters;
}

TripList& World::get_trip_list()
{
    return m_trip_list;
}

const TripList& World::get_trip_list() const
{
    return m_trip_list;
}

void World::use_migration_rules(bool param)
{
    m_use_migration_rules = param;
    // Set up global migration rules for all agents
    // check if a person has to go to the hospital, ICU or home due to quarantine/recovery
    if (m_use_migration_rules) {
        m_migration_rules = {
            std::make_pair(&get_buried, std::vector<LocationType>{LocationType::ICU, LocationType::Cemetery}),
            std::make_pair(&return_home_when_recovered,
                           std::vector<LocationType>{
                               LocationType::Home,
                               LocationType::Hospital}), //assumption: if there is an ICU, there is also an hospital
            std::make_pair(&go_to_hospital, std::vector<LocationType>{LocationType::Home, LocationType::Hospital}),
            std::make_pair(&go_to_icu, std::vector<LocationType>{LocationType::Hospital, LocationType::ICU}),
            std::make_pair(&go_to_school, std::vector<LocationType>{LocationType::School, LocationType::Home}),
            std::make_pair(&go_to_work, std::vector<LocationType>{LocationType::Home, LocationType::Work}),
            std::make_pair(&go_to_shop, std::vector<LocationType>{LocationType::Home, LocationType::BasicsShop}),
            std::make_pair(&go_to_event, std::vector<LocationType>{LocationType::Home, LocationType::SocialEvent}),
            std::make_pair(&go_to_quarantine, std::vector<LocationType>{LocationType::Home})};
    }
    else {
        m_migration_rules = {
            std::make_pair(&get_buried, std::vector<LocationType>{LocationType::ICU, LocationType::Cemetery}),
            std::make_pair(&return_home_when_recovered,
                           std::vector<LocationType>{
                               LocationType::Home,
                               LocationType::Hospital}), //assumption: if there is an ICU, there is also an hospital
            std::make_pair(&go_to_hospital, std::vector<LocationType>{LocationType::Home, LocationType::Hospital}),
            std::make_pair(&go_to_icu, std::vector<LocationType>{LocationType::Hospital, LocationType::ICU}),
            std::make_pair(&go_to_quarantine, std::vector<LocationType>{LocationType::Home})};
    }
}

bool World::use_migration_rules() const
{
    return m_use_migration_rules;
}

TestingStrategy& World::get_testing_strategy()
{
    return m_testing_strategy;
}

const TestingStrategy& World::get_testing_strategy() const
{
    return m_testing_strategy;
}

} // namespace abm
} // namespace mio
