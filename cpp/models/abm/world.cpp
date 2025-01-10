/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Majid Abedi, Elisabeth Kluth, Carlotta Gerstein, Martin J. Kuehn, David Kerkmann, Khoa Nguyen
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
#include "abm/location_type.h"
#include "abm/mask_type.h"
#include "abm/person.h"
#include "abm/location.h"
#include "abm/migration_rules.h"
#include "abm/infection.h"
#include "abm/vaccine.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/mioomp.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/stl_util.h"

namespace mio
{
namespace abm
{

LocationId World::add_location(LocationType type, uint32_t num_cells)
{
    LocationId id = {static_cast<uint32_t>(m_locations.size()), type};
    m_locations.emplace_back(std::make_unique<Location>(id, parameters.get_num_groups(), num_cells));
    m_has_locations[size_t(type)] = true;
    return id;
}

Person& World::add_person(const LocationId id, AgeGroup age)
{
    assert(age.get() < parameters.get_num_groups());
    uint32_t person_id = static_cast<uint32_t>(m_persons.size());
    m_persons.push_back(std::make_unique<Person>(m_rng, get_individualized_location(id), age, person_id));
    auto& person = *m_persons.back();
    person.set_assigned_location(m_cemetery_id);
    get_individualized_location(id).add_person(person);
    return person;
}

void World::evolve(TimePoint t, TimeSpan dt)
{
    begin_step(t, dt);
    log_info("ABM World interaction.");
    interaction(t, dt);
    log_info("ABM World migration.");
    migration(t, dt);
}

void World::interaction(TimePoint t, TimeSpan dt)
{
    PRAGMA_OMP(parallel for)
    for (auto i = size_t(0); i < m_persons.size(); ++i) {
        auto&& person     = m_persons[i];
        auto personal_rng = Person::RandomNumberGenerator(m_rng, *person);
        person->interact(personal_rng, t, dt, parameters);
    }
}

void World::migration(TimePoint t, TimeSpan dt)
{
    PRAGMA_OMP(parallel for)
    for (auto i = size_t(0); i < m_persons.size(); ++i) {
        auto&& person     = m_persons[i];
        auto personal_rng = Person::RandomNumberGenerator(m_rng, *person);

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
            (has_locations({LocationType::Cemetery}) && try_migration_rule(&get_buried)) ||
                (has_locations({LocationType::Home}) && try_migration_rule(&return_home_when_recovered)) ||
                (has_locations({LocationType::Hospital}) && try_migration_rule(&go_to_hospital)) ||
                (has_locations({LocationType::ICU}) && try_migration_rule(&go_to_icu)) ||
                (has_locations({LocationType::School, LocationType::Home}) && try_migration_rule(&go_to_school)) ||
                (has_locations({LocationType::Work, LocationType::Home}) && try_migration_rule(&go_to_work)) ||
                (has_locations({LocationType::BasicsShop, LocationType::Home}) && try_migration_rule(&go_to_shop)) ||
                (has_locations({LocationType::SocialEvent, LocationType::Home}) && try_migration_rule(&go_to_event)) ||
                (has_locations({LocationType::Home}) && try_migration_rule(&go_to_quarantine));
        }
        else {
            //no daily routine migration, just infection related
            (has_locations({LocationType::Cemetery}) && try_migration_rule(&get_buried)) ||
                (has_locations({LocationType::Home}) && try_migration_rule(&return_home_when_recovered)) ||
                (has_locations({LocationType::Hospital}) && try_migration_rule(&go_to_hospital)) ||
                (has_locations({LocationType::ICU}) && try_migration_rule(&go_to_icu)) ||
                (has_locations({LocationType::Home}) && try_migration_rule(&go_to_quarantine));
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
            auto personal_rng = Person::RandomNumberGenerator(m_rng, *person);
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

void World::begin_step(TimePoint t, TimeSpan dt)
{
    m_testing_strategy.update_activity_status(t);
    PRAGMA_OMP(parallel for)
    for (auto i = size_t(0); i < m_locations.size(); ++i) {
        auto&& location = m_locations[i];
        location->cache_exposure_rates(t, dt, parameters.get_num_groups());
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

const Location& World::find_location(LocationType type, const Person& person) const
{
    auto index = person.get_assigned_location_index(type);
    assert(index != INVALID_LOCATION_INDEX && "unexpected error.");
    return get_individualized_location({index, type});
}

Location& World::find_location(LocationType type, const Person& person)
{
    auto index = person.get_assigned_location_index(type);
    assert(index != INVALID_LOCATION_INDEX && "unexpected error.");
    return get_individualized_location({index, type});
}

size_t World::get_subpopulation_combined(TimePoint t, InfectionState s) const
{
    return std::accumulate(m_locations.begin(), m_locations.end(), (size_t)0,
                           [t, s](size_t running_sum, const std::unique_ptr<Location>& loc) {
                               return running_sum + loc->get_subpopulation(t, s);
                           });
}

size_t World::get_subpopulation_combined_per_location_type(TimePoint t, InfectionState s, LocationType type) const
{
    return std::accumulate(m_locations.begin(), m_locations.end(), (size_t)0,
                           [t, s, type](size_t running_sum, const std::unique_ptr<Location>& loc) {
                               return loc->get_type() == type ? running_sum + loc->get_subpopulation(t, s)
                                                              : running_sum;
                           });
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
