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
#include "abm/person.h"
#include "abm/location.h"
#include "abm/migration_rules.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/mioomp.h"
#include "memilio/utils/stl_util.h"

namespace mio
{
namespace abm
{

LocationId World::add_location(LocationType type, uint32_t num_cells)
{
    LocationId id = {static_cast<uint32_t>(m_locations.size()), type};
    m_locations.emplace_back(id, parameters.get_num_groups(), num_cells);
    m_has_locations[size_t(type)] = true;

    // mark caches for rebuild
    m_local_population_cache.invalidate();
    m_air_exposure_rates_cache.invalidate();
    m_contact_exposure_rates_cache.invalidate();
    m_exposure_rates_need_rebuild = true;

    return id;
}

PersonId World::add_person(const LocationId id, AgeGroup age)
{
    return add_person(Person(m_rng, id, age));
}

PersonId World::add_person(Person&& person)
{
    assert(person.get_location().index != INVALID_LOCATION_INDEX);
    assert(person.get_location().index < m_locations.size());
    assert(person.get_age().get() < parameters.get_num_groups());

    PersonId new_id = static_cast<PersonId>(m_persons.size());
    m_persons.emplace_back(person, new_id);
    auto& new_person = m_persons.back();
    new_person.set_assigned_location(m_cemetery_id);

    if (m_local_population_cache.is_valid()) {
        ++m_local_population_cache.write()[new_person.get_location().index];
    }
    return new_id;
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
    const auto num_persons = m_persons.size();
    PRAGMA_OMP(parallel for)
    for (auto i = size_t(0); i < num_persons; ++i) {
        interact(i, t, dt);
    }
}

void World::migration(TimePoint t, TimeSpan dt)
{
    const auto num_persons = m_persons.size();
    PRAGMA_OMP(parallel for)
    for (auto i = size_t(0); i < num_persons; ++i) {
        Person& person    = m_persons[i];
        auto personal_rng = PersonalRandomNumberGenerator(m_rng, person);

        auto try_migration_rule = [&](auto rule) -> bool {
            //run migration rule and check if migration can actually happen
            auto target_type                  = rule(personal_rng, person, t, dt, parameters);
            const Location& target_location   = get_location(find_location(target_type, i));
            const LocationId current_location = person.get_location();
            if (m_testing_strategy.run_strategy(personal_rng, person, target_location, t)) {
                if (target_location.get_id() != current_location &&
                    get_number_persons(target_location.get_id()) < target_location.get_capacity().persons) {
                    bool wears_mask = person.apply_mask_intervention(personal_rng, target_location);
                    if (wears_mask) {
                        migrate(i, target_location.get_id());
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
            auto& person      = get_person(trip.person_id);
            auto personal_rng = PersonalRandomNumberGenerator(m_rng, person);
            if (!person.is_in_quarantine(t, parameters) && person.get_infection_state(t) != InfectionState::Dead) {
                auto& target_location = get_location(trip.migration_destination);
                if (m_testing_strategy.run_strategy(personal_rng, person, target_location, t)) {
                    person.apply_mask_intervention(personal_rng, target_location);
                    migrate(person.get_person_id(), target_location.get_id(), trip.trip_mode);
                }
            }
            m_trip_list.increase_index();
        }
    }
    if (((t).days() < std::floor((t + dt).days()))) {
        m_trip_list.reset_index();
    }
}

void World::build_compute_local_population_cache() const
{
    PRAGMA_OMP(single)
    {
        const auto num_locations = m_locations.size();
        const auto num_persons   = m_persons.size();
        m_local_population_cache.write().resize(num_locations);
        PRAGMA_OMP(taskloop)
        for (size_t i = 0; i < num_locations; i++) {
            m_local_population_cache.write()[i] = 0.;
        } // implicit taskloop barrier
        PRAGMA_OMP(taskloop)
        for (size_t i = 0; i < num_persons; i++) {
            ++m_local_population_cache.write()[m_persons[i].get_location().index];
        } // implicit taskloop barrier
    } // implicit single barrier
}

void World::build_exposure_caches()
{
    PRAGMA_OMP(single)
    {
        const size_t num_locations = m_locations.size();
        m_air_exposure_rates_cache.write().resize(num_locations);
        m_contact_exposure_rates_cache.write().resize(num_locations);
        PRAGMA_OMP(taskloop)
        for (size_t i = 0; i < num_locations; i++) {
            m_air_exposure_rates_cache.write()[i].resize(
                {CellIndex(m_locations[i].get_cells().size()), VirusVariant::Count});
            m_contact_exposure_rates_cache.write()[i].resize({CellIndex(m_locations[i].get_cells().size()),
                                                              VirusVariant::Count,
                                                              AgeGroup(parameters.get_num_groups())});
        } // implicit taskloop barrier
        m_air_exposure_rates_cache.invalidate();
        m_contact_exposure_rates_cache.invalidate();
        m_exposure_rates_need_rebuild = false;
    } // implicit single barrier
}

void World::compute_exposure_caches(TimePoint t, TimeSpan dt)
{
    PRAGMA_OMP(single)
    {
        // if cache shape was changed (e.g. by add_location), rebuild it
        if (m_exposure_rates_need_rebuild) {
            build_exposure_caches();
        }
        // use these const values to help omp recognize that the for loops are bounded
        const auto num_locations = m_locations.size();
        const auto num_persons   = m_persons.size();

        // 1) reset all cached values
        // Note: we cannot easily reuse values, as they are time dependant (get_infection_state)
        PRAGMA_OMP(taskloop)
        for (size_t i = 0; i < num_locations; ++i) {
            const auto index         = i;
            auto& local_air_exposure = m_air_exposure_rates_cache.write()[index];
            std::for_each(local_air_exposure.begin(), local_air_exposure.end(), [](auto& r) {
                r = 0.0;
            });
            auto& local_contact_exposure = m_contact_exposure_rates_cache.write()[index];
            std::for_each(local_contact_exposure.begin(), local_contact_exposure.end(), [](auto& r) {
                r = 0.0;
            });
        } // implicit taskloop barrier
        // here is an implicit (and needed) barrier from parallel for

        // 2) add all contributions from each person
        PRAGMA_OMP(taskloop)
        for (size_t i = 0; i < num_persons; ++i) {
            const Person& person = m_persons[i];
            const auto location  = person.get_location().index;
            mio::abm::add_exposure_contribution(m_air_exposure_rates_cache.write()[location],
                                                m_contact_exposure_rates_cache.write()[location], person,
                                                get_location(person.get_person_id()), t, dt);
        } // implicit taskloop barrier
    } // implicit single barrier
}

void World::begin_step(TimePoint t, TimeSpan dt)
{
    m_testing_strategy.update_activity_status(t);

    if (!m_local_population_cache.is_valid()) {
        build_compute_local_population_cache();
        m_local_population_cache.validate();
    }
    compute_exposure_caches(t, dt);
    m_air_exposure_rates_cache.validate();
    m_contact_exposure_rates_cache.validate();
}

auto World::get_locations() const -> Range<std::pair<ConstLocationIterator, ConstLocationIterator>>
{
    return std::make_pair(m_locations.cbegin(), m_locations.cend());
}
auto World::get_locations() -> Range<std::pair<LocationIterator, LocationIterator>>
{
    return std::make_pair(m_locations.begin(), m_locations.end());
}

auto World::get_persons() const -> Range<std::pair<ConstPersonIterator, ConstPersonIterator>>
{
    return std::make_pair(m_persons.cbegin(), m_persons.cend());
}

auto World::get_persons() -> Range<std::pair<PersonIterator, PersonIterator>>
{
    return std::make_pair(m_persons.begin(), m_persons.end());
}

LocationId World::find_location(LocationType type, const PersonId person) const
{
    auto index = get_person(person).get_assigned_location_index(type);
    assert(index != INVALID_LOCATION_INDEX && "unexpected error.");
    return {index, type};
}

size_t World::get_subpopulation_combined(TimePoint t, InfectionState s) const
{
    return std::accumulate(m_locations.begin(), m_locations.end(), (size_t)0,
                           [t, s, this](size_t running_sum, const Location& loc) {
                               return running_sum + get_subpopulation(loc.get_id(), t, s);
                           });
}

size_t World::get_subpopulation_combined_per_location_type(TimePoint t, InfectionState s, LocationType type) const
{
    return std::accumulate(
        m_locations.begin(), m_locations.end(), (size_t)0, [t, s, type, this](size_t running_sum, const Location& loc) {
            return loc.get_type() == type ? running_sum + get_subpopulation(loc.get_id(), t, s) : running_sum;
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
