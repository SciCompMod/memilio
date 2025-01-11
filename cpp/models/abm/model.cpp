/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "abm/model.h"
#include "abm/location_id.h"
#include "abm/location_type.h"
#include "abm/intervention_type.h"
#include "abm/model_functions.h"
#include "abm/person.h"
#include "abm/location.h"
#include "abm/mobility_rules.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/mioomp.h"
#include "memilio/utils/stl_util.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>

namespace mio
{
namespace abm
{

LocationId Model::add_location(LocationType type, uint32_t num_cells)
{
    LocationId id{static_cast<uint32_t>(m_locations.size())};
    m_locations.emplace_back(type, id, parameters.get_num_groups(), num_cells);
    m_has_locations[size_t(type)] = true;

    // mark caches for rebuild
    m_is_local_population_cache_valid = false;
    m_are_exposure_caches_valid       = false;
    m_exposure_caches_need_rebuild    = true;

    return id;
}

PersonId Model::add_person(const LocationId id, AgeGroup age)
{
    return add_person(Person(m_rng, get_location(id).get_type(), id, age));
}

PersonId Model::add_person(Person&& person)
{
    assert(person.get_location() != LocationId::invalid_id() && "Added Person's location must be valid.");
    assert(person.get_location() < LocationId((uint32_t)m_locations.size()) &&
           "Added Person's location is not in Model.");
    assert(person.get_age() < (AgeGroup)parameters.get_num_groups() && "Added Person's AgeGroup is too large.");

    PersonId new_id{static_cast<uint32_t>(m_persons.size())};
    m_persons.emplace_back(person, new_id);
    auto& new_person = m_persons.back();
    new_person.set_assigned_location(LocationType::Cemetery, m_cemetery_id);

    if (m_is_local_population_cache_valid) {
        ++m_local_population_cache[new_person.get_location().get()];
    }
    return new_id;
}

void Model::evolve(TimePoint t, TimeSpan dt)
{
    begin_step(t, dt);
    log_info("ABM Model interaction.");
    interaction(t, dt);
    log_info("ABM Model mobility.");
    perform_mobility(t, dt);
}

void Model::interaction(TimePoint t, TimeSpan dt)
{
    const uint32_t num_persons = static_cast<uint32_t>(m_persons.size());
    PRAGMA_OMP(parallel for)
    for (uint32_t person_id = 0; person_id < num_persons; ++person_id) {
        interact(person_id, t, dt);
    }
}

void Model::perform_mobility(TimePoint t, TimeSpan dt)
{
    const uint32_t num_persons = static_cast<uint32_t>(m_persons.size());
    PRAGMA_OMP(parallel for)
    for (uint32_t person_id = 0; person_id < num_persons; ++person_id) {
        Person& person    = m_persons[person_id];
        auto personal_rng = PersonalRandomNumberGenerator(m_rng, person);

        auto try_mobility_rule = [&](auto rule) -> bool {
            // run mobility rule and check if change of location can actually happen
            auto target_type                  = rule(personal_rng, person, t, dt, parameters);
            const Location& target_location   = get_location(find_location(target_type, person_id));
            const LocationId current_location = person.get_location();

            // the Person cannot move if they do not wear mask as required at targeted location
            if (target_location.is_mask_required() && !person.is_compliant(personal_rng, InterventionType::Mask)) {
                return false;
            }
            // the Person cannot move if the capacity of targeted Location is reached
            if (target_location.get_id() == current_location ||
                get_number_persons(target_location.get_id()) >= target_location.get_capacity().persons) {
                return false;
            }
            // the Person cannot move if the performed TestingStrategy is positive
            if (!m_testing_strategy.run_strategy(personal_rng, person, target_location, t)) {
                return false;
            }
            // update worn mask to target location's requirements
            if (target_location.is_mask_required()) {
                // if the current MaskProtection level is lower than required, the Person changes mask
                if (parameters.get<MaskProtection>()[person.get_mask().get_type()] <
                    parameters.get<MaskProtection>()[target_location.get_required_mask()]) {
                    person.set_mask(target_location.get_required_mask(), t);
                }
            }
            else {
                person.set_mask(MaskType::None, t);
            }
            // all requirements are met, move to target location
            change_location(person_id, target_location.get_id());
            return true;
        };

        // run mobility rules one after the other if the corresponding location type exists
        // shortcutting of bool operators ensures the rules stop after the first rule is applied
        if (m_use_mobility_rules) {
            (has_locations({LocationType::Cemetery}) && try_mobility_rule(&get_buried)) ||
                (has_locations({LocationType::Home}) && try_mobility_rule(&return_home_when_recovered)) ||
                (has_locations({LocationType::Hospital}) && try_mobility_rule(&go_to_hospital)) ||
                (has_locations({LocationType::ICU}) && try_mobility_rule(&go_to_icu)) ||
                (has_locations({LocationType::School, LocationType::Home}) && try_mobility_rule(&go_to_school)) ||
                (has_locations({LocationType::Work, LocationType::Home}) && try_mobility_rule(&go_to_work)) ||
                (has_locations({LocationType::BasicsShop, LocationType::Home}) && try_mobility_rule(&go_to_shop)) ||
                (has_locations({LocationType::SocialEvent, LocationType::Home}) && try_mobility_rule(&go_to_event)) ||
                (has_locations({LocationType::Home}) && try_mobility_rule(&go_to_quarantine));
        }
        else {
            // no daily routine mobility, just infection related
            (has_locations({LocationType::Cemetery}) && try_mobility_rule(&get_buried)) ||
                (has_locations({LocationType::Home}) && try_mobility_rule(&return_home_when_recovered)) ||
                (has_locations({LocationType::Hospital}) && try_mobility_rule(&go_to_hospital)) ||
                (has_locations({LocationType::ICU}) && try_mobility_rule(&go_to_icu)) ||
                (has_locations({LocationType::Home}) && try_mobility_rule(&go_to_quarantine));
        }
    }

    // check if a person makes a trip
    bool weekend     = t.is_weekend();
    size_t num_trips = m_trip_list.num_trips(weekend);

    for (; m_trip_list.get_current_index() < num_trips &&
           m_trip_list.get_next_trip_time(weekend).seconds() < (t + dt).time_since_midnight().seconds();
         m_trip_list.increase_index()) {
        auto& trip        = m_trip_list.get_next_trip(weekend);
        auto& person      = get_person(trip.person_id);
        auto personal_rng = PersonalRandomNumberGenerator(m_rng, person);
        // skip the trip if the person is in quarantine or is dead
        if (person.is_in_quarantine(t, parameters) || person.get_infection_state(t) == InfectionState::Dead) {
            continue;
        }
        auto& target_location = get_location(trip.destination);
        // skip the trip if the Person wears mask as required at targeted location
        if (target_location.is_mask_required() && !person.is_compliant(personal_rng, InterventionType::Mask)) {
            continue;
        }
        // skip the trip if the performed TestingStrategy is positive
        if (!m_testing_strategy.run_strategy(personal_rng, person, target_location, t)) {
            continue;
        }
        // all requirements are met, move to target location
        change_location(person.get_id(), target_location.get_id(), trip.trip_mode);
        // update worn mask to target location's requirements
        if (target_location.is_mask_required()) {
            // if the current MaskProtection level is lower than required, the Person changes mask
            if (parameters.get<MaskProtection>()[person.get_mask().get_type()] <
                parameters.get<MaskProtection>()[target_location.get_required_mask()]) {
                person.set_mask(target_location.get_required_mask(), t);
            }
        }
        else {
            person.set_mask(MaskType::None, t);
        }
    }
    if (((t).days() < std::floor((t + dt).days()))) {
        m_trip_list.reset_index();
    }
}

void Model::build_compute_local_population_cache() const
{
    PRAGMA_OMP(single)
    {
        const size_t num_locations = m_locations.size();
        const size_t num_persons   = m_persons.size();
        m_local_population_cache.resize(num_locations);
        PRAGMA_OMP(taskloop)
        for (size_t i = 0; i < num_locations; i++) {
            m_local_population_cache[i] = 0;
        } // implicit taskloop barrier
        PRAGMA_OMP(taskloop)
        for (size_t i = 0; i < num_persons; i++) {
            ++m_local_population_cache[m_persons[i].get_location().get()];
        } // implicit taskloop barrier
    } // implicit single barrier
}

void Model::build_exposure_caches()
{
    PRAGMA_OMP(single)
    {
        const size_t num_locations = m_locations.size();
        m_air_exposure_rates_cache.resize(num_locations);
        m_contact_exposure_rates_cache.resize(num_locations);
        PRAGMA_OMP(taskloop)
        for (size_t i = 0; i < num_locations; i++) {
            m_air_exposure_rates_cache[i].resize({CellIndex(m_locations[i].get_cells().size()), VirusVariant::Count});
            m_contact_exposure_rates_cache[i].resize({CellIndex(m_locations[i].get_cells().size()), VirusVariant::Count,
                                                      AgeGroup(parameters.get_num_groups())});
        } // implicit taskloop barrier
        m_are_exposure_caches_valid    = false;
        m_exposure_caches_need_rebuild = false;
    } // implicit single barrier
}

void Model::compute_exposure_caches(TimePoint t, TimeSpan dt)
{
    PRAGMA_OMP(single)
    {
        // if cache shape was changed (e.g. by add_location), rebuild it
        if (m_exposure_caches_need_rebuild) {
            build_exposure_caches();
        }
        // use these const values to help omp recognize that the for loops are bounded
        const auto num_locations = m_locations.size();
        const auto num_persons   = m_persons.size();

        // 1) reset all cached values
        // Note: we cannot easily reuse values, as they are time dependent (get_infection_state)
        PRAGMA_OMP(taskloop)
        for (size_t i = 0; i < num_locations; ++i) {
            mio::abm::adjust_contact_rates(m_locations[i], parameters.get_num_groups());
            const auto index         = i;
            auto& local_air_exposure = m_air_exposure_rates_cache[index];
            std::for_each(local_air_exposure.begin(), local_air_exposure.end(), [](auto& r) {
                r = 0.0;
            });
            auto& local_contact_exposure = m_contact_exposure_rates_cache[index];
            std::for_each(local_contact_exposure.begin(), local_contact_exposure.end(), [](auto& r) {
                r = 0.0;
            });
        } // implicit taskloop barrier
        // here is an implicit (and needed) barrier from parallel for

        // 2) add all contributions from each person
        PRAGMA_OMP(taskloop)
        for (size_t i = 0; i < num_persons; ++i) {
            const Person& person = m_persons[i];
            const auto location  = person.get_location().get();
            mio::abm::add_exposure_contribution(m_air_exposure_rates_cache[location],
                                                m_contact_exposure_rates_cache[location], person,
                                                get_location(person.get_id()), t, dt);
        } // implicit taskloop barrier
        //normalize cached exposure rates
        for (size_t i = 0; i < num_locations; ++i) {
            for (auto age_group = AgeGroup(0); age_group < AgeGroup(parameters.get_num_groups()); ++age_group) {
                auto num_persons_in_location =
                    std::count_if(m_persons.begin(), m_persons.end(), [age_group](Person& p) {
                        return p.get_age() == age_group;
                    });
                if (num_persons_in_location > 0) {
                    for (auto& v : m_contact_exposure_rates_cache[i].slice(AgeGroup(age_group))) {
                        v = v / num_persons_in_location;
                    }
                }
            }
        }
    } // implicit single barrier
}

void Model::begin_step(TimePoint t, TimeSpan dt)
{
    m_testing_strategy.update_activity_status(t);

    if (!m_is_local_population_cache_valid) {
        build_compute_local_population_cache();
        m_is_local_population_cache_valid = true;
    }
    compute_exposure_caches(t, dt);
    m_are_exposure_caches_valid = true;
}

auto Model::get_locations() const -> Range<std::pair<ConstLocationIterator, ConstLocationIterator>>
{
    return std::make_pair(m_locations.cbegin(), m_locations.cend());
}
auto Model::get_locations() -> Range<std::pair<LocationIterator, LocationIterator>>
{
    return std::make_pair(m_locations.begin(), m_locations.end());
}

auto Model::get_persons() const -> Range<std::pair<ConstPersonIterator, ConstPersonIterator>>
{
    return std::make_pair(m_persons.cbegin(), m_persons.cend());
}

auto Model::get_persons() -> Range<std::pair<PersonIterator, PersonIterator>>
{
    return std::make_pair(m_persons.begin(), m_persons.end());
}

LocationId Model::find_location(LocationType type, const PersonId person) const
{
    auto location_id = get_person(person).get_assigned_location(type);
    assert(location_id != LocationId::invalid_id() && "The person has no assigned location of that type.");
    return location_id;
}

size_t Model::get_subpopulation_combined(TimePoint t, InfectionState s) const
{
    return std::accumulate(m_locations.begin(), m_locations.end(), (size_t)0,
                           [t, s, this](size_t running_sum, const Location& loc) {
                               return running_sum + get_subpopulation(loc.get_id(), t, s);
                           });
}

size_t Model::get_subpopulation_combined_per_location_type(TimePoint t, InfectionState s, LocationType type) const
{
    return std::accumulate(
        m_locations.begin(), m_locations.end(), (size_t)0, [t, s, type, this](size_t running_sum, const Location& loc) {
            return loc.get_type() == type ? running_sum + get_subpopulation(loc.get_id(), t, s) : running_sum;
        });
}

TripList& Model::get_trip_list()
{
    return m_trip_list;
}

const TripList& Model::get_trip_list() const
{
    return m_trip_list;
}

void Model::use_mobility_rules(bool param)
{
    m_use_mobility_rules = param;
}

bool Model::use_mobility_rules() const
{
    return m_use_mobility_rules;
}

TestingStrategy& Model::get_testing_strategy()
{
    return m_testing_strategy;
}

const TestingStrategy& Model::get_testing_strategy() const
{
    return m_testing_strategy;
}

} // namespace abm
} // namespace mio
