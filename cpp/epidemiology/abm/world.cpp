/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: Daniel Abele, Majid Abedi
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
#include "epidemiology/abm/world.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/abm/location.h"
#include "epidemiology/abm/migration_rules.h"
#include "epidemiology/utils/random_number_generator.h"
#include "epidemiology/utils/stl_util.h"

namespace epi
{

LocationId World::add_location(LocationType type)
{
    auto& locations = m_locations[(uint32_t)type];
    uint32_t index =static_cast<uint32_t> (locations.size());
    locations.emplace_back(Location(type, index));
    return {index, type};
}

Person& World::add_person(LocationId id, InfectionState state, AbmAgeGroup age)
{
    m_persons.push_back(std::make_unique<Person>(id, state, age, m_infection_parameters, m_testing_parameters));
    auto& person = *m_persons.back();
    get_location(person).add_person(person);
    return person;
}

void World::evolve(TimePoint t, TimeSpan dt)
{
    begin_step(t, dt);
    interaction(t, dt);
    migration(t, dt);
}

void World::interaction(TimePoint /*t*/, TimeSpan dt)
{
    for (auto&& person : m_persons) {
        auto& loc = get_location(*person);
        person->interact(dt, m_infection_parameters, loc, m_testing_parameters);
    }
}

void World::migration(TimePoint t, TimeSpan dt)
{
    using migration_rule   = LocationType (*)(const Person&, TimePoint, TimeSpan, const AbmMigrationParameters&);
    const std::pair<migration_rule, std::vector<LocationType>> rules [] = {
        std::make_pair(&return_home_when_recovered, std::vector<LocationType> {LocationType::Home, LocationType::Hospital}), //assumption: if there is an ICU, there is also an hospital
        std::make_pair(&go_to_hospital, std::vector<LocationType> {LocationType::Home, LocationType::Hospital}),
        std::make_pair(&go_to_icu, std::vector<LocationType> {LocationType::Hospital, LocationType::ICU}),
        std::make_pair(&go_to_school, std::vector<LocationType> {LocationType::School, LocationType::Home}),
        std::make_pair(&go_to_work, std::vector<LocationType> {LocationType::Home, LocationType::Work}),
        std::make_pair(&go_to_shop, std::vector<LocationType> {LocationType::Home, LocationType::BasicsShop}),
        std::make_pair(&go_to_event, std::vector<LocationType> {LocationType::Home, LocationType::SocialEvent})
    };
    for (auto&& person : m_persons) {
        for (auto rule : rules) {

            //check if transition rule can be applied
            const auto& locs = rule.second;
            bool nonempty = !locs.empty(); 
            nonempty = std::all_of(locs.begin(), locs.end(), [this](LocationType type) {return !m_locations[(uint32_t)type].empty();});
            if (nonempty){
                auto target_type = rule.first(*person, t, dt, m_migration_parameters);
                Location* target      = find_location(target_type, *person);
                if (target != & get_location(*person)) {
                    target->get_testing_scheme().run_scheme(*person, m_testing_parameters);
                    person->migrate_to(get_location(*person), *target);
                    break;
                }
            }
            
        }
    }
}

void World::begin_step(TimePoint /*t*/, TimeSpan dt)
{
    for (auto&& locations : m_locations) {
        for (auto& location : locations){
            location.begin_step(dt, m_infection_parameters);
        }
    }
}

auto World::get_locations() const -> Range<std::pair<std::vector<std::vector<Location>>::const_iterator, std::vector<std::vector<Location>>::const_iterator>>
{
    return std::make_pair(m_locations.begin(), m_locations.end());
}

auto World::get_persons() const -> Range<std::pair<ConstPersonIterator, ConstPersonIterator>>
{
    return std::make_pair(ConstPersonIterator(m_persons.begin()), ConstPersonIterator(m_persons.end()));
}

const Location& World::get_individualized_location(LocationId id) const
{
    return m_locations[(uint32_t)id.type][id.index];
}

Location& World::get_individualized_location(LocationId id)
{
    return m_locations[(uint32_t)id.type][id.index];
}

const Location& World::get_location(const Person& person) const
{
    return get_individualized_location(person.get_location_id());
}

Location* World::find_location(LocationType type, const Person& person)
{ 
    auto index = person.get_assigned_location_index(type);
    assert(index != INVALID_LOCATION_INDEX && "unexpected error.");
    return &get_individualized_location({index, type});    
}

Location& World::get_location(const Person& person)
{
    return get_individualized_location(person.get_location_id());
}

int World::get_subpopulation_combined(InfectionState s, LocationType type) const
{
    auto& locs = m_locations[(uint32_t)type];
    return std::accumulate(locs.begin(), locs.end(), 0, [&](int running_sum, const Location& loc)
           { return running_sum + loc.get_subpopulation(s); });
}

AbmMigrationParameters& World::get_migration_parameters(){
    return m_migration_parameters;
} 

const AbmMigrationParameters& World::get_migration_parameters() const{
    return m_migration_parameters;
} 

GlobalInfectionParameters& World::get_global_infection_parameters(){
    return m_infection_parameters;
}

const GlobalInfectionParameters& World::get_global_infection_parameters() const{
    return m_infection_parameters;
}

GlobalTestingParameters& World::get_global_testing_parameters(){
    return m_testing_parameters;
}

const GlobalTestingParameters& World::get_global_testing_parameters() const{
    return m_testing_parameters;
}

} // namespace epi
