#include "epidemiology/abm/world.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/abm/location.h"
#include "epidemiology/abm/migration_rules.h"
#include "epidemiology/utils/random_number_generator.h"
#include "epidemiology/utils/stl_util.h"

namespace epi
{

Location& World::add_location(LocationType type)
{
    m_locations.push_back(std::make_unique<Location>(type));
    return *m_locations.back();
}

Person& World::add_person(Location& location, InfectionState state, Index<AbmAgeGroup> age)
{
    m_persons.push_back(std::make_unique<Person>(location, state, age));
    auto& person = *m_persons.back();
    location.add_person(person);
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
        person->interact(dt, m_infection_parameters);
    }
}

void World::migration(TimePoint t, TimeSpan dt)
{
    using migration_rule   = LocationType (*)(const Person&, TimePoint, TimeSpan, const MigrationParameters&);
    migration_rule rules[] = {&return_home_when_recovered,
                              &go_to_hospital,
                              &go_to_icu,
                              &go_to_school,
                              &go_to_work,
                              &go_to_shop,
                              &go_to_event};
    for (auto& person : m_persons) {
        for (auto rule : rules) {
            auto destination_type = rule(*person, t, dt, m_migration_parameters);
            auto destination_node =
                std::find_if(m_locations.begin(), m_locations.end(), [destination_type](auto& location) {
                    return location->get_type() == destination_type;
                });
            if (destination_node != m_locations.end() && destination_node->get() != &person->get_location()) {
                person->migrate_to(**destination_node);
                break;
            }
        }
    }
}

void World::begin_step(TimePoint /*t*/, TimeSpan dt)
{
    for (auto&& location : m_locations) {
        location->begin_step(dt, m_infection_parameters);
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

} // namespace epi
