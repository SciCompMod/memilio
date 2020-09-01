#include "epidemiology/abm/world.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/abm/location.h"
#include "epidemiology/abm/random_number_generator.h"
#include "epidemiology/abm/migration_rules.h"
#include "epidemiology/utils/stl_util.h"

namespace epi
{

Location& World::add_location(LocationType type)
{
    m_locations.push_back(std::make_unique<Location>(type));
    return *m_locations.back();
}

Person& World::add_person(Location& location, InfectionState state, AbmAgeGroup age)
{
    m_persons.push_back(std::make_unique<Person>(location, state, age));
    auto& person = *m_persons.back();
    location.add_person(person);
    return person;
}

void World::evolve(double dt)
{
    begin_step(dt);
    interaction(dt);
    migration(dt);
}

void World::interaction(double dt)
{
    for (auto&& person : m_persons) {
        person->interact(dt, m_infection_parameters);
    }
}

void World::migration(double dt)
{
    using migration_rule   = LocationType (*)(const Person&, double, double, const MigrationParameters&);
    migration_rule rules[] = {&random_migration};
    for (auto& person : m_persons) {
        for (auto rule : rules) {
            auto target_type = rule(*person, 0, dt, m_migration_parameters);
            auto target      = std::find_if(m_locations.begin(), m_locations.end(), [target_type](auto& location) {
                return location->get_type() == target_type;
            });
            if (target != m_locations.end() && target->get() != &person->get_location()) {
                person->migrate_to(**target);
            }
        }
    }
}

void World::begin_step(double dt)
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