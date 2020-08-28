#include "epidemiology/abm/world.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/abm/location.h"
#include "epidemiology/abm/random_number_generator.h"
#include "epidemiology/stl_util.h"

namespace epi
{

Location& World::add_location(LocationType type)
{
    m_locations.push_back(std::make_unique<Location>(type));
    return *m_locations.back();
}

Person& World::add_person(Location& location, InfectionState state)
{
    m_persons.push_back(std::make_unique<Person>(location, state));
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
    //migrate a few times per day to a random different location
    auto u = std::uniform_real_distribution<double>()(thread_local_rng());
    if (u < 2 * dt) {
        for (auto&& person : m_persons) {
            auto idx_distribution    = std::uniform_int_distribution<size_t>(0, m_locations.size() - 2);
            auto random_location_idx = idx_distribution(thread_local_rng());
            //exclude the current location from the random selection
            if (contains(m_locations.begin(), m_locations.begin() + random_location_idx, [&person](auto& e) {
                    return e.get() == &person->get_location();
                })) {
                ++random_location_idx;
            }
            person->migrate_to(*m_locations[random_location_idx]);
        }
    }
    //TODO: more sensible migration probabilities
    //TODO: migration by complex rules
}

void World::begin_step(double dt)
{
    for (auto&& location : m_locations) {
        location->begin_step(dt);
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