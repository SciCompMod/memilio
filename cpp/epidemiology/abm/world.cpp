#include "epidemiology/abm/world.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/abm/node.h"

namespace epi
{

Node& World::add_node(NodeType type)
{
    m_nodes.emplace_back(type);
    return m_nodes.back();
}

Person& World::add_person(Node& node, State state)
{
    m_persons.emplace_back(node, state);
    auto& person = m_persons.back();
    node.add_person(person);
    return m_persons.back();
}

void World::evolve(double dt)
{
    interaction(dt);
    migration(dt);
    end_migration(dt);
}

void World::interaction(double dt)
{
    for (auto&& person : m_persons) {
        person.interact(dt);
    }
}

void World::migration(double dt)
{
    //TODO: variable migration probabilities
    //TODO: migration by complex rules
    std::uniform_int_distribution<size_t> random_node_idx(0, m_nodes.size() - 1);
    for (auto&& person : m_persons) {
        person.migrate_to(m_nodes[random_node_idx(m_rng)]);
    }
}

void World::end_migration(double dt)
{
    for (auto&& node : m_nodes) {
        node.end_migration(dt);
    }
}

} // namespace epi