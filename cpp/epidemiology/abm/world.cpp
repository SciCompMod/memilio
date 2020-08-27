#include "epidemiology/abm/world.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/abm/node.h"
#include "epidemiology/abm/random_number_generator.h"
#include "epidemiology/stl_util.h"

namespace epi
{

Node& World::add_node(NodeType type)
{
    m_nodes.push_back(std::make_unique<Node>(type));
    return *m_nodes.back();
}

Person& World::add_person(Node& node, InfectionState state)
{
    m_persons.push_back(std::make_unique<Person>(node, state));
    auto& person = *m_persons.back();
    node.add_person(person);
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
    //migrate a few times per day to a random different node
    auto u = std::uniform_real_distribution<double>()(thread_local_rng());
    if (u < 2 * dt) {
        for (auto&& person : m_persons) {
            auto idx_distribution = std::uniform_int_distribution<size_t>(0, m_nodes.size() - 2);
            auto random_node_idx  = idx_distribution(thread_local_rng());
            //exclude the current node from the random selection
            if (contains(m_nodes.begin(), m_nodes.begin() + random_node_idx, [&person](auto& e) {
                    return e.get() == &person->get_node();
                })) {
                ++random_node_idx;
            }
            person->migrate_to(*m_nodes[random_node_idx]);
        }
    }
    //TODO: more sensible migration probabilities
    //TODO: migration by complex rules
}

void World::begin_step(double dt)
{
    for (auto&& node : m_nodes) {
        node->begin_step(dt);
    }
}

auto World::get_nodes() const -> Range<std::pair<ConstNodeIterator, ConstNodeIterator>>
{
    return std::make_pair(ConstNodeIterator(m_nodes.begin()), ConstNodeIterator(m_nodes.end()));
}

auto World::get_persons() const -> Range<std::pair<ConstPersonIterator, ConstPersonIterator>>
{
    return std::make_pair(ConstPersonIterator(m_persons.begin()), ConstPersonIterator(m_persons.end()));
}

} // namespace epi