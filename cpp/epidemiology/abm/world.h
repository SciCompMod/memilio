#ifndef EPI_ABM_WORLD_H
#define EPI_ABM_WORLD_H

#include "epidemiology/abm/node.h"
#include "epidemiology/abm/person.h"
#include <vector>
#include <random>

namespace epi
{

class World
{
public:
    void evolve(double dt);
    Node& add_node(NodeType type);
    Person& add_person(Node& node, State state);
    void end_migration(double dt);
    const std::vector<Person>& get_persons() const;
    const std::vector<Node>& get_nodes() const;

private:
    void interaction(double dt);
    void migration(double dt);

    std::vector<Person> m_persons;
    std::vector<Node> m_nodes;
    std::mt19937_64 m_rng; //TODO: seed
};

} // namespace epi

#endif