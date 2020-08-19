#ifndef EPI_ABM_WORLD_H
#define EPI_ABM_WORLD_H

#include "epidemiology/abm/parameters.h"
#include "epidemiology/abm/node.h"
#include "epidemiology/abm/person.h"
#include <vector>
#include <memory>

namespace epi
{

class World
{
public:
    World(const GlobalInfectionParameters& params = {})
        : m_infection_parameters(params)
    {
    }

    //move-only for stable references of persons/nodes
    World(World&& other) = default;
    World& operator=(World&& other) = default;
    World(const World&) = delete;
    World& operator=(const World&) = delete;

    void evolve(double dt);
    Node& add_node(NodeType type);
    Person& add_person(Node& node, InfectionState state);
    void end_step(double dt);

private:
    void interaction(double dt);
    void migration(double dt);

    std::vector<std::unique_ptr<Person>> m_persons;
    std::vector<std::unique_ptr<Node>> m_nodes;
    GlobalInfectionParameters m_infection_parameters;
};

} // namespace epi

#endif