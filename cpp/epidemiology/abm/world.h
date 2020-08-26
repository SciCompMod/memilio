#ifndef EPI_ABM_WORLD_H
#define EPI_ABM_WORLD_H

#include "epidemiology/abm/parameters.h"
#include "epidemiology/abm/node.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/pointer_dereferencing_iterator.h"
#include "epidemiology/stl_util.h"
#include <vector>
#include <memory>

namespace epi
{

/**
 * The world of the simulation.
 * consists of Nodes (Locations) and Persons (Actors)
 */
class World
{
public:
    using NodeIterator        = PointerDereferencingIterator<std::vector<std::unique_ptr<Node>>::iterator>;
    using ConstNodeIterator   = PointerDereferencingIterator<std::vector<std::unique_ptr<Node>>::const_iterator>;
    using PersonIterator      = PointerDereferencingIterator<std::vector<std::unique_ptr<Person>>::iterator>;
    using ConstPersonIterator = PointerDereferencingIterator<std::vector<std::unique_ptr<Person>>::const_iterator>;

    World(const GlobalInfectionParameters& params = {})
        : m_infection_parameters(params)
    {
    }

    //move-only for stable references of persons/nodes
    World(World&& other) = default;
    World& operator=(World&& other) = default;
    World(const World&)             = delete;
    World& operator=(const World&) = delete;

    /** prepare the world for the next simulation step */
    void begin_step(double dt);
    /** evolve the world one discrete time step */
    void evolve(double dt);
    /** add a node to the world */
    Node& add_node(NodeType type);
    /** add a person to the world */
    Person& add_person(Node& node, InfectionState state);

    Range<std::pair<ConstNodeIterator, ConstNodeIterator>> get_nodes() const;
    Range<std::pair<ConstPersonIterator, ConstPersonIterator>> get_persons() const;

private:
    void interaction(double dt);
    void migration(double dt);

    std::vector<std::unique_ptr<Person>> m_persons;
    std::vector<std::unique_ptr<Node>> m_nodes;
    GlobalInfectionParameters m_infection_parameters;
};

} // namespace epi

#endif