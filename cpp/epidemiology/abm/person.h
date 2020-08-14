#ifndef EPI_ABM_PERSON_H
#define EPI_ABM_PERSON_H

#include "epidemiology/abm/state.h"
#include <functional>

namespace epi
{

class Node;

class Person
{
public:
    Person(Node& node, InfectionState state);
    void interact(double dt);
    void migrate_to(Node& node); //could also be migrate() with the decision made internally

    InfectionState get_infection_state() const
    {
        return m_state;
    }
    const Node& get_node() const
    {
        return m_node;
    }

private:
    Node& get_node()
    {
        return m_node;
    }

    InfectionState m_state;
    std::reference_wrapper<Node> m_node;
    //age, times, ...
};

} // namespace epi

#endif