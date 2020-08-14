#include "epidemiology/abm/node.h"
#include "epidemiology/abm/person.h"

namespace epi
{

Node::Node(NodeType type)
    : m_type(type)
    , m_subpopulations{}
{
}

InfectionState Node::next_infection_state(const Person& person, double dt) const
{
    //TODO: markov transition
    return person.get_infection_state();
}

void Node::end_migration(double dt)
{
    //TODO: markov coefficients update
}

void Node::add_person(Person& p)
{
    ++m_num_persons;
    InfectionState s = p.get_infection_state();
    change_subpopulation(s, +1);
}

void Node::remove_person(Person& p)
{
    --m_num_persons;
    InfectionState s = p.get_infection_state();
    change_subpopulation(s, -1);
}

void Node::change_subpopulation(InfectionState s, int delta)
{
    m_subpopulations[size_t(s)] += delta;
}

int Node::get_subpopulation(InfectionState s)
{
    return m_subpopulations[size_t(s)];
}
} // namespace epi
