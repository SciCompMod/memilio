#include "epidemiology/abm/node.h"
#include "epidemiology/abm/person.h"

namespace epi
{

Node::Node(NodeType type)
    : m_type(type)
{
}

State Node::next_state(const Person& person, double dt) const
{
    //TODO: markov transition
    return person.get_state();
}

void Node::end_migration(double dt)
{
    //TODO: markov coefficients update
}

void Node::add_person(Person& p)
{
    ++m_num_persons;
}

void Node::remove_person(Person& p)
{
    --m_num_persons;
}

} // namespace epi