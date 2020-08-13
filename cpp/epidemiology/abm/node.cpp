#include "epidemiology/abm/node.h"
#include "epidemiology/abm/person.h"

namespace epi
{

Node::Node(NodeType type)
    : m_type(type)
{
//initialize the subpopulations of node
//TODO: How to iterate over all enum values?
	set_subpopulation(epi::State::Susceptible,0);
	set_subpopulation(epi::State::Exposed,0);
	set_subpopulation(epi::State::Recovered,0);
	set_subpopulation(epi::State::Infected,0);

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
    State s = p.get_state();
    change_subpopulation(s, +1);
}

void Node::remove_person(Person& p)
{
    --m_num_persons;
    State s = p.get_state();
    change_subpopulation(s, -1);
}

void Node::change_subpopulation(State& s, int delta){
    subpopulations[s] += delta;
}

void Node::set_subpopulation(State s, int v){
    subpopulations[s] = v;
}

int Node::get_subpopulation(State& s){
    return subpopulations[s];
}
} // namespace epi
