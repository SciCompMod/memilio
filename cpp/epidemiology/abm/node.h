#ifndef EPI_ABM_NODE_H
#define EPI_ABM_NODE_H

#include "epidemiology/abm/state.h"
#include "epidemiology/abm/node_type.h"
#include <map>

namespace epi
{
class Person;

class Node
{
public:
    Node(NodeType type);
    State next_state(const Person& person, double dt) const;
    void add_person(Person& person);
    void remove_person(Person& person);
    void end_migration(double dt);
    void change_subpopulation(State& s, int delta);
    int get_subpopulation(State& s);
    void set_subpopulation(State s, int v);

private:
    NodeType m_type;
    int m_num_persons = 0;
    std::map<State, int> subpopulations;
};
} // namespace epi

#endif
