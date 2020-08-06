#ifndef EPI_ABM_NODE_H
#define EPI_ABM_NODE_H

#include "epidemiology/abm/state.h"
#include "epidemiology/abm/node_type.h"

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

private:
    NodeType m_type;
    int m_num_persons = 0;
    //markov parameters, ...
};
} // namespace epi

#endif