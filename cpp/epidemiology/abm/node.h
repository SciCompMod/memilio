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
    InfectionState next_infection_state(const Person& person, double dt) const;
    void add_person(Person& person);
    void remove_person(Person& person);
    void end_migration(double dt);
    int get_subpopulation(InfectionState s);

private:
    void change_subpopulation(InfectionState s, int delta);

private:
    NodeType m_type;
    int m_num_persons = 0;
    std::array<int, size_t(InfectionState::Count)> m_subpopulations;
};
} // namespace epi

#endif
