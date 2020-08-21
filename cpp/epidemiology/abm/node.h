#ifndef EPI_ABM_NODE_H
#define EPI_ABM_NODE_H

#include "epidemiology/abm/parameters.h"
#include "epidemiology/abm/state.h"
#include "epidemiology/abm/node_type.h"
#include <array>
#include <random>
#include <Eigen/Core>

namespace epi
{
class Person;

class Node
{
public:
    Node(NodeType type);
    InfectionState interact(const Person& person, double dt) const;
    void add_person(Person& person);
    void remove_person(Person& person);
    void changed_state(const Person& person, InfectionState old_state);
    void begin_step(double dt);
    int get_subpopulation(InfectionState s) const;
    
private:
    void change_subpopulation(InfectionState s, int delta);

private:
    NodeType m_type;
    int m_num_persons = 0;
    std::array<int, size_t(InfectionState::Count)> m_subpopulations;
    LocalInfectionParameters m_parameters;
    double m_cached_exposure_rate;
};
} // namespace epi

#endif
