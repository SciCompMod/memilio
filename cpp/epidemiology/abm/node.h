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

/**
 * all locations in the simulated world where persons gather.
 */
class Node
{
public:
    Node(NodeType type);
    /** a person interacts with the population at this node, may change infection state */
    InfectionState interact(const Person& person, double dt) const;
    /** add a person to this node */
    void add_person(Person& person);
    /** remove a person from this node */
    void remove_person(Person& person);
    /** one person in this node changed state */
    void changed_state(const Person& person, InfectionState old_state);
    /** prepare the node for the next simulation step */
    void begin_step(double dt);
    /** number of persons in one infection state */
    int get_subpopulation(InfectionState s) const;
    /** number of persons in all states (states are indices) */
    Eigen::Ref<const Eigen::VectorXi> get_subpopulations() const;

    LocalInfectionParameters& get_infection_parameters()
    {
        return m_parameters;
    }
    const LocalInfectionParameters& get_infection_parameters() const
    {
        return m_parameters;
    }

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
