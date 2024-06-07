#ifndef ABM_GRAPH_MOBILITY_H
#define ABM_GRAPH_MOBILITY_H

#include "abm/simulation.h"
#include "abm/time.h"
#include <utility>

namespace mio
{
/**
* @brief ABM simulation in one node of the abm graph model
*/
class ABMSimulationNode
{

public:
    using Sim = mio::abm::Simulation;

    template <class... Args, typename = std::enable_if_t<std::is_constructible<Sim, Args...>::value, void>>
    ABMSimulationNode(Args&&... args)
        : m_simulation(std::forward<Args>(args)...)
    {
    }

    /**
    *@brief get abm simulation in this node.
    */
    Sim& get_simulation()
    {
        return m_simulation;
    }
    const Sim& get_simulation() const
    {
        return m_simulation;
    }

    /**
    * @brief advances the simulation in this node by t+dt and logs information in History object(s)
    * @tparam History history object type(s)
    * @param[in] t Current time point
    * @param[in] dt Time span that shoulb be advanced
    * @param[in, out] history History object(s) storing simulation information
    */
    template <class... History>
    void evolve(mio::abm::TimePoint t, mio::abm::TimeSpan dt, History&... history)
    {
        m_simulation.advance(t + dt, history...);
    }

private:
    Sim m_simulation;
};

class ABMMobilityEdge
{

public:
    ABMMobilityEdge()
    {
    }

    void apply_migration(ABMSimulationNode& node_from, ABMSimulationNode& node_to);
};

void ABMMobilityEdge::apply_migration(ABMSimulationNode& node_from, ABMSimulationNode& node_to)
{
}

} // namespace mio

#endif // ABM_GRAPH_MOBILITY_H