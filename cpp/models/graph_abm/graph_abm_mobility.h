/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Julia Bicker
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#ifndef MIO_ABM_GRAPH_MOBILITY_H
#define MIO_ABM_GRAPH_MOBILITY_H

#include "abm/simulation.h"
#include "abm/time.h"
#include "abm/location_type.h"
#include "abm/parameters.h"
#include "abm/person.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/graph.h"
#include <cstddef>
#include <vector>

namespace mio
{
/**
* @brief Represents the ABM simulation in one node of the abm graph model.
*/
template <class... History>
class ABMSimulationNode
{

public:
    using Sim = abm::Simulation;

    template <class... Args, typename = std::enable_if_t<std::is_constructible<Sim, Args...>::value, void>>
    ABMSimulationNode(std::tuple<History...> history, Args&&... args)
        : m_history(history)
        , m_simulation(std::forward<Args>(args)...)
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
    void evolve(mio::abm::TimePoint t, mio::abm::TimeSpan dt)
    {
        m_simulation.advance(t + dt, std::get<0>(m_history));
    }

private:
    Sim m_simulation; ///< ABM Simulation of the node
    std::tuple<History...> m_history;
};

/**
 * @brief Parameters influencing the mobility between two abm graph nodes.
 */
class ABMMobilityParameters
{

public:
    using MobilityRuleType = abm::LocationType (*)(const abm::Person&, abm::TimePoint, const abm::Parameters&);

    /**
     * Constructor for initializing commuting persons
     * @param commuting_persons Vector holding commuting persons' ids
     */
    ABMMobilityParameters(const std::vector<uint32_t>& commuting_persons,
                          const std::vector<MobilityRuleType>& mobility_rules)
        : m_commuting_persons(commuting_persons)
        , m_mobility_rules(mobility_rules)
    {
    }

    /**
     * Equality comparison operators
     */
    bool operator==(const ABMMobilityParameters& other) const
    {
        return m_commuting_persons == other.m_commuting_persons;
    }
    bool operator!=(const ABMMobilityParameters& other) const
    {
        return m_commuting_persons != other.m_commuting_persons;
    }

    /**
      * Get/Set the commuting persons vector.
      * The vector represents the persons (by their id) that commute from one node to another 
      * according to mobility rules.
      */
    const std::vector<uint32_t>& get_commuting_persons() const
    {
        return m_commuting_persons;
    }

    std::vector<uint32_t>& get_commuting_persons()
    {
        return m_commuting_persons;
    }
    /**
     * @param[in] commuting_persons Vector with commuting person ids.
     */
    void set_commuting_persons(const std::vector<uint32_t>& commuting_persons)
    {
        m_commuting_persons = commuting_persons;
    }

    /**
     * Get/ the mobility rules.
     * The rules are applied to the persons in m_commuting_persons every time exchange betwen two nodes is triggered.
     */
    const std::vector<MobilityRuleType>& get_mobility_rules() const
    {
        return m_mobility_rules;
    }
    std::vector<MobilityRuleType>& get_mobility_rules()
    {
        return m_mobility_rules;
    }
    /**
      * @brief Add mobility rule to member vector.
      * @param[in] mobility_rule Rule to be added for mobility between nodes.
      */
    void add_mobility_rule(const MobilityRuleType& mobility_rule)
    {
        m_mobility_rules.push_back(mobility_rule);
    }

private:
    std::vector<uint32_t> m_commuting_persons; ///< Person ids that are commuting via an edge
    std::vector<MobilityRuleType> m_mobility_rules; ///< Rules for moving persons from one node to another
};

/**
 * Represents the mobility between two nodes.
 */
template <class... History>
class ABMMobilityEdge
{
    using MobilityRuleType = abm::LocationType (*)(const abm::Person&, abm::TimePoint, const abm::Parameters&);

public:
    /**
     * Creates edge with mobility parameters
     * @param params mobility parameters including people commuting via the edge and mobility rules
     */
    ABMMobilityEdge(const ABMMobilityParameters& params)
        : m_parameters(params)
    {
    }

    ABMMobilityEdge(const std::vector<uint32_t>& commuting_persons,
                    const std::vector<MobilityRuleType>& mobility_rules = {})
        : m_parameters(commuting_persons, mobility_rules)
    {
    }

    /**
     * @brief Get mobility paramters.
     */
    const ABMMobilityParameters& get_parameters() const
    {
        return m_parameters;
    }

    /**
     * @brief Exchanges persons via the edge. 
     * Commuters are given by the ABMMobilityParameters and exchanged via mobility rules also given ABMMobilityParameters.
     * @param[in] node_from Commuters home node
     * @param[in] node_to Node commuters (temporarily) move to
     * @param[in] t Echange time point
     */
    void apply_mobility(ABMSimulationNode<History...>& node_from, ABMSimulationNode<History...>& node_to,
                        abm::TimePoint t)
    {
        // iterate over all persons that could commute via the edge
        for (auto p : m_parameters.get_commuting_persons()) {
            // as all nodes have all person it doesn't matter which node's persons we take here
            auto& person           = node_from.get_simulation().get_world().get_persons()[p];
            auto& params           = node_from.get_simulation().get_world().parameters;
            auto& current_location = person.get_location();
            for (auto& rule : m_parameters.get_mobility_rules()) {
                auto target_type = rule(person, t, params);
                abm::Location& target_location =
                    node_from.get_simulation().get_world().find_location(target_type, person);
                assert((node_from.get_simulation().get_world().get_id() == target_location.get_world_id() ||
                        node_to.get_simulation().get_world().get_id() == target_location.get_world_id()) &&
                       "Wrong graph edge. Target location is no edge node.");
                if (target_location != current_location &&
                    target_location.get_number_persons() < target_location.get_capacity().persons) {
                    person.migrate_to(target_location);
                    // change activeness status for commuted person
                    node_to.get_simulation().get_world().change_activeness(p);
                    node_from.get_simulation().get_world().change_activeness(p);
                    // only one mobility rule per person can be applied
                    break;
                }
            }
        }
    }

private:
    ABMMobilityParameters m_parameters; ///< Mobility parameters
};

/**
 * @brief Edge functor for abm graph simulation. 
 * @see ABMMobilityEdge::apply_mobility
 * @param[in] t Time point the functor is applied.
 * @param[in] edge ABMMobilityEdge for which the functor is applied.
 * @param[in] node_from Edge start node.
 * @param[in] node_to Edge end node.
 */
template <class... History>
void apply_mobility(abm::TimePoint t, abm::TimeSpan /*dt*/, ABMMobilityEdge<History...>& edge,
                    ABMSimulationNode<History...>& node_from, ABMSimulationNode<History...>& node_to)
{
    edge.apply_mobility(node_from, node_to, t);
}

/**
 * @brief Node functor for abm graph simulation.
 * @see ABMSimulationNode::evolve
 * @param[in] t Time point the functor is applied.
 * @param[in] dt Time interval the node is evolved.
 * @param[in] node ABMSimulationNode to which the functor is applied.
 */
template <class... History>
void evolve_model(abm::TimePoint t, abm::TimeSpan dt, ABMSimulationNode<History...>& node)
{
    node.evolve(t, dt);
}

/**
 * @brief Creates an abm graph simulation.
 * Every dt time step for each edge the persons given in ABMMobilityParameters move from one node to the other 
 * according to mobility rules also given by ABMMobilityParameters.
 * @param[in] t0 Start time point of the simulation.
 * @param[in] dt Step between mobility on edges.
 * @param[in] graph Graph for simulation.
 */
template <class... History>
GraphSimulation<Graph<ABMSimulationNode<History...>, ABMMobilityEdge<History...>>, abm::TimePoint, abm::TimeSpan>
make_abm_graph_sim(abm::TimePoint t0, abm::TimeSpan dt,
                   Graph<ABMSimulationNode<History...>, ABMMobilityEdge<History...>>&& graph)
{
    return make_graph_sim(t0, dt, std::move(graph), &evolve_model<History...>, &apply_mobility<History...>);
}

} // namespace mio

#endif // MIO_ABM_GRAPH_MOBILITY_H