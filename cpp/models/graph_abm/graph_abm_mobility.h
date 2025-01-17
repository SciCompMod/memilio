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
#include "abm/person_id.h"
#include "abm/model_functions.h"
#include "graph_abm/graph_abmodel.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/graph.h"
#include "memilio/utils/compiler_diagnostics.h"
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

namespace mio
{
/**
* @brief Represents the ABM simulation in one node of the ABM graph model.
*/
template <class... History>
class ABMSimulationNode
{

public:
    using Sim = abm::Simulation<GraphABModel>;

    template <class... Args, typename = std::enable_if_t<std::is_constructible<Sim, Args...>::value, void>>
    ABMSimulationNode(std::tuple<History...>&& history, Args&&... args)
        : m_simulation(std::forward<Args>(args)...)
        , m_history(history)
    {
    }

    /**
    *@brief Get abm simulation in this node.
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
     * @brief Get history object(s) in this node.
     */
    std::tuple<History...>& get_history()
    {
        return m_history;
    }
    const std::tuple<History...>& get_history() const
    {
        return m_history;
    }

    /**
    * @brief advances the simulation in this node by t+dt and logs information in History object(s)
    * @tparam History history object type(s)
    * @param[in] t Current time point
    * @param[in] dt Time span that shoulb be advanced
    * @param[in, out] history History object(s) storing simulation information
    */
    void advance(mio::abm::TimePoint t, mio::abm::TimeSpan dt)
    {
        m_simulation.advance(t + dt, std::get<0>(m_history));
    }

private:
    Sim m_simulation; ///< ABM Simulation of the node
    std::tuple<History...> m_history;
};

/**
 * Represents the mobility between two nodes.
 */
template <class... History>
class ABMMobilityEdge
{

public:
    /**
     * @brief Exchanges persons via the edge. 
     * Commuters are given by the person buffer of node_from.
     * @param[in] node_from Commuters home node
     * @param[in] node_to Node commuters (temporarily) move to
     * @param[in] t Echange time point
     */
    void apply_mobility(ABMSimulationNode<History...>& node_from, ABMSimulationNode<History...>& node_to,
                        abm::TimePoint /*t*/)
    {
        auto& model_from        = node_from.get_simulation().get_model();
        auto& model_to          = node_to.get_simulation().get_model();
        auto& persons_to_change = model_from.get_person_buffer();
        //sort vector such that we start removing the persons from the bottom
        std::sort(persons_to_change.begin(), persons_to_change.end());
        //iterate over all persons that change from node_from
        for (int i = int(persons_to_change.size()) - 1; i >= 0; --i) {
            auto& person     = model_from.get_persons()[persons_to_change[i]];
            auto target_type = person.get_location_type();
            //check if Person uses this edge
            if (person.get_assigned_location_model_id(target_type) == model_to.get_id()) {
                auto target_id = person.get_assigned_location(target_type);
                //set correct location for person
                person.set_location(target_type, target_id, model_to.get_id());
                //add person to model_to
                model_to.add_person(std::move(person));
                //remove person from model_from
                model_from.remove_person(persons_to_change[i]);
                //correct indices in persons buffer from node_from
                //here it is required that the vector is sorted
                for (size_t j = i + 1; j < persons_to_change.size(); ++j) {
                    persons_to_change[j]--;
                }
                //delete current index from list
                persons_to_change.erase(persons_to_change.begin() + i);
            }
        }
    }
};

/**
 * @brief Edge functor for abm graph simulation. 
 * @see ABMMobilityEdge::apply_mobility
 * The attribute dt is required by the GraphSimulation class and therefore an input argument of the function.
 * However it is not used in ABMMobilityEdge::apply_mobility.
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
 * @see ABMSimulationNode::advance
 * @param[in] t Time point the functor is applied.
 * @param[in] dt Time interval the node is advanced.
 * @param[in] node ABMSimulationNode to which the functor is applied.
 */
template <class... History>
void advance_model(abm::TimePoint t, abm::TimeSpan dt, ABMSimulationNode<History...>& node)
{
    node.advance(t, dt);
}

/**
 * @brief Creates an abm graph simulation.
 * Every dt time step for each edge the persons that want to change to a location in another node 
 * are removed from the model in their former location's node and added to the model of the new location.
 * @param[in] t0 Start time point of the simulation.
 * @param[in] dt Step between mobility on edges.
 * @param[in] graph Graph for simulation.
 */
template <class... History>
GraphSimulation<Graph<ABMSimulationNode<History...>, ABMMobilityEdge<History...>>, abm::TimePoint, abm::TimeSpan,
                void (*)(mio::abm::TimePoint, mio::abm::TimeSpan, mio::ABMMobilityEdge<History...>&,
                         mio::ABMSimulationNode<History...>&, mio::ABMSimulationNode<History...>&),
                void (*)(mio::abm::TimePoint, mio::abm::TimeSpan, mio::ABMSimulationNode<History...>&)>
make_abm_graph_sim(abm::TimePoint t0, abm::TimeSpan dt,
                   Graph<ABMSimulationNode<History...>, ABMMobilityEdge<History...>>&& graph)
{
    return make_graph_sim(t0, dt, std::move(graph), &advance_model<History...>, &apply_mobility<History...>);
}

} // namespace mio

#endif // MIO_ABM_GRAPH_MOBILITY_H
