/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Julia Bicker, Daniel Abele, Martin J. Kuehn
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
#ifndef HYBRID_GRAPH_H
#define HYBRID_GRAPH_H

#include "memilio/mobility/graph.h"
#include "memilio/utils/stl_util.h"

namespace mio
{

/**
 * @brief generic hybrid graph structure
 */
template <class ABMGraph, class ODEGraph, class ABM, class MigrationEdgeHybrid, class SetupABMFromOdeFct>
class HybridGraph
{
public:
    using HybridEdge = Edge<MigrationEdgeHybrid>;

    HybridGraph(const ABMGraph& abm_graph, const ODEGraph& ode_graph, const SetupABMFromOdeFct& setup_func,
                std::vector<std::vector<double>> household_distributions = {})
        : m_abm_graph(abm_graph)
        , m_ode_graph(ode_graph)
        , m_abm_nodes_to_ode_nodes(setup_func(household_distributions))
    {
    }

    HybridGraph(ABMGraph&& abm_graph, ODEGraph&& ode_graph, const SetupABMFromOdeFct& setup_func,
                std::vector<std::vector<double>> household_distributions = {})
        : m_abm_graph(std::move(abm_graph))
        , m_ode_graph(std::move(ode_graph))
        , m_abm_nodes_to_ode_nodes(setup_func(household_distributions))
    {
    }

    ABMGraph& get_abm_graph() &
    {
        return m_abm_graph;
    }

    const ABMGraph& get_abm_graph() const&
    {
        return m_abm_graph;
    }

    ODEGraph& get_ode_graph() &
    {
        return m_ode_graph;
    }

    const ODEGraph& get_ode_graph() const&
    {
        return m_ode_graph;
    }

    auto hybrid_edges()
    {
        return make_range(begin(m_hybrid_edges), end(m_hybrid_edges));
    }

    auto hybrid_edges() const
    {
        return make_range(begin(m_hybrid_edges), end(m_hybrid_edges));
    }

    Node<typename ABMGraph::NodeProperty>& get_abm_node_from_hybrid_edge(Edge<MigrationEdgeHybrid>& edge)
    {
        auto start_idx = edge.start_node_idx;
        auto end_idx   = edge.end_node_idx;
        auto iter      = std::find_if(m_abm_graph.nodes().begin(), m_abm_graph.nodes().end(),
                                      [start_idx](const Node<typename ABMGraph::NodeProperty>& node) {
                                     return node.id == start_idx;
                                 });
        if (iter != m_abm_graph.nodes().end()) {
            return *iter;
        }
        else {
            iter = std::find_if(m_abm_graph.nodes().begin(), m_abm_graph.nodes().end(),
                                [end_idx](const Node<typename ABMGraph::NodeProperty>& node) {
                                    return node.id == end_idx;
                                });
            if (iter != m_abm_graph.nodes().end()) {
                return *iter;
            }
            else {
                mio::log_error("Hybrid edge is not in graph");
            }
        }
    }

    Node<typename ODEGraph::NodeProperty>& get_ode_node_from_hybrid_edge(Edge<MigrationEdgeHybrid>& edge)
    {
        auto start_idx = edge.start_node_idx;
        auto end_idx   = edge.end_node_idx;
        auto iter      = std::find_if(m_ode_graph.nodes().begin(), m_ode_graph.nodes().end(),
                                      [start_idx](const Node<typename ODEGraph::NodeProperty>& node) {
                                     return node.id == start_idx;
                                 });
        if (iter != m_ode_graph.nodes().end()) {
            return *iter;
        }
        else {
            iter = std::find_if(m_ode_graph.nodes().begin(), m_ode_graph.nodes().end(),
                                [end_idx](const Node<typename ODEGraph::NodeProperty>& node) {
                                    return node.id == end_idx;
                                });
            if (iter != m_ode_graph.nodes().end()) {
                return *iter;
            }
            else {
                mio::log_error("Hybrid edge is not in graph");
            }
        }
    }

private:
    ABMGraph m_abm_graph;
    ODEGraph m_ode_graph;
    std::vector<Node<typename ABMGraph::NodeProperty>> m_abm_nodes_to_ode_nodes;
    std::vector<HybridEdge> m_hybrid_edges;
};

} // namespace mio

#endif //HYBRID_GRAPH_H
