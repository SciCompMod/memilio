/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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
#ifndef EPI_MIGRATION_GRAPH_SIMULATION_H
#define EPI_MIGRATION_GRAPH_SIMULATION_H

#include "epidemiology/utils/graph.h"

namespace epi
{

/**
 * @brief abstract simulation on a graph with alternating node and edge actions
 */
template <class Graph>
class GraphSimulation
{
public:
    using node_function = std::function<void(double, double, typename Graph::NodeProperty&)>;

    using edge_function = std::function<void(double, double, typename Graph::EdgeProperty&,
                                             typename Graph::NodeProperty&, typename Graph::NodeProperty&)>;

    GraphSimulation(double t0, double dt, const Graph& g, const node_function& node_func,
                    const edge_function&& edge_func)
        : m_t(t0)
        , m_dt(dt)
        , m_graph(g)
        , m_node_func(node_func)
        , m_edge_func(edge_func)
    {
    }

    GraphSimulation(double t0, double dt, Graph&& g, const node_function& node_func, const edge_function&& edge_func)
        : m_t(t0)
        , m_dt(dt)
        , m_graph(std::move(g))
        , m_node_func(node_func)
        , m_edge_func(edge_func)
    {
    }

    void advance(double t_max = 1.0)
    {
        auto dt = m_dt;
        while (m_t < t_max) {
            if (m_t + dt > t_max) {
                dt = t_max - m_t;
            }

            for (auto& n : m_graph.nodes()) {
                if (n.id == 16069) {
                    std::cout << "arrived" << std::endl;
                }
                m_node_func(m_t, dt, n.property);
            }

            m_t += dt;

            for (auto& e : m_graph.edges()) {
                /*if (e.start_node_idx == 316 && e.end_node_idx == 259) {

                    std::cout << "From ID: " << e.start_node_idx << " To ID: " << e.end_node_idx << ", t: " << m_t
                              << std::endl;
                }*/
                m_edge_func(m_t, dt, e.property, m_graph.nodes()[e.start_node_idx].property,
                            m_graph.nodes()[e.end_node_idx].property);
            }
        }
    }

    double get_t() const
    {
        return m_t;
    }

    Graph& get_graph() &
    {
        return m_graph;
    }

    const Graph& get_graph() const&
    {
        return m_graph;
    }

    Graph&& get_graph() &&
    {
        return std::move(m_graph);
    }

private:
    double m_t;
    double m_dt;
    Graph m_graph;
    node_function m_node_func;
    edge_function m_edge_func;
};

template <class Graph, class NodeF, class EdgeF>
auto make_graph_sim(double t0, double dt, Graph&& g, NodeF&& node_func, EdgeF&& edge_func)
{
    return GraphSimulation<std::decay_t<Graph>>(t0, dt, std::forward<Graph>(g), std::forward<NodeF>(node_func),
                                                std::forward<EdgeF>(edge_func));
}

} // namespace epi
#endif //EPI_MIGRATION_GRAPH_SIMULATION_H
