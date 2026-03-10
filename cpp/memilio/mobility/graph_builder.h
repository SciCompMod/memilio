/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Kilian Volmer
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
#ifndef GRAPH_BUILDER_H
#define GRAPH_BUILDER_H

#include "memilio/mobility/graph.h"
#include "memilio/utils/logging.h"

#include <algorithm>

namespace mio
{

/**
 * @brief A builder class for constructing graphs.
 *
 * This class provides a interface for adding nodes and edges to a graph. It allows for efficient construction of large 
 * graphs by reserving space for nodes and edges in advance. The build method finalizes the graph by sorting edges and 
 * optionally removing duplicates.
 * The advantage over the :ref add_edge function of the Graph class is that edges are only sorted once during the build 
 * process, improving performance when adding many edges.
 * 
 * @tparam NodePropertyT Type of the node property.
 * @tparam EdgePropertyT Type of the edge property.
 */
template <class NodePropertyT, class EdgePropertyT>
class GraphBuilder
{
public:
    using NodeProperty = NodePropertyT;
    using EdgeProperty = EdgePropertyT;

    GraphBuilder() = default;
    GraphBuilder(const size_t num_nodes, const size_t num_edges)
    {
        m_nodes.reserve(num_nodes);
        m_edges.reserve(num_edges);
    }

    /**
     * @brief Add a node to the GraphBuilder.
     *
     * The property of the node is constructed from arguments.
     * @param id Id for the node.
     * @tparam args Additional arguments for node construction.
     */
    template <class... Args>
    void add_node(int id, Args&&... args)
    {
        m_nodes.emplace_back(id, std::forward<Args>(args)...);
    }

    /**
     * @brief Add an edge to the GraphBuilder.
     * 
     * @param start_node_idx Id of start node
     * @param end_node_idx Id of end node
     * @tparam args Additional arguments for edge construction
     */
    template <class... Args>
    void add_edge(size_t start_node_idx, size_t end_node_idx, Args&&... args)
    {
        assert(m_nodes.size() > start_node_idx && m_nodes.size() > end_node_idx);
        m_edges.emplace_back(start_node_idx, end_node_idx, std::forward<Args>(args)...);
    }

    /**
     * @brief Build the graph from the added nodes and edges.
     * 
     * Sorts the edges and optionally removes duplicate edges (same start and end node indices).
     * Without duplicate removal, multiple edges between the same nodes are allowed and the order of insertion is stable.
     * @param make_unique If true, duplicate edges are removed. The last added edge is kept!
     * @return Graph<NodePropertyT, EdgePropertyT> The constructed graph.
     */
    Graph<NodeProperty, EdgeProperty> build(bool make_unique = false) &&
    {
        sort_edges();
        if (make_unique) {
            remove_duplicate_edges();
        }
        Graph<NodeProperty, EdgeProperty> graph(std::move(m_nodes), std::move(m_edges));
        return graph;
    }

private:
    /**
     * @brief Sort the edge vector of a graph.
     * 
     * Sorts the edges first by start node index, then by end node index. We use stable_sort to keep the order of insertion
     * for edges with the same start and end node indices.
     */
    void sort_edges()
    {
        std::stable_sort(m_edges.begin(), m_edges.end(), [](auto&& e1, auto&& e2) {
            return e1.start_node_idx == e2.start_node_idx ? e1.end_node_idx < e2.end_node_idx
                                                          : e1.start_node_idx < e2.start_node_idx;
        });
    }

    /**
     * @brief Remove duplicate edges from a sorted edge vector.
     * 
     * Copies all the unique edges to a new vector and replaces the original edge vector with it. Unique means that
     * the start and end node indices are unique. Other edge properties are not checked and may get lost. Only the last 
     * edge in the vector is kept.
     */
    void remove_duplicate_edges()
    {
        std::vector<Edge<EdgePropertyT>> unique_edges;
        unique_edges.reserve(m_edges.size());
        bool duplicate_edges = false;
        auto curr_elem       = m_edges.begin();

        while (curr_elem != m_edges.end()) {
            auto next_elem = std::next(curr_elem);
            if (next_elem != m_edges.end() && curr_elem->start_node_idx == next_elem->start_node_idx &&
                curr_elem->end_node_idx == next_elem->end_node_idx) {
                duplicate_edges = true;
            }
            else if (next_elem != m_edges.end()) {
                std::copy(curr_elem, next_elem, std::back_inserter(unique_edges));
            }
            curr_elem = next_elem;
        }
        std::copy(std::prev(m_edges.end()), m_edges.end(), std::back_inserter(unique_edges));
        m_edges = std::move(unique_edges);
        if (duplicate_edges) {
            mio::log_warning("Removed duplicate edge(s)");
        }
    }

private:
    std::vector<Node<NodePropertyT>> m_nodes;
    std::vector<Edge<EdgePropertyT>> m_edges;
};

} // namespace mio

#endif // GRAPH_BUILDER_H
