#ifndef GRAPH_H
#define GRAPH_H

#include <epidemiology/stl_util.h>
#include <epidemiology/seir.h>
#include <epidemiology/secir.h>
#include <iostream>

namespace epi
{

struct OutEdgeBase {
    OutEdgeBase(size_t start)
        : start_node_idx(start)
    {
    }
    size_t start_node_idx;
};
struct InEdgeBase {
    InEdgeBase(size_t end)
        : end_node_idx(end)
    {
    }
    size_t end_node_idx;
};
struct EdgeBase : public OutEdgeBase, public InEdgeBase {
    EdgeBase(size_t start, size_t end)
        : OutEdgeBase(start)
        , InEdgeBase(end)
    {
    }
};

/**
 * @brief represents an edge of the graph
 */
template<class EdgePropertyT>
struct Edge : public EdgeBase {
    template <class... Args>
    Edge(size_t start, size_t end, Args&&... args)
        : EdgeBase{start, end}
        , property(std::forward<Args>(args)...)
    {
    }

    EdgePropertyT property;
};

/**
 * @brief comparison operator if edge property type is equality comparable
 */
template <class T>
std::enable_if_t<has_eq_op<T>::value, bool> operator==(const Edge<T>& e1, const Edge<T>& e2)
{
    return e1.start_node_idx == e2.start_node_idx && e1.end_node_idx == e2.end_node_idx && e1.property == e2.property;
}

/**
 * @brief out stream operator for edges if edge property type has stream operator defined
 */
template<class T>
std::enable_if_t<has_ostream_op<T>::value, std::ostream&> operator<<(std::ostream& os, const Edge<T>& e)
{
    os << e.start_node_idx << " > " << e.end_node_idx << " : " << e.property; 
    return os;
}

/**
 * @brief out stream operator for edges if edge property type does not have stream operator defined
 */
template<class T>
std::enable_if_t<!has_ostream_op<T>::value, std::ostream&> operator<<(std::ostream& os, const Edge<T>& e)
{
    os << e.start_node_idx << " > " << e.end_node_idx;
    return os;
}

/**
 * @brief generic graph structure
 */
template <class NodePropertyT, class EdgePropertyT>
class Graph
{
public:
    /**
     * @brief add a node to the graph. property of the node is constructed from arguments.
     */
    template <class... Args>
    NodePropertyT& add_node(Args&&... args)
    {
        m_nodes.emplace_back(std::forward<Args>(args)...);
        return m_nodes.back();
    }

    /**
     * @brief add an edge to the graph. property of the edge is constructed from arguments.
     */
    template <class... Args>
    Edge<EdgePropertyT>& add_edge(size_t start_node_idx, size_t end_node_idx, Args&&... args)
    {
        assert(m_nodes.size() > start_node_idx && m_nodes.size() > end_node_idx);
        return *insert_sorted_replace(
            m_edges, Edge<EdgePropertyT>(start_node_idx, end_node_idx, std::forward<Args>(args)...), [](auto&& e1, auto&& e2) {
                return e1.start_node_idx == e2.start_node_idx ? e1.end_node_idx < e2.end_node_idx
                                                              : e1.start_node_idx < e2.start_node_idx;
            });
    }

    /**
     * @brief range of nodes
     */
    auto nodes()
    {
        return make_range(begin(m_nodes), end(m_nodes));
    }

    /**
     * @brief range of nodes
     */
    auto nodes() const
    {
        return make_range(begin(m_nodes), end(m_nodes));
    };

    /**
     * @brief range of edges
     */
    auto edges()
    {
        return make_range(begin(m_edges), end(m_edges));
    }

    /**
     * @brief range of edges
     */
    auto edges() const
    {
        return make_range(begin(m_edges), end(m_edges));
    }

    /**
     * @brief range of edges going out from a specific node
     */
    auto out_edges(size_t node_idx)
    {
        return out_edges(begin(m_edges), end(m_edges), node_idx);
    }

    /**
     * @brief range of edges going out from a specific node
     */
    auto out_edges(size_t node_idx) const
    {
        return out_edges(begin(m_edges), end(m_edges), node_idx);
    }

private:
    template <typename Iter>
    static auto out_edges(Iter b, Iter e, size_t idx)
    {
        return make_range(std::equal_range(b, e, OutEdgeBase(idx), [](auto&& e1, auto&& e2) {
            return e1.start_node_idx < e2.start_node_idx;
        }));
    }

private:
    std::vector<NodePropertyT> m_nodes;
    std::vector<Edge<EdgePropertyT>> m_edges;
};

template <class Graph, class EdgeF, class NodeF>
void simulate_graph(double t0, double tmax, double dt, Graph&& g, NodeF&& node_func, EdgeF&& edge_func)
{
    auto t = t0;
    while (t < tmax) {
        auto dt_eff = std::min(dt, tmax - t);
        for (auto& n : g.nodes()) {
            node_func(t, dt_eff, n);
        }

        t += dt_eff;

        for (auto& e : g.edges()) {
            edge_func(t, dt_eff, e.property, g.nodes()[e.start_node_idx], g.nodes()[e.end_node_idx]);
        }
    }
}

template <class T>
std::enable_if_t<!has_ostream_op<T>::value, void> print_graph_object(std::ostream& os, size_t idx, const T& o)
{
    os << idx;
}

template <class T>
std::enable_if_t<has_ostream_op<T>::value, void> print_graph_object(std::ostream& os, size_t idx, const T& o)
{
    os << idx << " [" << o << "]";
}

template <class Graph>
void print_graph(std::ostream& os, const Graph& g)
{
    auto nodes = g.nodes();
    for (size_t i = 0; i < nodes.size(); ++i) {
        os << "NODE ";
        print_graph_object(os, i, nodes[i]);
        os << '\n';
    }

    auto edges = g.edges();
    for (size_t i = 0; i < edges.size(); ++i) {
        auto& e = edges[i];
        os << "EDGE ";
        print_graph_object(os, i, e.property);
        os << " FROM NODE ";
        print_graph_object(os, e.start_node_idx, nodes[e.start_node_idx]);
        os << " TO ";
        print_graph_object(os, e.end_node_idx, nodes[e.end_node_idx]);
        os << '\n';
    }
}

} // namespace epi

#endif //GRAPH_H