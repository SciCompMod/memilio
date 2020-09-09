#ifndef MIGRATION_H
#define MIGRATION_H

#include "epidemiology/utils/graph.h"

#include <Eigen/Core>
#include <vector>

namespace epi
{

/**
 * @brief abstract simulation on a graph with alternating node and edge actions
 */
template <class Graph, class NodeF, class EdgeF>
class GraphSimulation
{
public:
    GraphSimulation(double t0, double dt, const Graph& g, const NodeF& node_func, const EdgeF&& edge_func)
        : m_graph(g)
        , m_t(t0)
        , m_dt(dt)
        , m_node_func(node_func)
        , m_edge_func(edge_func)
    {
    }

    void advance(int n_steps = 1)
    {
        for (int i = 0; i < n_steps; i++) {
            for (auto& n : m_graph.nodes()) {
                m_node_func(m_t, m_dt, n);
            }

            m_t += m_dt;

            for (auto& e : m_graph.edges()) {
                m_edge_func(m_t, m_dt, e.property, m_graph.nodes()[e.start_node_idx], m_graph.nodes()[e.end_node_idx]);
            }
        }
    }

    double get_t() const
    {
        return m_t;
    }

    Graph& get_graph()
    {
        return m_graph;
    }

    const Graph& get_graph() const
    {
        return m_graph;
    }

private:
    double m_t;
    double m_dt;
    Graph m_graph;
    NodeF m_node_func;
    EdgeF m_edge_func;
};

template <class Graph, class NodeF, class EdgeF>
auto make_graph_sim(double t0, double dt, Graph&& g, NodeF&& node_func, EdgeF&& edge_func)
{
    return GraphSimulation<std::decay_t<Graph>, NodeF, EdgeF>(t0, dt, g, std::forward<NodeF>(node_func),
                                                              std::forward<EdgeF>(edge_func));
}

template <class Model>
class ModelNode
{
public:
    template <class... Args, typename = std::enable_if_t<std::is_constructible<Model, Args...>::value, void>>
    ModelNode(Args&&... args)
        : model(std::forward<Args>(args)...)
    {
    }
    Model model;
    Eigen::VectorXd last_state;
};

class MigrationEdge
{
public:
    MigrationEdge(const Eigen::VectorXd& coefficients)
        : coefficients(coefficients)
    {
    }
    //per group and compartment
    Eigen::VectorXd coefficients;
};

template <class Model>
void evolve_model(double t, double dt, ModelNode<Model>& model)
{
    model.model.advance(t + dt);
    model.last_state = model.model.get_result().get_last_value();
}

template <class Model>
void apply_migration(double t, double dt, MigrationEdge& migrationEdge, ModelNode<Model>& model1,
                     ModelNode<Model>& model2)
{
    auto migration = (dt * model1.last_state.array() * migrationEdge.coefficients.array()).matrix();
    model2.model.get_result().get_last_value() += migration;
    model1.model.get_result().get_last_value() -= migration;
}

template <typename Model>
auto make_migration_sim(double t0, double dt, const Graph<ModelNode<Model>, MigrationEdge>& graph)
{
    return make_graph_sim(t0, dt, graph, &evolve_model<Model>, &apply_migration<Model>);
}

} // namespace epi

#endif //MIGRATION_H
