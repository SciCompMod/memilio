#ifndef MIGRATION_H
#define MIGRATION_H

#include "epidemiology/utils/graph.h"

#include "epidemiology/utils/eigen.h"
#include <vector>

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

    void advance(double t_max = 1.0)
    {
        auto dt = m_dt;
        while (m_t < t_max) {
            if (m_t + dt > t_max) {
                dt = t_max - m_t;
            }

            for (auto& n : m_graph.nodes()) {
                m_node_func(m_t, dt, n);
            }

            m_t += dt;

            for (auto& e : m_graph.edges()) {
                m_edge_func(m_t, dt, e.property, m_graph.nodes()[e.start_node_idx], m_graph.nodes()[e.end_node_idx]);
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
    node_function m_node_func;
    edge_function m_edge_func;
};

template <class Graph, class NodeF, class EdgeF>
auto make_graph_sim(double t0, double dt, Graph&& g, NodeF&& node_func, EdgeF&& edge_func)
{
    return GraphSimulation<std::decay_t<Graph>>(t0, dt, g, std::forward<NodeF>(node_func),
                                                std::forward<EdgeF>(edge_func));
}

/**
 * represents the simulation in one node of the graph.
 */
template <class Model>
class ModelNode
{
public:
    template <class... Args, typename = std::enable_if_t<std::is_constructible<Model, Args...>::value, void>>
    ModelNode(Args&&... args)
        : model(std::forward<Args>(args)...)
    {
    }

    /**
     * get the result of the simulation in this node.
     */
    decltype(auto) get_result() const
    {
        return model.get_result();
    }

    /**
     * return the model
     */
    Model& get_model()
    {
        return model;
    }

    /**
     * return the model
     */
    Model const & get_model() const
    {
        return model;
    }
    
    /**
     * get the parameters of the simulation in this node.
     */
    decltype(auto) get_params() const
    {
        return model.get_params();
    }

    Model model;
    Eigen::VectorXd last_state;
};

class MigrationEdge
{
public:
    MigrationEdge(const Eigen::VectorXd& coeffs)
        : coefficients(coeffs)
    {
    }
    //per group and compartment
    Eigen::VectorXd coefficients;

    bool operator==(const MigrationEdge& other) const 
    {
        return coefficients == other.coefficients;
    }
};

template <class Model>
void evolve_model(double t, double dt, ModelNode<Model>& model)
{
    model.model.advance(t + dt);
    model.last_state = model.model.get_result().get_last_value();
}

template <class Model>
void apply_migration(double /*t*/, double dt, MigrationEdge& migrationEdge, ModelNode<Model>& model1,
                     ModelNode<Model>& model2)
{
    auto migration = (dt * model1.last_state.array() * migrationEdge.coefficients.array()).matrix();
    model2.model.get_result().get_last_value() += migration;
    model1.model.get_result().get_last_value() -= migration;
}

template <typename Model>
GraphSimulation<Graph<ModelNode<Model>, MigrationEdge>>
make_migration_sim(double t0, double dt, const Graph<ModelNode<Model>, MigrationEdge>& graph)
{
    return make_graph_sim(t0, dt, graph, &evolve_model<Model>, &apply_migration<Model>);
}

} // namespace epi

#endif //MIGRATION_H
