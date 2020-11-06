#ifndef MIGRATION_H
#define MIGRATION_H

#include "epidemiology/migration/graph_simulation.h"
#include "epidemiology/utils/time_series.h"
#include "epidemiology/utils/eigen.h"
#include "epidemiology/utils/eigen_util.h"
#include "epidemiology/secir/populations.h"
#include "epidemiology/utils/compiler_diagnostics.h"

#include <cassert>

namespace epi
{

class SecirParams;
class SeirParams;

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
        , m_last_state(model.get_result().get_last_value())
        , m_t0(model.get_result().get_last_time())
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
     * get the result of the simulation in this node.
     */
    decltype(auto) get_result()
    {
        return model.get_result();
    }

    /**
     * get the parameters of the simulation in this node.
     */
    decltype(auto) get_params() const
    {
        return model.get_params();
    }

    Eigen::Ref<const Eigen::VectorXd> get_last_state() const
    {
        return m_last_state;
    }

    double get_t0() const
    {
        return m_t0;
    }

    void evolve(double t, double dt)
    {
        model.advance(t + dt);
        m_last_state = model.get_result().get_last_value();
    }

    Model model;

private:
    Eigen::VectorXd m_last_state;
    double m_t0;
};

/** 
 * represents the migration between two nodes.
 */
class MigrationEdge
{
public:
    /**
     * create edge with coefficients.
     * @param coeffs % of people in each group and compartment that migrate in each time step.
     */
    MigrationEdge(const Eigen::VectorXd& coeffs)
        : m_coefficients(coeffs)
        , m_migrated(coeffs.rows())
        , m_return_times(0)
        , m_return_migrated(false)
    {
        assert((m_coefficients.array() >= 0).all() && "migration coefficients must be nonnegative");
    }

    bool operator==(const MigrationEdge& other) const
    {
        return m_coefficients == other.m_coefficients;
    }

    bool operator!=(const MigrationEdge& other) const
    {
        return m_coefficients != other.m_coefficients;
    }

    /**
     * get the migration coefficients.
     */
    Eigen::Ref<const Eigen::VectorXd> get_coefficients() const
    {
        return m_coefficients;
    }
    Eigen::Ref<Eigen::VectorXd> get_coefficients()
    {
        return m_coefficients;
    }

    /**
     * adjust number of migrated people when they return.
     * @param[inout] migrated number of people that migrated as input, number of people that return as output
     * @param params parameters of model in the node that the people migrated to.
     * @param total total population in the node that the people migrated to.
     * @param t time of migration
     * @param dt time between migration and return
     */
    void calculate_returns_ode(Eigen::Ref<TimeSeries<double>::Vector> migrated, const SecirParams& params,
                               Eigen::Ref<const TimeSeries<double>::Vector> total, double t, double dt);
    void calculate_returns_ode(Eigen::Ref<TimeSeries<double>::Vector> migrated, const SeirParams& params,
                               Eigen::Ref<const TimeSeries<double>::Vector> total, double t, double dt);

    /**
     * compute migration from node_from to node_to.
     * migration is based on coefficients.
     * migrants are added to the current state of node_to, subtracted from node_from.
     * on return, migrants (adjusted for infections) are subtracted from node_to, added to node_from.
     * @param t current time
     * @param dt last time step (fixed to 0.5 for migration model)
     * @param node_from node that people migrated from, return to
     * @param node_to node that people migrated to, return from
     */
    template <class Model>
    void apply_migration(double t, double dt, ModelNode<Model>& node_from, ModelNode<Model>& node_to);

private:
    Eigen::VectorXd m_coefficients; //one per group and compartment
    TimeSeries<double> m_migrated;
    TimeSeries<double> m_return_times;
    bool m_return_migrated;
};

template <class Model>
void MigrationEdge::apply_migration(double t, double dt, ModelNode<Model>& node_from, ModelNode<Model>& node_to)
{
    //returns
    for (Eigen::Index i = m_return_times.get_num_time_points() - 1; i >= 0; --i) {
        if (m_return_times.get_time(i) <= t) {
            auto v0 = find_value_reverse(node_to.get_result(), m_migrated.get_time(i), 1e-10, 1e-10);
            assert(v0 != node_to.get_result().rend() && "unexpected error.");
            calculate_returns_ode(m_migrated[i], node_to.get_params(), *v0, m_migrated.get_time(i), dt);
            node_from.get_result().get_last_value() += m_migrated[i];
            node_to.get_result().get_last_value() -= m_migrated[i];
            m_migrated.remove_time_point(i);
            m_return_times.remove_time_point(i);
        }
    }

    if (!m_return_migrated) {
        //normal daily migration
        m_migrated.add_time_point(t, (node_from.get_last_state().array() * m_coefficients.array()).matrix());
        m_return_times.add_time_point(t + dt);

        node_to.get_result().get_last_value() += m_migrated.get_last_value();
        node_from.get_result().get_last_value() -= m_migrated.get_last_value();
    }
    m_return_migrated = !m_return_migrated;
}

/**
 * edge functor for migration simulation.
 * @see ModelNode::evolve
 */
template <class Model>
void evolve_model(double t, double dt, ModelNode<Model>& node)
{
    node.evolve(t, dt);
}

/**
 * edge functor for migration simulation.
 * @see MigrationEdge::apply_migration
 */
template <class Model>
void apply_migration(double t, double dt, MigrationEdge& migrationEdge, ModelNode<Model>& node_from,
                     ModelNode<Model>& node_to)
{
    migrationEdge.apply_migration(t, dt, node_from, node_to);
}

/**
 * create a migration simulation.
 * After every second time step, for each edge a portion of the population corresponding to the coefficients of the edge
 * moves from one node to the other. In the next timestep, the migrated population return to their "home" node. 
 * Returns are adjusted based on the development in the target node. 
 * @param t0 start time of the simulation
 * @param dt time step between migrations
 * @param graph set up for migration simulation
 */
template <typename Model>
GraphSimulation<Graph<ModelNode<Model>, MigrationEdge>>
make_migration_sim(double t0, double dt, const Graph<ModelNode<Model>, MigrationEdge>& graph)
{
    return make_graph_sim(t0, dt, graph, &evolve_model<Model>, &apply_migration<Model>);
}

} // namespace epi

#endif //MIGRATION_H
