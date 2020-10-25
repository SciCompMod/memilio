#ifndef MIGRATION_H
#define MIGRATION_H

#include "epidemiology/migration/graph_simulation.h"
#include "epidemiology/utils/time_series.h"
#include "epidemiology/utils/eigen.h"
#include "epidemiology/utils/eigen_util.h"
#include "epidemiology/secir/populations.h"

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

class MigrationEdge
{
public:
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

    Eigen::Ref<const Eigen::VectorXd> get_coefficients() const
    {
        return m_coefficients;
    }

    Eigen::Ref<Eigen::VectorXd> get_coefficients()
    {
        return m_coefficients;
    }
    
    void calculate_returns_ode(Eigen::Ref<TimeSeries<double>::Vector> migrated, const SecirParams& params,
                               Eigen::Ref<const TimeSeries<double>::Vector> total, double t, double dt);

    void calculate_returns_ode(Eigen::Ref<TimeSeries<double>::Vector> migrated, const SeirParams& params,
                               Eigen::Ref<const TimeSeries<double>::Vector> total, double t, double dt);

    template <class Model>
    void apply_migration(double t, double dt, ModelNode<Model>& node1, ModelNode<Model>& node2);

private:
    Eigen::VectorXd m_coefficients; //one per group and compartment
    TimeSeries<double> m_migrated;
    TimeSeries<double> m_return_times;
    bool m_return_migrated;
};

template <class FP>
typename TimeSeries<FP>::const_reverse_iterator find_value_reverse(const TimeSeries<FP> ts, FP t_search)
{
    auto iter_t = find_if(ts.get_reverse_times().begin(), ts.get_reverse_times().end(), [=](auto t) {
        return std::abs(t - t_search) < 1e-10;
    });
    if (iter_t != ts.get_reverse_times().end()) {
        return ts.rbegin() + (iter_t - ts.get_reverse_times().begin());
    }
    return ts.rend();
}

template <class Model>
void MigrationEdge::apply_migration(double t, double dt, ModelNode<Model>& node1, ModelNode<Model>& node2)
{
    assert(dt == 0.5);

    //returns
    for (Eigen::Index i = m_return_times.get_num_time_points() - 1; i >= 0; --i) {
        if (m_return_times.get_time(i) <= t) {
            auto v0 = *find_value_reverse(node2.get_result(), m_migrated.get_time(i));
            calculate_returns_ode(m_migrated[i], node2.get_params(), v0, m_migrated.get_time(i), dt);
            node1.get_result().get_last_value() += m_migrated[i];
            node2.get_result().get_last_value() -= m_migrated[i];
            m_migrated.remove_time_point(i);
            m_return_times.remove_time_point(i);
        }
    }

    if (!m_return_migrated) {
        //normal daily migration
        const auto migration = (node1.get_last_state().array() * m_coefficients.array()).matrix();

        node2.get_result().get_last_value() += migration;
        node1.get_result().get_last_value() -= migration;

        //I don't understand why this doesn't work: m_migrated.add_time_point(t, migration);
        m_migrated.add_time_point(t, (node1.get_last_state().array() * m_coefficients.array()).matrix());
        m_return_times.add_time_point(t + 0.5);
    }
    m_return_migrated = !m_return_migrated;
}

template <class Model>
void evolve_model(double t, double dt, ModelNode<Model>& node)
{
    node.evolve(t, dt);
}

template <class Model>
void apply_migration(double t, double dt, MigrationEdge& migrationEdge, ModelNode<Model>& node1,
                     ModelNode<Model>& node2)
{
    migrationEdge.apply_migration(t, dt, node1, node2);
}

template <typename Model>
GraphSimulation<Graph<ModelNode<Model>, MigrationEdge>>
make_migration_sim(double t0, const Graph<ModelNode<Model>, MigrationEdge>& graph)
{
    return make_graph_sim(t0, 0.5, graph, &evolve_model<Model>, &apply_migration<Model>);
}

} // namespace epi

#endif //MIGRATION_H
