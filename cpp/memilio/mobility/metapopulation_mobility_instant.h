/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#ifndef METAPOPULATION_MOBILITY_INSTANT_H
#define METAPOPULATION_MOBILITY_INSTANT_H

#include "memilio/mobility/graph_simulation.h"
#include "memilio/utils/time_series.h"
#include "memilio/math/eigen.h"
#include "memilio/math/eigen_util.h"
#include "memilio/utils/metaprogramming.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/math/euler.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/dynamic_npis.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/date.h"

#include "boost/filesystem.hpp"

#include <cassert>

namespace mio
{

/**
 * represents the simulation in one node of the graph.
 */
template <class Sim>
class SimulationNode
{
public:
    template <class... Args, typename = std::enable_if_t<std::is_constructible<Sim, Args...>::value, void>>
    SimulationNode(Args&&... args)
        : m_simulation(std::forward<Args>(args)...)
        , m_last_state(m_simulation.get_result().get_last_value())
        , m_t0(m_simulation.get_result().get_last_time())
    {
    }

    /**
     * get the result of the simulation in this node.
     * @{
     */
    decltype(auto) get_result() const
    {
        return m_simulation.get_result();
    }
    decltype(auto) get_result()
    {
        return m_simulation.get_result();
    }
    /**@}*/

    /**
     * get the the simulation in this node.
     * @{
     */
    Sim& get_simulation()
    {
        return m_simulation;
    }
    const Sim& get_simulation() const
    {
        return m_simulation;
    }
    /**@}*/

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
        m_simulation.advance(t + dt);
        m_last_state = m_simulation.get_result().get_last_value();
    }

private:
    Sim m_simulation;
    Eigen::VectorXd m_last_state;
    double m_t0;
};

/**
 * time dependent migration coefficients.
 */
using MigrationCoefficients = DampingMatrixExpression<VectorDampings>;

/**
 * sum of time dependent migration coefficients.
 * differentiate between sources of migration.
 */
using MigrationCoefficientGroup = DampingMatrixExpressionGroup<MigrationCoefficients>;

/**
 * parameters that influence migration.
 */
class MigrationParameters
{
public:
    /**
     * constructor from migration coefficients.
     * @param coeffs migration coefficients
     */
    MigrationParameters(const MigrationCoefficientGroup& coeffs)
        : m_coefficients(coeffs)
    {
    }

    /**
     * constructor from migration coefficients.
     * @param coeffs migration coefficients
     */
    MigrationParameters(const Eigen::VectorXd& coeffs)
        : m_coefficients({MigrationCoefficients(coeffs)})
    {
    }

    /** 
     * equality comparison operators
     */
    //@{
    bool operator==(const MigrationParameters& other) const
    {
        return m_coefficients == other.m_coefficients;
    }
    bool operator!=(const MigrationParameters& other) const
    {
        return m_coefficients != other.m_coefficients;
    }
    //@}

    /**
     * Get/Setthe migration coefficients.
     * The coefficients represent the (time-dependent) percentage of people migrating 
     * from one node to another by age and infection compartment. 
     * @{
     */
    /**
     * @return the migration coefficients.
     */
    const MigrationCoefficientGroup& get_coefficients() const
    {
        return m_coefficients;
    }
    MigrationCoefficientGroup& get_coefficients()
    {
        return m_coefficients;
    }
    /**
     * @param coeffs the migration coefficients.
     */
    void set_coefficients(const MigrationCoefficientGroup& coeffs)
    {
        m_coefficients = coeffs;
    }
    /** @} */

    /**
     * Get/Set dynamic NPIs that are implemented when relative infections exceed thresholds.
     * This feature is optional. The simulation model needs to overload the get_infected_relative function.
     * @{
     */
    /**
     * @return dynamic NPIs for relative infections.
     */
    const DynamicNPIs& get_dynamic_npis_infected() const
    {
        return m_dynamic_npis;
    }
    DynamicNPIs& get_dynamic_npis_infected()
    {
        return m_dynamic_npis;
    }
    /**
     * @param v dynamic NPIs for relative infections.
     */
    void set_dynamic_npis_infected(const DynamicNPIs& v)
    {
        m_dynamic_npis = v;
    }
    /** @} */

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("MigrationParameters");
        obj.add_element("Coefficients", m_coefficients);
        obj.add_element("DynamicNPIs", m_dynamic_npis);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<MigrationParameters> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("MigrationParameters");
        auto c   = obj.expect_element("Coefficients", Tag<MigrationCoefficientGroup>{});
        auto d   = obj.expect_element("DynamicNPIs", Tag<DynamicNPIs>{});
        return apply(
            io,
            [](auto&& c_, auto&& d_) {
                MigrationParameters params(c_);
                params.set_dynamic_npis_infected(d_);
                return params;
            },
            c, d);
    }

private:
    MigrationCoefficientGroup m_coefficients; //one per group and compartment
    DynamicNPIs m_dynamic_npis;
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
    MigrationEdge(const MigrationParameters& params)
        : m_parameters(params)
        , m_migrated(params.get_coefficients().get_shape().rows())
        , m_return_times(0)
        , m_return_migrated(false)
    {
    }

    /**
     * create edge with coefficients.
     * @param coeffs % of people in each group and compartment that migrate in each time step.
     */
    MigrationEdge(const Eigen::VectorXd& coeffs)
        : m_parameters(coeffs)
        , m_migrated(coeffs.rows())
        , m_return_times(0)
        , m_return_migrated(false)
    {
    }

    /**
     * get the migration parameters.
     */
    const MigrationParameters& get_parameters() const
    {
        return m_parameters;
    }

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
    template <class Sim>
    void apply_migration(double t, double dt, SimulationNode<Sim>& node_from, SimulationNode<Sim>& node_to);

private:
    MigrationParameters m_parameters;
    TimeSeries<double> m_migrated;
    TimeSeries<double> m_return_times;
    bool m_return_migrated;
    double m_t_last_dynamic_npi_check               = -std::numeric_limits<double>::infinity();
    std::pair<double, SimulationTime> m_dynamic_npi = {-std::numeric_limits<double>::max(), SimulationTime(0)};
};

/**
 * adjust number of migrated people when they return according to the model.
 * E.g. during the time in the other node, some people who left as susceptible will return exposed.
 * Implemented for general compartmentmodel simulations, overload for your custom model if necessary
 * so that it can be found with argument-dependent lookup, i.e. in the same namespace as the model.
 * @param[inout] migrated number of people that migrated as input, number of people that return as output
 * @param params parameters of model in the node that the people migrated to.
 * @param total total population in the node that the people migrated to.
 * @param t time of migration
 * @param dt time between migration and return
 */
template <class Sim, class = std::enable_if_t<is_compartment_model_simulation<Sim>::value>>
void calculate_migration_returns(Eigen::Ref<TimeSeries<double>::Vector> migrated, const Sim& sim,
                                 Eigen::Ref<const TimeSeries<double>::Vector> total, double t, double dt)
{
    auto y0 = migrated.eval();
    auto y1 = migrated;
    EulerIntegratorCore().step(
        [&](auto&& y, auto&& t_, auto&& dydt) {
            sim.get_model().get_derivatives(total, y, t_, dydt);
        },
        y0, t, dt, y1);
}

/**
 * detect a get_infections_relative function for the Model type.
 */
template <class Sim>
using get_infections_relative_expr_t = decltype(get_infections_relative(
    std::declval<const Sim&>(), std::declval<double>(), std::declval<const Eigen::Ref<const Eigen::VectorXd>&>()));

/**
 * get the percantage of infected people of the total population in the node
 * If dynamic NPIs are enabled, there needs to be an overload of get_infections_relative(model, y)
 * for the Model type that can be found with argument-dependent lookup. Ideally define get_infections_relative 
 * in the same namespace as the Model type.
 * @param node a node of a migration graph.
 * @param y the current value of the simulation.
 * @param t the current simulation time
 */
template <class Sim,
          std::enable_if_t<!is_expression_valid<get_infections_relative_expr_t, Sim>::value, void*> = nullptr>
double get_infections_relative(const SimulationNode<Sim>& /*node*/, double /*t*/,
                               const Eigen::Ref<const Eigen::VectorXd>& /*y*/)
{
    assert(false && "Overload get_infections_relative for your own model/simulation if you want to use dynamic NPIs.");
    return 0;
}
template <class Sim, std::enable_if_t<is_expression_valid<get_infections_relative_expr_t, Sim>::value, void*> = nullptr>
double get_infections_relative(const SimulationNode<Sim>& node, double t, const Eigen::Ref<const Eigen::VectorXd>& y)
{
    return get_infections_relative(node.get_simulation(), t, y);
}

/**
 * detect a get_migration_factors function for the Model type.
 */
template <class Sim>
using get_migration_factors_expr_t = decltype(get_migration_factors(
    std::declval<const Sim&>(), std::declval<double>(), std::declval<const Eigen::Ref<const Eigen::VectorXd>&>()));

/**
 * Get an additional migration factor.
 * The absolute migration for each compartment is computed by c_i * y_i * f_i, wher c_i is the coefficient set in 
 * MigrationParameters, y_i is the current compartment population, f_i is the factor returned by this function.
 * This factor is optional, default 1.0. If you need to adjust migration in that way, overload get_migration_factors(model, t, y) 
 * for your Model type so that can be found with argument-dependent lookup.
 * @param node a node of a migration graph.
 * @param y the current value of the simulation.
 * @param t the current simulation time
 * @return a vector expression, same size as y, with the factor for each compartment.
 */
template <class Sim, std::enable_if_t<!is_expression_valid<get_migration_factors_expr_t, Sim>::value, void*> = nullptr>
auto get_migration_factors(const SimulationNode<Sim>& /*node*/, double /*t*/,
                           const Eigen::Ref<const Eigen::VectorXd>& y)
{
    return Eigen::VectorXd::Ones(y.rows());
}
template <class Sim, std::enable_if_t<is_expression_valid<get_migration_factors_expr_t, Sim>::value, void*> = nullptr>
auto get_migration_factors(const SimulationNode<Sim>& node, double t, const Eigen::Ref<const Eigen::VectorXd>& y)
{
    return get_migration_factors(node.get_simulation(), t, y);
}

/**
 * detect a get_migration_factors function for the Model type.
 */
template <class Sim>
using test_commuters_expr_t = decltype(test_commuters(
    std::declval<Sim&>(), std::declval<Eigen::Ref<const Eigen::VectorXd>&>(), std::declval<double>()));

/**
 * Test persons when migrating from their source node.
 * May transfer persons between compartments, e.g., if an infection was detected.
 * This feature is optional, default implementation does nothing.
 * In order to support this feature for your model, implement a test_commuters overload 
 * that can be found with argument-dependent lookup.
 * @param node a node of a migration graph.
 * @param migrated mutable reference to vector of persons per compartment that migrate.
 * @param t the current simulation time.
 */
template <class Sim, std::enable_if_t<!is_expression_valid<test_commuters_expr_t, Sim>::value, void*> = nullptr>
void test_commuters(SimulationNode<Sim>& /*node*/, Eigen::Ref<Eigen::VectorXd> /*migrated*/, double /*time*/)
{
}
template <class Sim, std::enable_if_t<is_expression_valid<test_commuters_expr_t, Sim>::value, void*> = nullptr>
void test_commuters(SimulationNode<Sim>& node, Eigen::Ref<Eigen::VectorXd> migrated, double time)
{
    return test_commuters(node.get_simulation(), migrated, time);
}

template <class Sim>
void MigrationEdge::apply_migration(double t, double dt, SimulationNode<Sim>& node_from, SimulationNode<Sim>& node_to)
{
    //check dynamic npis
    if (m_t_last_dynamic_npi_check == -std::numeric_limits<double>::infinity()) {
        m_t_last_dynamic_npi_check = node_from.get_t0();
    }

    auto& dyn_npis = m_parameters.get_dynamic_npis_infected();
    if (dyn_npis.get_thresholds().size() > 0 &&
        floating_point_greater_equal(t, m_t_last_dynamic_npi_check + dyn_npis.get_interval().get())) {
        auto inf_rel = get_infections_relative(node_from, t, node_from.get_last_state()) * dyn_npis.get_base_value();
        auto exceeded_threshold = dyn_npis.get_max_exceeded_threshold(inf_rel);
        if (exceeded_threshold != dyn_npis.get_thresholds().end() &&
            (exceeded_threshold->first > m_dynamic_npi.first ||
             t > double(m_dynamic_npi.second))) { //old NPI was weaker or is expired
            auto t_end    = SimulationTime(t + double(dyn_npis.get_duration()));
            m_dynamic_npi = std::make_pair(exceeded_threshold->first, t_end);
            implement_dynamic_npis(
                m_parameters.get_coefficients(), exceeded_threshold->second, SimulationTime(t), t_end, [this](auto& g) {
                    return make_migration_damping_vector(m_parameters.get_coefficients().get_shape(), g);
                });
        }
        m_t_last_dynamic_npi_check = t;
    }

    //returns
    for (Eigen::Index i = m_return_times.get_num_time_points() - 1; i >= 0; --i) {
        if (m_return_times.get_time(i) <= t) {
            auto v0 = find_value_reverse(node_to.get_result(), m_migrated.get_time(i), 1e-10, 1e-10);
            assert(v0 != node_to.get_result().rend() && "unexpected error.");
            calculate_migration_returns(m_migrated[i], node_to.get_simulation(), *v0, m_migrated.get_time(i), dt);

            //the lower-order return calculation may in rare cases produce negative compartments,
            //especially at the beginning of the simulation.
            //fix by subtracting the supernumerous returns from the biggest compartment of the age group.
            Eigen::VectorXd remaining_after_return = (node_to.get_result().get_last_value() - m_migrated[i]).eval();
            for (Eigen::Index j = 0; j < node_to.get_result().get_last_value().size(); ++j) {
                if (remaining_after_return(j) < 0) {
                    auto num_comparts = (Eigen::Index)Sim::Model::Compartments::Count;
                    auto group        = Eigen::Index(j / num_comparts);
                    auto compart      = j % num_comparts;
                    log(remaining_after_return(j) < -1e-3 ? LogLevel::warn : LogLevel::info,
                        "Underflow during migration returns at time {}, compartment {}, age group {}: {}", t, compart,
                        group, remaining_after_return(j));
                    Eigen::Index max_index;
                    slice(remaining_after_return, {group * num_comparts, num_comparts}).maxCoeff(&max_index);
                    log_info("Transferring to compartment {}", max_index);
                    max_index += group * num_comparts;
                    m_migrated[i](max_index) -= remaining_after_return(j);
                    m_migrated[i](j) += remaining_after_return(j);
                }
            }
            node_from.get_result().get_last_value() += m_migrated[i];
            node_to.get_result().get_last_value() -= m_migrated[i];
            m_migrated.remove_time_point(i);
            m_return_times.remove_time_point(i);
        }
    }

    if (!m_return_migrated && (m_parameters.get_coefficients().get_matrix_at(t).array() > 0.0).any()) {
        //normal daily migration
        m_migrated.add_time_point(
            t, (node_from.get_last_state().array() * m_parameters.get_coefficients().get_matrix_at(t).array() *
                get_migration_factors(node_from, t, node_from.get_last_state()).array())
                   .matrix());
        m_return_times.add_time_point(t + dt);

        test_commuters(node_from, m_migrated.get_last_value(), t);

        node_to.get_result().get_last_value() += m_migrated.get_last_value();
        node_from.get_result().get_last_value() -= m_migrated.get_last_value();
    }
    m_return_migrated = !m_return_migrated;
}

/**
 * edge functor for migration simulation.
 * @see SimulationNode::evolve
 */
template <class Sim>
void evolve_model(double t, double dt, SimulationNode<Sim>& node)
{
    node.evolve(t, dt);
}

/**
 * edge functor for migration simulation.
 * @see MigrationEdge::apply_migration
 */
template <class Sim>
void apply_migration(double t, double dt, MigrationEdge& migrationEdge, SimulationNode<Sim>& node_from,
                     SimulationNode<Sim>& node_to)
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
 * @{
 */
template <class Sim>
GraphSimulation<Graph<SimulationNode<Sim>, MigrationEdge>>
make_migration_sim(double t0, double dt, const Graph<SimulationNode<Sim>, MigrationEdge>& graph)
{
    return make_graph_sim(t0, dt, graph, &evolve_model<Sim>, &apply_migration<Sim>);
}

template <class Sim>
GraphSimulation<Graph<SimulationNode<Sim>, MigrationEdge>>
make_migration_sim(double t0, double dt, Graph<SimulationNode<Sim>, MigrationEdge>&& graph)
{
    return make_graph_sim(t0, dt, std::move(graph), &evolve_model<Sim>, &apply_migration<Sim>);
}

/** @} */

} // namespace mio

#endif //METAPOPULATION_MOBILITY_INSTANT_H
