/*
* Copyright (C) 2020-2025 MEmilio
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
template <typename FP, class Sim>
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

    Eigen::Ref<const Eigen::VectorX<FP>> get_last_state() const
    {
        return m_last_state;
    }

    FP get_t0() const
    {
        return m_t0;
    }

    void advance(FP t, FP dt)
    {
        m_simulation.advance(t + dt);
        m_last_state = m_simulation.get_result().get_last_value();
    }

private:
    Sim m_simulation;
    Eigen::VectorX<FP> m_last_state;
    FP m_t0;
};

/**
 * time dependent mobility coefficients.
 */
template <typename FP>
using MobilityCoefficients = DampingMatrixExpression<FP, VectorDampings<FP>>;

/**
 * sum of time dependent mobility coefficients.
 * differentiate between sources of mobility.
 */
template <typename FP>
using MobilityCoefficientGroup = DampingMatrixExpressionGroup<FP, MobilityCoefficients<FP>>;

/**
 * parameters that influence mobility.
 * @tparam FP the underlying floating point type, e.g., double
 */
template <typename FP>
class MobilityParameters
{
public:
    /**
     * constructor from mobility coefficients.
     * @param coeffs mobility coefficients
     */
    MobilityParameters(const MobilityCoefficientGroup<FP>& coeffs)
        : m_coefficients(coeffs)
        , m_saved_compartment_indices(0)
    {
    }

    /**
     * constructor from mobility coefficients.
     * @param coeffs mobility coefficients
     */
    MobilityParameters(const Eigen::VectorX<FP>& coeffs)
        : m_coefficients({MobilityCoefficients<FP>(coeffs)})
        , m_saved_compartment_indices(0)
    {
    }

    /**
     * @brief Constructor for initializing mobility parameters with coefficients from type `MobilityCoefficientGroup`
     * and specific save indices.
     *
     * @param[in] coeffs A group of mobility coefficients represented by a `MobilityCoefficientGroup` object, defining
     * how individuals move between nodes.
     * @param[in] save_indices A 2D vector of indices. The outer vector represents different sets of compartments.
     * Each inner vector represents a set of compartments whose data will be saved in the member `m_mobility_results`
     * of the `MobilityEdge` class during the simulation using the `add_mobility_result_time_point` function.
     */
    MobilityParameters(const MobilityCoefficientGroup<FP>& coeffs, const std::vector<std::vector<size_t>>& save_indices)
        : m_coefficients(coeffs)
        , m_saved_compartment_indices(save_indices)
    {
    }

    /**
     * @brief Constructor for initializing mobility parameters with coefficients from an Eigen Vector
     * and specific save indices.
     *
     * @param[in] coeffs An `Eigen::VectorX<FP>` containing mobility coefficients.
     * @param[in] save_indices A 2D vector of indices. The outer vector represents different sets of compartments.
     * Each inner vector represents a set of compartments whose data will be saved in the member `m_mobility_results`
     * of the `MobilityEdge` class during the simulation using the `add_mobility_result_time_point` function.
     */
    MobilityParameters(const Eigen::VectorX<FP>& coeffs, const std::vector<std::vector<size_t>>& save_indices)
        : m_coefficients({MobilityCoefficients<FP>(coeffs)})
        , m_saved_compartment_indices(save_indices)
    {
    }

    /**
     * equality comparison operators
     */
    //@{
    bool operator==(const MobilityParameters& other) const
    {
        return m_coefficients == other.m_coefficients;
    }
    bool operator!=(const MobilityParameters& other) const
    {
        return m_coefficients != other.m_coefficients;
    }
    //@}

    /**
     * Get/Setthe mobility coefficients.
     * The coefficients represent the (time-dependent) percentage of people moving
     * from one node to another by age and infection compartment.
     * @{
     */
    /**
     * @return the mobility coefficients.
     */
    const MobilityCoefficientGroup<FP>& get_coefficients() const
    {
        return m_coefficients;
    }
    MobilityCoefficientGroup<FP>& get_coefficients()
    {
        return m_coefficients;
    }
    /**
     * @param coeffs the mobility coefficients.
     */
    void set_coefficients(const MobilityCoefficientGroup<FP>& coeffs)
    {
        m_coefficients = coeffs;
    }

    /**
     * @brief Get the indices of compartments to be saved during mobility.
     *
     * This function returns a reference to the vector of `m_saved_compartment_indices`, which specifies the groups of
     * compartments that are saved in the member `m_mobility_results` of the `MobilityEdge` class during the simulation
     * using the `add_mobility_result_time_point` function.
     *
     * @return A reference to the 2D vector containing indices of compartments to be saved. The outer vector represents different sets of compartments.
     * Each inner vector represents a group of compartments defined by indices.
     */
    const auto& get_save_indices() const
    {
        return m_saved_compartment_indices;
    }

    /**
     * Get/Set dynamic NPIs that are implemented when relative infections exceed thresholds.
     * This feature is optional. The simulation model needs to overload the get_infected_relative function.
     * @{
     */
    /**
     * @return dynamic NPIs for relative infections.
     */
    const DynamicNPIs<FP>& get_dynamic_npis_infected() const
    {
        return m_dynamic_npis;
    }
    DynamicNPIs<FP>& get_dynamic_npis_infected()
    {
        return m_dynamic_npis;
    }
    /**
     * @param v dynamic NPIs for relative infections.
     */
    void set_dynamic_npis_infected(const DynamicNPIs<FP>& v)
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
        auto obj = io.create_object("MobilityParameters");
        obj.add_element("Coefficients", m_coefficients);
        obj.add_element("DynamicNPIs", m_dynamic_npis);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<MobilityParameters> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("MobilityParameters");
        auto c   = obj.expect_element("Coefficients", Tag<MobilityCoefficientGroup<FP>>{});
        auto d   = obj.expect_element("DynamicNPIs", Tag<DynamicNPIs<FP>>{});
        return apply(
            io,
            [](auto&& c_, auto&& d_) {
                MobilityParameters params(c_);
                params.set_dynamic_npis_infected(d_);
                return params;
            },
            c, d);
    }

private:
    MobilityCoefficientGroup<FP> m_coefficients; //one per group and compartment
    DynamicNPIs<FP> m_dynamic_npis;
    std::vector<std::vector<size_t>> m_saved_compartment_indices; // groups of indices from compartments to save
};

/**
 * represents the mobility between two nodes.
 */
template <typename FP>
class MobilityEdge
{
public:
    /**
     * @brief Create edge with coefficients.
     * @param coeffs % of people in each group and compartment that change node in each time step.
     */
    MobilityEdge(const MobilityParameters<FP>& params)
        : m_parameters(params)
        , m_mobile_population(params.get_coefficients().get_shape().rows())
        , m_return_times(0)
        , m_return_mobile_population(false)
        , m_saved_compartment_indices(params.get_save_indices())
        , m_mobility_results(m_saved_compartment_indices.size() + 1)
    {
    }

    /**
     * @brief Create edge with coefficients.
     *
     * @param[in] coeffs An `Eigen::VectorX<FP>` representing the percentage of people in each group and compartment
     * that change nodes in each time step.
     */
    MobilityEdge(const Eigen::VectorX<FP>& coeffs)
        : m_parameters(coeffs)
        , m_mobile_population(coeffs.rows())
        , m_return_times(0)
        , m_return_mobile_population(false)
        , m_saved_compartment_indices(0)
        , m_mobility_results(m_saved_compartment_indices.size() + 1)
    {
    }

    /**
     * @brief Create edge with coefficients as MobilityParameters object and a 2D vector of indices which determine which compartments are saved.
     *
     * @param[in] params A `MobilityParameters` object representing the percentage of people in each group and compartment
     * that change nodes in each time step.
     * @param[in] save_indices A 2D vector of indices. The outer vector represents different sets of compartments.
     * Each inner vector represents a group of indices to be saved.
     */
    MobilityEdge(const MobilityParameters<FP>& params, const std::vector<std::vector<size_t>>& save_indices)
        : m_parameters(params)
        , m_mobile_population(params.get_coefficients().get_shape().rows())
        , m_return_times(0)
        , m_return_mobile_population(false)
        , m_saved_compartment_indices(save_indices)
        , m_mobility_results(m_saved_compartment_indices.size() + 1)
    {
    }

    /**
     * @brief Create edge with coefficients and a 2D vector of indices which determine which compartments are saved.
     *
     * @param[in] coeffs An `Eigen::VectorX<FP>` representing the percentage of people in each group and compartment that migrate
     * in each time step.
     * @param[in] save_indices A 2D vector of indices. The outer vector represents different sets of compartments, while each
     * inner vector represents a group of indices for compartments to be saved.
     */
    MobilityEdge(const Eigen::VectorX<FP>& coeffs, const std::vector<std::vector<size_t>>& save_indices)
        : m_parameters(coeffs)
        , m_mobile_population(coeffs.rows())
        , m_return_times(0)
        , m_return_mobile_population(false)
        , m_saved_compartment_indices(save_indices)
        , m_mobility_results(m_saved_compartment_indices.size() + 1)
    {
    }

    /**
     * get the mobility parameters.
     */
    const MobilityParameters<FP>& get_parameters() const
    {
        return m_parameters;
    }

    /**
     * @brief Get the count of commuters in selected compartments, along with the total number of commuters.
     *
     * @return A reference to the TimeSeries object representing the mobility results.
     * @{
     */
    TimeSeries<FP>& get_mobility_results()
    {
        return m_mobility_results;
    }
    const TimeSeries<FP>& get_mobility_results() const
    {
        return m_mobility_results;
    }
    /** @} */

    /**
     * compute mobility from node_from to node_to.
     * mobility is based on coefficients.
     * The mobile population is added to the current state of node_to, subtracted from node_from.
     * on return, the mobile population (adjusted for infections) is subtracted from node_to, added to node_from.
     * @param t current time
     * @param dt last time step (fixed to 0.5 for mobility model)
     * @param node_from node that people changed from, return to
     * @param node_to node that people changed to, return from
     */
    template <class Sim>
    void apply_mobility(FP t, FP dt, SimulationNode<FP, Sim>& node_from, SimulationNode<FP, Sim>& node_to);

private:
    MobilityParameters<FP> m_parameters;
    TimeSeries<FP> m_mobile_population;
    TimeSeries<FP> m_return_times;
    bool m_return_mobile_population;
    FP m_t_last_dynamic_npi_check                   = -std::numeric_limits<FP>::infinity();
    std::pair<FP, SimulationTime<FP>> m_dynamic_npi = {-std::numeric_limits<FP>::max(), SimulationTime<FP>(0)};
    std::vector<std::vector<size_t>> m_saved_compartment_indices; // groups of indices from compartments to save
    TimeSeries<FP> m_mobility_results; // save results from edges + entry for the total number of commuters

    /**
     * @brief Computes a condensed version of `m_mobile_population` and stores it in `m_mobility_results`.
     *
     * The `m_mobility_results` then only contains commuters with infection states `InfectedNoSymptoms` and
     * `InfectedSymptoms`. Additionally, the total number of commuters is stored in the last entry of `m_mobility_results`.
     *
     * @param[in] t The current time.
     */
    void add_mobility_result_time_point(const FP t);
};

template <typename FP>
void MobilityEdge<FP>::add_mobility_result_time_point(const FP t)
{
    const size_t save_indices_size = this->m_saved_compartment_indices.size();
    if (save_indices_size > 0) {
        const auto& last_value = m_mobile_population.get_last_value();

        Eigen::VectorX<FP> condensed_values = Eigen::VectorX<FP>::Zero(save_indices_size + 1);

        // sum up the values of m_saved_compartment_indices for each group (e.g. Age groups)
        for (size_t i = 0; i < save_indices_size; ++i) {
            FP sum = 0.0;
            for (auto index : this->m_saved_compartment_indices[i]) {
                sum += last_value[index];
            }
            condensed_values[i] = sum;
        }

        // the last value is the sum of commuters
        condensed_values[save_indices_size] = m_mobile_population.get_last_value().sum();

        // Move the condensed values to the m_mobility_results time series
        m_mobility_results.add_time_point(t, std::move(condensed_values));
    }
}

/**
 * adjust number of people that changed node when they return according to the model.
 * E.g. during the time in the other node, some people who left as susceptible will return exposed.
 * Implemented for general compartmentmodel simulations, overload for your custom model if necessary
 * so that it can be found with argument-dependent lookup, i.e. in the same namespace as the model.
 * @param[inout] mobile_population number of people that changed node as input, number of people that return as output
 * @param params parameters of model in the node that the people changed to.
 * @param total total population in the node that the people changed to.
 * @param t time of mobility
 * @param dt time between mobility and return
 */
template <typename FP, class Sim, class = std::enable_if_t<is_compartment_model_simulation<FP, Sim>::value>>
void calculate_mobility_returns(Eigen::Ref<typename TimeSeries<FP>::Vector> mobile_population, const Sim& sim,
                                Eigen::Ref<const typename TimeSeries<FP>::Vector> total, FP t, FP dt)
{
    auto y0 = mobile_population.eval();
    auto y1 = mobile_population;
    EulerIntegratorCore<FP>().step(
        [&](auto&& y, auto&& t_, auto&& dydt) {
            sim.get_model().get_derivatives(total, y, t_, dydt);
        },
        y0, t, dt, y1);
}

/**
 * detect a get_infections_relative function for the Model type.
 */
template <typename FP, class Sim>
using get_infections_relative_expr_t = decltype(get_infections_relative(
    std::declval<const Sim&>(), std::declval<FP>(), std::declval<const Eigen::Ref<const Eigen::VectorX<FP>>&>()));

/**
 * get the percantage of infected people of the total population in the node
 * If dynamic NPIs are enabled, there needs to be an overload of get_infections_relative(model, y)
 * for the Model type that can be found with argument-dependent lookup. Ideally define get_infections_relative
 * in the same namespace as the Model type.
 * @param node a node of a mobility graph.
 * @param y the current value of the simulation.
 * @param t the current simulation time
 */
template <typename FP, class Sim,
          std::enable_if_t<!is_expression_valid<get_infections_relative_expr_t, Sim>::value, void*> = nullptr>
FP get_infections_relative(const SimulationNode<FP, Sim>& /*node*/, FP /*t*/,
                           const Eigen::Ref<const Eigen::VectorX<FP>>& /*y*/)
{
    assert(false && "Overload get_infections_relative for your own model/simulation if you want to use dynamic NPIs.");
    return 0;
}
template <typename FP, class Sim,
          std::enable_if_t<is_expression_valid<get_infections_relative_expr_t, Sim>::value, void*> = nullptr>
FP get_infections_relative(const SimulationNode<FP, Sim>& node, FP t, const Eigen::Ref<const Eigen::VectorX<FP>>& y)
{
    return get_infections_relative<FP, Sim>(node.get_simulation(), t, y);
}

/**
 * detect a get_mobility_factors function for the Model type.
 */
template <typename FP, class Sim>
using get_mobility_factors_expr_t = decltype(get_mobility_factors(
    std::declval<const Sim&>(), std::declval<FP>(), std::declval<const Eigen::Ref<const Eigen::VectorX<FP>>&>()));

/**
 * Get an additional mobility factor.
 * The absolute mobility for each compartment is computed by c_i * y_i * f_i, wher c_i is the coefficient set in
 * MobilityParameters, y_i is the current compartment population, f_i is the factor returned by this function.
 * This factor is optional, default 1.0. If you need to adjust mobility in that way, overload get_mobility_factors(model, t, y)
 * for your Model type so that can be found with argument-dependent lookup.
 * @param node a node of a mobility graph.
 * @param y the current value of the simulation.
 * @param t the current simulation time
 * @return a vector expression, same size as y, with the factor for each compartment.
 */
template <typename FP, class Sim,
          std::enable_if_t<!is_expression_valid<get_mobility_factors_expr_t, Sim>::value, void*> = nullptr>
auto get_mobility_factors(const SimulationNode<FP, Sim>& /*node*/, FP /*t*/,
                          const Eigen::Ref<const Eigen::VectorX<FP>>& y)
{
    return Eigen::VectorX<FP>::Ones(y.rows());
}
template <typename FP, class Sim,
          std::enable_if_t<is_expression_valid<get_mobility_factors_expr_t, Sim>::value, void*> = nullptr>
auto get_mobility_factors(const SimulationNode<FP, Sim>& node, FP t, const Eigen::Ref<const Eigen::VectorX<FP>>& y)
{
    return get_mobility_factors<FP, Sim>(node.get_simulation(), t, y);
}

/**
 * detect a get_mobility_factors function for the Model type.
 */
template <typename FP, class Sim>
using test_commuters_expr_t = decltype(test_commuters(
    std::declval<Sim&>(), std::declval<Eigen::Ref<const Eigen::VectorX<FP>>&>(), std::declval<FP>()));

/**
 * Test persons when moving from their source node.
 * May transfer persons between compartments, e.g., if an infection was detected.
 * This feature is optional, default implementation does nothing.
 * In order to support this feature for your model, implement a test_commuters overload
 * that can be found with argument-dependent lookup.
 * @param node a node of a mobility graph.
 * @param mobile_population mutable reference to vector of persons per compartment that change nodes.
 * @param t the current simulation time.
 */
template <typename FP, class Sim,
          std::enable_if_t<!is_expression_valid<test_commuters_expr_t, Sim>::value, void*> = nullptr>
void test_commuters(SimulationNode<FP, Sim>& /*node*/, Eigen::Ref<Eigen::VectorX<FP>> /*mobile_population*/,
                    FP /*time*/)
{
}
template <typename FP, class Sim,
          std::enable_if_t<is_expression_valid<test_commuters_expr_t, Sim>::value, void*> = nullptr>
void test_commuters(SimulationNode<FP, Sim>& node, Eigen::Ref<Eigen::VectorX<FP>> mobile_population, FP time)
{
    return test_commuters<FP, Sim>(node.get_simulation(), mobile_population, time);
}

template <typename FP>
template <class Sim>
void mio::MobilityEdge<FP>::apply_mobility(FP t, FP dt, SimulationNode<FP, Sim>& node_from,
                                           SimulationNode<FP, Sim>& node_to)
{
    //check dynamic npis
    if (m_t_last_dynamic_npi_check == -std::numeric_limits<FP>::infinity()) {
        m_t_last_dynamic_npi_check = node_from.get_t0();
    }

    auto& dyn_npis = m_parameters.get_dynamic_npis_infected();
    if (dyn_npis.get_thresholds().size() > 0 &&
        floating_point_greater_equal<FP>(t, m_t_last_dynamic_npi_check + dyn_npis.get_interval().get())) {
        auto inf_rel =
            get_infections_relative<FP, Sim>(node_from, t, node_from.get_last_state()) * dyn_npis.get_base_value();
        auto exceeded_threshold = dyn_npis.get_max_exceeded_threshold(inf_rel);
        if (exceeded_threshold != dyn_npis.get_thresholds().end() &&
            (exceeded_threshold->first > m_dynamic_npi.first ||
             t > FP(m_dynamic_npi.second))) { //old NPI was weaker or is expired
            auto t_end    = SimulationTime<FP>(t + dyn_npis.get_duration().get());
            m_dynamic_npi = std::make_pair(exceeded_threshold->first, t_end);
            implement_dynamic_npis<FP>(m_parameters.get_coefficients(), exceeded_threshold->second,
                                       SimulationTime<FP>(t), t_end, [this](auto& g) {
                                           return make_mobility_damping_vector<FP>(
                                               m_parameters.get_coefficients().get_shape(), g);
                                       });
        }
        m_t_last_dynamic_npi_check = t;
    }

    //returns
    for (Eigen::Index i = m_return_times.get_num_time_points() - 1; i >= 0; --i) {
        if (m_return_times.get_time(i) <= t) {
            auto v0 = find_value_reverse<FP>(node_to.get_result(), m_mobile_population.get_time(i), 1e-10, 1e-10);
            assert(v0 != node_to.get_result().rend() && "unexpected error.");
            calculate_mobility_returns<FP, Sim>(m_mobile_population[i], node_to.get_simulation(), *v0,
                                                m_mobile_population.get_time(i), dt);

            //the lower-order return calculation may in rare cases produce negative compartments,
            //especially at the beginning of the simulation.
            //fix by subtracting the supernumerous returns from the biggest compartment of the age group.
            Eigen::VectorX<FP> remaining_after_return =
                (node_to.get_result().get_last_value() - m_mobile_population[i]).eval();
            for (Eigen::Index j = 0; j < node_to.get_result().get_last_value().size(); ++j) {
                if (remaining_after_return(j) < 0) {
                    auto num_comparts = (Eigen::Index)Sim::Model::Compartments::Count;
                    auto group        = Eigen::Index(j / num_comparts);
                    auto compart      = j % num_comparts;
                    log(remaining_after_return(j) < -1e-3 ? LogLevel::warn : LogLevel::info,
                        "Underflow during mobility returns at time {}, compartment {}, age group {}: {}", t, compart,
                        group, remaining_after_return(j));
                    Eigen::Index max_index;
                    slice(remaining_after_return, {group * num_comparts, num_comparts}).maxCoeff(&max_index);
                    log_info("Transferring to compartment {}", max_index);
                    max_index += group * num_comparts;
                    m_mobile_population[i](max_index) -= remaining_after_return(j);
                    m_mobile_population[i](j) += remaining_after_return(j);
                }
            }
            node_from.get_result().get_last_value() += m_mobile_population[i];
            node_to.get_result().get_last_value() -= m_mobile_population[i];
            add_mobility_result_time_point(t);
            m_mobile_population.remove_time_point(i);
            m_return_times.remove_time_point(i);
        }
    }

    if (!m_return_mobile_population &&
        (m_parameters.get_coefficients().get_matrix_at(SimulationTime<FP>(t)).array() > 0.0).any()) {
        //normal daily mobility
        m_mobile_population.add_time_point(
            t, (node_from.get_last_state().array() *
                m_parameters.get_coefficients().get_matrix_at(SimulationTime<FP>(t)).array() *
                get_mobility_factors<FP>(node_from, t, node_from.get_last_state()).array())
                   .matrix());
        m_return_times.add_time_point(t + dt);

        test_commuters<FP>(node_from, m_mobile_population.get_last_value(), t);

        node_to.get_result().get_last_value() += m_mobile_population.get_last_value();
        node_from.get_result().get_last_value() -= m_mobile_population.get_last_value();

        add_mobility_result_time_point(t);
    }
    m_return_mobile_population = !m_return_mobile_population;
}

/**
 * edge functor for mobility-based simulation.
 * @see SimulationNode::advance
 */
template <typename FP, class Sim>
void advance_model(FP t, FP dt, SimulationNode<FP, Sim>& node)
{
    node.advance(t, dt);
}

/**
 * edge functor for mobility-based simulation.
 * @see MobilityEdge::apply_mobility
 */
template <typename FP, class Sim>
void apply_mobility(FP t, FP dt, MobilityEdge<FP>& mobilityEdge, SimulationNode<FP, Sim>& node_from,
                    SimulationNode<FP, Sim>& node_to)
{
    mobilityEdge.apply_mobility(t, dt, node_from, node_to);
}

/**
 * create a mobility-based simulation.
 * After every second time step, for each edge a portion of the population corresponding to the coefficients of the edge
 * changes from one node to the other. In the next timestep, the mobile population returns to their "home" node.
 * Returns are adjusted based on the development in the target node.
 * @param t0 start time of the simulation
 * @param dt time step between mobility
 * @param graph set up for mobility-based simulation
 * @{
 */
template <typename FP, class Sim>
GraphSimulation<FP, Graph<SimulationNode<FP, Sim>, MobilityEdge<FP>>, FP, FP,
                void (*)(FP, FP, mio::MobilityEdge<FP>&, mio::SimulationNode<FP, Sim>&, mio::SimulationNode<FP, Sim>&),
                void (*)(FP, FP, mio::SimulationNode<FP, Sim>&)>
make_mobility_sim(FP t0, FP dt, const Graph<SimulationNode<FP, Sim>, MobilityEdge<FP>>& graph)
{
    return make_graph_sim<FP>(
        t0, dt, graph, static_cast<void (*)(FP, FP, SimulationNode<FP, Sim>&)>(&advance_model<FP, Sim>),
        static_cast<void (*)(FP, FP, MobilityEdge<FP>&, SimulationNode<FP, Sim>&, SimulationNode<FP, Sim>&)>(
            &apply_mobility<FP>));
}

template <typename FP, class Sim>
GraphSimulation<FP, Graph<SimulationNode<FP, Sim>, MobilityEdge<FP>>, FP, FP,
                void (*)(FP, FP, mio::MobilityEdge<FP>&, mio::SimulationNode<FP, Sim>&, mio::SimulationNode<FP, Sim>&),
                void (*)(FP, FP, mio::SimulationNode<FP, Sim>&)>
make_mobility_sim(FP t0, FP dt, Graph<SimulationNode<FP, Sim>, MobilityEdge<FP>>&& graph)
{
    return make_graph_sim<FP>(
        t0, dt, std::move(graph), static_cast<void (*)(FP, FP, SimulationNode<FP, Sim>&)>(&advance_model<FP, Sim>),
        static_cast<void (*)(FP, FP, MobilityEdge<FP>&, SimulationNode<FP, Sim>&, SimulationNode<FP, Sim>&)>(
            &apply_mobility<FP>));
}

/** @} */

/**
 * Create a graph simulation without mobility.
 *
 * Note that we set the time step of the graph simulation to infinity since we do not require any exchange between the
 * nodes. Hence, in each node, the simulation runs until tmax when advancing the simulation without interruption.
 *
 * @param t0 Start time of the simulation.
 * @param graph Set up for graph-based simulation.
 * @{
 */
template <typename FP, class Sim>
auto make_no_mobility_sim(FP t0, Graph<SimulationNode<FP, Sim>, MobilityEdge<FP>>& graph)
{
    using GraphSim = GraphSimulation<FP, Graph<SimulationNode<FP, Sim>, MobilityEdge<FP>>, FP, FP,
                                     void (*)(FP, FP, mio::MobilityEdge<FP>&, mio::SimulationNode<FP, Sim>&,
                                              mio::SimulationNode<FP, Sim>&),
                                     void (*)(FP, FP, mio::SimulationNode<FP, Sim>&)>;
    return GraphSim(t0, std::numeric_limits<FP>::infinity(), graph, &advance_model<FP, Sim>,
                    [](FP, FP, MobilityEdge<FP>&, SimulationNode<FP, Sim>&, SimulationNode<FP, Sim>&) {});
}

template <typename FP, class Sim>
auto make_no_mobility_sim(FP t0, Graph<SimulationNode<FP, Sim>, MobilityEdge<FP>>&& graph)
{
    using GraphSim = GraphSimulation<FP, Graph<SimulationNode<FP, Sim>, MobilityEdge<FP>>, FP, FP,
                                     void (*)(FP, FP, mio::MobilityEdge<FP>&, mio::SimulationNode<FP, Sim>&,
                                              mio::SimulationNode<FP, Sim>&),
                                     void (*)(FP, FP, mio::SimulationNode<FP, Sim>&)>;
    return GraphSim(t0, std::numeric_limits<FP>::infinity(), std::move(graph), &advance_model<FP, Sim>,
                    [](FP, FP, MobilityEdge<FP>&, SimulationNode<FP, Sim>&, SimulationNode<FP, Sim>&) {});
}

/** @} */

} // namespace mio

#endif //METAPOPULATION_MOBILITY_INSTANT_H
