/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Daniel Abele, Henrik Zunker
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
#ifndef MIO_MOBILITY_GRAPH_SIMULATION_H
#define MIO_MOBILITY_GRAPH_SIMULATION_H

#include "memilio/mobility/graph.h"
#include "memilio/mobility/trade.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/compartments/feedback_simulation.h"
#include "memilio/geography/regions.h"
#include "models/smm/model.h"

namespace mio
{

/**
 * @brief abstract simulation on a graph with alternating node and edge actions
 */
template <class GraphT, class Timepoint, class Timespan, class edge_f, class node_f>
class GraphSimulationBase
{
public:
    using Graph = GraphT;

    using node_function = node_f;
    using edge_function = edge_f;

    GraphSimulationBase(Timepoint t0, Timespan dt, const Graph& g, const node_function& node_func,
                        const edge_function&& edge_func)
        : m_t(t0)
        , m_dt(dt)
        , m_graph(g)
        , m_node_func(node_func)
        , m_edge_func(edge_func)
    {
    }

    GraphSimulationBase(Timepoint t0, Timespan dt, Graph&& g, const node_function& node_func,
                        const edge_function&& edge_func)
        : m_t(t0)
        , m_dt(dt)
        , m_graph(std::move(g))
        , m_node_func(node_func)
        , m_edge_func(edge_func)
    {
    }

    Timepoint get_t() const
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

protected:
    Timepoint m_t;
    Timespan m_dt;
    Graph m_graph;
    node_function m_node_func;
    edge_function m_edge_func;
};

template <typename FP, class Graph, class Timepoint, class Timespan,
          class edge_f = void (*)(Timepoint, Timespan, typename Graph::EdgeProperty&, typename Graph::NodeProperty&,
                                  typename Graph::NodeProperty&),
          class node_f = void (*)(Timepoint, Timespan, typename Graph::NodeProperty&)>
class GraphSimulation : public GraphSimulationBase<Graph, Timepoint, Timespan, edge_f, node_f>
{
    using Base = GraphSimulationBase<Graph, Timepoint, Timespan, edge_f, node_f>;
    using Base::Base;

public:
    void advance(Timepoint t_max = 1.0)
    {
        auto dt = Base::m_dt;
        while (Base::m_t < t_max) {
            if (Base::m_t + dt > t_max) {
                dt = t_max - Base::m_t;
            }

            for (auto& n : Base::m_graph.nodes()) {
                Base::m_node_func(Base::m_t, dt, n.property);
            }

            Base::m_t += dt;

            for (auto& e : Base::m_graph.edges()) {
                Base::m_edge_func(Base::m_t, dt, e.property, Base::m_graph.nodes()[e.start_node_idx].property,
                                  Base::m_graph.nodes()[e.end_node_idx].property);
            }
        }
    }
};

template <typename FP, class Graph>
class JollyGraphSimulation
    : public GraphSimulationBase<Graph, FP, FP,
                                 std::function<void(FP, FP, typename Graph::EdgeProperty&,
                                                    typename Graph::NodeProperty&, typename Graph::NodeProperty&)>,
                                 std::function<void(FP, FP, typename Graph::NodeProperty&)>>
{
    using Base = GraphSimulationBase<Graph, FP, FP,
                                     std::function<void(FP, FP, typename Graph::EdgeProperty&,
                                                        typename Graph::NodeProperty&, typename Graph::NodeProperty&)>,
                                     std::function<void(FP, FP, typename Graph::NodeProperty&)>>;

    using node_function = typename Base::node_function;
    using edge_function = typename Base::edge_function;

public:
    JollyGraphSimulation(FP t0, FP dt, const Graph& g, const node_function& node_func, const edge_function&& edge_func)
        : Base(t0, dt, g, node_func, std::move(edge_func))
        , m_rates(Base::m_graph.edges().size() * Base::m_graph.edges()[0].property.get_parameters().size())
        , _trade(Base::m_graph)
        , infection_indices()
    {
    }

    JollyGraphSimulation(FP t0, FP dt, Graph&& g, const node_function& node_func, const edge_function&& edge_func)
        : Base(t0, dt, std::forward<Graph>(g), node_func, std::move(edge_func))
        , m_rates(Base::m_graph.edges().size() * Base::m_graph.edges()[0].property.get_parameters().size())
        , _trade(Base::m_graph)
        , infection_indices()
    {
    }

    void advance(FP t_max = 1.0)
    {

        using Model        = std::decay_t<decltype(Base::m_graph.nodes()[0].property.get_simulation().get_model())>;
        using AdoptionRate = mio::smm::AdoptionRates<ScalarType, typename Model::Status, mio::regions::Region>;

        auto dt              = Base::m_dt;
        auto next_trade_time = _trade.next_trade_time();
        while (Base::m_t < t_max) {
            if (Base::m_t + dt > t_max) {
                dt = t_max - Base::m_t;
            }
            if (next_trade_time <= Base::m_t) {
                _trade.execute_trades(Base::m_t);
                next_trade_time = _trade.next_trade_time();
            }

            for (auto& n : Base::m_graph.nodes()) {
                Base::m_node_func(Base::m_t, dt, n.property);
            }
            const auto current_infections = sum_current_infections();
            mio::log_debug("Current infections at time {}: {}", Base::m_t, current_infections);
            const auto culling_time = 6. / current_infections;
            for (auto& n : Base::m_graph.nodes()) {
                for (auto& index : culling_indices) {
                    n.property.get_simulation().get_model().parameters.template get<AdoptionRate>()[index].factor =
                        culling_time;
                }
            }

            Base::m_t += dt;

            get_values(m_rates);
            size_t j = 0;
            for (auto& e : Base::m_graph.edges()) {
                auto edge_rates = e.property.get_parameters().get_coefficients();
                for (const auto& [index, rate] : edge_rates) {
                    if (m_rates[j] > 0) {
                        Base::m_edge_func(index, m_rates[j], e.property,
                                          Base::m_graph.nodes()[e.start_node_idx].property,
                                          Base::m_graph.nodes()[e.end_node_idx].property);
                    }
                    j++;
                }
            }
        }
    }

    auto sum_nodes()
    {
        const auto size = Base::m_graph.nodes()[0].property.get_result().get_last_value().size();
        std::vector<double> results(size, 0.0);
        for (auto& n : Base::m_graph.nodes()) {
            const auto& node_result = n.property.get_result().get_last_value();
            assert(n.property.get_result().get_last_value().size() == size);
            for (int i = 0; i < size; i++) {
                results[i] += node_result[i];
            }
        }
        return results;
    }

    auto statistics_per_timestep()
    {
        const auto size = Base::m_graph.nodes()[0].property.get_result().get_last_value().size();
        std::vector<double> timepoints;
        // Collect all exchange timepoints => All simulations are stopped and write results at those

        for (auto& n : Base::m_graph.edges()) {
            auto local_timepoints = n.property.get_mobility_results().get_time_points();
            for (auto t : local_timepoints) {
                if (std::find(timepoints.begin(), timepoints.end(), t) == timepoints.end()) {
                    timepoints.push_back(t);
                }
            }
        }
        std::sort(timepoints.begin(), timepoints.end());
        auto results = TimeSeries<ScalarType>::zero(timepoints.size(), size, timepoints);
        for (auto& n : Base::m_graph.nodes()) {
            assert(n.property.get_result().get_last_value().size() == size);
            auto node_timeseries = n.property.get_result();
            for (Eigen::Index time_index = 0; time_index < node_timeseries.get_num_time_points(); ++time_index) {
                auto time  = node_timeseries.get_time(time_index);
                auto index = results.get_index_of_time(time);
                if (index == Eigen::Index(-1)) {
                    continue;
                }
                for (int i = 0; i < size; i++) {
                    results.get_value(index)[i] += node_timeseries.get_value(time_index)[i];
                }
            }
        }
        return results;
    }

    auto return_all_time_series()
    {
        std::vector<TimeSeries<ScalarType>> all_time_series;
        for (auto& n : Base::m_graph.nodes()) {
            all_time_series.push_back(n.property.get_result());
        }
        return all_time_series;
    }

    RandomNumberGenerator& get_rng()
    {
        return m_rng;
    }

    void get_values(std::vector<FP>& rates)
    {
        size_t j           = 0;
        auto& distribution = mio::PoissonDistribution<int>::get_instance();
        for (auto& e : Base::m_graph.edges()) {
            auto edge_rates   = e.property.get_parameters().get_coefficients();
            const auto result = Base::m_graph.nodes()[e.start_node_idx].property.get_result().get_last_value();
            for (const auto& [index, rate] : edge_rates) {
                const auto compartment_value = result[index];
                rates[j]                     = (compartment_value < 1.0) ? 0.0 : distribution(m_rng, rate);
                j++;
            }
        }
    }

    void set_infection_indices(const std::vector<size_t>& indices)
    {
        infection_indices = indices;
    }

    void set_culling_indices(const std::vector<size_t>& indices)
    {
        culling_indices = indices;
    }

private:
    FP sum_current_infections()
    {
        FP sum_inf = 0;
        for (auto& n : Base::m_graph.nodes()) {
            const auto& node_result = n.property.get_result().get_last_value();
            for (const auto& index : infection_indices) {
                sum_inf += node_result[index];
            }
        }
        return sum_inf;
    }

    std::vector<FP> m_rates;
    RandomNumberGenerator m_rng;
    mio::Trade<FP, Graph> _trade;
    std::vector<size_t> infection_indices;
    std::vector<size_t> culling_indices;
};

template <typename FP, class Graph>
class GraphSimulationStochastic
    : public GraphSimulationBase<Graph, FP, FP,
                                 std::function<void(typename Graph::EdgeProperty&, size_t,
                                                    typename Graph::NodeProperty&, typename Graph::NodeProperty&)>,
                                 std::function<void(FP, FP, typename Graph::NodeProperty&)>>
{
    using Base = GraphSimulationBase<Graph, FP, FP,
                                     std::function<void(typename Graph::EdgeProperty&, size_t,
                                                        typename Graph::NodeProperty&, typename Graph::NodeProperty&)>,
                                     std::function<void(FP, FP, typename Graph::NodeProperty&)>>;

    using node_function = typename Base::node_function;
    using edge_function = typename Base::edge_function;

public:
    GraphSimulationStochastic(FP t0, FP dt, const Graph& g, const node_function& node_func,
                              const edge_function&& edge_func)
        : Base(t0, dt, g, node_func, std::move(edge_func))
        , m_rates(Base::m_graph.edges().size() *
                  Base::m_graph.edges()[0].property.get_parameters().get_coefficients().get_shape().rows())
        , _trade(Base::m_graph)
    {
    }

    GraphSimulationStochastic(FP t0, FP dt, Graph&& g, const node_function& node_func, const edge_function&& edge_func)
        : Base(t0, dt, std::forward<Graph>(g), node_func, std::move(edge_func))
        , m_rates(Base::m_graph.edges().size() *
                  Base::m_graph.edges()[0].property.get_parameters().get_coefficients().get_shape().rows())
        , _trade(Base::m_graph)
    {
    }

    void advance(FP t_max)
    {
        //draw normalized waiting time
        ScalarType normalized_waiting_time = ExponentialDistribution<ScalarType>::get_instance()(m_rng, 1.0);
        std::vector<ScalarType> dt_cand(Base::m_graph.nodes().size());
        ScalarType cumulative_rate = 0; //cumulative transition rate
        size_t parameters_per_edge =
            size_t(Base::m_graph.edges()[0].property.get_parameters().get_coefficients().get_shape().rows());
        std::vector<ScalarType> transition_rates(parameters_per_edge * Base::m_graph.edges().size());
        auto next_trade_time = _trade.next_trade_time();
        while (Base::m_t < t_max) {
            Base::m_dt = std::min({Base::m_dt, t_max - Base::m_t, next_trade_time - Base::m_t});
            //calculate current transition rates and cumulative rate
            cumulative_rate = get_cumulative_transition_rate();
            if (cumulative_rate * Base::m_dt >
                normalized_waiting_time) { //at least one transition event during current time step
                do {
                    //evaluate rates
                    get_rates(m_rates);
                    //draw transition event
                    size_t event = mio::DiscreteDistribution<size_t>::get_instance()(m_rng, m_rates);
                    //edge that performs transition event
                    auto& event_edge = Base::m_graph.edges()[event / parameters_per_edge];
                    //index for compartment and age group moving
                    auto flat_index = event % parameters_per_edge;

                    //advance nodes until t + (waiting_time / cumulative_rate)
                    // #pragma omp parallel for
                    for (size_t node_iter = 0; node_iter < Base::m_graph.nodes().size(); ++node_iter) {
                        auto& node = Base::m_graph.nodes()[node_iter];
                        Base::m_node_func(Base::m_t, normalized_waiting_time / cumulative_rate, node.property);
                    }

                    //advance time
                    Base::m_t += normalized_waiting_time / cumulative_rate;

                    //reduce remaining time of current time step
                    Base::m_dt -= normalized_waiting_time / cumulative_rate;

                    //perform transition
                    Base::m_edge_func(event_edge.property, flat_index,
                                      Base::m_graph.nodes()[event_edge.start_node_idx].property,
                                      Base::m_graph.nodes()[event_edge.end_node_idx].property);

                    //calculate new cumulative rate
                    cumulative_rate = get_cumulative_transition_rate();

                    //draw new normalized waiting time
                    normalized_waiting_time = ExponentialDistribution<ScalarType>::get_instance()(m_rng, 1.0);

                } while (cumulative_rate * Base::m_dt > normalized_waiting_time);
            }
            else { //no transition event in current time step
                normalized_waiting_time -= cumulative_rate * Base::m_dt; //reduce waiting time by current time step
            }

            //advance nodes until t+dt
            for (size_t node_iter = 0; node_iter < Base::m_graph.nodes().size(); ++node_iter) {
                auto& node = Base::m_graph.nodes()[node_iter];
                Base::m_node_func(Base::m_t, Base::m_dt, node.property);
                //get new dt of each node
                dt_cand[node_iter] = node.property.get_simulation().get_dt();
            }
            if (next_trade_time <= Base::m_t) {
                _trade.execute_trades(Base::m_t);
                next_trade_time = _trade.next_trade_time();
            }

            //advance time
            Base::m_t += Base::m_dt;
            //new dt ist the minimal dt of all nodes
            Base::m_dt = *std::min_element(dt_cand.begin(), dt_cand.end());
        }
    }

    RandomNumberGenerator& get_rng()
    {
        return m_rng;
    }

private:
    FP get_cumulative_transition_rate()
    {
        //compute current cumulative transition rate
        FP cumulative_transition_rate = 0.0;
        for (auto& e : Base::m_graph.edges()) {
            cumulative_transition_rate +=
                e.property.get_transition_rates(Base::m_graph.nodes()[e.start_node_idx].property).sum();
        }
        return cumulative_transition_rate;
    }

    void get_rates(std::vector<FP>& rates)
    {
        size_t j = 0;
        for (auto& e : Base::m_graph.edges()) {
            auto edge_rates   = e.property.get_transition_rates(Base::m_graph.nodes()[e.start_node_idx].property);
            const auto result = Base::m_graph.nodes()[e.start_node_idx].property.get_result().get_last_value();
            for (Eigen::Index i = 0; i < edge_rates.size(); ++i) {
                const auto compartment_value = result[i]; // Move last value out?
                rates[j]                     = (compartment_value < 1.0) ? 0.0 : edge_rates(i);

                j++;
            }
        }
    }

    std::vector<FP> m_rates;
    RandomNumberGenerator m_rng;
    mio::Trade<FP, Graph> _trade;
};

template <typename FP, typename Timepoint, class Timespan, class Graph, class NodeF, class EdgeF>
auto make_graph_sim(Timepoint t0, Timespan dt, Graph&& g, NodeF&& node_func, EdgeF&& edge_func)
{
    return GraphSimulation<FP, std::decay_t<Graph>, Timepoint, Timespan, EdgeF, NodeF>(
        t0, dt, std::forward<Graph>(g), std::forward<NodeF>(node_func), std::forward<EdgeF>(edge_func));
}

template <typename FP, class Graph, class NodeF, class EdgeF>
auto make_graph_sim_stochastic(FP t0, FP dt, Graph&& g, NodeF&& node_func, EdgeF&& edge_func)
{
    return GraphSimulationStochastic<FP, std::decay_t<Graph>>(
        t0, dt, std::forward<Graph>(g), std::forward<NodeF>(node_func), std::forward<EdgeF>(edge_func));
}

template <typename FP, class Graph, class NodeF, class EdgeF>
auto make_graph_sim_jolly(FP t0, FP dt, Graph&& g, NodeF&& node_func, EdgeF&& edge_func)
{
    return JollyGraphSimulation<FP, std::decay_t<Graph>>(t0, dt, std::forward<Graph>(g), std::forward<NodeF>(node_func),
                                                         std::forward<EdgeF>(edge_func));
}

// FeedbackGraphSimulation is only allowed to be used with local FeedbackSimulation.
// Therefore, we use type traits to check if the type is a specialization of FeedbackSimulation
template <class T>
struct is_feedback_simulation : std::false_type {
};

template <typename FP, typename Sim, typename ContactPatterns>
struct is_feedback_simulation<FeedbackSimulation<FP, Sim, ContactPatterns>> : std::true_type {
};

template <typename FP, class Graph>
class FeedbackGraphSimulation
{
public:
    FeedbackGraphSimulation(Graph&& g, FP t0, FP dt)
        : m_graph(std::move(g))
        , m_t(t0)
        , m_dt(dt)
        , m_initialized(false)
        , m_global_icu_occupancy(0)
    {
        using SimT = decltype(m_graph.nodes()[0].property.get_simulation());
        static_assert(is_feedback_simulation<std::decay_t<SimT>>::value,
                      "Graph node simulation must be a FeedbackSimulation.");
    }

    void advance(FP t_max)
    {
        // Initialize global and regional ICU occupancy if not done yet
        if (!m_initialized) {
            calculate_global_icu_occupancy();
            calculate_regional_icu_occupancy();
            m_initialized = true;
        }

        while (m_t < t_max) {
            FP dt_eff = std::min(m_dt, t_max - m_t);

            // assign the regional and global ICU occupancy to each node
            distribute_icu_data();

            // apply feedback and advance simulation for each node
            for (auto& node : m_graph.nodes()) {
                node.property.get_simulation().advance(m_t + dt_eff, m_dt);
            }

            // apply mobility between nodes
            for (auto& edge : m_graph.edges()) {
                auto& node_from = m_graph.nodes()[edge.start_node_idx];
                auto& node_to   = m_graph.nodes()[edge.end_node_idx];
                edge.property.apply_mobility(m_t, m_dt, node_from.property, node_to.property);
            }

            // update ICU occupancy for each node
            update_global_icu_occupancy();
            update_regional_icu_occupancy();

            m_t += dt_eff;
        }
    }

    Graph& get_graph()
    {
        return m_graph;
    }

private:
    /**
     * @brief Calculates the global ICU occupancy based on the ICU history of all nodes.
     * The result is stored in m_global_icu_occupancy.
     */
    void calculate_global_icu_occupancy()
    {
        auto& first_node_sim   = m_graph.nodes()[0].property.get_simulation();
        auto& first_node_icu   = first_node_sim.get_parameters().template get<ICUOccupancyHistory<FP>>();
        m_global_icu_occupancy = mio::TimeSeries<FP>(first_node_icu.get_num_elements());
        m_global_icu_occupancy.add_time_point(0.0, Eigen::VectorXd::Zero(first_node_icu.get_num_elements()));

        FP total_population = 0;
        for (Eigen::Index i = 0; i < first_node_icu.get_num_time_points(); ++i) {
            Eigen::VectorXd sum = Eigen::VectorXd::Zero(first_node_icu.get_num_elements());
            for (auto& node : m_graph.nodes()) {
                auto& sim         = node.property.get_simulation();
                auto& icu_history = sim.get_parameters().template get<ICUOccupancyHistory<FP>>();
                sum += icu_history.get_value(i) * node.property.get_simulation().get_model().populations.get_total();
                total_population += node.property.get_simulation().get_model().populations.get_total();
            }
            m_global_icu_occupancy.add_time_point(first_node_icu.get_time(i), sum / total_population);
        }
    }

    /**
     * @brief Calculates the regional ICU occupancy for each region.
     * The calculation is based on the ICU history of the nodes in each region.
     * The results are stored in m_regional_icu_occupancy.
     */
    void calculate_regional_icu_occupancy()
    {
        std::unordered_map<int, FP> regional_population;
        for (auto& node : m_graph.nodes()) {
            auto region_id = mio::regions::get_state_id(node.id).get();
            regional_population[region_id] += node.property.get_simulation().get_model().populations.get_total();
        }

        auto& first_node_sim         = m_graph.nodes()[0].property.get_simulation();
        auto& first_node_icu         = first_node_sim.get_parameters().template get<ICUOccupancyHistory<FP>>();
        Eigen::Index num_time_points = first_node_icu.get_num_time_points();
        Eigen::Index num_elements    = first_node_icu.get_num_elements();

        // For each region: vector of summed ICU values for all time points
        std::unordered_map<int, std::vector<Eigen::VectorXd>> regional_sums;
        for (auto const& [region_id, _] : regional_population) {
            regional_sums[region_id] =
                std::vector<Eigen::VectorXd>(num_time_points, Eigen::VectorXd::Zero(num_elements));
        }

        // Sum up
        for (auto& node : m_graph.nodes()) {
            auto region_id    = mio::regions::get_state_id(node.id).get();
            auto& sim         = node.property.get_simulation();
            auto& icu_history = sim.get_parameters().template get<ICUOccupancyHistory<FP>>();
            FP pop            = node.property.get_simulation().get_model().populations.get_total();
            for (Eigen::Index i = 0; i < num_time_points; ++i) {
                regional_sums[region_id][i] += icu_history.get_value(i) * pop;
            }
        }

        // Initialize TimeSeries and insert values
        m_regional_icu_occupancy.clear();
        for (auto const& [region_id, sum_vec] : regional_sums) {
            mio::TimeSeries<FP> ts(num_elements);
            for (Eigen::Index i = 0; i < num_time_points; ++i) {
                ts.add_time_point(first_node_icu.get_time(i), sum_vec[i] / regional_population[region_id]);
            }
            m_regional_icu_occupancy.emplace(region_id, std::move(ts));
        }
    }

    /**
     * @brief Updates the global ICU occupancy with the latest values from all nodes.
     * A new time point is added to m_global_icu_occupancy.
     */
    void update_global_icu_occupancy()
    {
        FP total_population = 0;
        Eigen::VectorXd sum = Eigen::VectorXd::Zero(m_global_icu_occupancy.get_num_elements());
        for (auto& node : m_graph.nodes()) {
            auto& sim         = node.property.get_simulation();
            auto& icu_history = sim.get_parameters().template get<ICUOccupancyHistory<FP>>();
            sum += icu_history.get_last_value() * node.property.get_simulation().get_model().populations.get_total();
            total_population += node.property.get_simulation().get_model().populations.get_total();
        }
        m_global_icu_occupancy.add_time_point(m_t, sum / total_population);
    }

    /**
     * @brief Updates the regional ICU occupancy for each region with the latest values.
     * A new time point is added to each TimeSeries in m_regional_icu_occupancy.
     */
    void update_regional_icu_occupancy()
    {

        std::unordered_map<int, FP> regional_population;
        for (auto& [region_id, regional_data] : m_regional_icu_occupancy) {
            Eigen::VectorXd sum = Eigen::VectorXd::Zero(regional_data.get_num_elements());
            for (auto& node : m_graph.nodes()) {
                if (mio::regions::get_state_id(node.id).get() == region_id) {
                    auto& sim         = node.property.get_simulation();
                    auto& icu_history = sim.get_parameters().template get<ICUOccupancyHistory<FP>>();
                    sum += icu_history.get_last_value() *
                           node.property.get_simulation().get_model().populations.get_total();
                    regional_population[region_id] +=
                        node.property.get_simulation().get_model().populations.get_total();
                }
            }
            regional_data.add_time_point(m_t, sum / regional_population[region_id]);
        }
    }

    /**
     * @brief Distributes the calculated global and regional ICU occupancy data to each node.
     * This makes the data available for feedback mechanisms within the node-local simulations.
     */
    void distribute_icu_data()
    {
        for (auto& node : m_graph.nodes()) {
            auto region_id = mio::regions::get_state_id(node.id).get();
            auto& sim      = node.property.get_simulation();
            sim.set_regional_icu_occupancy(m_regional_icu_occupancy.at(region_id));
            sim.set_global_icu_occupancy(m_global_icu_occupancy);
        }
    }

    Graph m_graph;
    FP m_t;
    FP m_dt;
    bool m_initialized;
    mio::TimeSeries<FP> m_global_icu_occupancy;
    std::unordered_map<int, mio::TimeSeries<FP>> m_regional_icu_occupancy;
};

} // namespace mio
#endif //MIO_MOBILITY_GRAPH_SIMULATION_H
