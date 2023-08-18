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
#ifndef MIO_MOBILITY_GRAPH_SIMULATION_H
#define MIO_MOBILITY_GRAPH_SIMULATION_H

#include "memilio/compartments/simulation.h"
#include "memilio/config.h"
#include "memilio/mobility/graph.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/random_number_generator.h"
#include <cmath>
#include <cstddef>
#include <vector>

namespace mio
{

/**
 * @brief abstract simulation on a graph with alternating node and edge actions
 */
template <class Graph,
          class edge_f = std::function<void(double, double, typename Graph::EdgeProperty&,
                                            typename Graph::NodeProperty&, typename Graph::NodeProperty&, int)>>
class GraphSimulationBase
{
public:
    using node_function = std::function<void(double, double, typename Graph::NodeProperty&)>;

    using edge_function = edge_f;

    GraphSimulationBase(double t0, double dt, const Graph& g, const node_function& node_func,
                        const edge_function&& edge_func)
        : m_t(t0)
        , m_dt(dt)
        , m_graph(g)
        , m_node_func(node_func)
        , m_edge_func(edge_func)
    {
    }

    GraphSimulationBase(double t0, double dt, Graph&& g, const node_function& node_func,
                        const edge_function&& edge_func)
        : m_t(t0)
        , m_dt(dt)
        , m_graph(std::move(g))
        , m_node_func(node_func)
        , m_edge_func(edge_func)
    {
    }

    inline double round_second_decimal(double x)
    {
        return std::round(x * 100.0) / 100.0;
    }

    void advance(double t_max = 1.0)
    {
        // ScalarType traveltime_max = std::accumulate(
        //     m_graph.edges().begin(), m_graph.edges().end(), round_second_decimal(m_graph.edges()[0].traveltime),
        //     [&](ScalarType current_max, const auto& e) {
        //         return std::max(current_max, round_second_decimal(e.traveltime));
        //     });

        // ScalarType duration_max = std::accumulate(
        //     m_graph.nodes().begin(), m_graph.nodes().end(), round_second_decimal(m_graph.nodes()[0].stay_duration),
        //     [&](ScalarType current_max, const auto& n) {
        //         return std::max(current_max, round_second_decimal(n.stay_duration));
        //     });

        ScalarType dt_first_mobility =
            std::accumulate(m_graph.edges().begin(), m_graph.edges().end(), std::numeric_limits<ScalarType>::max(),
                            [&](ScalarType current_min, const auto& e) {
                                auto traveltime_per_region = round_second_decimal(e.traveltime / e.path.size());
                                if (traveltime_per_region < 0.01)
                                    traveltime_per_region = 0.01;
                                auto start_mobility =
                                    round_second_decimal(1 - 2 * (traveltime_per_region * e.path.size()) -
                                                         m_graph.nodes()[e.end_node_idx].stay_duration);
                                if (start_mobility < 0) {
                                    start_mobility = 0.;
                                }
                                return std::min(current_min, start_mobility);
                            });

        // ScalarType min_dt =
        //     std::accumulate(m_graph.edges().begin(), m_graph.edges().end(), std::numeric_limits<ScalarType>::max(),
        //                     [this](ScalarType minVal, const auto& e) {
        //                         ScalarType dt = e.traveltime / (e.path.size());
        //                         // round to third decimal
        //                         return std::min(minVal, std::trunc(dt * 1000.) / 1000.);
        //                     });

        // set population to zero in mobility nodes before starting
        for (auto& n : m_graph.nodes()) {
            n.node_pt.get_result().get_last_value().setZero();
            n.node_pt.get_simulation().get_model().populations.set_total(0);
        }

        auto min_dt    = 0.01;
        double t_begin = m_t - 1.;

        // calc schedule for each edge
        precompute_schedule();
        const double epsilon = std::numeric_limits<double>::epsilon();
        while (m_t - epsilon < t_max) {

            t_begin += 1;
            if (m_t + dt_first_mobility > t_max) {
                dt_first_mobility = t_max - m_t;
                for (auto& n : m_graph.nodes()) {
                    m_node_func(m_t, dt_first_mobility, n.property);
                }
                break;
            }

            size_t indx_schedule = 0;
            while (t_begin + 1 > m_t + 1e-10) {
                size_t indx_edge = 0;
                for (auto& e : m_graph.edges()) {
                    // find start of mobility
                    size_t first_mobility = std::distance(mobility_schedule_edges[indx_edge].begin(),
                                                          std::find(mobility_schedule_edges[indx_edge].begin(),
                                                                    mobility_schedule_edges[indx_edge].end(), true));
                    if (indx_schedule >= first_mobility) {
                        auto& node_from = mobility_schedule_edges[indx_edge][indx_schedule]
                                              ? m_graph.nodes()[schedule_edges[indx_edge][indx_schedule]].node_pt
                                              : m_graph.nodes()[schedule_edges[indx_edge][indx_schedule]].property;

                        auto& node_to = mobility_schedule_edges[indx_edge][indx_schedule + 1]
                                            ? m_graph.nodes()[schedule_edges[indx_edge][indx_schedule + 1]].node_pt
                                            : m_graph.nodes()[schedule_edges[indx_edge][indx_schedule + 1]].property;

                        if (indx_schedule == first_mobility)
                            m_edge_func(m_t, min_dt, e.property, node_from, node_to, 0);
                        else if (indx_schedule > first_mobility &&
                                 indx_schedule < mobility_schedule_edges[indx_edge].size() - 1)
                            m_edge_func(m_t, min_dt, e.property, node_from, node_to, 1);
                        else
                            m_edge_func(m_t, min_dt, e.property, node_from,
                                        m_graph.nodes()[schedule_edges[indx_edge][indx_schedule]].property, 2);
                    }
                    indx_edge++;
                }

                // integrate nodes after mobility updates
                for (auto& n : m_graph.nodes()) {
                    m_node_func(m_t, min_dt, n.property);
                    m_node_func(m_t, min_dt, n.node_pt);
                }
                indx_schedule++;
                m_t += min_dt;
            }
            // set each compartment zero for all mobility nodes since we only model daily mobility
            for (auto& n : m_graph.nodes()) {
                n.node_pt.get_result().get_last_value().setZero();
            }
            std::cout << "aktuell bei t = " << m_t << "\n";
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

protected:
    double m_t;
    double m_dt;
    Graph m_graph;
    node_function m_node_func;
    edge_function m_edge_func;

private:
    std::vector<std::vector<size_t>> schedule_edges;
    std::vector<std::vector<bool>> mobility_schedule_edges;

    void precompute_schedule()
    {
        schedule_edges.reserve(m_graph.edges().size());
        mobility_schedule_edges.reserve(m_graph.edges().size());

        const double epsilon = std::numeric_limits<double>::epsilon();

        for (auto& e : m_graph.edges()) {
            // 100 since we round to second decimal
            std::vector<size_t> schedule(100, 0.);
            std::vector<bool> is_mobility_node(100, false);

            double traveltime_per_region = round_second_decimal(e.traveltime / e.path.size());

            if (traveltime_per_region < 0.01)
                traveltime_per_region = 0.01;

            double start_mobility = round_second_decimal(0 + 1 - 2 * (traveltime_per_region * e.path.size()) -
                                                         m_graph.nodes()[e.end_node_idx].stay_duration);
            if (start_mobility < 0.0) {
                start_mobility = 0.;
            }

            double arrive_at = start_mobility + traveltime_per_region * e.path.size();

            // all values true between start mob und arrive
            std::fill(is_mobility_node.begin() + static_cast<int>((start_mobility + epsilon) * 100),
                      is_mobility_node.begin() + static_cast<int>((arrive_at + epsilon) * 100), true);
            int count_true       = std::count(is_mobility_node.begin(), is_mobility_node.end(), true);
            size_t current_index = static_cast<int>((start_mobility + epsilon) * 100);
            for (size_t county : e.path) {
                std::fill(schedule.begin() + current_index,
                          schedule.begin() + current_index +
                              static_cast<size_t>((traveltime_per_region + epsilon) * 100),
                          county);
                current_index += static_cast<size_t>((traveltime_per_region + epsilon) * 100);
            }

            // stay in destination county
            std::fill(schedule.begin() + current_index, schedule.end() - count_true, e.path[e.path.size() - 1]);

            // revert trip for return.
            std::fill(is_mobility_node.end() - count_true, is_mobility_node.end(), true);

            auto first_true = std::find(is_mobility_node.begin(), is_mobility_node.end(), true);
            auto last_true  = std::find(is_mobility_node.rbegin(), is_mobility_node.rend(), true);

            // If all values are false, we get an error
            if (first_true != is_mobility_node.end() && last_true != is_mobility_node.rend()) {
                int first_index = std::distance(is_mobility_node.begin(), first_true);
                std::vector<size_t> path_reversed(schedule.begin() + first_index,
                                                  schedule.begin() + first_index + count_true);
                std::reverse(path_reversed.begin(), path_reversed.end());
                std::copy(path_reversed.begin(), path_reversed.end(), schedule.end() - count_true);
                schedule_edges.push_back(schedule);
                mobility_schedule_edges.push_back(is_mobility_node);
            }
            else {
                log_error("Error in creating schedule.");
            }
        }
    }
}; // namespace mio

template <class Graph>
class GraphSimulation : public GraphSimulationBase<Graph>
{
    using GraphSimulationBase<Graph>::GraphSimulationBase;
};

template <class Graph>
class GraphSimulationStochastic
    : public GraphSimulationBase<Graph,
                                 std::function<void(typename Graph::EdgeProperty&, size_t,
                                                    typename Graph::NodeProperty&, typename Graph::NodeProperty&)>>
{
    using Base =
        GraphSimulationBase<Graph, std::function<void(typename Graph::EdgeProperty&, size_t,
                                                      typename Graph::NodeProperty&, typename Graph::NodeProperty&)>>;
    using node_function = typename Base::node_function;
    using edge_function = typename Base::edge_function;

public:
    GraphSimulationStochastic(double t0, double dt, const Graph& g, const node_function& node_func,
                              const edge_function&& edge_func)
        : Base(t0, dt, g, node_func, std::move(edge_func))
        , m_rates(Base::m_graph.edges().size() *
                  Base::m_graph.edges()[0].property.get_parameters().get_coefficients().get_shape().rows())
    {
    }

    GraphSimulationStochastic(double t0, double dt, Graph&& g, const node_function& node_func,
                              const edge_function&& edge_func)
        : Base(t0, dt, std::forward<Graph>(g), node_func, std::move(edge_func))
        , m_rates(Base::m_graph.edges().size() *
                  Base::m_graph.edges()[0].property.get_parameters().get_coefficients().get_shape().rows())
    {
    }

    void advance(double t_max)
    {
        //draw normalized waiting time
        ScalarType normalized_waiting_time = ExponentialDistribution<ScalarType>::get_instance()(1.0);
        std::vector<ScalarType> dt_cand(Base::m_graph.nodes().size());
        ScalarType cumulative_rate = 0; //cumulative transition rate
        size_t parameters_per_edge =
            size_t(Base::m_graph.edges()[0].property.get_parameters().get_coefficients().get_shape().rows());
        std::vector<ScalarType> transition_rates(parameters_per_edge * Base::m_graph.edges().size());
        while (Base::m_t < t_max) {
            Base::m_dt = std::min({Base::m_dt, t_max - Base::m_t});
            //calculate current transition rates and cumulative rate
            cumulative_rate = get_cumulative_transition_rate();
            if (cumulative_rate * Base::m_dt >
                normalized_waiting_time) { //at least one transition event during current time step
                do {
                    //evaluate rates
                    get_rates(m_rates);
                    //draw transition event
                    size_t event = mio::DiscreteDistribution<size_t>::get_instance()(m_rates);
                    //edge that performs transition event
                    auto& event_edge = Base::m_graph.edges()[event / parameters_per_edge];
                    //index for compartment and age group migrating
                    auto flat_index = event % parameters_per_edge;

                    //advance nodes until t + (waiting_time / cumulative_rate)
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
                    normalized_waiting_time = ExponentialDistribution<ScalarType>::get_instance()(1.0);

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

            //advance time
            Base::m_t += Base::m_dt;
            //new dt ist the minimal dt of all nodes
            Base::m_dt = *std::min_element(dt_cand.begin(), dt_cand.end());
        }
    }

private:
    ScalarType get_cumulative_transition_rate()
    {
        //compute current cumulative transition rate
        ScalarType cumulative_transition_rate = 0;
        for (auto& e : Base::m_graph.edges()) {
            cumulative_transition_rate +=
                e.property.get_transition_rates(Base::m_graph.nodes()[e.start_node_idx].property).sum();
        }
        return cumulative_transition_rate;
    }

    void get_rates(std::vector<ScalarType>& rates)
    {
        size_t j = 0;
        for (auto& e : Base::m_graph.edges()) {
            auto edge_rates = e.property.get_transition_rates(Base::m_graph.nodes()[e.start_node_idx].property);
            for (Eigen::Index i = 0; i < edge_rates.size(); ++i) {
                const auto compartment_value =
                    Base::m_graph.nodes()[e.start_node_idx].property.get_result().get_last_value()[i];
                rates[j] = (compartment_value < 1.) ? 0. : edge_rates(i);

                j++;
            }
        }
    }

    std::vector<ScalarType> m_rates;
};

template <class Graph, class NodeF, class EdgeF>
auto make_graph_sim(double t0, double dt, Graph&& g, NodeF&& node_func, EdgeF&& edge_func)
{
    return GraphSimulation<std::decay_t<Graph>>(t0, dt, std::forward<Graph>(g), std::forward<NodeF>(node_func),
                                                std::forward<EdgeF>(edge_func));
}

template <class Graph, class NodeF, class EdgeF>
auto make_graph_sim_stochastic(double t0, double dt, Graph&& g, NodeF&& node_func, EdgeF&& edge_func)
{
    return GraphSimulationStochastic<std::decay_t<Graph>>(
        t0, dt, std::forward<Graph>(g), std::forward<NodeF>(node_func), std::forward<EdgeF>(edge_func));
}

} // namespace mio
#endif //MIO_MOBILITY_GRAPH_SIMULATION_H
