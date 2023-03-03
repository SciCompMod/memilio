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
#ifndef EPI_MOBILITY_GRAPH_SIMULATION_H
#define EPI_MOBILITY_GRAPH_SIMULATION_H

#include "memilio/mobility/graph.h"
#include "memilio/utils/random_number_generator.h"

namespace mio
{

/**
 * @brief abstract simulation on a graph with alternating node and edge actions
 */
template <class Graph, class edge_f = std::function<void(double, double, typename Graph::EdgeProperty&,
                                                         typename Graph::NodeProperty&, typename Graph::NodeProperty&)>>
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

    // GraphSimulation(double t0, double dt, Graph&& g, const node_function& node_func)
    //     : m_t(t0)
    //     , m_dt(dt)
    //     , m_graph(std::move(g))
    //     , m_node_func(node_func)
    // {
    // }

    // GraphSimulation(double t0, double dt, Graph& g, const node_function& node_func)
    //     : m_t(t0)
    //     , m_dt(dt)
    //     , m_graph(g)
    //     , m_node_func(node_func)
    // {
    // }

    void advance(double t_max = 1.0)
    {
        auto dt = m_dt;
        while (m_t < t_max) {
            if (m_t + dt > t_max) {
                dt = t_max - m_t;
            }

            for (auto& n : m_graph.nodes()) {
                m_node_func(m_t, dt, n.property);
            }

            m_t += dt;

            for (auto& e : m_graph.edges()) {
                m_edge_func(m_t, dt, e.property, m_graph.nodes()[e.start_node_idx].property,
                            m_graph.nodes()[e.end_node_idx].property);
            }
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
};

template <class Graph>
class GraphSimulation : public GraphSimulationBase<Graph>
{
    using GraphSimulationBase<Graph>::GraphSimulationBase;
};

template <class NodeProperty, class EdgeProperty>
class GraphSimulationStochastic
    : public GraphSimulationBase<
          mio::Graph<NodeProperty, EdgeProperty>,
          std::function<void(typename EdgeProperty&, size_t, typename NodeProperty&, typename NodeProperty&)>>
{
    using GraphSimulationBase<mio::Graph<NodeProperty, EdgeProperty>>::GraphSimulationBase;

public:
    void advance(double t_max = 1.0)
    {
        //draw normalized waiting time
        ScalarType normalized_waiting_time = ExponentialDistribution<ScalarType>::get_instance(1.0);
        std::vector<ScalarType> dt_cand(m_graph.nodes().size());
        ScalarType cumulative_rate = get_cumulative_transition_rate();
        while (m_t < t_max) {
            m_dt              = std::min({m_dt, tmax - m_t});
            int node_iterator = 0;
            if (cumulative_rate * m_dt >
                normalized_waiting_time) { //at least one transition event during current time step
                double remaining_frac =
                    1 - normalized_waiting_time / (cumulative_rate * m_dt); //remaining fraction of current time step

                int parameters_per_edge = m_graph.edges()[0].property.get_parameters().get_coefficients().get_shape();

                //draw transition event
                size_t event = mio::DiscreteDistribution<size_t>::get_instance()(get_rates());
                //edge that performs transition event
                auto& event_edge = m_graph.edges()[int(event / parameters_per_edge)];
                //index for compartment and age group migrating
                auto flat_index = event % parameters_per_edge;

                //advance nodes until t + (waiting_time / cumulative_rate)
                for (auto& n : m_graph.nodes()) {
                    m_node_func(m_t, normalized_waiting_time / cumulative_rate, n.property);
                    //get new dt of each node
                    dt_cand[node_iterator] = n.property.get_simulation().get_integrator().get_dt();
                }
                //advance time
                m_t += normalized_waiting_time / cumulative_rate;
                //calculate new cumulative rate
                cumulative_rate = get_cumulative_transition_rate();
                //perform transition
                m_edge_func(event_edge.property, flat_index, m_graph.nodes()[event_edge.start_node_idx].property,
                            m_graph.nodes()[event_edge.end_node_idx].property);
                bool another_transition_event = true;
                while (another_transition_event) {
                    //draw new waiting time
                    normalized_waiting_time = ExponentialDistribution<ScalarType>::get_instance(1.0);
                    if (remaining_frac * cumulative_rate * m_dt >
                        normalized_waiting_time) { //another transition event takes place in current time step
                        remaining_frac -= normalized_waiting_time / (cumulative_rate * m_dt);
                        //draw transition event
                        size_t event = mio::DiscreteDistribution<size_t>::get_instance()(get_rates());
                        //edge that performs transition event
                        auto& event_edge = m_graph.edges()[int(event / parameters_per_edge)];
                        //index for compartment and age group migrating
                        auto flat_index = event % parameters_per_edge;

                        //advance nodes until t + (waiting_time / cumulative_rate)
                        for (auto& n : m_graph.nodes()) {
                            m_node_func(m_t, normalized_waiting_time / cumulative_rate, n.property);
                            //get new dt of each node
                            dt_cand[node_iterator] = n.property.get_simulation().get_integrator().get_dt();
                        }
                        //advance time
                        m_t += normalized_waiting_time / cumulative_rate;
                        //calculate new cumulative rate
                        cumulative_rate = get_cumulative_transition_rate();
                        //perform transition
                        m_edge_func(event_edge.property, flat_index,
                                    m_graph.nodes()[event_edge.start_node_idx].property,
                                    m_graph.nodes()[event_edge.end_node_idx].property);
                    }
                    else { //no other transition event takes place in current time step
                        normalized_waiting_time -= reamining_frac * cumulative_rate * m_dt;
                        another_transition_event = false;
                    }
                }
            }
            else { //no transition event in current time step
                normalized_waiting_time -= cumulative_rate * m_dt; //reduce waiting time by current time step
            }

            //advance nodes until t+dt
            for (auto& n : m_graph.nodes()) {
                m_node_func(m_t, m_dt, n.property);
                //get new dt of each node
                dt_cand[node_iterator] = n.property.get_simulation().get_integrator().get_dt();
            }
            //advance time
            m_t += m_dt;
            //calculate new cumulative rate
            cumulative_rate = get_cumulative_transition_rate();
            //new dt ist the minimal dt of all nodes
            m_dt = *std::min_element(dt_cand.begin(), dt_cand.end());
        }
    }

private:
    ScalarType get_cumulative_transition_rate()
    {
        //compute current cumulative transition rate
        ScalarType cumulative_transition_rate = 0;
        for (auto& e : m_graph.edges()) {
            cumulative_transition_rate += e.get_transition_rates(m_graph.nodes()[e.start_node_idx].property).sum();
        }
        return cumulative_transition_rate;
    }

    std::vector<ScalarType> get_rates()
    {
        std::vector<ScalarType> rates;
        for (auto& e : m_graph.edges()) {
            std::vector<ScalarType> edge_rates = e.get_transition_rates(m_graph.nodes()[e.start_node_idx].property);
            for (auto& coeff : edge_rates) {
                rates.emplace_back(coeff);
            }
        }
        return rates;
    }
};

template <class Graph, class NodeF, class EdgeF>
auto make_graph_sim(double t0, double dt, Graph&& g, NodeF&& node_func, EdgeF&& edge_func)
{
    return GraphSimulation<std::decay_t<Graph>>(t0, dt, std::forward<Graph>(g), std::forward<NodeF>(node_func),
                                                std::forward<EdgeF>(edge_func));
}

} // namespace mio
#endif //EPI_MOBILITY_GRAPH_SIMULATION_H
