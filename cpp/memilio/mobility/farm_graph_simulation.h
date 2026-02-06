/* 
* Copyright (C) 2020-2025 MEmilio
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

#ifndef MIO_FARM_GRAPH_SIMULATION_H
#define MIO_FARM_GRAPH_SIMULATION_H

#include "memilio/config.h"
#include "memilio/mobility/farm_simulation.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/trade.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/random_number_generator.h"
#include <cmath>
#include <map>
#include <numeric>
#include <queue>
#include "memilio/geography/regions.h"
#include "models/smm/parameters.h"

namespace mio
{
/**
 * @brief Variation of a graph simulation for farm nodes.
 * 
 * @tparam Graph 
 * @tparam Timepoint 
 * @tparam Timespan 
 * @tparam (*)(Timepoint, Timespan, typename Graph::EdgeProperty&, typename Graph::NodeProperty&,
 * typename Graph::NodeProperty&, mio::RandomNumberGenerator&) 
 * @tparam (*)(Timepoint, Timespan, typename Graph::NodeProperty&) 
 */
template <class Graph, class Timepoint = double, class Timespan = double,
          class edge_f = void (*)(Timepoint, Timespan, typename Graph::EdgeProperty&, typename Graph::NodeProperty&,
                                  typename Graph::NodeProperty&, mio::RandomNumberGenerator&),
          class node_f = void (*)(Timepoint, Timespan, typename Graph::NodeProperty&)>
class FarmSimulation : public GraphSimulationBase<Graph, Timepoint, Timespan, edge_f, node_f>
{
    using Base = GraphSimulationBase<Graph, Timepoint, Timespan, edge_f, node_f>;
    using Base::GraphSimulationBase;
    using Node = Graph::NodeProperty;

public:
    void advance(Timepoint t_max = 1.0)
    {
        mio::timing::AutoTimer<"Graph Simulation Advance"> timer;
        // auto dt = std::min(m_trade.next_trade_time() - Base::m_t, 1.0);
        auto dt     = Base::m_dt;
        auto& graph = Base::get_graph();
        while (Base::m_t < t_max) {
            dt = std::min(dt, t_max - Base::m_t);
            seed_infections(graph, Base::m_t);
            for (auto& n : graph.nodes()) {
                if (n.property.get_result().get_last_value()[2] > 0) {
                    n.property.set_infection_status(true);
                }
            }
            update_force_of_infection(graph);
            for (auto& n : graph.nodes()) {
                Base::m_node_func(Base::m_t, dt, n.property);
            }
            if (!m_standstill) {
                for (auto& n : graph.nodes()) {
                    update_population(n.property);
                }
            }
            for (auto& n : graph.nodes()) {
                if (std::accumulate(n.property.get_result().get_last_value().begin(),
                                    n.property.get_result().get_last_value().end(), 0) > 0) {
                    if (!n.property.is_suspicious()) {
                        check_for_cases(n.property, Base::m_t);
                    }
                    else {
                        if (n.property.get_date_confirmation() < Base::m_t + dt && !n.property.is_quarantined()) {
                            cull_node(n.id);
                            for (size_t index = 0; index < n.property.get_regional_neighbors()[1].size(); ++index) {
                                auto& neighbor = graph.nodes()[n.property.get_regional_neighbors()[1][index]];
                                neighbor.property.set_reg_zone_day(Base::m_t);
                            }
                            for (size_t index = 0; index < n.property.get_regional_neighbors()[0].size(); ++index) {
                                prev_cull_node(n.property.get_regional_neighbors()[0][index]);
                            }
                        }
                    }
                }
            }
            if (m_first_detection > Base::m_t) {
                for (auto& n : graph.nodes()) {
                    if (n.property.get_date_confirmation() <= Base::m_t) {
                        m_standstill      = true;
                        m_first_detection = Base::m_t;
                        break;
                    }
                }
            }
            else if (Base::m_t - m_first_detection >= 2) {
                m_standstill = false;
            }
            cull(dt);
            // dt = std::min(m_trade.next_trade_time() - Base::m_t, 1.0);
            Base::m_t += dt;
        }
    }

private:
    void seed_infections(Graph& graph, Timepoint t)
    {
        if (t == m_infection_dates[0]) {
            graph.nodes()[7658].property.get_result().get_last_value()[2] = 100;
            graph.nodes()[7658].property.get_result().get_last_value()[0] -= 100;
        }
        if (t == m_infection_dates[1]) {
            graph.nodes()[4117].property.get_result().get_last_value()[2] = 100;
            graph.nodes()[4117].property.get_result().get_last_value()[0] -= 100;
        }
        if (t == m_infection_dates[2]) {
            graph.nodes()[5369].property.get_result().get_last_value()[2] = 100;
            graph.nodes()[5369].property.get_result().get_last_value()[0] -= 100;
        }
    }
    /**
    * @brief Update the force of infection for all farms based on their infected neighbors.
    * 
    * @param graph 
    */
    void update_force_of_infection(Graph& graph)
    {
        using Model        = std::decay_t<decltype(Base::m_graph.nodes()[0].property.get_simulation().get_model())>;
        using AdoptionRate = mio::smm::AdoptionRates<ScalarType, typename Model::Status, mio::regions::Region>;
        for (auto& node : graph.nodes()) {
            if (std::accumulate(node.property.get_result().get_last_value().begin(),
                                node.property.get_result().get_last_value().end(), 0) == 0) {
                continue;
            }
            auto infection_pressure = 0.0;
            auto neighbors          = node.property.get_regional_neighbors();
            for (size_t index = 0; index < neighbors[2].size(); ++index) {
                auto& neighbour = graph.nodes()[neighbors[2][index]];
                if (neighbour.id != node.id && neighbour.property.get_infection_status()) {
                    infection_pressure +=
                        m_foi_inner_factor[neighbour.property.get_type()] *
                        (m_h0 /
                         (1 +
                          pow(node.property.get_location().distance(neighbour.property.get_location()).meters() / m_r0,
                              m_alpha)));
                }
            }
            node.property.get_simulation().get_model().parameters.template get<AdoptionRate>()[3].factor =
                m_foi_outer_factor[node.property.get_type()] * infection_pressure;
        }
    }

    void update_population(Node& node)
    {
        // ducks, layer
        if (node.get_type() < 3) {
            if (std::accumulate(node.get_result().get_last_value().begin(), node.get_result().get_last_value().end(),
                                0) > 0) {
                if (Base::m_t - node.get_population_date() + 1 > m_duration[node.get_type()]) {
                    for (auto& val : node.get_result().get_last_value()) {
                        val = 0;
                    }
                    node.set_slaughter_date(Base::m_t);
                    node.set_infection_status(false);
                }
            }
            else if (Base::m_t - node.get_slaughter_date() + 1 > m_downtime[node.get_type()] &&
                     !node.is_quarantined()) {
                node.get_result().get_last_value()[0] = 0.95 * node.get_capacity();
                node.set_population_date(Base::m_t);
            }
        }
        // broiler 2
        else if (node.get_type() == 4) {
            if (std::accumulate(node.get_result().get_last_value().begin(), node.get_result().get_last_value().end(),
                                0) > 0) {
                if (Base::m_t - node.get_population_date() + 1 > m_duration[4]) {
                    for (auto& val : node.get_result().get_last_value()) {
                        val = 0;
                    }
                    node.set_slaughter_date(Base::m_t);
                    node.set_infection_status(false);
                }
            }
        }
        // broiler 1
        else if (std::accumulate(node.get_result().get_last_value().begin(), node.get_result().get_last_value().end(),
                                 0) > 0) {
            if (Base::m_t - node.get_population_date() + 1 > m_duration[3]) {
                // If farm in regulation zone, only slaughter is allowed
                if (node.get_reg_zone_day() < Base::m_t && node.get_reg_zone_day() > Base::m_t - 28) {
                    for (auto& val : node.get_result().get_last_value()) {
                        val = 0;
                    }
                    node.set_slaughter_date(Base::m_t);
                    node.set_infection_status(false);
                }
                //...
            }
            else if (Base::m_t - node.get_slaughter_date() + 1 >= m_downtime[3] && !node.is_quarantined()) {
                node.get_result().get_last_value()[0] = 0.95 * node.get_capacity();
                node.set_population_date(Base::m_t);
            }
        }
    }

    void check_for_cases(Node& node, Timepoint t)
    {
        const auto& current_state = node.get_result().get_last_value();
        if (current_state[3] / std::accumulate(current_state.begin(), current_state.end(), 0.0) >
            m_suspicion_threshold) {
            node.set_date_suspicion(t);
            // if (mio::UniformDistribution<ScalarType>::get_instance()(m_rng, 0.0, 1.0) <
            //     1 - std::pow((1 - m_sensitivity), 20)) {
            node.set_date_confirmation(t + 2.0);
            // }
        }
    }

    /**
     * @brief Go through the culling queues and cull as many animals as possible within dt.
     * 
     * First the (reactive) culling queue is culled, then the preventive culling queue if it is past the 20th day (should be 2026-01-01).
     * Calls @ref cull_queue for both queues.
     * @param dt Time span for culling.
     */
    void cull(Timespan dt)
    {
        auto capacity = m_culling_capacity_per_day * dt;
        capacity      = cull_queue(m_culling_queue, capacity);
        if (Base::m_t > 20) {
            capacity = cull_queue(m_prev_culling_queue, capacity);
        }
    }

    /**
     * @brief Go through the culling queue and cull as many animals as possible with capacity.
     * 
     * Animals are killed in the order of compartments if not all animals of a farm can be culled within dt. 
     * @param queue Queue to cull from.
     * @param capacity Number of animals that can be culled.
     */
    ScalarType cull_queue(std::queue<std::pair<size_t, ScalarType>> queue, ScalarType capacity)
    {
        while (!queue.empty() && capacity > 0) {
            auto [node_id, day]     = queue.front();
            auto animals            = Base::m_graph.nodes()[node_id].property.get_result().get_last_value();
            auto num_living_animals = std::accumulate(animals.begin(), animals.end() - 1, 0);
            if (num_living_animals <= capacity) {
                mio::log_debug("Culling {} animals at node {} starting on day {}.", num_living_animals, node_id,
                               Base::m_t);
                for (auto index = 0; index < animals.size(); index++) {
                    animals[index] = 0;
                }
                queue.pop();
                capacity -= num_living_animals;
                Base::m_graph.nodes()[node_id].property.set_infection_status(false);
            }
            else {
                while (capacity > 0) {
                    mio::log_debug("Culling {} animals of node {} on day {}.", capacity, node_id, Base::m_t);
                    for (auto index = 0; index < animals.size(); index++) {
                        if (animals[index] < capacity) {
                            capacity -= animals[index];
                            animals[index] = 0;
                        }
                        else {
                            animals[index] -= capacity;
                            capacity = 0;
                            break;
                        }
                    }
                }
            }
        }
        return capacity;
    }

    /**
     * @brief Add a node to the culling queue and set it to quarantined.
     * 
     * @param node_id 
     */
    void cull_node(size_t node_id)
    {
        m_culling_queue.push(std::make_pair(node_id, Base::m_t));
        Base::m_graph.nodes()[node_id].property.set_quarantined(true);
    }

    /**
     * @brief Add a node to the preventive culling queue and set it to quarantined.
     * 
     * @param node_id 
     */
    void prev_cull_node(size_t node_id)
    {
        m_prev_culling_queue.push(std::make_pair(node_id, Base::m_t));
        Base::m_graph.nodes()[node_id].property.set_quarantined(true);
    }

    void vaccinate_node(size_t node_id)
    {
        m_vaccination_queue.push(std::make_pair(node_id, Base::m_t));
    }

    /**
     * get the mobility parameters.
     */
    // const Trade<ScalarType, Graph>& get_trade() const
    // {
    //     return m_trade;
    // }
public:
    auto sum_exchanges()
    {
        const auto size = Base::m_graph.edges()[0].property.get_mobility_results().get_num_elements();
        std::vector<double> results(size, 0.0);
        for (auto& n : Base::m_graph.edges()) {
            assert(n.property.get_mobility_results().get_num_elements() == size);
            for (auto result = n.property.get_mobility_results().begin();
                 result != n.property.get_mobility_results().end(); ++result) {
                for (int i = 0; i < size; i++) {
                    results[i] += (*result)[i];
                }
            }
        }
        return results;
    }

    auto exchanges_per_timestep()
    {
        const auto size = Base::m_graph.edges()[0].property.get_mobility_results().get_num_elements();
        std::vector<double> timepoints;
        // Collect all exchange timepoints

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
        for (auto& n : Base::m_graph.edges()) {
            assert(n.property.get_mobility_results().get_num_elements() == size);
            auto edge_result = n.property.get_mobility_results();
            // Add exchange data to big TimeSeries
            for (Eigen::Index time_index = 0; time_index < edge_result.get_num_time_points(); ++time_index) {
                auto time  = edge_result.get_time(time_index);
                auto index = results.get_index_of_time(time);
                if (index == Eigen::Index(-1)) {
                    continue;
                }
                for (int i = 0; i < size; i++) {
                    results.get_value(index)[i] += edge_result.get_value(time_index)[i];
                }
            }
        }
        return results;
    }

    auto sum_nodes()
    {
        const auto size = Base::m_graph.nodes()[0].property.get_result().get_num_elements();
        std::vector<double> results(size, 0.0);
        for (auto& n : Base::m_graph.nodes()) {
            assert(n.property.get_result().get_num_elements() == size);
            for (int i = 0; i < size; i++) {
                results[i] += n.property.get_result().get_last_value()[i];
            }
        }
        return results;
    }
    /* Is this correct or should it be:
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
*/

    auto statistics_per_timestep()
    {
        const auto size = Base::m_graph.nodes()[0].property.get_result().get_num_elements();
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
            assert(n.property.get_result().get_num_elements() == size);
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

    auto statistics_per_timestep(std::vector<size_t> node_indices)
    {
        assert(node_indices.size() > 0);
        const auto size = Base::m_graph.nodes()[node_indices[0]].property.get_result().get_num_elements();
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
        for (size_t node_index : node_indices) {
            auto node_timeseries = Base::m_graph.nodes()[node_index].property.get_result();
            assert(node_timeseries.get_num_elements() == size);
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

    /**
     * @brief Get the confirmation dates
     * 
     * Lists for all notes the date of confirmation value (-1 if no confirmation)
     * @return std::vector<int> 
     */
    std::vector<int> get_confirmation_dates()
    {
        std::vector<int> result;
        result.reserve(Base::m_graph.nodes().size());
        for (auto& n : Base::m_graph.nodes()) {
            result.push_back(n.property.get_date_confirmation());
        }
        return result;
    }

    // void apply_interventions()
    // {
    //     // Set quarantine
    //     // for (auto& n : Base::m_graph.nodes()) {
    //     //     if (n.property.get_result().get_last_value()[2] > 40) {
    //     //         mio::log_debug("Node {} is quarantined at time {}.", n.id, Base::m_t);
    //     //         n.property.set_quarantined(true);
    //     //     }
    //     //     else {
    //     //         mio::log_debug("Node {} is not quarantined at time {}.", n.id, Base::m_t);
    //     //         n.property.set_quarantined(false);
    //     //     }
    //     // }
    //     auto infectious_compartments = {2};
    //     for (auto& n : Base::m_graph.nodes()) {
    //         auto total_infections = 0;
    //         for (auto comp : infectious_compartments) {
    //             total_infections += n.property.get_result().get_last_value()[comp];
    //         }
    //         if (total_infections > 1) {
    //             for (size_t index = 0; index < n.property.get_regional_neighbors()[0].size(); ++index) {
    //                 const auto node_id = n.property.get_regional_neighbors()[0][index];
    //                 using Model =
    //                     std::decay_t<decltype(Base::m_graph.nodes()[node_id].property.get_simulation().get_model())>;
    //                 using AdoptionRate =
    //                     mio::smm::AdoptionRates<ScalarType, typename Model::Status, mio::regions::Region>;
    //                 Base::m_graph.nodes()[node_id]
    //                     .property.get_simulation()
    //                     .get_model()
    //                     .parameters.template get<AdoptionRate>()[1]
    //                     .factor += infectionrisk;
    //             }
    //         }
    //         if (total_infections > 4) {
    //             if (m_first_detection > Base::m_t) {
    //                 m_first_detection = Base::m_t;
    //             }
    //             mio::log_debug("Node {} is culled at time {} because there are {} total infections.", n.id, Base::m_t,
    //                            total_infections);
    //             cull_node(n.id);
    //             for (size_t index = 0; index < n.property.get_regional_neighbors()[0].size(); ++index) {
    //                 const auto node_id = n.property.get_regional_neighbors()[0][index];
    //                 vaccinate(node_id);
    //             }
    //         }
    //     }
    // }

    /**
     * @brief Get the random number generator object
     * 
     * @return RandomNumberGenerator&
     */
    RandomNumberGenerator& get_rng()
    {
        return m_rng;
    }
    /**
     * @brief Set the free parameters for the simulation.
     * 
     * @param suspicion_threshold 
     * @param sensitivity 
     * @param h0 
     * @param r0 
     * @param alpha
     * @param infection_dates
     * @param foi_inner_factors
     * @param foi_outer_factors
     */
    void set_parameters(ScalarType suspicion_threshold, ScalarType sensitivity, ScalarType h0, ScalarType r0,
                        ScalarType alpha, std::vector<Timepoint> infection_dates,
                        std::vector<ScalarType> foi_inner_factors, std::vector<ScalarType> foi_outer_factors)
    {
        m_suspicion_threshold = suspicion_threshold;
        m_sensitivity         = sensitivity;
        m_h0                  = h0;
        m_r0                  = r0;
        m_alpha               = alpha;
        m_infection_dates     = infection_dates;
        m_foi_inner_factor    = foi_inner_factors;
        m_foi_outer_factor    = foi_outer_factors;
    }

private:
    RandomNumberGenerator m_rng;
    std::queue<std::pair<size_t, ScalarType>> m_culling_queue;
    std::queue<std::pair<size_t, ScalarType>> m_prev_culling_queue;
    ScalarType m_culling_capacity_per_day = 2000;
    std::queue<std::pair<size_t, ScalarType>> m_vaccination_queue;
    ScalarType m_vaccination_capacity_per_day = 500;
    ScalarType m_first_detection              = std::numeric_limits<ScalarType>::max();
    // Trade<ScalarType, Graph> m_trade;
    bool m_standstill                          = false;
    ScalarType m_suspicion_threshold           = 0.05;
    ScalarType m_sensitivity                   = 0.95;
    std::vector<ScalarType> m_duration         = {21.0, 21.0, 400.0, 21.0, 28.0};
    std::vector<ScalarType> m_downtime         = {21.0, 21.0, 21.0, 21.0, 21.0};
    std::vector<ScalarType> m_foi_inner_factor = {1.0, 1.0, 1.0, 1.0, 1.0};
    std::vector<ScalarType> m_foi_outer_factor = {1.0, 1.0, 1.0, 1.0, 1.0};
    ScalarType m_h0                            = 0.0002;
    ScalarType m_r0                            = 4000;
    ScalarType m_alpha                         = 10;
    std::vector<Timepoint> m_infection_dates   = {0, 2, 2};
};

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
FarmSimulation<Graph<FarmNode<FP, Sim>, MobilityEdgeDirected<FP>>>
make_farm_sim(FP t0, FP dt, const Graph<FarmNode<FP, Sim>, MobilityEdgeDirected<FP>>& graph)
{
    using GraphSim = FarmSimulation<Graph<FarmNode<FP, Sim>, MobilityEdgeDirected<FP>>, FP, FP,
                                    void (*)(FP, FP, mio::MobilityEdgeDirected<FP>&, mio::FarmNode<FP, Sim>&,
                                             mio::FarmNode<FP, Sim>&, mio::RandomNumberGenerator&),
                                    void (*)(FP, FP, mio::FarmNode<FP, Sim>&)>;
    return GraphSim(t0, dt, graph, &advance_model<FP, Sim>, &apply_timed_mobility<FP, Sim>);
}

template <typename FP, class Sim>
FarmSimulation<Graph<FarmNode<FP, Sim>, MobilityEdgeDirected<FP>>>
make_farm_sim(FP t0, FP dt, Graph<FarmNode<FP, Sim>, MobilityEdgeDirected<FP>>&& graph)
{
    using GraphSim = FarmSimulation<Graph<FarmNode<FP, Sim>, MobilityEdgeDirected<FP>>, FP, FP,
                                    void (*)(FP, FP, mio::MobilityEdgeDirected<FP>&, mio::FarmNode<FP, Sim>&,
                                             mio::FarmNode<FP, Sim>&, mio::RandomNumberGenerator&),
                                    void (*)(FP, FP, mio::FarmNode<FP, Sim>&)>;
    return GraphSim(t0, dt, std::move(graph), &advance_model<FP, Sim>, &apply_timed_mobility<FP, Sim>);
}

} // namespace mio

#endif // MIO_FARM_GRAPH_SIMULATION_H