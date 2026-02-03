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
#include "memilio/mobility/graph_simulation.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/random_number_generator.h"
#include <queue>
#include "memilio/compartments/feedback_simulation.h"
#include "memilio/geography/regions.h"
#include "models/smm/parameters.h"

namespace mio
{

template <class Graph, class Timepoint = double, class Timespan = double,
          class edge_f = void (*)(Timepoint, Timespan, typename Graph::EdgeProperty&, typename Graph::NodeProperty&,
                                  typename Graph::NodeProperty&, mio::RandomNumberGenerator&),
          class node_f = void (*)(Timepoint, Timespan, typename Graph::NodeProperty&)>
class FarmSimulation : public GraphSimulationBase<Graph, Timepoint, Timespan, edge_f, node_f>
{
    using Base = GraphSimulationBase<Graph, Timepoint, Timespan, edge_f, node_f>;
    using Base::GraphSimulationBase;
    using Node = std::decay_t<decltype(Base::m_graph.nodes()[0].property)>;

public:
    void advance(Timepoint t_max = 1.0)
    {
        mio::timing::AutoTimer<"Graph Simulation Advance"> timer;
        auto dt = std::min(m_trade.next_event_time() - Base::m_t, 1.0);
        while (Base::m_t < t_max) {
            dt = std::min(dt, t_max - Base::m_t);
            for (auto& n : Base::m_graph::nodes()) {
                update_force_of_infection(n);
            }
            for (auto& n : Base::m_graph::nodes()) {
                Base::m_node_func(Base::m_t, dt, n.property);
            }
            if (!m_standstill) {
                for (auto& n : Base::m_graph::nodes()) {
                    update_population(n);
                }
            }
            for (auto& n : Base::m_graph::nodes()) {
                if (sum(n.property.get_last_result() > 0)) {
                    if (!n.property.is_suspicious()) {
                        check_for_cases(n, Base::m_t);
                    }
                    else {
                        if (n.property.get_date_confirmation() < Base::m_t + dt) {
                            //...
                        }
                    }
                }
                dt = std::min(m_trade.next_event_time() - Base::m_t, 1.0);
            }
        }

        void update_force_of_infection(Node & node)
        {
            mio::unused(node);
            //...
        }

        void update_population(Node & node)
        {
            mio::unused(node);
            //...
        }

        void check_for_cases(Node & node)
        {
            mio::unused(node);
            //...
        }

        /**
     * @brief Go through the culling queue and cull as many animals as possible within dt.
     * 
     * @param dt Time span for culling.
     */
        void cull(ScalarType dt)
        {
            auto capacity = m_culling_capacity_per_day * dt;
            while (!m_culling_queue.empty() && capacity > 0) {
                auto [node_id, day] = m_culling_queue.front();
                auto animals        = Base::m_graph.nodes()[node_id].property.get_result().get_last_value();
                auto num_animals    = std::accumulate(animals.begin(), animals.end() - 1, 0);
                if (num_animals <= capacity) {
                    mio::log_debug("Culling {} animals at node {} starting on day {} and going on for {} days.",
                                   num_animals, node_id, Base::m_t, dt);
                    for (auto index = 0; index < (animals.size() - 1); index++) {
                        animals[index] = 0;
                    }
                    animals[animals.size() - 1] += num_animals;
                    m_culling_queue.pop();
                    capacity -= num_animals;
                }
                else {
                    while (capacity > 0) {
                        mio::log_debug("Culling {} animals of node {} on day {} for {} days.", capacity, node_id,
                                       Base::m_t, dt);
                        for (auto index = 0; index < (animals.size() - 1); index++) {
                            if (animals[index] < capacity) {
                                capacity -= animals[index];
                                animals[animals.size() - 1] += animals[index];
                                animals[index] = 0;
                            }
                            else {
                                animals[animals.size() - 1] += capacity;
                                animals[index] -= capacity;
                                capacity = 0;
                                break;
                            }
                        }
                    }
                }
            }
        }

        /**
     * get the mobility parameters.
     */
        const Trade<ScalarType, Graph>& get_trade() const
        {
            return m_trade;
        }

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

        void apply_interventions()
        {
            // Set quarantine
            // for (auto& n : Base::m_graph.nodes()) {
            //     if (n.property.get_result().get_last_value()[2] > 40) {
            //         mio::log_debug("Node {} is quarantined at time {}.", n.id, Base::m_t);
            //         n.property.set_quarantined(true);
            //     }
            //     else {
            //         mio::log_debug("Node {} is not quarantined at time {}.", n.id, Base::m_t);
            //         n.property.set_quarantined(false);
            //     }
            // }
            auto infectious_compartments = {2};
            for (auto& n : Base::m_graph.nodes()) {
                auto total_infections = 0;
                for (auto comp : infectious_compartments) {
                    total_infections += n.property.get_result().get_last_value()[comp];
                }
                if (total_infections > 1) {
                    for (size_t index = 0; index < n.property.get_regional_neighbors()[0].size(); ++index) {
                        const auto node_id = n.property.get_regional_neighbors()[0][index];
                        using Model        = std::decay_t<
                                   decltype(Base::m_graph.nodes()[node_id].property.get_simulation().get_model())>;
                        using AdoptionRate =
                            mio::smm::AdoptionRates<ScalarType, typename Model::Status, mio::regions::Region>;
                        Base::m_graph.nodes()[node_id]
                            .property.get_simulation()
                            .get_model()
                            .parameters.template get<AdoptionRate>()[1]
                            .factor += infectionrisk;
                    }
                }
                if (total_infections > 4) {
                    if (m_first_detection > Base::m_t) {
                        m_first_detection = Base::m_t;
                    }
                    mio::log_debug("Node {} is culled at time {} because there are {} total infections.", n.id,
                                   Base::m_t, total_infections);
                    cull_node(n.id);
                    for (size_t index = 0; index < n.property.get_regional_neighbors()[0].size(); ++index) {
                        const auto node_id = n.property.get_regional_neighbors()[0][index];
                        vaccinate(node_id);
                    }
                }
            }
        }

        void cull_node(size_t node_id)
        {
            m_culling_queue.push(std::make_pair(node_id, Base::m_t));
            Base::m_graph.nodes()[node_id].property.set_quarantined(true);
        }

        void vaccinate_node(size_t node_id)
        {
            m_vaccination_queue.push(std::make_pair(node_id, Base::m_t));
        }

        RandomNumberGenerator& get_rng()
        {
            return m_rng;
        }

        ScalarType infectionrisk = 0.0;

    private:
        RandomNumberGenerator m_rng;
        std::queue<std::pair<size_t, ScalarType>> m_culling_queue;
        ScalarType m_culling_capacity_per_day = 2000;
        std::queue<std::pair<size_t, ScalarType>> m_vaccination_queue;
        ScalarType m_vaccination_capacity_per_day = 500;
        ScalarType m_first_detection              = std::numeric_limits<ScalarType>::max();
        Trade<ScalarType, Graph> m_trade;
        bool m_standstill = false;
    };

} // namespace mio

#endif // MIO_FARM_GRAPH_SIMULATION_H