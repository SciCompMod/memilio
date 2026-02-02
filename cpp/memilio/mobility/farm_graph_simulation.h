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

class MobilityParametersTimed
{

public:
    /**
     * @brief Construct a new Mobility Parameters Timed object
     * 
     */
    MobilityParametersTimed()
        : _exchanges{} {};
    MobilityParametersTimed(std::filebuf& input_data)
        : _exchanges{}
    {
        insert_input_data(input_data);
    };
    MobilityParametersTimed(double time, double number, size_t from, size_t to)
        : _exchanges{}
    {
        _exchanges.push(ExchangeData(time, number, from, to));
    };

    void add_exchange(double time, double number, size_t from, size_t to)
    {
        _exchanges.push(ExchangeData{time, number, from, to});
    }

    /**
    * @brief Return the number of exchanged items in the next exchange event.
    * 
    * @return auto 
    */
    auto next_event_number()
    {
        return _exchanges.top().number;
    };
    /**
     * @brief Return the timepoint of the next exchangeclass EdgePropertyT event.
     * 
     * @return auto 
     */
    auto next_event_time() const
    {
        if (_exchanges.empty()) {
            return std::numeric_limits<double>::max();
        }
        return _exchanges.top().time;
    };
    /**
     * @brief Return the destination node id of the next exchange
     * 
     * @return auto 
     */
    auto next_event_edge_id()
    {
        return _exchanges.top().edge_id;
    }
    auto next_event_to_id()
    {
        return _exchanges.top().to;
    }
    auto next_event_from_id()
    {
        return _exchanges.top().from;
    }
    /**
     * @brief Return a const reference to the next event
     * 
     * @return auto 
     */
    auto next_event()
    {
        return _exchanges.top();
    }
    /**
     * @brief Delete the next event from the heap
     * 
     * @return auto 
     */
    auto pop_next_event()
    {
        return _exchanges.pop();
    }
    /**
     * @brief Return the ExchangeData for the next exchange event and delete it from the list.
     * 
     * @return auto 
     */
    auto process_next_event()
    {
        auto next_event = _exchanges.top();
        _exchanges.pop();
        return next_event;
    };

    auto size() const
    {
        return _exchanges.size();
    }

private:
    void insert_input_data(std::filebuf& input_data)
    {
        //...
        mio::unused(input_data);
    };

    /**
     * @brief Stores Timepoint and number of exchanged items for an exchange process.
     * 
     * @param time Timepoint of the exchange process
     * @param number Number of exchanged items
     */
    struct ExchangeData {
        double time;
        double number;
        size_t from;
        size_t to;
        int edge_id{};
    };

    struct CompareExchangeData {
        bool operator()(const ExchangeData& left, const ExchangeData& right)
        {
            return left.time > right.time;
        };
    };

private:
    std::priority_queue<ExchangeData, std::vector<ExchangeData>, CompareExchangeData> _exchanges;
};

template <class Graph, class Timepoint = double, class Timespan = double,
          class edge_f = void (*)(Timepoint, Timespan, typename Graph::EdgeProperty&, typename Graph::NodeProperty&,
                                  typename Graph::NodeProperty&, mio::RandomNumberGenerator&),
          class node_f = void (*)(Timepoint, Timespan, typename Graph::NodeProperty&)>
class FarmSimulation : public GraphSimulationBase<Graph, Timepoint, Timespan, edge_f, node_f>
{
    using Base = GraphSimulationBase<Graph, Timepoint, Timespan, edge_f, node_f>;
    using Base::GraphSimulationBase;

public:
    void advance(Timepoint t_max = 1.0)
    {
        mio::timing::AutoTimer<"Graph Simulation Advance"> timer;
        auto dt = std::min(m_parameters.next_event_time() - Base::m_t, 1.0);
        while (Base::m_t < t_max) {
            mio::log_debug("Time: {}", Base::m_t);
            if (Base::m_t + dt > t_max) {
                dt = t_max - Base::m_t;
            }

            // cull(dt);
            // vaccinate(dt);

            for (auto& n : Base::m_graph.nodes()) {
                Base::m_node_func(Base::m_t, dt, n.property);
            }

            Base::m_t += dt;

            while (m_parameters.next_event_time() == Base::m_t) {
                auto next_event = m_parameters.process_next_event();
                if (Base::m_graph.nodes()[next_event.from].property.is_quarantined() ||
                    Base::m_graph.nodes()[next_event.to].property.is_quarantined()) {
                    mio::log_debug("Mobility from node {} to node {} at time {} skipped due to quarantine.",
                                   next_event.from, next_event.to, Base::m_t);
                    continue;
                }
                auto& e = Base::m_graph.get_edge(next_event.from, next_event.to);
                mio::log_debug("{}, {}", e.start_node_idx, e.end_node_idx);
                Base::m_edge_func(Base::m_t, next_event.number, e.property,
                                  Base::m_graph.nodes()[e.start_node_idx].property,
                                  Base::m_graph.nodes()[e.end_node_idx].property, m_rng);
            }
            dt = std::min(m_parameters.next_event_time() - Base::m_t, 1.0);
            apply_interventions();
        }
    }

    /**
     * @brief Go through the culling queue and cull as many animals as possible within dt.
     * 
     * @param dt Time span for culling.
     */
    void cull(ScalarType dt)
    {
        auto capacity = culling_capacity_per_day * dt;
        while (!culling_queue.empty() && capacity > 0) {
            auto [node_id, day] = culling_queue.front();
            auto animals        = Base::m_graph.nodes()[node_id].property.get_result().get_last_value();
            auto num_animals    = std::accumulate(animals.begin(), animals.end() - 1, 0);
            if (num_animals <= capacity) {
                mio::log_debug("Culling {} animals at node {} starting on day {} and going on for {} days.",
                               num_animals, node_id, Base::m_t, dt);
                for (auto index = 0; index < (animals.size() - 1); index++) {
                    animals[index] = 0;
                }
                animals[animals.size() - 1] += num_animals;
                culling_queue.pop();
                capacity -= num_animals;
            }
            else {
                while (capacity > 0) {
                    mio::log_debug("Culling {} animals of node {} on day {} for {} days.", capacity, node_id, Base::m_t,
                                   dt);
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
     * @brief Go through the vaccination queue and vaccinate as many animals as possible within dt.
     * 
     * @param dt Time span for vaccination.
     */
    void vaccinate(ScalarType dt)
    {
        auto capacity = vaccination_capacity_per_day * dt;
        while (!vaccination_queue.empty() && capacity > 0) {
            auto [node_id, day] = vaccination_queue.front();
            auto animals        = Base::m_graph.nodes()[node_id].property.get_result().get_last_value();
            auto num_animals    = std::accumulate(animals.begin(), animals.end() - 1, 0);
            if (num_animals <= capacity) {
                mio::log_debug("Vaccinating {} animals at node {} starting on day {} and going on for {} days.",
                               num_animals, node_id, Base::m_t, dt);
                for (auto index = 0; index < (animals.size() - 2); index++) {
                    animals[index] = 0;
                }
                animals[animals.size() - 2] += num_animals;
                vaccination_queue.pop();
                capacity -= num_animals;
            }
            else {
                while (capacity > 0) {
                    for (auto index = 0; index < (animals.size() - 1); index++) {
                        if (animals[index] < capacity) {
                            capacity -= animals[index];
                            animals[animals.size() - 2] += animals[index];
                            animals[index] = 0;
                        }
                        else {
                            animals[animals.size() - 2] += capacity;
                            animals[index] -= capacity;
                            capacity = 0;
                            break;
                        }
                    }
                }
            }
        }
    }

    void add_exchange(double time, double number, size_t from, size_t to)
    {
        m_parameters.add_exchange(time, number, from, to);
    }

    /**
         * get the mobility parameters.
         */
    const MobilityParametersTimed& get_parameters() const
    {
        return m_parameters;
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
        auto infectious_compartments = {1, 2, 3, 4};
        for (auto& n : Base::m_graph.nodes()) {
            auto total_infections = 0;
            for (auto comp : infectious_compartments) {
                total_infections += n.property.get_result().get_last_value()[comp];
            }
            if (total_infections > 1) {
                mio::log_debug("Node {} spreads higher risk category at time {} because there are {} total infections.",
                               n.id, Base::m_t, total_infections);
                mio::log_debug("This will affect {} other farms.", n.property.get_regional_neighbors()[0].size());
                for (size_t index = 0; index < n.property.get_regional_neighbors()[0].size(); ++index) {
                    const auto node_id = n.property.get_regional_neighbors()[0][index];
                    mio::log_debug("Neighbor node id: {}", node_id);
                    auto neighbour_property = Base::m_graph.nodes()[node_id].property;
                    using Model             = std::decay_t<decltype(neighbour_property.get_simulation().get_model())>;
                    using AdoptionRate =
                        mio::smm::AdoptionRates<ScalarType, typename Model::Status, mio::regions::Region>;
                    neighbour_property.get_simulation().get_model().parameters.template get<AdoptionRate>()[1].factor +=
                        infectionrisk;
                }
            }
            if (total_infections > 4) {
                if (first_detection > Base::m_t) {
                    first_detection = Base::m_t;
                }
                mio::log_debug("Node {} is culled at time {} because there are {} total infections.", n.id, Base::m_t,
                               total_infections);
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
        culling_queue.push(std::make_pair(node_id, Base::m_t));
        Base::m_graph.nodes()[node_id].property.set_quarantined(true);
        mio::log_debug("It is day {} and we already want to cull {} farms.", Base::m_t, culling_queue.size());
    }

    void vaccinate_node(size_t node_id)
    {
        vaccination_queue.push(std::make_pair(node_id, Base::m_t));
    }

    RandomNumberGenerator& get_rng()
    {
        return m_rng;
    }

    ScalarType infectionrisk = 0.0;

private:
    mio::MobilityParametersTimed m_parameters;
    RandomNumberGenerator m_rng;
    std::queue<std::pair<size_t, ScalarType>> culling_queue;
    ScalarType culling_capacity_per_day = 2000;
    std::queue<std::pair<size_t, ScalarType>> vaccination_queue;
    ScalarType vaccination_capacity_per_day = 500;
    ScalarType first_detection              = std::numeric_limits<ScalarType>::max();
};

} // namespace mio

#endif // MIO_FARM_GRAPH_SIMULATION_H