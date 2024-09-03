/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Henrik Zunker
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
#ifndef METAPOPULATION_MOBILITY_DETAILED_H
#define METAPOPULATION_MOBILITY_DETAILED_H

#include "memilio/compartments/parameter_studies.h"
#include "memilio/epidemiology/simulation_day.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/logging.h"
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
#include "memilio/mobility/graph.h"
#include "memilio/io/mobility_io.h"
#include "memilio/math/floating_point.h"

#include "boost/filesystem.hpp"

#include <cassert>
#include <cstddef>
#include <string>

namespace mio
{

template <typename Sim>
struct ExtendedNodeProperty {
    Sim base_sim;
    Sim mobility_sim;
    double stay_duration;

    ExtendedNodeProperty(Sim sim1, Sim sim2, double t)
        : base_sim(sim1)
        , mobility_sim(sim2)
        , stay_duration(t)
    {
    }
};

template <typename FP = double>
class ExtendedMobilityEdge : public MobilityEdge<FP>
{
public:
    double travel_time;
    std::vector<int> path;

    ExtendedMobilityEdge(const MobilityParameters<FP>& params, double tt, std::vector<int> p)
        : MobilityEdge<FP>(params)
        , travel_time(tt)
        , path(p)
    {
    }

    ExtendedMobilityEdge(const Eigen::VectorXd& coeffs, double tt, std::vector<int> p)
        : MobilityEdge<FP>(coeffs)
        , travel_time(tt)
        , path(p)
    {
    }

    auto& get_migrated()
    {
        return this->m_mobile_population;
    }
    auto& get_return_times()
    {
        return this->m_return_times;
    }
    auto& get_parameters() const
    {
        return this->m_parameters;
    }
};

template <class Sim>
using ExtendedGraph = Graph<ExtendedNodeProperty<Sim>, ExtendedMobilityEdge<double>>;

// Default implementation when get_mobility_factors is not defined for Sim
template <class Sim, std::enable_if_t<!is_expression_valid<get_mobility_factors_expr_t, Sim>::value, void*> = nullptr>
auto get_mobility_factors(const Sim& /*sim*/, double /*t*/, const Eigen::Ref<const Eigen::VectorXd>& y)
{
    return Eigen::VectorXd::Ones(y.rows());
}

template <typename FP>
void check_negative_values_vec(Eigen::Ref<Vector<FP>> vec, const size_t num_age_groups, FP tolerance = -1e-7,
                               const size_t max_iterations = 100)
{
    // before moving the commuters, we need to look for negative values in vec and correct them.
    const size_t num_comparts = vec.size() / num_age_groups;
    size_t iteration_count    = 0;

    // check for negative values in vec
    while (vec.minCoeff() < tolerance) {
        Eigen::Index min_index;
        const FP min_value = vec.minCoeff(&min_index);

        const auto curr_age_group = min_index / num_comparts;
        const auto indx_begin     = curr_age_group * num_comparts;
        const auto indx_end       = indx_begin + num_comparts;

        auto max_group_indx = indx_begin;
        for (auto i = indx_begin; i < indx_end; ++i) {
            if (vec(i) > vec(max_group_indx)) {
                max_group_indx = i;
            }
        }

        // Correct the negative value
        vec(min_index) = 0;
        vec(max_group_indx) += min_value;

        // check if the number of iterations is exceeded
        if (iteration_count > max_iterations) {
            std::vector<FP> vec_std(vec.data(), vec.data() + vec.size());
            mio::unused(vec_std);
            log_error("Number of iterations exceeded in check_negative_values_vec.");
            std::exit(1);
        }
        iteration_count++;
    }
}

template <typename FP, class Sim>
Eigen::Index find_time_index(Sim& simulation, FP t, bool create_new_tp)
{
    auto& results = simulation.get_result();
    auto& flows   = simulation.get_flows();
    if (results.get_num_time_points() != flows.get_num_time_points()) {
        log_error("Number of time points in results " + std::to_string(results.get_num_time_points()) + " and flows " +
                  std::to_string(flows.get_num_time_points()) + " do not match in find_time_index.");
    }
    Eigen::Index t_indx = results.get_num_time_points() - 1;
    for (; t_indx >= 0; --t_indx) {
        if (std::abs(results.get_time(t_indx) - t) < 1e-10) {
            break;
        }
    }

    if (t_indx < 0 && !create_new_tp) {
        log_error("Time point " + std::to_string(t) + " not found in find_time_index. Lates time point is " +
                  std::to_string(results.get_last_time()));
    }

    if (t_indx < 0 && results.get_last_time() < t && create_new_tp) {
        // if we allow to create a new time point, we initialize the compartments with zero
        // the flows are accumulated. Therefore, we can just copy the last value.
        Eigen::VectorXd results_vec = Eigen::VectorXd::Zero(results.get_last_value().size());
        results.add_time_point(t, results_vec);
        flows.add_time_point(t, flows.get_last_value());
        t_indx = results.get_num_time_points() - 1;
    }

    return t_indx;
}

template <typename FP, class Sim>
class MobilityFunctions
{
public:
    void init_mobility(FP t, ExtendedMobilityEdge<FP>& edge, Sim& from_sim, Sim& to_sim)
    {
        const auto t_indx_start_mobility_sim_from = find_time_index(from_sim, t, false);

        // initialize the number of commuters at the start of the mobility
        auto results_from = from_sim.get_result();
        edge.get_migrated().add_time_point(
            t, (results_from.get_value(t_indx_start_mobility_sim_from).array() *
                edge.get_parameters().get_coefficients().get_matrix_at(t).array() *
                get_mobility_factors(from_sim, t, results_from.get_value(t_indx_start_mobility_sim_from)).array())
                   .matrix());
        edge.get_return_times().add_time_point(t);

        // move them to the starting mobility model
        // if the simulation we are adding the commuters to is not having the same time point as the current time point,
        // we need to add a new time point to the simulation.
        const auto t_indx_start_mobility_sim_to = find_time_index(to_sim, t, true);
        to_sim.get_result().get_value(t_indx_start_mobility_sim_to) += edge.get_migrated().get_last_value();
        from_sim.get_result().get_last_value() -= edge.get_migrated().get_last_value();
    }

    void move_migrated(FP t, ExtendedMobilityEdge<FP>& edge, Sim& from_sim, Sim& to_sim)
    {
        // When moving from one regional entity/model to another, we need to update the local population.
        // check_negative_values_vec needs to be called once since its checks for negative values and corrects them.
        const size_t num_age_groups = static_cast<size_t>(from_sim.get_model().parameters.get_num_groups());
        check_negative_values_vec(edge.get_migrated().get_last_value(), num_age_groups);
        const auto t_indx_sim_to_arrival = find_time_index(to_sim, t, true);
        from_sim.get_result().get_last_value() -= edge.get_migrated().get_last_value();
        to_sim.get_result().get_value(t_indx_sim_to_arrival) += edge.get_migrated().get_last_value();

        // check each result for negative values and correct them if necessary
        check_negative_values_vec(from_sim.get_result().get_last_value(), num_age_groups);
        check_negative_values_vec(to_sim.get_result().get_value(t_indx_sim_to_arrival), num_age_groups);
    }

    void update_commuters(FP t, FP dt, ExtendedMobilityEdge<FP>& edge, Sim& sim, bool is_mobility_model)
    {
        const auto t_indx_start_mobility_sim = find_time_index(sim, t, true);
        Eigen::VectorXd flows                = Eigen::VectorXd::Zero(sim.get_flows().get_last_value().size());
        update_status_migrated(edge.get_migrated().get_last_value(), sim,
                               sim.get_result().get_value(t_indx_start_mobility_sim), t, dt, flows);

        // if the simulation is holding a mobility model, we need to update the mobility model as well
        // Therefore, we check if the time point already exists in the mobility model and create a new one if necessary.
        // Next, we build the population in the mobility model based on the commuters
        if (is_mobility_model) {
            const auto t_indx_mobility_model = find_time_index(sim, t + dt, true);
            if (t_indx_mobility_model != sim.get_result().get_num_time_points() - 1) {
                log_error("Time point " + std::to_string(t + dt) +
                          " not the lastest in update_commuters. Latest time point is " +
                          std::to_string(sim.get_result().get_last_time()));
            }
            sim.get_result().get_value(t_indx_mobility_model) += edge.get_migrated().get_last_value();
            sim.get_flows().get_value(t_indx_mobility_model) += flows;
        }
    }

    void delete_migrated(ExtendedMobilityEdge<FP>& edge)
    {
        for (Eigen::Index i = edge.get_return_times().get_num_time_points() - 1; i >= 0; --i) {
            edge.get_migrated().remove_time_point(i);
            edge.get_return_times().remove_time_point(i);
        }
    }
};

class ScheduleManager
{
public:
    struct Schedule {
        std::vector<std::vector<size_t>> schedule_edges;
        std::vector<std::vector<bool>> mobility_schedule_edges;
        std::vector<std::vector<size_t>> mobility_int_schedule;
        std::vector<std::vector<size_t>> local_int_schedule;
        std::vector<std::vector<size_t>> edges_mobility;
        std::vector<std::vector<size_t>> nodes_mobility;
        std::vector<std::vector<size_t>> nodes_mobility_m;
        std::vector<size_t> first_mobility;
    };

    ScheduleManager(size_t t, double eps = 1e-10)
        : timesteps(t)
        , epsilon(eps)
    {
    }

    template <typename Graph>
    Schedule compute_schedule(const Graph& graph)
    {
        Schedule schedule;
        calculate_edge_schedule(graph, schedule);
        calculate_node_schedule(graph, schedule);
        return schedule;
    }

    /**
     * @brief Calculates the edge schedule for a given graph.
     * 
     * This function computes the vectors schedule_edges and mobility_schedule_edges for each edge in the graph.
     * The schedule_eges describes for each time step in which node_id the individuals are located.
     * The mobility_schedule_edges describes for each time step if the individuals are in a mobility node or not.
     * 
     * @tparam Graph The type of the graph.
     * @param graph The graph object containing the nodes and edges.
     * @param schedule The Schedule object to store the computed schedules.
     */
    template <typename Graph>
    void calculate_edge_schedule(const Graph& graph, Schedule& schedule)
    {
        schedule.schedule_edges.reserve(graph.edges().size());
        schedule.mobility_schedule_edges.reserve(graph.edges().size());

        // calculate the schedule for each edge
        for (auto& e : graph.edges()) {
            std::vector<size_t> tmp_schedule(timesteps, 0);
            std::vector<bool> is_mobility_node(timesteps, false);

            // Calculate travel time per region, ensuring a minimum value of 0.01
            const double traveltime_per_region =
                std::max(0.01, round_nth_decimal(e.property.travel_time / e.property.path.size(), 2));

            // Calculate the start time for mobility, ensuring it greater or equal to 0.01
            const double start_mobility =
                std::max(0.01, round_nth_decimal(1 - 2 * traveltime_per_region * e.property.path.size() -
                                                     graph.nodes()[e.end_node_idx].property.stay_duration,
                                                 2));

            // Calculate the arrival time at the destination node
            const double arrive_at = start_mobility + traveltime_per_region * e.property.path.size();

            // Lambda to fill the schedule vector with the node index during the trip to the destination.
            auto fill_schedule = [&](size_t start_idx, size_t end_idx, size_t value) {
                std::fill(tmp_schedule.begin() + start_idx, tmp_schedule.begin() + end_idx, value);
            };

            // Lambda to fill the schedule for the mobility models with a bool during the trip to the destination.
            auto fill_mobility = [&](size_t start_idx, size_t end_idx, bool value) {
                std::fill(is_mobility_node.begin() + start_idx, is_mobility_node.begin() + end_idx, value);
            };

            // Indices for schedule filling
            const size_t start_idx    = static_cast<size_t>((start_mobility + epsilon) * 100);
            const size_t arrive_idx   = static_cast<size_t>((arrive_at + epsilon) * 100);
            const size_t stay_end_idx = timesteps - (arrive_idx - start_idx);

            // Fill the schedule up to the start of mobility with the start node index
            fill_schedule(0, start_idx, e.start_node_idx);

            // Mark the mobility period in the mobility node vector
            fill_mobility(start_idx, arrive_idx, true);

            // Fill the schedule for the path during the mobility period
            size_t current_index = start_idx;
            for (size_t county : e.property.path) {
                size_t next_index = current_index + static_cast<size_t>((traveltime_per_region + epsilon) * 100);
                fill_schedule(current_index, next_index, county);
                current_index = next_index;
            }

            // Fill the remaining schedule after mobility with the end node index
            fill_schedule(current_index, stay_end_idx, e.property.path.back());

            // Mark the return mobility period in the mobility node vector
            fill_mobility(stay_end_idx, timesteps, true);

            // Find the first and last true values in the mobility node vector
            auto first_true = std::find(is_mobility_node.begin(), is_mobility_node.end(), true);
            auto last_true  = std::find(is_mobility_node.rbegin(), is_mobility_node.rend(), true);

            // Ensure there is at least one true value
            if (first_true != is_mobility_node.end() && last_true != is_mobility_node.rend()) {
                size_t first_index = std::distance(is_mobility_node.begin(), first_true);
                size_t count_true  = arrive_idx - start_idx;

                // Create a reversed path segment for the return trip
                std::vector<size_t> path_reversed(tmp_schedule.begin() + first_index,
                                                  tmp_schedule.begin() + first_index + count_true);
                std::reverse(path_reversed.begin(), path_reversed.end());

                // Copy the reversed path segment to the end of the schedule for the return trip
                std::copy(path_reversed.begin(), path_reversed.end(), tmp_schedule.end() - count_true);

                // Add the schedule and mobility node vectors to their respective containers
                schedule.schedule_edges.push_back(std::move(tmp_schedule));
                schedule.mobility_schedule_edges.push_back(std::move(is_mobility_node));
            }
            else {
                log_error("Error in creating schedule.");
            }
        }
    }

    /**
     * @brief Calculates the node schedule for a given graph.
     * 
     * This function computes the vectors local_int_schedule, mobility_int_schedule, edges_mobility, nodes_mobility, 
     * nodes_mobility_m and first_mobility.
     * The local_int_schedule vector describes the time steps where we need to synchronize the integrator in the local models.
     * The mobility_int_schedule vector describes the time steps where we need to synchronize the integrator in the mobility models.
     * The edges_mobility vector describes all edges where mobility takes place for each time step.
     * The nodes_mobility vector describes all local models where mobility takes place for each time step.
     * The nodes_mobility_m vector describes all mobility models where mobility takes place at the beginning of each time step.
     * The first_mobility vector describes the first time step where mobility takes place for each edge.
     * 
     * @tparam Graph The type of the graph object.
     * @param graph The graph object containing the nodes and edges.
     * @param schedule The Schedule object to store the computed schedules.
     */
    template <typename Graph>
    void calculate_node_schedule(const Graph& graph, Schedule& schedule)
    {

        // iterate over nodes to create the integration schedule for each node. The integration schedule is necessary
        // to determine the correct time step for the integrator in the nodes.
        for (size_t indx_node = 0; indx_node < graph.nodes().size(); ++indx_node) {
            // Local node initialization with starting at t=0
            std::vector<size_t> order_local_node = {0};
            std::vector<size_t> indx_edges;

            // Find edges starting from the current node and not in mobility at t=0
            auto find_edges = [&](size_t t, bool mobility) {
                std::vector<size_t> edges;
                for (size_t indx_edge = 0; indx_edge < schedule.schedule_edges.size(); ++indx_edge) {
                    if (schedule.schedule_edges[indx_edge][t] == indx_node &&
                        schedule.mobility_schedule_edges[indx_edge][t] == mobility) {
                        edges.push_back(indx_edge);
                    }
                }
                return edges;
            };

            indx_edges = find_edges(0, false);

            // Iterate through each timestep to identify changes in local node schedule
            for (size_t indx_t = 1; indx_t < timesteps; ++indx_t) {
                auto indx_edges_new = find_edges(indx_t, false);

                if (indx_edges_new.size() != indx_edges.size() ||
                    !std::equal(indx_edges.begin(), indx_edges.end(), indx_edges_new.begin())) {
                    order_local_node.push_back(indx_t);
                    indx_edges = indx_edges_new;
                }
            }

            // Ensure the last timestep is included
            if (order_local_node.back() != timesteps - 1) {
                order_local_node.push_back(timesteps - 1);
            }
            schedule.local_int_schedule.push_back(order_local_node);

            // Mobility node initialization
            std::vector<size_t> order_mobility_node;
            std::vector<size_t> indx_edges_mobility;

            indx_edges_mobility = find_edges(0, true);

            // Iterate through each timestep to identify changes in mobility node schedule
            for (size_t indx_t = 1; indx_t < timesteps; ++indx_t) {
                auto indx_edges_mobility_new = find_edges(indx_t, true);

                if (indx_edges_mobility_new.size() != indx_edges_mobility.size() ||
                    !std::equal(indx_edges_mobility.begin(), indx_edges_mobility.end(),
                                indx_edges_mobility_new.begin())) {
                    order_mobility_node.push_back(indx_t);
                    indx_edges_mobility = indx_edges_mobility_new;
                }
            }

            schedule.mobility_int_schedule.push_back(order_mobility_node);
        }

        // Reserve space for first_mobility vector and initialize it
        schedule.first_mobility.reserve(graph.edges().size());
        for (size_t indx_edge = 0; indx_edge < graph.edges().size(); ++indx_edge) {
            // find the first time step where mobility takes place. If there is no mobility, it is set the size of the schedule.
            size_t index_time = 0;
            for (; index_time < schedule.mobility_schedule_edges[indx_edge].size(); ++index_time) {
                if (schedule.mobility_schedule_edges[indx_edge][index_time]) {
                    break;
                }
            }
            schedule.first_mobility.push_back(index_time);
        }

        // Reserve space for mobility-related vectors
        schedule.edges_mobility.reserve(timesteps);
        schedule.nodes_mobility.reserve(timesteps);
        schedule.nodes_mobility_m.reserve(timesteps);

        // Handle the case where indx_current = 0 separately
        std::vector<size_t> temp_edges_mobility;
        for (size_t indx_edge = 0; indx_edge < graph.edges().size(); ++indx_edge) {
            if (schedule.first_mobility[indx_edge] == 0) {
                temp_edges_mobility.push_back(indx_edge);
            }
        }
        schedule.edges_mobility.push_back(std::move(temp_edges_mobility));

        // Initialize nodes_mobility with all node indices for t=0
        std::vector<size_t> temp_nodes_mobility(graph.nodes().size());
        std::iota(temp_nodes_mobility.begin(), temp_nodes_mobility.end(), 0);
        schedule.nodes_mobility.emplace_back(std::move(temp_nodes_mobility));

        // Initialize nodes_mobility_m with nodes that have mobility activities at t=0
        std::vector<size_t> temp_nodes_mobility_m;
        for (size_t node_indx = 0; node_indx < graph.nodes().size(); ++node_indx) {
            if (std::binary_search(schedule.mobility_int_schedule[node_indx].begin(),
                                   schedule.mobility_int_schedule[node_indx].end(), 0)) {
                temp_nodes_mobility_m.push_back(node_indx);
            }
        }
        schedule.nodes_mobility_m.push_back(temp_nodes_mobility_m);

        for (size_t indx_current = 1; indx_current < timesteps; ++indx_current) {
            std::vector<size_t> temp_edge_mobility;
            for (size_t indx_edge = 0; indx_edge < graph.edges().size(); ++indx_edge) {
                size_t current_node_indx = schedule.schedule_edges[indx_edge][indx_current];
                if (indx_current >= schedule.first_mobility[indx_edge]) {
                    if (schedule.mobility_schedule_edges[indx_edge][indx_current] &&
                        std::binary_search(schedule.mobility_int_schedule[current_node_indx].begin(),
                                           schedule.mobility_int_schedule[current_node_indx].end(), indx_current)) {
                        temp_edge_mobility.push_back(indx_edge);
                    }
                    else if (!schedule.mobility_schedule_edges[indx_edge][indx_current] &&
                             std::binary_search(schedule.local_int_schedule[current_node_indx].begin(),
                                                schedule.local_int_schedule[current_node_indx].end(), indx_current)) {
                        temp_edge_mobility.push_back(indx_edge);
                    }
                }
            }
            schedule.edges_mobility.push_back(temp_edge_mobility);

            // Clear and fill temp_nodes_mobility and temp_nodes_mobility_m for current timestep
            temp_nodes_mobility.clear();
            temp_nodes_mobility_m.clear();
            for (size_t node_indx = 0; node_indx < graph.nodes().size(); ++node_indx) {
                if (std::binary_search(schedule.local_int_schedule[node_indx].begin(),
                                       schedule.local_int_schedule[node_indx].end(), indx_current)) {
                    temp_nodes_mobility.push_back(node_indx);
                }

                if (std::binary_search(schedule.mobility_int_schedule[node_indx].begin(),
                                       schedule.mobility_int_schedule[node_indx].end(), indx_current)) {
                    temp_nodes_mobility_m.push_back(node_indx);
                }
            }
            schedule.nodes_mobility.push_back(temp_nodes_mobility);
            schedule.nodes_mobility_m.push_back(temp_nodes_mobility_m);
        }
    }

    size_t timesteps;
    double epsilon;
};

template <typename Graph, typename MobilityFunctions>
class GraphSimulationExtended : public GraphSimulationBase<Graph>
{
public:
    using node_function = std::function<void(double, double, typename Graph::NodeProperty&)>;
    using edge_function =
        std::function<void(double, double, typename Graph::EdgeProperty&, typename Graph::NodeProperty&,
                           typename Graph::NodeProperty&, MobilityFunctions&)>;

    GraphSimulationExtended(double t0, double dt, const Graph& g, const node_function& node_func,
                            MobilityFunctions modes)
        : GraphSimulationBase<Graph>(t0, dt, g, node_func, {})
        , m_mobility_functions(modes)
    {
        ScheduleManager schedule_manager(100); // Assuming 100 timesteps
        schedules = schedule_manager.compute_schedule(this->m_graph);
    }

    GraphSimulationExtended(double t0, double dt, Graph&& g, const node_function& node_func, MobilityFunctions modes)
        : GraphSimulationBase<Graph>(t0, dt, std::move(g), node_func, {})
        , m_mobility_functions(modes)
    {
        ScheduleManager schedule_manager(100); // Assuming 100 timesteps
        schedules = schedule_manager.compute_schedule(this->m_graph);
    }

    void advance(double t_max = 1.0)
    {
        ScalarType dt_first_mobility =
            std::accumulate(this->m_graph.edges().begin(), this->m_graph.edges().end(),
                            std::numeric_limits<ScalarType>::max(), [&](ScalarType current_min, const auto& e) {
                                auto traveltime_per_region =
                                    round_nth_decimal(e.property.travel_time / e.property.path.size(), 2);
                                if (traveltime_per_region < 0.01)
                                    traveltime_per_region = 0.01;
                                auto start_mobility =
                                    round_nth_decimal(1 - 2 * (traveltime_per_region * e.property.path.size()) -
                                                          this->m_graph.nodes()[e.end_node_idx].property.stay_duration,
                                                      2);
                                if (start_mobility < 0) {
                                    start_mobility = 0.;
                                }
                                return std::min(current_min, start_mobility);
                            });

        // set population to zero in mobility nodes before starting
        for (auto& n : this->m_graph.nodes()) {
            n.property.mobility_sim.get_result().get_last_value().setZero();
        }

        auto min_dt    = 0.01;
        double t_begin = this->m_t - 1.;

        while (this->m_t < t_max - epsilon) {
            t_begin += 1;
            if (this->m_t + dt_first_mobility > t_max) {
                dt_first_mobility = t_max - this->m_t;
                for (auto& n : this->m_graph.nodes()) {
                    n.property.base_sim.advance(this->m_t + dt_first_mobility);
                }
                break;
            }

            size_t indx_schedule = 0;
            while (t_begin + 1 > this->m_t + 1e-10) {
                // the graph simulation is structured in 3 steps:
                // 1. move indivudals if necessary
                // 2. integrate the local nodes
                // 3. update the status of the commuter and use the obtained values to update the mobility models.
                move_commuters(indx_schedule);
                advance_local_nodes(indx_schedule);
                update_status_commuters(indx_schedule);

                // in the last time step we have to move the individuals back to their local model and
                // delete the commuters time series
                handle_last_time_step(indx_schedule);

                indx_schedule++;
                this->m_t += min_dt;
            }
            // At the end of the day. we set each compartment zero for all mobility nodes since we have to estimate
            // the state of the indivuals moving between different counties.
            // Therefore there can be differences with the states given by the integrator used for the mobility node.
            for (auto& n : this->m_graph.nodes()) {
                n.property.mobility_sim.get_result().get_last_value().setZero();
            }

            // to save memory, we interpolate each time series after every day if t > 20
            if (this->m_t > 20) {
                for (auto& n : this->m_graph.nodes()) {
                    // base sim
                    auto& result_node_local          = n.property.base_sim.get_result();
                    auto interpolated_result         = interpolate_simulation_result(result_node_local);
                    auto& flow_node_local            = n.property.base_sim.get_flows();
                    auto interpolated_flows          = interpolate_simulation_result(flow_node_local);
                    n.property.base_sim.get_result() = interpolated_result;
                    n.property.base_sim.get_flows()  = interpolated_flows;

                    // mobility sim
                    auto& result_node_mobility           = n.property.mobility_sim.get_result();
                    interpolated_result                  = interpolate_simulation_result(result_node_mobility);
                    auto& flow_node_mobility             = n.property.mobility_sim.get_flows();
                    interpolated_flows                   = interpolate_simulation_result(flow_node_mobility);
                    n.property.mobility_sim.get_result() = interpolated_result;
                    n.property.mobility_sim.get_flows()  = interpolated_flows;
                }
            }
        }
    }

private:
    MobilityFunctions m_mobility_functions;
    ScheduleManager::Schedule schedules;
    const double epsilon = 1e-10;

    ScalarType calculate_next_dt(size_t edge_indx, size_t indx_schedule)
    {
        auto current_node_indx = schedules.schedule_edges[edge_indx][indx_schedule];
        bool in_mobility_node  = schedules.mobility_schedule_edges[edge_indx][indx_schedule];

        // determine dt, which is equal to the last integration/synchronization point in the current node
        auto integrator_schedule_row = schedules.local_int_schedule[current_node_indx];
        if (in_mobility_node)
            integrator_schedule_row = schedules.mobility_int_schedule[current_node_indx];
        // search the index of indx_schedule in the integrator schedule
        const size_t indx_current = std::distance(
            integrator_schedule_row.begin(),
            std::lower_bound(integrator_schedule_row.begin(), integrator_schedule_row.end(), indx_schedule));

        if (integrator_schedule_row[indx_current] != indx_schedule)
            throw std::runtime_error("Error in schedule.");

        ScalarType dt_mobility;
        if (indx_schedule == 99 || indx_current == integrator_schedule_row.size() - 1) {
            // if we are at the last iteration, we choose the minimal time step
            dt_mobility = (100 - indx_schedule) * 0.01;
        }
        else {
            // else, we calculate the next time step based on the next integration point
            dt_mobility = round_nth_decimal((static_cast<double>(integrator_schedule_row[indx_current + 1]) -
                                             static_cast<double>(integrator_schedule_row[indx_current])) /
                                                    100 +
                                                epsilon,
                                            2);
        }

        return dt_mobility;
    }

    void move_commuters(size_t indx_schedule)
    {
        for (const auto& edge_indx : schedules.edges_mobility[indx_schedule]) {
            auto& e = this->m_graph.edges()[edge_indx];
            // start mobility by initializing the number of commuters and move to initial mobility model
            if (indx_schedule == schedules.first_mobility[edge_indx]) {
                auto& node_from =
                    this->m_graph.nodes()[schedules.schedule_edges[edge_indx][indx_schedule - 1]].property.base_sim;
                auto& node_to =
                    this->m_graph.nodes()[schedules.schedule_edges[edge_indx][indx_schedule]].property.mobility_sim;
                m_mobility_functions.init_mobility(this->m_t, e.property, node_from, node_to);
            }
            else if (indx_schedule > schedules.first_mobility[edge_indx]) {
                // send the individuals to the next node
                if ((schedules.schedule_edges[edge_indx][indx_schedule] !=
                     schedules.schedule_edges[edge_indx][indx_schedule - 1]) ||
                    (schedules.mobility_schedule_edges[edge_indx][indx_schedule] !=
                     schedules.mobility_schedule_edges[edge_indx][indx_schedule - 1])) {
                    auto& node_from =
                        schedules.mobility_schedule_edges[edge_indx][indx_schedule - 1]
                            ? this->m_graph.nodes()[schedules.schedule_edges[edge_indx][indx_schedule - 1]]
                                  .property.mobility_sim
                            : this->m_graph.nodes()[schedules.schedule_edges[edge_indx][indx_schedule - 1]]
                                  .property.base_sim;

                    auto& node_to = schedules.mobility_schedule_edges[edge_indx][indx_schedule]
                                        ? this->m_graph.nodes()[schedules.schedule_edges[edge_indx][indx_schedule]]
                                              .property.mobility_sim
                                        : this->m_graph.nodes()[schedules.schedule_edges[edge_indx][indx_schedule]]
                                              .property.base_sim;

                    m_mobility_functions.move_migrated(this->m_t, e.property, node_from, node_to);
                }
            }
        }
    }

    void update_status_commuters(size_t indx_schedule, const size_t max_num_contacts = 20)
    {
        for (const auto& edge_indx : schedules.edges_mobility[indx_schedule]) {
            auto& e      = this->m_graph.edges()[edge_indx];
            auto next_dt = calculate_next_dt(edge_indx, indx_schedule);
            auto& node_to =
                schedules.mobility_schedule_edges[edge_indx][indx_schedule]
                    ? this->m_graph.nodes()[schedules.schedule_edges[edge_indx][indx_schedule]].property.mobility_sim
                    : this->m_graph.nodes()[schedules.schedule_edges[edge_indx][indx_schedule]].property.base_sim;

            // get current contact  and scale it but only if mobility model
            auto contact_pattern_curr = get_contact_pattern(node_to.get_model());
            auto contacts_copy        = contact_pattern_curr;
            if (schedules.mobility_schedule_edges[edge_indx][indx_schedule]) {
                auto& contact_matrix          = contact_pattern_curr.get_cont_freq_mat();
                Eigen::MatrixXd scaled_matrix = contact_matrix[0].get_baseline().eval() / e.property.travel_time;
                // check if there a values greater max_num_contacts in the contact matrix. if higher, set to max_num_contacts
                for (auto i = 0; i < scaled_matrix.rows(); ++i) {
                    for (auto j = 0; j < scaled_matrix.cols(); ++j) {
                        if (scaled_matrix(i, j) > max_num_contacts) {
                            scaled_matrix(i, j) = max_num_contacts;
                        }
                    }
                }

                contact_matrix[0].get_baseline() = scaled_matrix;
                set_contact_pattern(node_to.get_model(), contact_pattern_curr);
            }

            m_mobility_functions.update_commuters(this->m_t, next_dt, e.property, node_to,
                                                  schedules.mobility_schedule_edges[edge_indx][indx_schedule]);

            // reset contact pattern after estimating the state of the commuters
            if (schedules.mobility_schedule_edges[edge_indx][indx_schedule])
                set_contact_pattern(node_to.get_model(), contacts_copy);
        }
    }

    void advance_local_nodes(size_t indx_schedule)
    {
        for (const auto& n_indx : schedules.nodes_mobility[indx_schedule]) {
            auto& n = this->m_graph.nodes()[n_indx];
            const size_t indx_current =
                std::distance(schedules.local_int_schedule[n_indx].begin(),
                              std::lower_bound(schedules.local_int_schedule[n_indx].begin(),
                                               schedules.local_int_schedule[n_indx].end(), indx_schedule));

            const size_t val_next = (indx_current == schedules.local_int_schedule[n_indx].size() - 1)
                                        ? 100
                                        : schedules.local_int_schedule[n_indx][indx_current + 1];
            const ScalarType next_dt =
                round_nth_decimal((static_cast<double>(val_next) - indx_schedule) / 100 + epsilon, 2);
            n.property.base_sim.advance(this->m_t + next_dt);
            // check if last value contains negative values or nan values
            if (n.property.base_sim.get_result().get_last_value().minCoeff() < -1e-7 ||
                std::isnan(n.property.base_sim.get_result().get_last_value().sum())) {
                auto last_value_as_vec =
                    std::vector<ScalarType>(n.property.base_sim.get_result().get_last_value().data(),
                                            n.property.base_sim.get_result().get_last_value().data() +
                                                n.property.base_sim.get_result().get_last_value().size());
                log_error("Negative Value " + std::to_string(n_indx) + " at time " + std::to_string(indx_schedule) +
                          " in local node detected.");
            }
        }
    }

    void handle_last_time_step(int indx_schedule)
    {
        if (indx_schedule == 99) {
            auto edge_index = 0;
            for (auto& e : this->m_graph.edges()) {
                auto& node_from =
                    this->m_graph.nodes()[schedules.schedule_edges[edge_index][indx_schedule]].property.mobility_sim;
                auto& node_to =
                    this->m_graph.nodes()[schedules.schedule_edges[edge_index][indx_schedule]].property.base_sim;

                if (schedules.schedule_edges[edge_index][indx_schedule] != e.start_node_idx)
                    log_error("Last node is not the start node in edge " + std::to_string(edge_index) + " at time " +
                              std::to_string(indx_schedule));
                // move the individuals back to the local model
                m_mobility_functions.move_migrated(this->m_t + 0.01, e.property, node_from, node_to);
                m_mobility_functions.delete_migrated(e.property);
                edge_index++;
            }
        }
    }
};

template <typename FP, class Sim>
GraphSimulationExtended<Graph<ExtendedNodeProperty<Sim>, ExtendedMobilityEdge<FP>>, MobilityFunctions<FP, Sim>>
make_mobility_sim(FP t0, FP dt, Graph<ExtendedNodeProperty<Sim>, ExtendedMobilityEdge<FP>>&& graph)
{
    auto migration_modes = MobilityFunctions<FP, Sim>();
    return GraphSimulationExtended<Graph<ExtendedNodeProperty<Sim>, ExtendedMobilityEdge<FP>>,
                                   MobilityFunctions<FP, Sim>>(t0, dt, std::move(graph), {}, migration_modes);
}

/**
 * @brief number of migrated people when they return according to the model.
 * E.g. during the time in the other node, some people who left as susceptible will return exposed.
 * Implemented for general compartmentmodel simulations, overload for your custom model if necessary
 * so that it can be found with argument-dependent lookup, i.e. in the same namespace as the model.
 * @param migrated number of people that migrated as input, number of people that return as output
 * @param sim Simulation that is used for the migration
 * @param integrator Integrator that is used for the estimation. Has to be a one-step scheme.
 * @param total total population in the node that the people migrated to.
 * @param t time of migration
 * @param dt timestep
 */

template <typename FP, class Sim, class = std::enable_if_t<is_compartment_model_simulation<FP, Sim>::value>>
void update_status_migrated(Eigen::Ref<typename TimeSeries<FP>::Vector> migrated, Sim& sim,
                            Eigen::Ref<const typename TimeSeries<FP>::Vector> total, FP t, FP dt,
                            Eigen::VectorXd& flows)
{
    auto y0 = migrated.eval();
    auto y1 = migrated.setZero();
    mio::EulerIntegratorCore<FP>().step(
        [&](auto&& y, auto&& t_, auto&& dydt) {
            sim.get_model().get_derivatives(total, y, t_, dydt);
        },
        y0, t, dt, y1);
    flows = sim.get_model().get_flow_values();
    flows *= dt;
}
} // namespace mio

#endif //METAPOPULATION_MOBILITY_DETAILED_H