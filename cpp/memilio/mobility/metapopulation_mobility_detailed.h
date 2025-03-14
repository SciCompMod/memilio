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

/**
 * @brief Aggregates two simulation models and a stay duration for a node.
 *
 * @tparam Sim Type representing the simulation or model.
 */
template <typename Sim>
struct ExtendedNodeProperty {
    Sim base_sim; ///< The base simulation model representing local dynamics.
    Sim mobility_sim; ///< The simulation model used for mobility-related processes.
    double stay_duration; ///< Duration individuals stay in the node.

    /**
     * @brief Constructs an ExtendedNodeProperty with the given simulations and stay duration.
     * @param sim1 The base simulation model.
     * @param sim2 The mobility simulation model.
     * @param t The stay duration.
     */
    ExtendedNodeProperty(Sim sim1, Sim sim2, double t)
        : base_sim(sim1)
        , mobility_sim(sim2)
        , stay_duration(t)
    {
    }
};

/**
 * @brief Extends the basic MobilityEdge with travel time and a detailed travel path.
 *
 * @tparam FP Floating-point type used (default is double).
 */
template <typename FP = double>
class ExtendedMobilityEdge : public MobilityEdge<FP>
{
public:
    double travel_time; ///< The travel time along this edge.
    std::vector<int> path; ///< A vector representing the travel path (e.g., sequence of node IDs).

    /**
     * @brief Constructs an ExtendedMobilityEdge using mobility parameters.
     * 
     * @param params The mobility parameters used to initialize the base MobilityEdge.
     * @param tt The travel time.
     * @param p A vector representing the travel path.
     */
    ExtendedMobilityEdge(const MobilityParameters<FP>& params, double tt, std::vector<int> p)
        : MobilityEdge<FP>(params)
        , travel_time(tt)
        , path(p)
    {
    }

    /**
     * @brief Constructs an ExtendedMobilityEdge using a vector of coefficients.
     *
     * @param coeffs A vector of mobility coefficients.
     * @param tt The travel time.
     * @param p A vector representing the travel path.
     */
    ExtendedMobilityEdge(const Eigen::VectorXd& coeffs, double tt, std::vector<int> p)
        : MobilityEdge<FP>(coeffs)
        , travel_time(tt)
        , path(p)
    {
    }

    /**
     * @brief Returns a reference to the mobile population that migrated along this edge.
     * 
     * @return A reference to the migrated mobile population.
     */
    auto& get_migrated()
    {
        return this->m_mobile_population;
    }

    /**
     * @brief Returns a reference to the return times associated with this edge.
     * 
     * @return A reference to the vector containing the return times.
     */
    auto& get_return_times()
    {
        return this->m_return_times;
    }

    /**
     * @brief Returns a const reference to the mobility parameters associated with this edge.
     * 
     * @return A const reference to the mobility parameters.
     */
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

/**
 * @brief Checks for negative values in an Eigen vector and corrects them.
 * 
 * The vector is assumed to be partitioned into groups corresponding to different age groups.
 * Each group has a size of (vec.size() / num_age_groups). If any element is below the specified 
 * tolerance, it is set to zero and its (negative) value is added to the maximum element in the same group.
 * This adjustment is repeated until no element is below the tolerance or the maximum iteration count is exceeded.
 *
 * @tparam FP Floating-point type.
 * @param vec An Eigen vector reference containing values to be checked and corrected.
 * @param num_age_groups The number of age groups used to determine the group size.
 * @param tolerance The threshold for negative values. Elements below this value are corrected (default is -1e-7).
 * @param max_iterations Maximum number of iterations allowed for corrections (default is 100).
 */
template <typename FP>
void check_negative_values_vec(Eigen::Ref<Eigen::VectorX<FP>> vec, const size_t num_age_groups, FP tolerance = -1e-7,
                               const size_t max_iterations = 100)
{
    // Determine the number of compartments per age group.
    const size_t num_comparts = vec.size() / num_age_groups;
    size_t iteration_count    = 0;

    // Loop until no negative values (below tolerance) remain in the vector.
    while (vec.minCoeff() < tolerance) {
        Eigen::Index min_index;
        // Get the minimum value in the vector and its index.
        const FP min_value = vec.minCoeff(&min_index);

        // Determine the current age group based on the index.
        const auto curr_age_group = min_index / num_comparts;
        // Compute the range of indices for the current age group.
        const auto indx_begin = curr_age_group * num_comparts;
        const auto indx_end   = indx_begin + num_comparts;

        // Find the index of the maximum value within the same age group.
        auto max_group_indx = indx_begin;
        for (auto i = indx_begin; i < indx_end; ++i) {
            if (vec(i) > vec(max_group_indx)) {
                max_group_indx = i;
            }
        }

        // Correct the negative value by setting it to zero and
        // adding its (negative) value to the maximum value of the group.
        vec(min_index) = 0;
        vec(max_group_indx) += min_value;

        // Check if the iteration count exceeds the allowed maximum.
        if (iteration_count > max_iterations) {
            log_error("Number of iterations exceeded in check_negative_values_vec.");
            std::exit(1);
        }
        iteration_count++;
    }
}

/**
 * @brief Finds the time index corresponding to a given time point in the simulation.
 * 
 * This function searches for a time point (t) in the simulation's results.
 * It compares the provided time with each time point stored in the simulation (starting from the latest)
 * until it finds one that differs from t by less than 1e-10. If the time point is not found and
 * create_new_tp is true, a new time point is added to the results and flows. Otherwise, an error is logged.
 *
 * @tparam FP Floating-point type used.
 * @tparam Sim Simulation type that has access to results and flows.
 * @param simulation The simulation containing the results and flows.
 * @param t The time point to search for.
 * @param create_new_tp Flag indicating whether to create a new time point if t is not found.
 * @return Eigen::Index The index corresponding to the time point.
 */
template <typename FP, class Sim>
Eigen::Index find_time_index(Sim& simulation, FP t, bool create_new_tp)
{
    // Get references to simulation results and flows.
    auto& results = simulation.get_result();
    auto& flows   = simulation.get_flows();

    // Ensure that the number of time points in results and flows match.
    if (results.get_num_time_points() != flows.get_num_time_points()) {
        log_error("Number of time points in results (" + std::to_string(results.get_num_time_points()) +
                  ") and flows (" + std::to_string(flows.get_num_time_points()) + ") do not match in find_time_index.");
    }

    // Start from the last time point index.
    Eigen::Index t_indx = results.get_num_time_points() - 1;
    // Iterate backwards until a time point is found that matches t within a small tolerance.
    for (; t_indx >= 0; --t_indx) {
        if (std::abs(results.get_time(t_indx) - t) < 1e-10) {
            break;
        }
    }

    // If no matching time point is found and create_new_tp is false, log an error.
    if (t_indx < 0 && !create_new_tp) {
        log_error("Time point " + std::to_string(t) + " not found in find_time_index. Latest time point is " +
                  std::to_string(results.get_last_time()));
    }

    // If t is greater than the last time point and create_new_tp is enabled, create a new time point.
    if (t_indx < 0 && results.get_last_time() < t && create_new_tp) {
        // Initialize a new result vector with zeros.
        Eigen::VectorXd results_vec = Eigen::VectorXd::Zero(results.get_last_value().size());
        results.add_time_point(t, results_vec);
        // For flows, copy the last value (since flows are accumulated).
        flows.add_time_point(t, flows.get_last_value());
        t_indx = results.get_num_time_points() - 1;
    }

    return t_indx;
}

/**
 * @brief Provides mobility-related functions for simulation models.

 * @tparam FP Floating-point type used.
 * @tparam Sim Simulation type that has access to results, flows, and model.
 */
template <typename FP, class Sim>
class MobilityFunctions
{
public:
    /**
     * @brief Initializes mobility for a given edge at time t.
     *
     * This function finds the starting time point in the source simulation, calculates the number
     * of commuters to be moved based on the simulation results and mobility coefficients, and adds
     * a corresponding time point to the return times of the edge. It then moves the commuters from the
     * source simulation to the target simulation.
     *
     * @param t The current time.
     * @param edge The extended mobility edge to be initialized.
     * @param sim_source The source simulation from which commuters are taken.
     * @param sim_target The target simulation to which commuters are added.
     */
    void init_mobility(FP t, ExtendedMobilityEdge<FP>& edge, Sim& sim_source, Sim& sim_target)
    {
        // Find the time index in the source simulation without creating a new time point.
        const auto t_indx_start_mobility_sim_from = find_time_index(sim_source, t, false);

        // Initialize the number of commuters at the start of mobility by computing
        // a weighted product of the result, mobility coefficients and mobility factors.
        auto results_from = sim_source.get_result();
        edge.get_migrated().add_time_point(
            t, (results_from.get_value(t_indx_start_mobility_sim_from).array() *
                edge.get_parameters().get_coefficients().get_matrix_at(t).array() *
                get_mobility_factors(sim_source, t, results_from.get_value(t_indx_start_mobility_sim_from)).array())
                   .matrix());
        // Save the time point for the return times.
        edge.get_return_times().add_time_point(t);

        // Move commuters to the target simulation.
        // If the target simulation does not have the time point, create one.
        const auto t_indx_start_mobility_sim_to = find_time_index(sim_target, t, true);
        sim_target.get_result().get_value(t_indx_start_mobility_sim_to) += edge.get_migrated().get_last_value();
        sim_source.get_result().get_last_value() -= edge.get_migrated().get_last_value();
    }

    /**
     * @brief Moves the migrated population from one simulation to another.
     *
     * This function transfers commuters from the source simulation to the target simulation,
     * updating the respective result vectors. It also calls check_negative_values_vec to ensure
     * that no negative values remain in the population vectors after moving.
     *
     * @param t The current time.
     * @param edge The extended mobility edge holding the migrated population.
     * @param sim_source The simulation from which commuters are removed.
     * @param sim_target The simulation to which commuters are added.
     */
    void move_migrated(FP t, ExtendedMobilityEdge<FP>& edge, Sim& sim_source, Sim& sim_target)
    {
        // Determine the number of age groups from the simulation's model parameters.
        const size_t num_age_groups = static_cast<size_t>(sim_source.get_model().parameters.get_num_groups());
        // Correct any negative values in the migrated population vector.
        check_negative_values_vec(edge.get_migrated().get_last_value(), num_age_groups);

        // Find the arrival time index in the target simulation, creating a new time point if needed.
        const auto t_indx_sim_to_arrival = find_time_index(sim_target, t, true);
        // Remove the migrated population from the source simulation.
        sim_source.get_result().get_last_value() -= edge.get_migrated().get_last_value();
        // Add the migrated population to the target simulation.
        sim_target.get_result().get_value(t_indx_sim_to_arrival) += edge.get_migrated().get_last_value();

        // Re-check both simulations for negative values after the move.
        check_negative_values_vec(sim_source.get_result().get_last_value(), num_age_groups);
        check_negative_values_vec(sim_target.get_result().get_value(t_indx_sim_to_arrival), num_age_groups);
    }

    /**
     * @brief Updates the status of commuters (migrated population) in the simulation.
     *
     * @param t The current time.
     * @param dt The time increment for updating commuters.
     * @param edge The extended mobility edge holding the migrated population.
     * @param sim The simulation whose mobility model is being updated.
     * @param is_mobility_model Flag indicating whether the simulation is a mobility model.
     */
    void update_commuters(FP t, FP dt, ExtendedMobilityEdge<FP>& edge, Sim& sim, bool is_mobility_model)
    {
        // Find or create the time point for the current time in the simulation.
        const auto t_indx_start_mobility_sim = find_time_index(sim, t, true);
        // Initialize flows vector with zeros.
        Eigen::VectorXd flows = Eigen::VectorXd::Zero(sim.get_flows().get_last_value().size());
        // Update the status of migrated commuters and update flows accordingly.
        update_status_migrated(edge.get_migrated().get_last_value(), sim,
                               sim.get_result().get_value(t_indx_start_mobility_sim), t, dt, flows);

        // If the simulation holds a mobility model, update it by creating a new time point if necessary.
        if (is_mobility_model) {
            const auto t_indx_mobility_model = find_time_index(sim, t + dt, true);
            if (t_indx_mobility_model != sim.get_result().get_num_time_points() - 1) {
                log_error("Time point " + std::to_string(t + dt) +
                          " not the latest in update_commuters. Latest time point is " +
                          std::to_string(sim.get_result().get_last_time()));
            }
            sim.get_result().get_value(t_indx_mobility_model) += edge.get_migrated().get_last_value();
            sim.get_flows().get_value(t_indx_mobility_model) += flows;
        }
    }

    /**
     * @brief Deletes all time points from the migrated and return times of an edge.
     *
     * @param edge The extended mobility edge whose time points are to be removed.
     */
    void delete_migrated(ExtendedMobilityEdge<FP>& edge)
    {
        // Remove time points in reverse order to avoid index shifting issues.
        for (Eigen::Index i = edge.get_return_times().get_num_time_points() - 1; i >= 0; --i) {
            edge.get_migrated().remove_time_point(i);
            edge.get_return_times().remove_time_point(i);
        }
    }
};
/**
 * @brief Manages the computation of schedules for a graph's edges and nodes.
 */
class ScheduleManager
{
public:
    /**
     * @brief Holds all schedule-related vectors.
     * Exlanation of the vectors:   
     * - "schedule_edges": For each edge, this vector contains the node index where an individual is located at each timestep for a day.
     * - "mobility_schedule_edges": For each edge and each timestep, a boolean flag indicating if the individual is in a mobility model.
     * - "mobility_int_schedule": For each node, a vector of timesteps at which the mobility model needs to be synchronized.
     * - "local_int_schedule": For each node, a vector of timesteps at which the local model needs to be synchronized.
     * - "edges_mobility": For each timestep, a list of edge indices where mobility is happening.
     * - "nodes_mobility": For each timestep, a list of node indices where the local model is exchanging individuals.
     * - "nodes_mobility_m": For each timestep, a list of node indices where the mobility model is exchanging individuals
     * - "first_mobility": For each edge, the first timestep at which mobility occurs.
     */
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

    /**
     * @brief Construct a new ScheduleManager object.
     * @param n_t Total number of timesteps.
     * @param eps A epsilon as tolerance.
     */
    ScheduleManager(size_t n_t, double eps = 1e-10)
        : n_timesteps(n_t)
        , epsilon(eps)
    {
    }

    /**
     * @brief Compute the complete schedule for the given graph.
     *
     * This function calculates both the edge schedule and the node schedule.
     * @tparam Graph Type of the graph.
     * @param graph The graph object containing nodes and edges.
     * @return Schedule The computed schedule.
     */
    template <typename Graph>
    Schedule compute_schedule(const Graph& graph)
    {
        Schedule schedule;
        calculate_edge_schedule(graph, schedule);
        calculate_node_schedule(graph, schedule);
        return schedule;
    }

private:
    size_t n_timesteps; ///< Total number of timesteps.
    double epsilon; ///< Tolerance

    /**
     * @brief Fills a range within a vector with a specified value.
     * @tparam T The type of elements in the vector.
     * @param vec The vector to fill.
     * @param start The starting index.
     * @param end The ending index.
     * @param value The value to fill with.
     */
    template <typename T>
    void fill_vector_range(std::vector<T>& vec, size_t start, size_t end, const T& value)
    {
        if (start < vec.size() && end <= vec.size() && start < end) {
            std::fill(vec.begin() + start, vec.begin() + end, value);
        }
    }

    /**
     * @brief Get indices of edges while iterating over all edges in the provided schedule and selects those edges
     * for which the node index at timestep t is equal to node_Idx, and the corresponding mobility flag matches the mobility parameter.
     *
     * @param schedule The precomputed schedule containing the node positions and mobility flags for each edge.
     * @param node_Idx The index of the node to check against the edge schedules.
     * @param t The timestep at which to evaluate the edge positions and mobility states.
     * @param mobility The mobility flag to match; set to true to select edges in mobility mode, false otherwise.
     * @return std::vector<size_t> A vector containing the indices of edges that satisfy both conditions.
     */
    std::vector<size_t> find_edges_at_time(const Schedule& schedule, size_t node_idx, size_t t, bool mobility) const
    {
        std::vector<size_t> edges;
        // Loop over all edges in the schedule.
        for (size_t edge_idx = 0; edge_idx < schedule.schedule_edges.size(); ++edge_idx) {
            // Find edges located at node_idx at timestep t with the specified mobility flag.
            if (schedule.schedule_edges[edge_idx][t] == node_idx &&
                schedule.mobility_schedule_edges[edge_idx][t] == mobility) {
                edges.push_back(edge_idx);
            }
        }
        return edges;
    }

    /**
     * @brief Calculates the schedule for all edges in the graph.
     *
     * @tparam Graph Type of the graph.
     * @param graph The graph object.
     * @param schedule The schedule object to populate.
     */
    template <typename Graph>
    void calculate_edge_schedule(const Graph& graph, Schedule& schedule)
    {
        schedule.schedule_edges.reserve(graph.edges().size());
        schedule.mobility_schedule_edges.reserve(graph.edges().size());

        for (auto& e : graph.edges()) {
            std::vector<size_t> edge_schedule(n_timesteps, 0);
            std::vector<bool> mobility_flag(n_timesteps, false);

            // Calculate travel time per region (minimum 0.01)
            const double travel_time_per_region =
                std::max(0.01, round_nth_decimal(e.property.travel_time / e.property.path.size(), 2));

            // Calculate start time for mobility (minimum 0.01)
            const double start_mobility =
                std::max(0.01, round_nth_decimal(1 - 2 * travel_time_per_region * e.property.path.size() -
                                                     graph.nodes()[e.end_node_idx].property.stay_duration,
                                                 2));

            // Calculate arrival time at destination
            const double arrival_time = start_mobility + travel_time_per_region * e.property.path.size();

            // Convert times to indices
            const size_t start_index   = static_cast<size_t>((start_mobility + epsilon) * 100);
            const size_t arrival_index = static_cast<size_t>((arrival_time + epsilon) * 100);
            const size_t end_index     = n_timesteps - (arrival_index - start_index);

            // Fill schedule before mobility with start node index.
            fill_vector_range(edge_schedule, 0, start_index, e.start_node_idx);

            // Fill mobility flag during mobility period.
            fill_vector_range(mobility_flag, start_index, arrival_index, true);

            // Fill schedule during mobility along the given path.
            size_t current_index = start_index;
            for (size_t county : e.property.path) {
                size_t next_index = current_index + static_cast<size_t>((travel_time_per_region + epsilon) * 100);
                fill_vector_range(edge_schedule, current_index, next_index, county);
                current_index = next_index;
            }

            // Fill remaining schedule after mobility with destination node.
            fill_vector_range(edge_schedule, current_index, end_index, e.property.path.back());

            // Fill return mobility period after staying.
            fill_vector_range(mobility_flag, end_index, n_timesteps, true);

            // Reverse the mobility for the return trip.
            auto first_true = std::find(mobility_flag.begin(), mobility_flag.end(), true);
            auto last_true  = std::find(mobility_flag.rbegin(), mobility_flag.rend(), true);
            if (first_true != mobility_flag.end() && last_true != mobility_flag.rend()) {
                size_t first_index_found = std::distance(mobility_flag.begin(), first_true);
                size_t mobility_duration = arrival_index - start_index;
                std::vector<size_t> reversed_path(edge_schedule.begin() + first_index_found,
                                                  edge_schedule.begin() + first_index_found + mobility_duration);
                std::reverse(reversed_path.begin(), reversed_path.end());
                std::copy(reversed_path.begin(), reversed_path.end(), edge_schedule.end() - mobility_duration);

                // Store the computed schedules.
                schedule.schedule_edges.push_back(std::move(edge_schedule));
                schedule.mobility_schedule_edges.push_back(std::move(mobility_flag));
            }
            else {
                log_error("Error in creating schedule for an edge.");
            }
        }
    }

    /**
     * @brief Calculates the node schedule based on edge schedules.
     *
     * @tparam Graph Type of the graph.
     * @param graph The graph object.
     * @param schedule The schedule object to populate.
     */
    template <typename Graph>
    void calculate_node_schedule(const Graph& graph, Schedule& schedule)
    {
        // Compute integration schedules per node.
        for (size_t node_idx = 0; node_idx < graph.nodes().size(); ++node_idx) {
            // Local integration schedule.
            std::vector<size_t> local_schedule{0}; // Always start at t=0.
            auto current_edges = find_edges_at_time(schedule, node_idx, 0, false);
            for (size_t t = 1; t < n_timesteps; ++t) {
                auto new_edges = find_edges_at_time(schedule, node_idx, t, false);
                if (new_edges.size() != current_edges.size() ||
                    !std::equal(current_edges.begin(), current_edges.end(), new_edges.begin())) {
                    local_schedule.push_back(t);
                    current_edges = new_edges;
                }
            }
            if (local_schedule.back() != n_timesteps - 1) {
                local_schedule.push_back(n_timesteps - 1);
            }
            schedule.local_int_schedule.push_back(local_schedule);

            // Mobility integration schedule.
            std::vector<size_t> mobility_schedule;
            auto current_mob_edges = find_edges_at_time(schedule, node_idx, 0, true);
            for (size_t t = 1; t < n_timesteps; ++t) {
                auto new_mob_edges = find_edges_at_time(schedule, node_idx, t, true);
                if (new_mob_edges.size() != current_mob_edges.size() ||
                    !std::equal(current_mob_edges.begin(), current_mob_edges.end(), new_mob_edges.begin())) {
                    mobility_schedule.push_back(t);
                    current_mob_edges = new_mob_edges;
                }
            }
            schedule.mobility_int_schedule.push_back(mobility_schedule);
        }

        // Determine the first timestep with mobility for each edge.
        schedule.first_mobility.reserve(graph.edges().size());
        for (size_t edge_idx = 0; edge_idx < graph.edges().size(); ++edge_idx) {
            size_t t = 0;
            for (; t < schedule.mobility_schedule_edges[edge_idx].size(); ++t) {
                if (schedule.mobility_schedule_edges[edge_idx][t]) {
                    break;
                }
            }
            schedule.first_mobility.push_back(t);
        }

        // Initialize mobility-related vectors per timestep.
        schedule.edges_mobility.reserve(n_timesteps);
        schedule.nodes_mobility.reserve(n_timesteps);
        schedule.nodes_mobility_m.reserve(n_timesteps);

        // At t = 0: collect edges with mobility starting at 0.
        std::vector<size_t> initial_edges;
        for (size_t edge_idx = 0; edge_idx < graph.edges().size(); ++edge_idx) {
            if (schedule.first_mobility[edge_idx] == 0) {
                initial_edges.push_back(edge_idx);
            }
        }
        schedule.edges_mobility.push_back(std::move(initial_edges));

        // At t = 0: initialize with all nodes.
        std::vector<size_t> initial_nodes(graph.nodes().size());
        std::iota(initial_nodes.begin(), initial_nodes.end(), 0);
        schedule.nodes_mobility.push_back(std::move(initial_nodes));

        // At t = 0: nodes with mobility activity (if present in mobility_int_schedule).
        std::vector<size_t> initial_nodes_mob_m;
        for (size_t node_idx = 0; node_idx < graph.nodes().size(); ++node_idx) {
            if (std::binary_search(schedule.mobility_int_schedule[node_idx].begin(),
                                   schedule.mobility_int_schedule[node_idx].end(), 0)) {
                initial_nodes_mob_m.push_back(node_idx);
            }
        }
        schedule.nodes_mobility_m.push_back(std::move(initial_nodes_mob_m));

        // For each subsequent timestep, update mobility edge and node lists.
        for (size_t t = 1; t < n_timesteps; ++t) {
            // Identify mobility-active edges at timestep t.
            std::vector<size_t> edges_at_t;
            for (size_t edge_idx = 0; edge_idx < graph.edges().size(); ++edge_idx) {
                size_t current_node = schedule.schedule_edges[edge_idx][t];
                if (t >= schedule.first_mobility[edge_idx]) {
                    if (schedule.mobility_schedule_edges[edge_idx][t] &&
                        std::binary_search(schedule.mobility_int_schedule[current_node].begin(),
                                           schedule.mobility_int_schedule[current_node].end(), t)) {
                        edges_at_t.push_back(edge_idx);
                    }
                    else if (!schedule.mobility_schedule_edges[edge_idx][t] &&
                             std::binary_search(schedule.local_int_schedule[current_node].begin(),
                                                schedule.local_int_schedule[current_node].end(), t)) {
                        edges_at_t.push_back(edge_idx);
                    }
                }
            }
            schedule.edges_mobility.push_back(edges_at_t);

            // Identify nodes with local and mobility integration at timestep t.
            std::vector<size_t> nodes_local;
            std::vector<size_t> nodes_mob_m;
            for (size_t node_idx = 0; node_idx < graph.nodes().size(); ++node_idx) {
                if (std::binary_search(schedule.local_int_schedule[node_idx].begin(),
                                       schedule.local_int_schedule[node_idx].end(), t)) {
                    nodes_local.push_back(node_idx);
                }
                if (std::binary_search(schedule.mobility_int_schedule[node_idx].begin(),
                                       schedule.mobility_int_schedule[node_idx].end(), t)) {
                    nodes_mob_m.push_back(node_idx);
                }
            }
            schedule.nodes_mobility.push_back(nodes_local);
            schedule.nodes_mobility_m.push_back(nodes_mob_m);
        }
    }
};

template <typename Graph, typename MobilityFunctions>
class GraphSimulationExtended
    : public GraphSimulationBase<Graph, double, double,
                                 std::function<void(typename Graph::EdgeProperty&, size_t,
                                                    typename Graph::NodeProperty&, typename Graph::NodeProperty&)>,
                                 std::function<void(double, double, typename Graph::NodeProperty&)>>
{
public:
    using node_function = std::function<void(double, double, typename Graph::NodeProperty&)>;
    using edge_function =
        std::function<void(double, double, typename Graph::EdgeProperty&, typename Graph::NodeProperty&,
                           typename Graph::NodeProperty&, MobilityFunctions&)>;

    GraphSimulationExtended(double t0, double dt, Graph& g, const node_function& node_func, MobilityFunctions modes)
        : GraphSimulationBase<Graph, double, double,
                              std::function<void(typename Graph::EdgeProperty&, size_t, typename Graph::NodeProperty&,
                                                 typename Graph::NodeProperty&)>,
                              std::function<void(double, double, typename Graph::NodeProperty&)>>(t0, dt, std::move(g),
                                                                                                  node_func, {})
        , m_mobility_functions(modes)
    {
        ScheduleManager schedule_manager(100); // Assuming 100 timesteps
        schedules = schedule_manager.compute_schedule(this->m_graph);
    }

    GraphSimulationExtended(double t0, double dt, Graph&& g, const node_function& node_func, MobilityFunctions modes)
        : GraphSimulationBase<Graph, double, double,
                              std::function<void(typename Graph::EdgeProperty&, size_t, typename Graph::NodeProperty&,
                                                 typename Graph::NodeProperty&)>,
                              std::function<void(double, double, typename Graph::NodeProperty&)>>(
              t0, dt, std::forward<Graph>(g), node_func, {})
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

    void update_status_commuters(size_t indx_schedule, const double max_num_contacts = 20.)
    {
        for (const auto& edge_indx : schedules.edges_mobility[indx_schedule]) {
            auto& e      = this->m_graph.edges()[edge_indx];
            auto next_dt = calculate_next_dt(edge_indx, indx_schedule);
            auto& node_to =
                schedules.mobility_schedule_edges[edge_indx][indx_schedule]
                    ? this->m_graph.nodes()[schedules.schedule_edges[edge_indx][indx_schedule]].property.mobility_sim
                    : this->m_graph.nodes()[schedules.schedule_edges[edge_indx][indx_schedule]].property.base_sim;

            // get current contact  and scale it but only if mobility model
            auto contact_pattern_curr = node_to.get_model().get_contact_pattern();
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
                node_to.get_model().set_contact_pattern(contact_pattern_curr);
            }

            m_mobility_functions.update_commuters(this->m_t, next_dt, e.property, node_to,
                                                  schedules.mobility_schedule_edges[edge_indx][indx_schedule]);

            // reset contact pattern after estimating the state of the commuters
            if (schedules.mobility_schedule_edges[edge_indx][indx_schedule])
                node_to.get_model().set_contact_pattern(contacts_copy);
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

    void handle_last_time_step(size_t indx_schedule)
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
