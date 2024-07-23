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
class ExtendedMigrationEdge : public MigrationEdge<FP>
{
public:
    double travel_time;
    std::vector<int> path;

    ExtendedMigrationEdge(const MigrationParameters<FP>& params, double tt, std::vector<int> p)
        : MigrationEdge<FP>(params)
        , travel_time(tt)
        , path(p)
    {
    }

    ExtendedMigrationEdge(const Eigen::VectorXd& coeffs, double tt, std::vector<int> p)
        : MigrationEdge<FP>(coeffs)
        , travel_time(tt)
        , path(p)
    {
    }

    auto& get_migrated()
    {
        return this->m_migrated;
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
using ExtendedGraph = Graph<ExtendedNodeProperty<Sim>, ExtendedMigrationEdge<double>>;

// Default implementation when get_migration_factors is not defined for Sim
template <class Sim, std::enable_if_t<!is_expression_valid<get_migration_factors_expr_t, Sim>::value, void*> = nullptr>
auto get_migration_factors(const Sim& /*sim*/, double /*t*/, const Eigen::Ref<const Eigen::VectorXd>& y)
{
    return Eigen::VectorXd::Ones(y.rows());
}

template <typename FP>
void move_migrated(Eigen::Ref<Vector<FP>> migrated, Eigen::Ref<Vector<FP>> results_from,
                   Eigen::Ref<Vector<FP>> results_to)
{
    // before moving the commuters, we need to look for negative values in migrated and correct them.
    const auto group        = 6;
    const auto num_comparts = results_to.size() / group;

    // check for negative values in migrated
    for (Eigen::Index j = 0; j < migrated.size(); ++j) {
        if (migrated(j) < -1e-8) {
            std::cout << "Negative Value in migration detected. With value" << migrated(j) << "\n";
            auto curr_age_group = int(j / num_comparts);
            auto indx_begin     = curr_age_group * num_comparts;
            auto indx_end       = (curr_age_group + 1) * num_comparts;
            // calculate max index in indx boundaries
            Eigen::Index max_index = indx_begin;
            for (Eigen::Index i = indx_begin; i < indx_end; ++i) {
                if (migrated(i) > migrated(max_index)) {
                    max_index = i;
                }
            }
            // we assume that the solution from migrated is bettter because there is contact with other nodes
            migrated(max_index) = migrated(max_index) + migrated(j);
            migrated(j)         = migrated(j) - migrated(j);
        }
    }

    // calc sum of migrated and from
    auto sum_migrated = migrated.sum();
    auto sum_from     = results_from.sum();
    if (std::abs(sum_migrated - sum_from) < 1e-8) {
        results_from = migrated;
    }
    else {
        Eigen::VectorXd remaining_after_return = (results_from - migrated).eval();
        // auto remaining_after_return_as_vector  = std::vector<double>(
        //     remaining_after_return.data(), remaining_after_return.data() + remaining_after_return.size());
        for (Eigen::Index j = 0; j < results_to.size(); ++j) {
            if (remaining_after_return(j) < -1e-8) {
                auto curr_age_group = int(j / num_comparts);
                auto indx_begin     = curr_age_group * num_comparts;
                auto indx_end       = (curr_age_group + 1) * num_comparts;
                // calculate max index in indx boundaries
                Eigen::Index max_index = indx_begin;
                for (Eigen::Index i = indx_begin; i < indx_end; ++i) {
                    if (remaining_after_return(i) > remaining_after_return(max_index)) {
                        max_index = i;
                    }
                }

                // its possible that the max value in the boundaries is not enough to fill the negative value.
                // Therefore we have to find multiple max values
                while (remaining_after_return(max_index) + remaining_after_return(j) < -1e-10) {

                    // calculate sum between indx_begin and indx_end
                    double result_from_sum_group      = 0;
                    double result_migration_sum_group = 0;
                    for (Eigen::Index i = indx_begin; i < indx_end; ++i) {
                        result_from_sum_group += results_from(i);
                        result_migration_sum_group += migrated(i);
                    }
                    auto diff = result_from_sum_group - result_migration_sum_group;
                    if (diff < -1e-8) {
                        std::cout << "Sum of results_from is smaller than sum of migrated. Diff is "
                                  << result_from_sum_group - result_migration_sum_group << "\n";
                        // transfer values from migrated to results_from
                        for (Eigen::Index i = indx_begin; i < indx_end; ++i) {
                            results_from(i) = migrated(i);
                        }
                    }

                    results_from(j)         = results_from(j) + remaining_after_return(max_index);
                    results_from(max_index) = results_from(max_index) - remaining_after_return(max_index);
                    remaining_after_return  = (results_from - migrated).eval();

                    max_index = indx_begin;
                    for (Eigen::Index i = indx_begin; i < indx_end; ++i) {
                        if (remaining_after_return(i) > remaining_after_return(max_index)) {
                            max_index = i;
                        }
                    }
                    if (max_index == indx_begin && remaining_after_return(max_index) == 0) {
                        std::cout << "Fixing negative Value in migration not possible."
                                  << "\n";
                    }
                }

                // we assume that the solution from migrated is bettter because there is contact with other nodes
                results_from(j)         = results_from(j) - remaining_after_return(j);
                results_from(max_index) = results_from(max_index) + remaining_after_return(j);
                remaining_after_return  = (results_from - migrated).eval();
            }
        }
    }

    results_from -= migrated;
    results_to += migrated;
}

template <typename FP, class Sim>
class MobilityFunctions
{
public:
    void init_mobility(FP t, FP dt, ExtendedMigrationEdge<FP>& edge, Sim& from_sim, Sim& to_sim)
    {
        edge.get_migrated().add_time_point(
            t, (from_sim.get_result().get_last_value().array() *
                edge.get_parameters().get_coefficients().get_matrix_at(t).array() *
                get_migration_factors(from_sim, t, from_sim.get_result().get_last_value()).array())
                   .matrix());
        edge.get_return_times().add_time_point(t + dt);
        move_migrated(edge.get_migrated().get_last_value(), from_sim.get_result().get_last_value(),
                      to_sim.get_result().get_last_value());
    }

    void update_and_move(FP t, FP dt, ExtendedMigrationEdge<FP>& edge, Sim& from_sim, Sim& to_sim)
    {
        auto& integrator_node = from_sim.get_integrator();
        update_status_migrated(edge.get_migrated().get_last_value(), from_sim, integrator_node,
                               from_sim.get_result().get_last_value(), t, dt);
        move_migrated(edge.get_migrated().get_last_value(), from_sim.get_result().get_last_value(),
                      to_sim.get_result().get_last_value());
    }

    void update_only(FP t, FP dt, ExtendedMigrationEdge<FP>& edge, Sim& from_sim)
    {
        auto& integrator_node = from_sim.get_integrator();
        update_status_migrated(edge.get_migrated().get_last_value(), from_sim, integrator_node,
                               from_sim.get_result().get_last_value(), t, dt);
        move_migrated(edge.get_migrated().get_last_value(), from_sim.get_result().get_last_value(),
                      from_sim.get_result().get_last_value());
    }

    void move_and_delete(ExtendedMigrationEdge<FP>& edge, Sim& from_sim, Sim& to_sim)
    {
        move_migrated(edge.get_migrated().get_last_value(), from_sim.get_result().get_last_value(),
                      to_sim.get_result().get_last_value());

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

            // Calculate the start time for mobility, ensuring it is not negative
            const double start_mobility =
                std::max(0.0, round_nth_decimal(1 - 2 * traveltime_per_region * e.property.path.size() -
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
            size_t start_idx    = static_cast<size_t>((start_mobility + epsilon) * 100);
            size_t arrive_idx   = static_cast<size_t>((arrive_at + epsilon) * 100);
            size_t stay_end_idx = timesteps - (arrive_idx - start_idx);

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
                size_t count_true  = std::count(is_mobility_node.begin(), is_mobility_node.end(), true);

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

            for (auto& node : this->m_graph.nodes()) {
                node.property.mobility_sim.set_integrator(std::make_shared<mio::EulerIntegratorCore<double>>());
            }

            size_t indx_schedule = 0;
            while (t_begin + 1 > this->m_t + 1e-10) {
                advance_edges(indx_schedule);

                // first we integrate the nodes in time. Afterwards the update on the edges is done.
                // We start with the edges since the values for t0 are given.
                advance_local_nodes(indx_schedule);
                advance_mobility_nodes(indx_schedule);

                indx_schedule++;
                this->m_t += min_dt;
            }
            // At the end of the day. we set each compartment zero for all mobility nodes since we have to estimate
            // the state of the indivuals moving between different counties.
            // Therefore there can be differences with the states given by the integrator used for the mobility node.
            for (auto& n : this->m_graph.nodes()) {
                n.property.mobility_sim.get_result().get_last_value().setZero();
            }
        }
    }

private:
    MobilityFunctions m_mobility_functions;
    ScheduleManager::Schedule schedules;
    const double epsilon = 1e-10;

    void advance_edges(size_t indx_schedule)
    {
        for (const auto& edge_indx : schedules.edges_mobility[indx_schedule]) {
            auto& e = this->m_graph.edges()[edge_indx];
            // first mobility activity
            if (indx_schedule == schedules.first_mobility[edge_indx]) {
                auto& node_from =
                    this->m_graph.nodes()[schedules.schedule_edges[edge_indx][indx_schedule - 1]].property.base_sim;
                auto& node_to =
                    this->m_graph.nodes()[schedules.schedule_edges[edge_indx][indx_schedule]].property.mobility_sim;
                // m_edge_func(m_t, 0.0, e.property, node_from, node_to, 0);
                m_mobility_functions.init_mobility(this->m_t, 0.0, e.property, node_from, node_to);
            }
            // next mobility activity
            else if (indx_schedule > schedules.first_mobility[edge_indx]) {
                auto current_node_indx = schedules.schedule_edges[edge_indx][indx_schedule];
                bool in_mobility_node  = schedules.mobility_schedule_edges[edge_indx][indx_schedule];

                // determine dt, which is equal to the last integration/syncronization point in the current node
                auto integrator_schedule_row = schedules.local_int_schedule[current_node_indx];
                if (in_mobility_node)
                    integrator_schedule_row = schedules.mobility_int_schedule[current_node_indx];
                // search the index of indx_schedule in the integrator schedule
                const size_t indx_current = std::distance(
                    integrator_schedule_row.begin(),
                    std::lower_bound(integrator_schedule_row.begin(), integrator_schedule_row.end(), indx_schedule));

                if (integrator_schedule_row[indx_current] != indx_schedule)
                    std::cout << "Error in schedule."
                              << "\n";

                ScalarType dt_mobility;
                if (indx_current == 0) {
                    dt_mobility = round_nth_decimal(e.property.travel_time / e.property.path.size(), 2);
                    if (dt_mobility < 0.01)
                        dt_mobility = 0.01;
                }
                else {
                    dt_mobility = round_nth_decimal((static_cast<double>(integrator_schedule_row[indx_current]) -
                                                     static_cast<double>(integrator_schedule_row[indx_current - 1])) /
                                                            100 +
                                                        epsilon,
                                                    2);
                }

                // We have two cases. Either, we want to send the individuals to the next node, or we just want
                // to update their state since a syncronization step is necessary in the current node.
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

                    if (indx_schedule < schedules.mobility_schedule_edges[edge_indx].size() - 1) {
                        // m_edge_func(m_t, dt_mobility, e.property, node_from, node_to, 1);
                        m_mobility_functions.update_and_move(this->m_t, dt_mobility, e.property, node_from, node_to);
                    }
                }
                else {
                    auto& node_from =
                        schedules.mobility_schedule_edges[edge_indx][indx_schedule - 1]
                            ? this->m_graph.nodes()[schedules.schedule_edges[edge_indx][indx_schedule - 1]]
                                  .property.mobility_sim
                            : this->m_graph.nodes()[schedules.schedule_edges[edge_indx][indx_schedule - 1]]
                                  .property.base_sim;

                    assert(node_from.get_result().get_last_value() ==
                           (schedules.mobility_schedule_edges[edge_indx][indx_schedule]
                                ? this->m_graph.nodes()[schedules.schedule_edges[edge_indx][indx_schedule]]
                                      .property.mobility_sim
                                : this->m_graph.nodes()[schedules.schedule_edges[edge_indx][indx_schedule]]
                                      .property.base_sim)
                               .get_result()
                               .get_last_value());
                    // m_edge_func(m_t, dt_mobility, e.property, node_from, node_to, 3);
                    m_mobility_functions.update_only(this->m_t, dt_mobility, e.property, node_from);
                }
            }
        }
        // in the last time step we have to move the individuals back to their local model
        if (indx_schedule == 99) {
            auto edge_index = 0;
            for (auto& e : this->m_graph.edges()) {
                auto current_node_indx = schedules.schedule_edges[edge_index][indx_schedule];
                bool in_mobility_node  = schedules.mobility_schedule_edges[edge_index][indx_schedule];

                // determine dt, which is equal to the last integration/syncronization point in the current node
                auto integrator_schedule_row = schedules.local_int_schedule[current_node_indx];
                if (in_mobility_node)
                    integrator_schedule_row = schedules.mobility_int_schedule[current_node_indx];
                // search the index of indx_schedule in the integrator schedule
                const size_t indx_current = std::distance(
                    integrator_schedule_row.begin(),
                    std::lower_bound(integrator_schedule_row.begin(), integrator_schedule_row.end(), indx_schedule));

                ScalarType dt_mobility;
                if (indx_current == 0) {
                    dt_mobility = round_nth_decimal(e.property.travel_time / e.property.path.size(), 2);
                    if (dt_mobility < 0.01)
                        dt_mobility = 0.01;
                }
                else {
                    dt_mobility = round_nth_decimal((static_cast<double>(integrator_schedule_row[indx_current]) -
                                                     static_cast<double>(integrator_schedule_row[indx_current - 1])) /
                                                            100 +
                                                        epsilon,
                                                    2);
                }

                auto& node_from =
                    this->m_graph.nodes()[schedules.schedule_edges[edge_index][indx_schedule]].property.mobility_sim;

                auto& node_to =
                    this->m_graph.nodes()[schedules.schedule_edges[edge_index][indx_schedule]].property.base_sim;

                m_mobility_functions.move_and_delete(e.property, node_from, node_to);
                edge_index++;
            }
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
            // m_node_func(this->m_t, next_dt, n.property.base_sim);
        }
    }

    void advance_mobility_nodes(size_t indx_schedule)
    {
        for (const size_t& n_indx : schedules.nodes_mobility_m[indx_schedule]) {
            auto& n = this->m_graph.nodes()[n_indx];
            // determine in which index of mobility_int_schedule we are
            const size_t indx_current =
                std::distance(schedules.mobility_int_schedule[n_indx].begin(),
                              std::lower_bound(schedules.mobility_int_schedule[n_indx].begin(),
                                               schedules.mobility_int_schedule[n_indx].end(), indx_schedule));
            // determine the next time step
            const size_t val_next = (indx_current == schedules.mobility_int_schedule[n_indx].size() - 1)
                                        ? 100
                                        : schedules.mobility_int_schedule[n_indx][indx_current + 1];
            const ScalarType next_dt =
                round_nth_decimal((static_cast<double>(val_next) - indx_schedule) / 100 + epsilon, 2);

            // get all time points from the last integration step
            auto& last_time_point = n.property.mobility_sim.get_result().get_time(
                n.property.mobility_sim.get_result().get_num_time_points() - 1);
            // if the last time point is not within an interval of 1-e10 from t, then set the last time point to m_t
            if (std::fabs(last_time_point - this->m_t) > 1e-10) {
                n.property.mobility_sim.get_result().get_last_time() = this->m_t;
            }
            // only advance in time if there are individuals in the mobility model
            if (n.property.mobility_sim.get_result().get_last_value().sum() > 1e-8) {
                n.property.mobility_sim.advance(this->m_t + next_dt);
            }
            Eigen::Index indx_min;
            while (n.property.mobility_sim.get_result().get_last_value().minCoeff(&indx_min) < -1e-7) {
                Eigen::Index indx_max;
                n.property.mobility_sim.get_result().get_last_value().maxCoeff(&indx_max);
                n.property.mobility_sim.get_result().get_last_value()[indx_max] -=
                    n.property.mobility_sim.get_result().get_last_value()[indx_min];
                n.property.mobility_sim.get_result().get_last_value()[indx_min] = 0;
            }
        }
    }
};

template <typename FP, class Sim>
GraphSimulationExtended<Graph<ExtendedNodeProperty<Sim>, ExtendedMigrationEdge<FP>>, MobilityFunctions<FP, Sim>>
make_migration_sim(FP t0, FP dt, Graph<ExtendedNodeProperty<Sim>, ExtendedMigrationEdge<FP>>&& graph)
{
    auto migration_modes = MobilityFunctions<FP, Sim>();
    return GraphSimulationExtended<Graph<ExtendedNodeProperty<Sim>, ExtendedMigrationEdge<FP>>,
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
                            mio::IntegratorCore<FP>& integrator,
                            Eigen::Ref<const typename TimeSeries<FP>::Vector> total, FP t, FP dt)
{
    auto y0 = migrated.eval();
    auto y1 = migrated;
    integrator.step(
        [&](auto&& y, auto&& t_, auto&& dydt) {
            sim.get_model().get_derivatives(total, y, t_, dydt);
        },
        y0, t, dt, y1);
    auto flows_model = sim.get_model().get_flow_values();
    flows_model *= dt;
    sim.get_model().set_flow_values(flows_model);
}
} // namespace mio

#endif //METAPOPULATION_MOBILITY_DETAILED_H