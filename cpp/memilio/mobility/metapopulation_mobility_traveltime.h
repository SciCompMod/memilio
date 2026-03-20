/* 
* Copyright (C) 2020-2026 MEmilio
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
#ifndef MIO_MOBILITY_TRAVELTIME_H
#define MIO_MOBILITY_TRAVELTIME_H

#include "memilio/mobility/graph.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/math/eigen.h"
#include "memilio/math/euler.h"
#include "memilio/math/floating_point.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <limits>
#include <numeric>
#include <vector>

namespace mio
{

/**
 * @brief Property for a graph node in the travel-time-aware mobility model.
 *
 * Each node holds two SimulationNode instances:
 * - `local_sim`:    Standard local ODE simulation
 * - `mobility_sim`: ODE simulation representing individuals currently in transit
 *                   through this node's mobility model
 *
 * Additionally, `stay_duration` defines how long inbound travelers remain at
 * this node before returning home. On a day-scale of [0,1], this is the fraction
 * of the day spent at the destination.
 *
 * @tparam FP Floating-point type (e.g., double).
 * @tparam Sim Simulation type (e.g., mio::osecir::Simulation<double>).
 */
template <typename FP, class Sim>
struct TravelTimeNodeProperty {

    /**
     * @brief Constructor if both nodes share the same model and t0.
     *
     * @param model  The epidemiological model.
     * @param t0     Start time.
     * @param stay   Stay duration in [0, 1).
     */
    TravelTimeNodeProperty(const typename Sim::Model& model, FP t0, FP stay)
        : local_sim(model, t0)
        , mobility_sim(model, t0)
        , stay_duration(stay)
    {
    }

    /**
     * @brief Construct from two simulation arguments and a stay duration.
     *
     * @param local_args   Arguments forwarded to local_sim.
     * @param mobility_args Arguments forwarded to mobility_sim.
     * @param stay         Stay duration in [0, 1).
     */
    template <class... LocalArgs, class... MobArgs>
    TravelTimeNodeProperty(std::tuple<LocalArgs...> local_args, std::tuple<MobArgs...> mobility_args, FP stay)
        : local_sim(std::make_from_tuple<SimulationNode<FP, Sim>>(std::move(local_args)))
        , mobility_sim(std::make_from_tuple<SimulationNode<FP, Sim>>(std::move(mobility_args)))
        , stay_duration(stay)
    {
    }

    SimulationNode<FP, Sim> local_sim; ///< Local population dynamics.
    SimulationNode<FP, Sim> mobility_sim; ///< Transit population dynamics.
    FP stay_duration; ///< Fraction of day spent at destination, in [0, 1).
};

/**
 * @brief Edge property for the travel-time-aware mobility model.
 *
 * Extends the standard MobilityParameters with:
 * - `travel_time`: Total travel time for one direction (as fraction of a day).
 * - `path`:        Ordered list of node indices traversed during the trip
 *                  (includes origin and destination, i.e., [from, r1, r2, ..., to]).
 *
 * The travel time is divided equally among all units in the path.
 *
 * @tparam FP Floating-point type.
 */
template <typename FP>
struct TravelTimeEdge {

    /**
     * @brief Construct from mobility parameters, travel time, and path.
     * @param params      Mobility coefficients.
     * @param tt          One-way travel time (fraction of day).
     * @param route       Ordered list of node indices (origin, intermediates, destination).
     */
    TravelTimeEdge(const MobilityParameters<FP>& params, FP tt, std::vector<size_t> route)
        : parameters(params)
        , travel_time(tt)
        , path(std::move(route))
        , mobile_population(params.get_coefficients().get_shape().rows())
        , return_times(0)
    {
    }

    /**
     * @brief Construct from coefficient vector, travel time, and path.
     * @param coeffs      Mobility coefficient vector.
     * @param tt          One-way travel time (fraction of day).
     * @param route       Ordered list of node indices.
     */
    TravelTimeEdge(const Eigen::VectorX<FP>& coeffs, FP tt, std::vector<size_t> route)
        : parameters(coeffs)
        , travel_time(tt)
        , path(std::move(route))
        , mobile_population(coeffs.rows())
        , return_times(0)
    {
    }

    MobilityParameters<FP> parameters; ///< Mobility coefficients and NPIs.
    FP travel_time; ///< One-way travel time as fraction of day.
    std::vector<size_t> path; ///< Node indices along the route (origin -> destination).

    // Commuter state tracking (one entry per commuter group, for advance function).
    TimeSeries<FP> mobile_population; ///< Compartment counts for people currently in transit.
    TimeSeries<FP> return_times; ///< Times at which the mobile population returns.
};

/**
 * @brief Correct small negative compartment values caused by the auxiliary Euler
 * step to estimate the infection states of travelers.
 *
 * The vector is partitioned into `num_age_groups` equally sized blocks.
 * Any value below `tolerance` is zeroed out and its deficit is added to the
 * largest compartment in the same age-group block. This preserves the total
 * population within each age group.
 *
 * @tparam FP           Floating-point type.
 * @param vec           Compartment vector to correct (in-place).
 * @param num_age_groups Number of age groups.
 * @param tolerance     Threshold below which a value is considered negative.
 * @param max_iter      Maximum number of correction iterations.
 */

template <typename FP>
void correct_negative_compartments(Eigen::Ref<Eigen::VectorX<FP>> vec, size_t num_age_groups, FP tolerance = -1e-7,
                                   size_t max_iter = 100)
{
    const size_t n_comparts = static_cast<size_t>(vec.size()) / num_age_groups;
    for (size_t iter = 0; iter < max_iter; ++iter) {
        if (vec.minCoeff() >= tolerance) {
            return;
        }
        Eigen::Index min_idx;
        const FP min_val = vec.minCoeff(&min_idx);
        const auto grp   = static_cast<size_t>(min_idx) / n_comparts;
        const auto beg   = static_cast<Eigen::Index>(grp * n_comparts);
        const auto end   = beg + static_cast<Eigen::Index>(n_comparts);
        Eigen::Index max_idx;
        vec.segment(beg, end - beg).maxCoeff(&max_idx);
        max_idx += beg;
        vec(max_idx) += min_val;
        vec(min_idx) = FP{0};
    }
    if (vec.minCoeff() < tolerance) {
        log_error("correct_negative_compartments: could not correct all negative values after {} iterations.",
                  max_iter);
    }
}

/**
 * @brief Daily mobility schedule computed once from the graph topology.
 *
 * The day is discretised into `n_steps` equal time slices of width `1/n_steps`.
 * For every edge the schedule records:
 * - `node_at_step[edge][step]`:     which node index holds the travelers at that step.
 * - `in_mobility[edge][step]`:      true, if the travelers are in a mobility model in that step.
 * - `first_mobility_step[edge]`:    first step at which mobility starts.
 *
 * Derived per-node schedules tell the integrator when it must stop and
 * synchronise with the mobility exchange:
 * - `local_breakpoints[node]`:      steps at which local_sim must be synchronized.
 * - `mobility_breakpoints[node]`:   steps at which mobility_sim must be synchronized.
 */
struct TravelTimeSchedule {
    std::vector<std::vector<size_t>> node_at_step; ///< [edge][step] -> node index
    std::vector<std::vector<bool>> in_mobility; ///< [edge][step] -> mobility flag
    std::vector<size_t> first_mobility_step; ///< [edge] -> first step with mobility
    std::vector<std::vector<size_t>> local_breakpoints; ///< [node] -> sorted step indices
    std::vector<std::vector<size_t>> mobility_breakpoints; ///< [node] -> sorted step indices
    /// For each step: which edges exchange travelers
    std::vector<std::vector<size_t>> edges_at_step;
    /// For each step: which nodes advance their local_sim
    std::vector<std::vector<size_t>> local_nodes_at_step;
    /// For each step: which nodes advance their mobility_sim
    std::vector<std::vector<size_t>> mobility_nodes_at_step;
};

/**
 * @brief Compute the TravelTimeSchedule for a given graph.
 *
 * @tparam FP           Floating-point type.
 * @tparam NodeProp     TravelTimeNodeProperty type.
 * @tparam EdgeProp     TravelTimeEdge type.
 * @param graph         The simulation graph.
 * @param n_steps       Number of time steps per day (default: 100).
 * @param eps           Floating-point tolerance for time comparisons.
 * @return              Fully populated TravelTimeSchedule.
 */
template <typename FP, class NodeProp, class EdgeProp>
TravelTimeSchedule compute_traveltime_schedule(const Graph<NodeProp, EdgeProp>& graph, size_t n_steps = 100,
                                               double eps = 1e-10)
{
    const size_t n_edges = graph.edges().size();
    const size_t n_nodes = graph.nodes().size();
    const double dt_step = 1.0 / static_cast<double>(n_steps);

    TravelTimeSchedule sched;
    sched.node_at_step.resize(n_edges);
    sched.in_mobility.resize(n_edges);
    sched.first_mobility_step.resize(n_edges, n_steps);

    // -- 1. Build per-edge schedules --------------------------------------
    for (size_t ei = 0; ei < n_edges; ++ei) {
        const auto& e          = graph.edges()[ei];
        const auto& edge       = e.property;
        const size_t dest_node = e.end_node_idx;
        const size_t src_node  = e.start_node_idx;
        const FP stay          = graph.nodes()[dest_node].property.stay_duration;
        const size_t n_patches = edge.path.size(); // includes start and end nodes
        // Travel time per patch (split evenly).
        const double tt_per_patch =
            std::max(dt_step, std::round(edge.travel_time / static_cast<double>(n_patches) / dt_step) * dt_step);
        const double total_travel = tt_per_patch * static_cast<double>(n_patches);

        // Start of outbound travel: everything before is local.
        const double t_depart    = std::max(0.0, 1.0 - 2.0 * total_travel - static_cast<double>(stay));
        const size_t step_depart = static_cast<size_t>((t_depart + eps) / dt_step);
        // Arrival at destination.
        const size_t steps_travel = static_cast<size_t>((total_travel + eps) / dt_step);
        const size_t step_arrive  = step_depart + steps_travel;
        // Departure from destination for return trip.
        const size_t step_leave_dst = n_steps - steps_travel;

        sched.node_at_step[ei].resize(n_steps, src_node);
        sched.in_mobility[ei].resize(n_steps, false);

        // Outbound travel: step_depart … step_arrive-1
        sched.first_mobility_step[ei] = step_depart;
        size_t current_step           = step_depart;
        for (size_t patch = 0; patch < n_patches; ++patch) {
            const size_t patch_steps = static_cast<size_t>((tt_per_patch + eps) / dt_step);
            for (size_t s = 0; s < patch_steps && (current_step + s) < n_steps; ++s) {
                sched.node_at_step[ei][current_step + s] = edge.path[patch];
                sched.in_mobility[ei][current_step + s]  = true;
            }
            current_step += patch_steps;
        }

        // Stay at destination: step_arrive … step_leave_dst-1
        for (size_t s = step_arrive; s < step_leave_dst && s < n_steps; ++s) {
            sched.node_at_step[ei][s] = dest_node;
            sched.in_mobility[ei][s]  = false;
        }

        // Return travel: step_leave_dst … n_steps-1 (reversed path)
        current_step = step_leave_dst;
        for (size_t patch = n_patches; patch > 0 && current_step < n_steps; --patch) {
            const size_t patch_steps = static_cast<size_t>((tt_per_patch + eps) / dt_step);
            for (size_t s = 0; s < patch_steps && (current_step + s) < n_steps; ++s) {
                sched.node_at_step[ei][current_step + s] = edge.path[patch - 1];
                sched.in_mobility[ei][current_step + s]  = true;
            }
            current_step += patch_steps;
        }
    }

    // -- 2. Derive per-node breakpoint schedules ---------------------------
    sched.local_breakpoints.resize(n_nodes);
    sched.mobility_breakpoints.resize(n_nodes);

    for (size_t ni = 0; ni < n_nodes; ++ni) {
        sched.local_breakpoints[ni].push_back(0);
        sched.mobility_breakpoints[ni].push_back(0);

        auto prev_local_edges = std::vector<size_t>{};
        auto prev_mob_edges   = std::vector<size_t>{};
        for (size_t t = 1; t < n_steps; ++t) {
            std::vector<size_t> local_e, mob_e;
            for (size_t ei = 0; ei < n_edges; ++ei) {
                if (sched.node_at_step[ei][t] == ni) {
                    if (sched.in_mobility[ei][t])
                        mob_e.push_back(ei);
                    else
                        local_e.push_back(ei);
                }
            }
            if (local_e != prev_local_edges) {
                sched.local_breakpoints[ni].push_back(t);
                prev_local_edges = local_e;
            }
            if (mob_e != prev_mob_edges) {
                sched.mobility_breakpoints[ni].push_back(t);
                prev_mob_edges = mob_e;
            }
        }
        // Ensure final step is always present.
        if (sched.local_breakpoints[ni].back() != n_steps - 1)
            sched.local_breakpoints[ni].push_back(n_steps - 1);
        if (sched.mobility_breakpoints[ni].back() != n_steps - 1)
            sched.mobility_breakpoints[ni].push_back(n_steps - 1);
    }

    // -- 3. Build per-step lookup tables ----------------------------------
    sched.edges_at_step.resize(n_steps);
    sched.local_nodes_at_step.resize(n_steps);
    sched.mobility_nodes_at_step.resize(n_steps);

    for (size_t t = 0; t < n_steps; ++t) {
        // Edges active at this step (first mobility step or node/mode change).
        for (size_t ei = 0; ei < n_edges; ++ei) {
            if (t < sched.first_mobility_step[ei])
                continue;
            bool add = (t == sched.first_mobility_step[ei]);
            if (!add && t > 0) {
                add = (sched.node_at_step[ei][t] != sched.node_at_step[ei][t - 1]) ||
                      (sched.in_mobility[ei][t] != sched.in_mobility[ei][t - 1]);
            }
            if (add)
                sched.edges_at_step[t].push_back(ei);
        }
        // Nodes whose local_sim breakpoint falls at t.
        for (size_t ni = 0; ni < n_nodes; ++ni) {
            if (std::binary_search(sched.local_breakpoints[ni].begin(), sched.local_breakpoints[ni].end(), t))
                sched.local_nodes_at_step[t].push_back(ni);
            if (std::binary_search(sched.mobility_breakpoints[ni].begin(), sched.mobility_breakpoints[ni].end(), t))
                sched.mobility_nodes_at_step[t].push_back(ni);
        }
    }

    return sched;
}

/**
 * @brief Graph simulation with travel-time-aware mobility.
 *
 * This class implements the travel-time-aware metapopulation model from:
 *   H. Zunker et al., "Novel travel time aware metapopulation models ...",
 *   (2024), https://doi.org/10.1371/journal.pcbi.1012630
 *
 * Key differences to `metapopulation_mobility_instant.h`:
 * - Realistic travel and stay durations
 * - Transmission on travel (even with other travelers from different origins/destinations)
 * - Explicit trip paths with intermediate nodes (e.g., for transit hubs)
 * - More complex day schedule with multiple synchronisation points
 *
 * Day structure (normalised to [0, 1]):
 * 0 --- local --> t_depart --- outbound transit --> t_arrive --- stay --> t_return --- inbound transit --> 1
 *
 * The day is discretised into `n_steps`. The schedule is precomputed
 * once by  compute_traveltime_schedule before the simulation starts.
 *
 *
 * @tparam FP       Floating-point type (e.g., double).
 * @tparam GraphT   Graph type (Graph<TravelTimeNodeProperty<FP,Sim>, TravelTimeEdge<FP>>).
 */
template <typename FP, class GraphT>
class GraphSimulationTravelTime
{
public:
    using Graph    = GraphT;
    using NodeProp = typename GraphT::NodeProperty;
    using EdgeProp = typename GraphT::EdgeProperty;

    /**
     * @brief Construct the simulation.
     *
     * @param t0      Simulation start time.
     * @param dt      Time step.
     * @param graph   The simulation graph.
     * @param n_steps Number of sub-steps per day (default 100).
     */
    GraphSimulationTravelTime(FP t0, FP dt, GraphT graph, size_t n_steps = 100)
        : m_t(t0)
        , m_dt(dt)
        , m_graph(std::move(graph))
        , m_n_steps(n_steps)
        , m_schedule(compute_traveltime_schedule<FP>(m_graph, n_steps))
    {
        // Reset mobility_sim compartments to zero as transit nodes always start empty.
        for (auto& n : m_graph.nodes()) {
            n.property.mobility_sim.get_result().get_last_value().setZero();
        }
    }

    FP get_t() const
    {
        return m_t;
    }

    GraphT& get_graph()
    {
        return m_graph;
    }

    const GraphT& get_graph() const
    {
        return m_graph;
    }

    /**
     * @brief Advance the simulation up to t_max.
     *
     * Processes each day in three phases per sub-step:
     * 1. Move commuters along active edges (init or transfer to next node).
     * 2. Advance local nodes up to the next breakpoint.
     * 3. Update commuter infection states in the transit nodes.
     *
     * At the end of each day, transit nodes are emptied (they should already be empty up to a certain tolerance)
     * and for longer simulations (more than interpolation_threshold (default 20) days) results are interpolated to restrict memory usage to one time point per day.
     *
     * @param t_max Maximum simulation time.
     * @param interpolation_threshold Time after which results are interpolated to daily resolution.
     */
    void advance(FP t_max, FP interpolation_threshold = FP{20})
    {
        const FP sub_dt = FP{1} / static_cast<FP>(m_n_steps);

        // Earliest commuter departure across all edges.
        FP t_first_departure = compute_first_departure();

        while (m_t < t_max - FP{1e-10}) {
            if (m_t + t_first_departure > t_max) {
                // If t_max is reached before the first commuter departure of the next day,
                // no mobility exchange occurs. Only advance local sims and stop.
                for (auto& n : m_graph.nodes()) {
                    n.property.local_sim.get_simulation().advance(t_max);
                }
                m_t = t_max;
                break;
            }

            for (size_t step = 0; step < m_n_steps; ++step) {
                const FP t_sub = m_t + static_cast<FP>(step) * sub_dt;

                // Phase 1: Move commuters along edges.
                move_commuters(step, t_sub);

                // Phase 2: Advance local simulations.
                advance_local_sims(step, t_sub, sub_dt);

                // Phase 3: Update infection states of commuters in transit.
                update_commuter_states(step, t_sub, sub_dt);
            }

            // End-of-day: return remaining transit individuals, reset mobility nodes.
            end_of_day(m_t + FP{1});

            m_t += FP{1};

            // Interpolate to daily resolution after interpolation_threshold days to limit memory usage.
            if (m_t > interpolation_threshold) {
                interpolate_results();
            }
        }
    }

private:
    FP m_t;
    FP m_dt;
    GraphT m_graph;
    size_t m_n_steps;
    TravelTimeSchedule m_schedule;

    FP compute_first_departure() const
    {
        FP earliest = FP{1};
        for (size_t ei = 0; ei < m_graph.edges().size(); ++ei) {
            const FP t_dep = static_cast<FP>(m_schedule.first_mobility_step[ei]) / m_n_steps;
            earliest       = std::min(earliest, t_dep);
        }
        return earliest;
    }

    void move_commuters(size_t step, FP t_sub)
    {
        for (size_t ei : m_schedule.edges_at_step[step]) {
            if (step == m_schedule.first_mobility_step[ei]) {
                // First mobility step: extract commuters from local_sim of origin.
                init_mobility(ei, t_sub);
            }
            else {
                // Subsequent steps: move commuters from previous to current transit node.
                transfer_commuters(ei, step, t_sub);
            }
        }
    }

    /// Extract commuters from local_sim at the origin node and place them in the
    /// first transit node's mobility_sim.
    void init_mobility(size_t ei, FP t_sub)
    {
        auto& e         = m_graph.edges()[ei];
        auto& edge      = e.property;
        auto& node_from = m_graph.nodes()[e.start_node_idx].property.local_sim;

        const auto& last_val = node_from.get_result().get_last_value();
        const auto& coeffs   = edge.parameters.get_coefficients().get_matrix_at(SimulationTime<FP>(t_sub));
        const auto factors   = get_mobility_factors<FP>(node_from, t_sub, last_val);

        Eigen::VectorX<FP> travel_subgroup = (last_val.array() * coeffs.array() * factors.array()).matrix();

        // Clamp negatives
        travel_subgroup = travel_subgroup.cwiseMax(FP{0});

        // Optional: test commuters
        test_commuters<FP>(node_from, travel_subgroup, t_sub);

        edge.mobile_population.add_time_point(t_sub, travel_subgroup);
        edge.return_times.add_time_point(t_sub);

        // Move from origin local_sim.
        Eigen::Ref<Eigen::VectorX<FP>> origin_val = node_from.get_result().get_last_value();
        origin_val -= travel_subgroup;

        // Place in first transit node's mobility_sim.
        const size_t first_transit = edge.path.front();
        auto& mob_sim              = m_graph.nodes()[first_transit].property.mobility_sim;
        mob_sim.get_result().get_last_value() += travel_subgroup;
    }

    /// Move commuters from the previous node to the current node in the schedule.
    void transfer_commuters(size_t ei, size_t step, FP /*t_sub*/)
    {
        auto& e    = m_graph.edges()[ei];
        auto& edge = e.property;

        if (edge.mobile_population.get_num_time_points() == 0)
            return;

        const size_t prev_node = m_schedule.node_at_step[ei][step - 1];
        const size_t curr_node = m_schedule.node_at_step[ei][step];
        const bool prev_mob    = m_schedule.in_mobility[ei][step - 1];
        const bool curr_mob    = m_schedule.in_mobility[ei][step];

        if (prev_node == curr_node && prev_mob == curr_mob)
            return;

        Eigen::Ref<Eigen::VectorX<FP>> travel_subgroup = edge.mobile_population.get_last_value();
        size_t num_age                                 = static_cast<size_t>(
            m_graph.nodes()[prev_node].property.local_sim.get_simulation().get_model().parameters.get_num_groups());
        correct_negative_compartments<FP>(travel_subgroup, num_age);

        // Remove from previous location.
        auto& sim_from =
            prev_mob ? m_graph.nodes()[prev_node].property.mobility_sim : m_graph.nodes()[prev_node].property.local_sim;
        sim_from.get_result().get_last_value() -= travel_subgroup;

        // Add to new location.
        auto& sim_to =
            curr_mob ? m_graph.nodes()[curr_node].property.mobility_sim : m_graph.nodes()[curr_node].property.local_sim;
        sim_to.get_result().get_last_value() += travel_subgroup;

        correct_negative_compartments<FP>(sim_from.get_result().get_last_value(), num_age);
        correct_negative_compartments<FP>(sim_to.get_result().get_last_value(), num_age);
    }

    void advance_local_sims(size_t step, FP t_sub, FP sub_dt)
    {
        for (size_t ni : m_schedule.local_nodes_at_step[step]) {
            // Determine the next breakpoint to find dt for this node.
            const FP node_dt = get_local_dt(ni, step, sub_dt);
            m_graph.nodes()[ni].property.local_sim.get_simulation().advance(t_sub + node_dt);
        }
    }

    FP get_local_dt(size_t ni, size_t step, FP sub_dt) const
    {
        const auto& bps = m_schedule.local_breakpoints[ni];
        const auto it   = std::lower_bound(bps.begin(), bps.end(), step);
        if (it == bps.end() || std::next(it) == bps.end()) {
            // Last breakpoint: advance to end of day.
            return (static_cast<FP>(m_n_steps) - static_cast<FP>(step)) * sub_dt;
        }
        return static_cast<FP>(*std::next(it) - step) * sub_dt;
    }

    void update_commuter_states(size_t step, FP t_sub, FP sub_dt)
    {
        for (size_t ei : m_schedule.edges_at_step[step]) {
            auto& e    = m_graph.edges()[ei];
            auto& edge = e.property;

            if (edge.mobile_population.get_num_time_points() == 0)
                continue;
            if (!m_schedule.in_mobility[ei][step])
                continue; // If not in transit, skip.

            const size_t curr_node = m_schedule.node_at_step[ei][step];
            auto& mob_node         = m_graph.nodes()[curr_node].property.mobility_sim;
            auto& model            = mob_node.get_simulation().get_model();

            // Auxiliary Euler step for commuter sub-population state estimation
            FP mob_dt                                      = get_mobility_dt(ei, step, sub_dt);
            Eigen::Ref<Eigen::VectorX<FP>> travel_subgroup = edge.mobile_population.get_last_value();
            const Eigen::VectorX<FP> total_mob             = mob_node.get_result().get_last_value();

            Eigen::VectorX<FP> y0 = travel_subgroup.eval();
            Eigen::VectorX<FP> y1 = Eigen::VectorX<FP>::Zero(y0.size());
            FP t_euler            = t_sub;
            auto deriv_fn = [&](Eigen::Ref<const Eigen::VectorX<FP>> y, FP /*t*/, Eigen::Ref<Eigen::VectorX<FP>> dydt) {
                model.get_derivatives(total_mob, y, t_sub, dydt);
            };
            DerivFunction<FP> f = deriv_fn;
            EulerIntegratorCore<FP>().step(f, y0, t_euler, mob_dt, y1);
            travel_subgroup = y1;

            // The auxiliary Euler heuristic can is prone to overshooting, especially when a subpopulation share is
            // high and the dynamics are fast. To prevent negative compartment values, we apply a correction step.
            size_t num_age = static_cast<size_t>(model.parameters.get_num_groups());
            correct_negative_compartments<FP>(travel_subgroup, num_age);
        }
    }

    FP get_mobility_dt(size_t ei, size_t step, FP sub_dt) const
    {
        const size_t curr_node = m_schedule.node_at_step[ei][step];
        const auto& bps        = m_schedule.mobility_breakpoints[curr_node];
        const auto it          = std::lower_bound(bps.begin(), bps.end(), step);
        if (it == bps.end() || std::next(it) == bps.end()) {
            return (static_cast<FP>(m_n_steps) - static_cast<FP>(step)) * sub_dt;
        }
        return static_cast<FP>(*std::next(it) - step) * sub_dt;
    }

    void end_of_day(FP t_end)
    {
        // Return all still-mobile commuters to their home local_sim.
        for (size_t ei = 0; ei < m_graph.edges().size(); ++ei) {
            auto& e    = m_graph.edges()[ei];
            auto& edge = e.property;

            if (edge.mobile_population.get_num_time_points() == 0)
                continue;

            // Final node in schedule is always the origin (home).
            const size_t home_node = e.start_node_idx;
            const size_t last_step = m_n_steps - 1;
            const bool last_mob    = m_schedule.in_mobility[ei][last_step];
            auto& sim_last = last_mob ? m_graph.nodes()[m_schedule.node_at_step[ei][last_step]].property.mobility_sim
                                      : m_graph.nodes()[m_schedule.node_at_step[ei][last_step]].property.local_sim;

            Eigen::Ref<Eigen::VectorX<FP>> travel_subgroup = edge.mobile_population.get_last_value();
            size_t num_age                                 = static_cast<size_t>(
                m_graph.nodes()[home_node].property.local_sim.get_simulation().get_model().parameters.get_num_groups());
            correct_negative_compartments<FP>(travel_subgroup, num_age);

            sim_last.get_result().get_last_value() -= travel_subgroup;
            m_graph.nodes()[home_node].property.local_sim.get_result().get_last_value() += travel_subgroup;

            // Clear edge commuter tracking.
            for (Eigen::Index i = edge.mobile_population.get_num_time_points() - 1; i >= 0; --i) {
                edge.mobile_population.remove_time_point(i);
                edge.return_times.remove_time_point(i);
            }
        }

        // Reset mobility_sim to zero: transit nodes start empty each day.
        for (auto& n : m_graph.nodes()) {
            n.property.mobility_sim.get_result().get_last_value().setZero();
        }

        mio::unused(t_end);
    }

    /// Keep only one time point per day in local_sim results (daily resolution).
    void interpolate_results()
    {
        for (auto& n : m_graph.nodes()) {
            auto& res = n.property.local_sim.get_simulation().get_result();
            if (res.get_num_time_points() > 2) {
                const FP t_last   = res.get_last_time();
                const auto v_last = res.get_last_value().eval();
                while (res.get_num_time_points() > 1) {
                    res.remove_last_time_point();
                }
                res.add_time_point(t_last, v_last);
            }
        }
    }
}; // class GraphSimulationTravelTime

/**
 * @brief Create a travel-time-aware graph simulation.
 *
 * @tparam FP       Floating-point type.
 * @tparam NodeProp TravelTimeNodeProperty specialisation.
 * @tparam EdgeProp TravelTimeEdge specialisation.
 * @param t0        Simulation start time.
 * @param dt        Day-level time step.
 * @param graph     Simulation graph.
 * @param n_steps   Number of day sub-steps (default 100).
 * @return          GraphSimulationTravelTime instance.
 */
template <typename FP, class NodeProp, class EdgeProp>
GraphSimulationTravelTime<FP, Graph<NodeProp, EdgeProp>>
make_traveltime_sim(FP t0, FP dt, Graph<NodeProp, EdgeProp> graph, size_t n_steps = 100)
{
    return GraphSimulationTravelTime<FP, Graph<NodeProp, EdgeProp>>(t0, dt, std::move(graph), n_steps);
}

} // namespace mio

#endif // MIO_MOBILITY_TRAVELTIME_H
