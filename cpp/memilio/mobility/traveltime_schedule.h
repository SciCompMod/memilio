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
#ifndef MIO_MOBILITY_TRAVELTIME_SCHEDULE_H
#define MIO_MOBILITY_TRAVELTIME_SCHEDULE_H

#include "memilio/math/floating_point.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <vector>

namespace mio
{

/**
 * @brief Precomputed mobility schedule for the travel-time-aware metapopulation model.
 *
 * Discretises one normalized time period [0, 1) into n_steps equal slices of width
 * `1/n_steps`, and records for every graph edge and every sub-step:
 * - which node index currently holds the edge's mobile sub-population, and
 * - whether that sub-population is in a `mobility_sim` (in transit) or a
 *   `local_sim` (at home or at the destination).
 *
 * Derived per-node breakpoint lists tell the ODE integrator when it must stop and
 * synchronise state with the commuter exchange:
 * - `local_breakpoints(ni)`:    sub-steps at which `local_sim` of node ni must synchronise.
 * - `mobility_breakpoints(ni)`: sub-steps at which `mobility_sim` of node ni must synchronise.
 *
 * Per-step tables to accelerate the inner advance loop:
 * - `edges_at(step)`:          edges whose commuters change node or mode at this step.
 * - `local_nodes_at(step)`:    nodes whose `local_sim` breakpoint falls at this step.
 * - `mobility_nodes_at(step)`: nodes whose `mobility_sim` breakpoint falls at this step.
 *
 * The schedule is computed once and then reused throughout the simulation.  
 * Use the constructor to build it and all data is read-only afterwards.
 */
class TravelTimeSchedule
{
public:
    /**
     * @brief Compute the schedule from a graph.
     *
     * @tparam FP     Floating-point type of the graph's node/edge properties.
     * @tparam GraphT Graph type, expected to be
     *                `Graph<TravelTimeNodeProperty<FP,Sim>, TravelTimeEdge<FP>>`.
     * @param graph   The simulation graph.
     * @param n_steps Number of sub-steps per normalized time period (default: 120).
     *                120 corresponds to 12-minute intervals when the period is one day
     *                and is evenly divisible by common commute-time fractions such as
     *                0.05, 0.1, and 0.25, so typical trip durations map to exact integer
     *                step counts without rounding.
     * @param eps     Floating-point tolerance used when converting continuous time
     *                fractions to step indices to avoid off-by-one errors due to
     *                rounding.
     */
    template <typename FP, class GraphT>
    TravelTimeSchedule(const GraphT& graph, size_t n_steps = 120, FP eps = Limits<FP>::zero_tolerance());

    /**
     * @brief Node index at which edge mobile sub-population are at sub-step step.
     *
     * Before mobility starts (step < first_mobility_step(ei)) the edge's commuters are
     * implicitly at the origin. Calling this function for those steps still returns the
     * source node.
     */
    size_t node_at(size_t ei, size_t step) const
    {
        assert(ei < m_n_edges && step < m_n_steps);
        return m_node_at_step[ei * m_n_steps + step];
    }

    /**
     * @brief Whether edge commuters are in a `mobility_sim` at sub-step step.
     *
     * Returns `true` during in-transit phases (outbound and return travel) and `false`
     * during the stay at the destination or while still at the origin.
     */
    bool in_mobility_at(size_t ei, size_t step) const
    {
        assert(ei < m_n_edges && step < m_n_steps);
        return m_in_mobility[ei * m_n_steps + step];
    }

    /**
     * @brief First sub-step at which edge commuters begin their outbound trip.
     *
     * All sub-steps before this index belong to the local phase at the origin.
     * Equals n_steps if the edge is never active (e.g. zero travel time or zero
     * mobility coefficient).
     */
    size_t first_mobility_step(size_t ei) const
    {
        assert(ei < m_n_edges);
        return m_first_mobility_step[ei];
    }

    /**
     * @brief Sorted list of sub-steps at which node local_sim must synchronize.
     *
     * A breakpoint is inserted whenever the set of edges whose commuters are present in
     * the local population of this node changes.  The integrator must stop at each
     * breakpoint to add or remove commuters before continuing.
     */
    const std::vector<size_t>& local_breakpoints(size_t ni) const
    {
        assert(ni < m_n_nodes);
        return m_local_breakpoints[ni];
    }

    /**
     * @brief Sorted list of sub-steps at which node mobility_sim must synchronize.
     *
     * Analogous to `local_breakpoints` but for the transit population simulation.
     */
    const std::vector<size_t>& mobility_breakpoints(size_t ni) const
    {
        assert(ni < m_n_nodes);
        return m_mobility_breakpoints[ni];
    }

    /**
     * @brief Edge indices whose commuters change their current node or mobility mode at
     *        sub-step step (including the first departure step).
     *
     * Used by the advance loop to efficiently process only the edges that require a
     * commuter transfer at each step.
     */
    const std::vector<size_t>& edges_at(size_t step) const
    {
        assert(step < m_n_steps);
        return m_edges_at_step[step];
    }

    /**
     * @brief Indices of nodes whose `local_sim` breakpoint falls at sub-step step.
     *
     * Allows the advance loop to advance only the affected nodes at each sub-step.
     */
    const std::vector<size_t>& local_nodes_at(size_t step) const
    {
        assert(step < m_n_steps);
        return m_local_nodes_at_step[step];
    }

    /**
     * @brief Indices of nodes whose `mobility_sim` breakpoint falls at sub-step step.
     */
    const std::vector<size_t>& mobility_nodes_at(size_t step) const
    {
        assert(step < m_n_steps);
        return m_mobility_nodes_at_step[step];
    }

    size_t num_steps() const
    {
        return m_n_steps;
    }
    size_t num_edges() const
    {
        return m_n_edges;
    }
    size_t num_nodes() const
    {
        return m_n_nodes;
    }

private:
    size_t m_n_edges;
    size_t m_n_nodes;
    size_t m_n_steps;

    /// Flat per-edge/per-step node index: `m_node_at_step[ei * n_steps + step]` returns the node index at which
    /// the edge's commuters are located at this step.
    std::vector<size_t> m_node_at_step;

    /// Flat per-edge/per-step in-transit flag: `m_in_mobility[ei * n_steps + step]` returns whether
    /// the edge's commuters are in transit at this step.
    std::vector<bool> m_in_mobility;

    std::vector<size_t> m_first_mobility_step; ///< first step with commuter movement

    std::vector<std::vector<size_t>> m_local_breakpoints; ///< sorted step indices
    std::vector<std::vector<size_t>> m_mobility_breakpoints; ///< sorted step indices

    std::vector<std::vector<size_t>> m_edges_at_step; ///< edge indices at each step
    std::vector<std::vector<size_t>> m_local_nodes_at_step; ///< node indices at each step
    std::vector<std::vector<size_t>> m_mobility_nodes_at_step; ///< node indices at each step
};

template <typename FP, class GraphT>
TravelTimeSchedule::TravelTimeSchedule(const GraphT& graph, size_t n_steps, FP eps)
    : m_n_edges(graph.edges().size())
    , m_n_nodes(graph.nodes().size())
    , m_n_steps(n_steps)
    , m_node_at_step(m_n_edges * n_steps, 0)
    , m_in_mobility(m_n_edges * n_steps, false)
    , m_first_mobility_step(m_n_edges, n_steps)
    , m_local_breakpoints(m_n_nodes)
    , m_mobility_breakpoints(m_n_nodes)
    , m_edges_at_step(n_steps)
    , m_local_nodes_at_step(n_steps)
    , m_mobility_nodes_at_step(n_steps)
{
    const double dt_step = 1.0 / static_cast<double>(n_steps);
    const double eps_d   = static_cast<double>(eps);

    auto all_edges = graph.edges();

    // 1. Build per-edge schedules: initialise default node to source, then overwrite the mobility and stay phases
    for (size_t ei = 0; ei < m_n_edges; ++ei) {
        const size_t src = all_edges[ei].start_node_idx;
        std::fill(m_node_at_step.begin() + static_cast<std::ptrdiff_t>(ei * n_steps),
                  m_node_at_step.begin() + static_cast<std::ptrdiff_t>((ei + 1) * n_steps), src);
    }

    for (size_t ei = 0; ei < m_n_edges; ++ei) {
        const auto& e          = all_edges[ei];
        const auto& edge       = e.property;
        const size_t dest_node = e.end_node_idx;
        const size_t n_patches = edge.path.size(); // includes origin and destination nodes

        // Travel time per path patch, snapped to the nearest step boundary.
        const double tt_per_patch = std::max(
            dt_step,
            std::round(static_cast<double>(edge.travel_time) / static_cast<double>(n_patches) / dt_step) * dt_step);
        const double total_travel = tt_per_patch * static_cast<double>(n_patches);

        // The normalized time period is partitioned as:
        //   [0, t_depart)              local at origin
        //   [t_depart, t_arrive)       outbound transit
        //   [t_arrive, t_leave_dst)    stay at destination
        //   [t_leave_dst, 1)           return transit
        const double t_depart = std::max(0.0, 1.0 - 2.0 * total_travel -
                                                  static_cast<double>(graph.nodes()[dest_node].property.stay_duration));

        const size_t step_depart    = static_cast<size_t>((t_depart + eps_d) / dt_step);
        const size_t steps_travel   = static_cast<size_t>((total_travel + eps_d) / dt_step);
        const size_t step_arrive    = step_depart + steps_travel;
        const size_t step_leave_dst = n_steps - steps_travel;

        m_first_mobility_step[ei] = step_depart;

        // Outbound travel: step_depart … step_arrive-1
        size_t current_step = step_depart;
        for (size_t patch = 0; patch < n_patches; ++patch) {
            const size_t patch_steps = static_cast<size_t>((tt_per_patch + eps_d) / dt_step);
            for (size_t s = 0; s < patch_steps && (current_step + s) < n_steps; ++s) {
                m_node_at_step[ei * n_steps + current_step + s] = edge.path[patch];
                m_in_mobility[ei * n_steps + current_step + s]  = true;
            }
            current_step += patch_steps;
        }

        // Stay at destination: step_arrive … step_leave_dst-1
        for (size_t s = step_arrive; s < step_leave_dst && s < n_steps; ++s) {
            m_node_at_step[ei * n_steps + s] = dest_node;
            m_in_mobility[ei * n_steps + s]  = false;
        }

        // Return travel: step_leave_dst … n_steps-1 (reversed path)
        current_step = step_leave_dst;
        for (size_t patch = n_patches; patch > 0 && current_step < n_steps; --patch) {
            const size_t patch_steps = static_cast<size_t>((tt_per_patch + eps_d) / dt_step);
            for (size_t s = 0; s < patch_steps && (current_step + s) < n_steps; ++s) {
                m_node_at_step[ei * n_steps + current_step + s] = edge.path[patch - 1];
                m_in_mobility[ei * n_steps + current_step + s]  = true;
            }
            current_step += patch_steps;
        }
    }

    // 2. Derive per-node breakpoint schedules
    for (size_t ni = 0; ni < m_n_nodes; ++ni) {
        m_local_breakpoints[ni].push_back(0);
        m_mobility_breakpoints[ni].push_back(0);

        auto prev_local_edges = std::vector<size_t>{};
        auto prev_mob_edges   = std::vector<size_t>{};
        for (size_t t = 1; t < n_steps; ++t) {
            std::vector<size_t> local_e, mob_e;
            for (size_t ei = 0; ei < m_n_edges; ++ei) {
                if (m_node_at_step[ei * n_steps + t] == ni) {
                    if (m_in_mobility[ei * n_steps + t])
                        mob_e.push_back(ei);
                    else
                        local_e.push_back(ei);
                }
            }
            if (local_e != prev_local_edges) {
                m_local_breakpoints[ni].push_back(t);
                prev_local_edges = local_e;
            }
            if (mob_e != prev_mob_edges) {
                m_mobility_breakpoints[ni].push_back(t);
                prev_mob_edges = mob_e;
            }
        }
        // Ensure the final sub-step is always present as a breakpoint.
        if (m_local_breakpoints[ni].back() != n_steps - 1)
            m_local_breakpoints[ni].push_back(n_steps - 1);
        if (m_mobility_breakpoints[ni].back() != n_steps - 1)
            m_mobility_breakpoints[ni].push_back(n_steps - 1);
    }

    // 3. Build per-step lookup tables
    for (size_t t = 0; t < n_steps; ++t) {
        // Edges whose commuters start moving or change node/mode at this step.
        for (size_t ei = 0; ei < m_n_edges; ++ei) {
            if (t < m_first_mobility_step[ei])
                continue;
            bool add = (t == m_first_mobility_step[ei]);
            if (!add && t > 0) {
                add = (m_node_at_step[ei * n_steps + t] != m_node_at_step[ei * n_steps + t - 1]) ||
                      (m_in_mobility[ei * n_steps + t] != m_in_mobility[ei * n_steps + t - 1]);
            }
            if (add)
                m_edges_at_step[t].push_back(ei);
        }
        // Nodes with a breakpoint at this step.
        for (size_t ni = 0; ni < m_n_nodes; ++ni) {
            if (std::binary_search(m_local_breakpoints[ni].begin(), m_local_breakpoints[ni].end(), t))
                m_local_nodes_at_step[t].push_back(ni);
            if (std::binary_search(m_mobility_breakpoints[ni].begin(), m_mobility_breakpoints[ni].end(), t))
                m_mobility_nodes_at_step[t].push_back(ni);
        }
    }
}

} // namespace mio

#endif // MIO_MOBILITY_TRAVELTIME_SCHEDULE_H
