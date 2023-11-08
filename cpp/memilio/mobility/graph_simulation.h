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
#include <algorithm>
#include <cassert>
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

        // set population to zero in mobility nodes before starting
        for (auto& n : m_graph.nodes()) {
            n.node_pt.get_result().get_last_value().setZero();
        }

        auto min_dt    = 0.01;
        double t_begin = m_t - 1.;

        // calc schedule for each edge
        precompute_schedule();
        const double epsilon = std::numeric_limits<double>::epsilon();
        while (m_t - epsilon < t_max) {
            // auto start = std::chrono::system_clock::now();

            t_begin += 1;
            if (m_t + dt_first_mobility > t_max) {
                dt_first_mobility = t_max - m_t;
                for (auto& n : m_graph.nodes()) {
                    m_node_func(m_t, dt_first_mobility, n.property);
                }
                break;
            }

            for (auto& node : m_graph.nodes()) {
                node.node_pt.get_simulation().set_integrator(std::make_shared<mio::EulerIntegratorCore>());
            }

            size_t indx_schedule = 0;
            while (t_begin + 1 > m_t + 1e-10) {
                for (const auto& edge_indx : edges_mobility[indx_schedule]) {
                    auto& e = m_graph.edges()[edge_indx];
                    // first mobility activity
                    if (indx_schedule == first_mobility[edge_indx]) {
                        auto& node_from = m_graph.nodes()[schedule_edges[edge_indx][indx_schedule - 1]].property;
                        auto& node_to   = m_graph.nodes()[schedule_edges[edge_indx][indx_schedule]].node_pt;
                        m_edge_func(m_t, 0.0, e.property, node_from, node_to, 0);
                    }
                    // next mobility activity
                    else if (indx_schedule > first_mobility[edge_indx]) {
                        auto current_node_indx = schedule_edges[edge_indx][indx_schedule];
                        bool in_mobility_node  = mobility_schedule_edges[edge_indx][indx_schedule];

                        // determine dt, which is equal to the last integration/syncronization point in the current node
                        auto integrator_schedule_row = ln_int_schedule[current_node_indx];
                        if (in_mobility_node)
                            integrator_schedule_row = mb_int_schedule[current_node_indx];
                        // search the index of indx_schedule in the integrator schedule
                        const size_t indx_current =
                            std::distance(integrator_schedule_row.begin(),
                                          std::lower_bound(integrator_schedule_row.begin(),
                                                           integrator_schedule_row.end(), indx_schedule));

                        if (integrator_schedule_row[indx_current] != indx_schedule)
                            std::cout << "Error in schedule."
                                      << "\n";

                        ScalarType dt_mobility;
                        if (indx_current == 0) {
                            dt_mobility = round_second_decimal(e.traveltime / e.path.size());
                            if (dt_mobility < 0.01)
                                dt_mobility = 0.01;
                        }
                        else {
                            dt_mobility =
                                round_second_decimal((static_cast<double>(integrator_schedule_row[indx_current]) -
                                                      static_cast<double>(integrator_schedule_row[indx_current - 1])) /
                                                         100 +
                                                     epsilon);
                        }

                        // We have two cases. Either, we want to send the individuals to the next node, or we just want
                        // to update their state since a syncronization step is necessary in the current node.
                        if ((schedule_edges[edge_indx][indx_schedule] !=
                             schedule_edges[edge_indx][indx_schedule - 1]) ||
                            (mobility_schedule_edges[edge_indx][indx_schedule] !=
                             mobility_schedule_edges[edge_indx][indx_schedule - 1])) {
                            auto& node_from =
                                mobility_schedule_edges[edge_indx][indx_schedule - 1]
                                    ? m_graph.nodes()[schedule_edges[edge_indx][indx_schedule - 1]].node_pt
                                    : m_graph.nodes()[schedule_edges[edge_indx][indx_schedule - 1]].property;

                            auto& node_to = mobility_schedule_edges[edge_indx][indx_schedule]
                                                ? m_graph.nodes()[schedule_edges[edge_indx][indx_schedule]].node_pt
                                                : m_graph.nodes()[schedule_edges[edge_indx][indx_schedule]].property;

                            if (indx_schedule < mobility_schedule_edges[edge_indx].size() - 1) {
                                m_edge_func(m_t, dt_mobility, e.property, node_from, node_to, 1);
                            }
                            else
                                // the last time step is handled differently since we have to close the timeseries
                                m_edge_func(m_t, dt_mobility, e.property, node_from,
                                            m_graph.nodes()[schedule_edges[edge_indx][indx_schedule]].property, 2);
                        }
                        else {
                            auto& node_from =
                                mobility_schedule_edges[edge_indx][indx_schedule - 1]
                                    ? m_graph.nodes()[schedule_edges[edge_indx][indx_schedule - 1]].node_pt
                                    : m_graph.nodes()[schedule_edges[edge_indx][indx_schedule - 1]].property;

                            auto& node_to = mobility_schedule_edges[edge_indx][indx_schedule]
                                                ? m_graph.nodes()[schedule_edges[edge_indx][indx_schedule]].node_pt
                                                : m_graph.nodes()[schedule_edges[edge_indx][indx_schedule]].property;

                            assert(node_from.get_result().get_last_value() == node_to.get_result().get_last_value());
                            m_edge_func(m_t, dt_mobility, e.property, node_from, node_to, 3);
                        }
                    }
                }

                // first we integrate the nodes in time. Afterwards the update on the edges is done.
                // We start with the edges since the values for t0 are given.
                // iterate over all local nodes and integrate them to the syncronization point
                for (const auto& n_indx : nodes_mobility[indx_schedule]) {
                    auto& n = m_graph.nodes()[n_indx];
                    const size_t indx_current =
                        std::distance(ln_int_schedule[n_indx].begin(),
                                      std::lower_bound(ln_int_schedule[n_indx].begin(), ln_int_schedule[n_indx].end(),
                                                       indx_schedule));

                    const size_t val_next = (indx_current == ln_int_schedule[n_indx].size() - 1)
                                                ? 100
                                                : ln_int_schedule[n_indx][indx_current + 1];
                    const ScalarType next_dt =
                        round_second_decimal((static_cast<double>(val_next) - indx_schedule) / 100 + epsilon);
                    m_node_func(m_t, next_dt, n.property);
                }

                // iterate over all mobility nodes and integrate them to the syncronization point
                for (const size_t& n_indx : nodes_mobility_m[indx_schedule]) {
                    auto& n = m_graph.nodes()[n_indx];
                    const size_t indx_current =
                        std::distance(mb_int_schedule[n_indx].begin(),
                                      std::lower_bound(mb_int_schedule[n_indx].begin(), mb_int_schedule[n_indx].end(),
                                                       indx_schedule));
                    const size_t val_next = (indx_current == mb_int_schedule[n_indx].size() - 1)
                                                ? 100
                                                : mb_int_schedule[n_indx][indx_current + 1];
                    const ScalarType next_dt =
                        round_second_decimal((static_cast<double>(val_next) - indx_schedule) / 100 + epsilon);

                    // get all time points from the last integration step
                    auto& last_time_point =
                        n.node_pt.get_result().get_time(n.node_pt.get_result().get_num_time_points() - 1);
                    // wenn last_time_point nicht innerhalb eines intervalls von 1-e10 von t liegt, dann setzte den letzten Zeitpunkt auf m_t
                    if (std::fabs(last_time_point - m_t) > 1e-10) {
                        n.node_pt.get_result().get_last_time() = m_t;
                    }
                    std::cout << "population in mobility node with id " << n.id << " at time " << m_t << " is "
                              << n.node_pt.get_result().get_last_value().sum() << "\n";
                    Eigen::Index indx_min;
                    if (n.node_pt.get_result().get_last_value().minCoeff(&indx_min) < 0) {
                        std::cout << "value schon vorher negativ " << n.node_pt.get_result().get_last_value()[indx_min]
                                  << "\n";
                    }
                    if (n.id == 7233 && m_t == 0.16) {
                        std::cout << "ja"
                                  << "\n";
                    }
                    std::cout << "Before num time points " << n.node_pt.get_result().get_num_time_points() << "\n";
                    m_node_func(m_t, next_dt, n.node_pt);
                    while (n.node_pt.get_result().get_last_value().minCoeff() < 0) {
                        std::cout << "id = " << n.id << "\n";
                        Eigen::Index indx_max;
                        n.node_pt.get_result().get_last_value().maxCoeff(&indx_max);
                        n.node_pt.get_result().get_last_value()[indx_max] -=
                            n.node_pt.get_result().get_last_value()[indx_min];
                        n.node_pt.get_result().get_last_value()[indx_min] = 0;
                        std::cout << "After num time points " << n.node_pt.get_result().get_num_time_points() << "\n";

                        std::cout << "second last results  = "
                                  << n.node_pt.get_result().get_value(n.node_pt.get_result().get_num_time_points() - 2)
                                  << "\n";
                        std::cout << "last results  = " << n.node_pt.get_result().get_last_value() << "\n";
                        std::cout
                            << "sum second last results  = "
                            << n.node_pt.get_result().get_value(n.node_pt.get_result().get_num_time_points() - 2).sum()
                            << "\n";
                        std::cout << "sum last results  = " << n.node_pt.get_result().get_last_value().sum() << "\n";
                        std::exit(1);
                    }
                }
                indx_schedule++;
                m_t += min_dt;
            }
            // set each compartment zero for all mobility nodes since we only model daily mobility
            for (auto& n : m_graph.nodes()) {
                n.node_pt.get_result().get_last_value().setZero();
            }
            // std::cout << "aktuell bei t = " << m_t << "\n";
            // // messe die zeit, wie lange eine iteration bis zu dieser stelle dauert
            // auto end                                      = std::chrono::system_clock::now();
            // std::chrono::duration<double> elapsed_seconds = end - start;
            // std::cout << "Time (min) needed per Iteration is " << elapsed_seconds.count() / 60 << "min\n";
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
    // describes the schedule for each edge, i.e. which node is visited at which time step
    std::vector<std::vector<size_t>> schedule_edges;

    // defines if people from a edge are currently in a mobility node. This is necessary to determine the correct
    // position in the schedule. Otherwise we could add a specific identifier in schedule_edges.
    std::vector<std::vector<bool>> mobility_schedule_edges;

    // describes the syncronization steps which are necessary for the integrator in the edges and nodes.
    std::vector<std::vector<size_t>> mb_int_schedule;
    std::vector<std::vector<size_t>> ln_int_schedule;

    // defines the edges on which mobility takes place for each time step
    std::vector<std::vector<size_t>> edges_mobility;

    // same for the nodes
    std::vector<std::vector<size_t>> nodes_mobility;
    std::vector<std::vector<size_t>> nodes_mobility_m;

    // first mobility activites per edge
    std::vector<size_t> first_mobility;

    void precompute_schedule()
    {
        const size_t timesteps = 100;
        schedule_edges.reserve(m_graph.edges().size());
        mobility_schedule_edges.reserve(m_graph.edges().size());

        const double epsilon = std::numeric_limits<double>::epsilon();

        for (auto& e : m_graph.edges()) {
            // 100 since we round to second decimal
            std::vector<size_t> schedule(timesteps, 0.);
            std::vector<bool> is_mobility_node(timesteps, false);

            double traveltime_per_region = round_second_decimal(e.traveltime / e.path.size());

            if (traveltime_per_region < 0.01)
                traveltime_per_region = 0.01;

            double start_mobility = round_second_decimal(0 + 1 - 2 * (traveltime_per_region * e.path.size()) -
                                                         m_graph.nodes()[e.end_node_idx].stay_duration);
            if (start_mobility < 0.0) {
                start_mobility = 0.;
            }

            double arrive_at = start_mobility + traveltime_per_region * e.path.size();

            std::fill(schedule.begin(), schedule.begin() + static_cast<int>((start_mobility + epsilon) * 100),
                      e.start_node_idx);

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

        // iterate over nodes
        size_t indx_node = 0;
        for (auto& n : m_graph.nodes()) {
            // local node and automatical add zero, since we want to start at t=0 and start with integrating all nodes to the next dt
            std::vector<size_t> order_local_node = {0};
            std::vector<size_t> indx_edges;
            for (size_t indx_edge = 0; indx_edge < schedule_edges.size(); ++indx_edge) {
                if (schedule_edges[indx_edge][0] == indx_node && !mobility_schedule_edges[indx_edge][0]) {
                    indx_edges.push_back(indx_edge);
                }
            }

            for (size_t indx_t = 1; indx_t < timesteps; ++indx_t) {
                std::vector<size_t> indx_edges_new;
                for (size_t indx_edge = 0; indx_edge < schedule_edges.size(); ++indx_edge) {
                    if (schedule_edges[indx_edge][indx_t] == indx_node && !mobility_schedule_edges[indx_edge][indx_t]) {
                        indx_edges_new.push_back(indx_edge);
                    }
                }

                if (indx_edges_new.size() != indx_edges.size() ||
                    !std::equal(indx_edges.begin(), indx_edges.end(), indx_edges_new.begin())) {
                    order_local_node.push_back(indx_t);
                    indx_edges = indx_edges_new;
                }
            }
            if (order_local_node[order_local_node.size() - 1] != timesteps - 1)
                order_local_node.push_back(timesteps - 1);
            ln_int_schedule.push_back(order_local_node);

            // mobility node
            std::vector<size_t> order_mobility_node;
            std::vector<size_t> indx_edges_mobility;
            for (size_t indx_edge = 0; indx_edge < schedule_edges.size(); ++indx_edge) {
                if (schedule_edges[indx_edge][0] == indx_node && mobility_schedule_edges[indx_edge][0]) {
                    indx_edges_mobility.push_back(indx_edge);
                }
            }
            for (size_t indx_t = 1; indx_t < timesteps; ++indx_t) {
                std::vector<size_t> indx_edges_mobility_new;
                for (size_t indx_edge = 0; indx_edge < schedule_edges.size(); ++indx_edge) {
                    if (schedule_edges[indx_edge][indx_t] == indx_node && mobility_schedule_edges[indx_edge][indx_t]) {
                        indx_edges_mobility_new.push_back(indx_edge);
                    }
                }

                if (indx_edges_mobility_new.size() != indx_edges_mobility.size() ||
                    !std::equal(indx_edges_mobility.begin(), indx_edges_mobility.end(),
                                indx_edges_mobility_new.begin())) {
                    order_mobility_node.push_back(indx_t);
                    indx_edges_mobility = indx_edges_mobility_new;
                }
            }
            mb_int_schedule.push_back(order_mobility_node);
            indx_node++;
        }

        auto indx_edge = 0;
        first_mobility.reserve(m_graph.edges().size());
        for (auto& e : m_graph.edges()) {
            this->first_mobility[indx_edge] = std::distance(
                mobility_schedule_edges[indx_edge].begin(),
                std::find(mobility_schedule_edges[indx_edge].begin(), mobility_schedule_edges[indx_edge].end(), true));
            indx_edge++;
        }

        edges_mobility.reserve(timesteps);
        nodes_mobility.reserve(timesteps);
        nodes_mobility_m.reserve(timesteps);

        // we handle indx_current = 0 separately since we want have added them always in the
        // intregration schedule so this would lead to wrong results
        // here we iterate over first mobility activities and add the edges indx when there is a start at t = 0
        std::vector<size_t> temp_edges_mobility;
        indx_edge = 0;
        for (auto& start_time : first_mobility) {
            if (start_time == 0) {
                temp_edges_mobility.push_back(indx_edge);
            }
            indx_edge++;
        }
        edges_mobility.push_back(std::move(temp_edges_mobility));

        // same for nodes
        std::vector<size_t> temp_nodes_mobility(m_graph.nodes().size());
        std::iota(temp_nodes_mobility.begin(), temp_nodes_mobility.end(), 0);
        nodes_mobility.emplace_back(std::move(temp_nodes_mobility));

        std::vector<size_t> temp_nodes_mobility_m;
        size_t node_indx = 0;
        for (auto& n : m_graph.nodes()) {
            if (std::binary_search(mb_int_schedule[node_indx].begin(), mb_int_schedule[node_indx].end(), 0)) {
                temp_nodes_mobility_m.push_back(node_indx);
                node_indx++;
            }
        }
        nodes_mobility_m.push_back(temp_nodes_mobility_m);

        for (int indx_current = 1; indx_current < timesteps; ++indx_current) {
            std::vector<size_t> temp_edge_mobility;
            indx_edge = 0;
            for (auto& e : m_graph.edges()) {
                auto current_node_indx = schedule_edges[indx_edge][indx_current];
                if (indx_current >= first_mobility[indx_edge]) {
                    if (mobility_schedule_edges[indx_edge][indx_current] &&
                        std::binary_search(mb_int_schedule[current_node_indx].begin(),
                                           mb_int_schedule[current_node_indx].end(), indx_current))
                        temp_edge_mobility.push_back(indx_edge);
                    else if (!mobility_schedule_edges[indx_edge][indx_current] &&
                             std::binary_search(ln_int_schedule[current_node_indx].begin(),
                                                ln_int_schedule[current_node_indx].end(), indx_current))
                        temp_edge_mobility.push_back(indx_edge);
                }
                indx_edge++;
            }
            edges_mobility.push_back(temp_edge_mobility);

            // reset temp_nodes_mobility
            temp_nodes_mobility.clear();
            temp_nodes_mobility_m.clear();
            node_indx = 0;
            for (auto& n : m_graph.nodes()) {
                if (std::binary_search(ln_int_schedule[node_indx].begin(), ln_int_schedule[node_indx].end(),
                                       indx_current)) {
                    temp_nodes_mobility.push_back(node_indx);
                }

                if (std::binary_search(mb_int_schedule[node_indx].begin(), mb_int_schedule[node_indx].end(),
                                       indx_current)) {
                    temp_nodes_mobility_m.push_back(node_indx);
                }
                node_indx++;
            }
            nodes_mobility.push_back(temp_nodes_mobility);
            nodes_mobility_m.push_back(temp_nodes_mobility_m);
        }

        // Finally, we want to count the number of interactions, when there a at least two edges within the same mobility node
        // Count the number of intersections for each time points
        // std::vector<size_t> intersections_mobility(timesteps);
        // for (size_t indx_current = 0; indx_current < timesteps; ++indx_current) {
        //     std::vector<size_t> num_groups_in_mobility_node(m_graph.nodes().size());
        //     size_t inndx_edge = 0;
        //     for (auto& e : m_graph.nodes()) {
        //         if (mobility_schedule_edges[inndx_edge][indx_current]) {
        //             auto curr_node = schedule_edges[inndx_edge][indx_current];
        //             num_groups_in_mobility_node[curr_node]++;
        //         }
        //     }
        //     // zähle die anzahl an einträgen größer 2 in num_groups_in_mobility_node
        //     intersections_mobility[indx_current] =
        //         std::count_if(num_groups_in_mobility_node.begin(), num_groups_in_mobility_node.end(), [](size_t i) {
        //             return i > 1;
        //         });
        //     num_groups_in_mobility_node.clear();

        //     // TODO: SPEICHERE das in eine TXT datei im schönen Format
        // }
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
        ScalarType normalized_waiting_time = ExponentialDistribution<ScalarType>::get_instance()(m_rng, 1.0);
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
                    size_t event = mio::DiscreteDistribution<size_t>::get_instance()(m_rng, m_rates);
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
                    normalized_waiting_time = ExponentialDistribution<ScalarType>::get_instance()(m_rng, 1.0);

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

    RandomNumberGenerator& get_rng()
    {
        return m_rng;
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
    RandomNumberGenerator m_rng;
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
