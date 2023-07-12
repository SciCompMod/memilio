/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Julia Bicker, Daniel Abele
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

#ifndef MIO_MOBILITY_HYBRID_GRAPH_SIMULATION_H
#define MIO_MOBILITY_HYBRID_GRAPH_SIMULATION_H

#include "memilio/mobility/hybrid_graph_simulation.h"
#include <functional>

namespace mio
{
template <class Timepoint, class Timespan, class HybridGraph>
class HybridGraphSimulation
{
    using node_function_abm = std::function<void(Timepoint, Timespan, typename HybridGraph::NodeProperty&)>;
    using node_function_ode = std::function<void(double, double, typename HybridGraph::ODEGraph::NodeProperty&)>;

    using edge_function_abm =
        std::function<void(Timepoint, Timespan, typename HybridGraph::ABMGraph::EdgeProperty&,
    typename HybridGraph::ABMGraph::NodeProperty&,
    typename HybridGraph::ABMGraph::NodeProperty&)>;
    using edge_function_ode = std::function<void(double, double, typename HybridGraph::ODEGraph::EdgeProperty&,
        typename HybridGraph::ODEGraph::NodeProperty&,
        typename HybridGraph::ODEGraph::NodeProperty&)>;
    using edge_function_hybrid =
        std::function<void(Timepoint, Timespan, typename HybridGraph::MigrationEdgeHybrid&,
    typename HybridGraph::ABMGraph::NodeProperty&,
    typename HybridGraph::ODEGraph::NodeProperty&)>;

    using time_conversion_function = std::function<double(Timepoint)>;
    HybridGraphSimulation(Timepoint t0, Timespan dt, const HybridGraph& g, const node_function_abm& node_func_abm,
                          const node_function_ode& node_func_ode, const edge_function_abm& edge_func_abm,
                          const edge_function_ode& edge_func_ode, const edge_function_hybrid& edge_func_hybrid,
                          const time_conversion_function& time_conversion_func)
        : m_t_abm(t0)
        , m_dt_abm(dt)
        , m_t_ode(time_conversion_func(t0))
        , m_dt_ode(time_conversion_func(Timepoint(0) + dt))
        , m_graph(g)
        , m_node_func_abm(node_func_abm)
        , m_node_func_ode(node_func_ode)
        , m_edge_func_abm(edge_func_abm)
        , m_edge_func_ode(edge_func_ode)
        , m_edge_func_hybrid(edge_func_hybrid)
    {
    }

    HybridGraphSimulation(Timepoint t0, Timespan dt, HybridGraph&& g, const node_function_abm& node_func_abm,
                          const node_function_ode& node_func_ode, const edge_function_abm& edge_func_abm,
                          const edge_function_ode& edge_func_ode, const edge_function_hybrid& edge_func_hybrid)
        : m_t_abm(t0)
        , m_dt_abm(dt)
        , m_t_ode(time_conversion_func(t0))
        , m_dt_ode(time_conversion_func(Timepoint(0) + dt))
        , m_graph(std::move(g))
        , m_node_func_abm(node_func_abm)
        , m_node_func_ode(node_func_ode)
        , m_edge_func_abm(edge_func_abm)
        , m_edge_func_ode(edge_func_ode)
        , m_edge_func_hybrid(edge_func_hybrid)
    {
    }

    // void advance(Timepoint tmax)
    // {
    //     auto dt = m_t;
    // }

private:
    Timepoint m_t_abm;
    Timespan m_dt_abm;
    double m_t_ode;
    double m_dt_ode;
    HybridGraph m_graph;
    node_function_abm m_node_func_abm;
    node_function_ode m_node_func_ode;
    edge_function_abm m_edge_func_abm;
    edge_function_ode m_edge_func_ode;
    edge_function_hybrid m_edge_func_hybrid;
};

} // namespace mio
#endif //MIO_MOBILITY_HYBRID_GRAPH_SIMULATION_H
