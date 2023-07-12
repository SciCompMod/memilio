/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Julia Bicker, Daniel Abele, Martin J. Kuehn
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
#ifndef HYBRID_GRAPH_H
#define HYBRID_GRAPH_H

#include "memilio/mobility/graph.h"

namespace mio
{

/**
 * @brief generic hybrid graph structure
 */
template <class ABMGraph, class ODEGraph, class ABM, class MigrationEdgeHybrid, class SetupABMFromOdeFct>
class HybridGraph
{
public:
    HybridGraph(const ABMGraph& abm_graph, const ODEGraph& ode_graph, const SetupABMFromOdeFct& setup_func,
                std::vector<std::vector<double>> household_distributions = {})
        : m_abm_graph(abm_graph)
        , m_ode_graph(ode_graph)
        , m_abm_nodes_to_ode_nodes(setup_func(household_distributions))
    {
    }

    HybridGraph(ABMGraph&& abm_graph, ODEGraph&& ode_graph, const SetupABMFromOdeFct& setup_func,
                std::vector<std::vector<double>> household_distributions = {})
        : m_abm_graph(std::move(abm_graph))
        , m_ode_graph(std::move(ode_graph))
        , m_abm_nodes_to_ode_nodes(setup_func(household_distributions))
    {
    }

private:
    ABMGraph m_abm_graph;
    ODEGraph m_ode_graph;
    std::vector<Node<ABM>> m_abm_nodes_to_ode_nodes;
    std::vector<Edge<MigrationEdgeHybrid>> m_hybrid_edges;
};

} // namespace mio

#endif //HYBRID_GRAPH_H
