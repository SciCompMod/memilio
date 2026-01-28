/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Kilian Volmer
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
#ifndef METAPOPULATION_MOBILITY_JOLLY_H
#define METAPOPULATION_MOBILITY_JOLLY_H

#include "graph_simulation.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_stochastic.h"
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
#include "memilio/config.h"

#include "boost/filesystem.hpp"
#include "metapopulation_mobility_instant.h"

#include <cassert>

namespace mio
{
/**
 * represents the mobility between two nodes.
 */
template <typename FP>
class JollyEdge
{
public:
    /**
     * create edge with coefficients.
     * @param coeffs mobility rate for each group and compartment
     */
    JollyEdge(const MobilityParametersStochastic<FP>& params)
        : m_parameters(params)
    {
    }

    /**
     * create edge with coefficients.
     * @param coeffs mobility rate for each group and compartment
     */
    JollyEdge(const Eigen::VectorX<FP>& coeffs)
        : m_parameters(coeffs)
    {
    }

    /**
     * get the mobility parameters.
     */
    const MobilityParametersStochastic<FP>& get_parameters() const
    {
        return m_parameters;
    }

    /**
     * get the cumulative transition rate of the edge.
    */
    template <class Sim>
    Eigen::VectorX<FP> get_transition_rates(SimulationNode<FP, Sim>& node_from)
    {
        Eigen::VectorX<FP> transitionRates(node_from.get_last_state().size());
        for (Eigen::Index i = 0; i < node_from.get_last_state().size(); ++i) {
            transitionRates[i] =
                node_from.get_last_state()(i) * m_parameters.get_coefficients().get_baseline()[(size_t)i];
        }
        return transitionRates;
    }

    /**
     * compute mobility from node_from to node_to for a given event
     * @param[in] event index specifying which compartment and age group change nodes
     * @param node_from node that people changed from
     * @param node_to node that people changed to
     */
    template <class Sim>
    void apply_mobility(size_t event, FP batch_size, SimulationNode<FP, Sim>& node_from,
                        SimulationNode<FP, Sim>& node_to);

private:
    MobilityParametersStochastic<FP> m_parameters;
};

template <typename FP>
template <class Sim>
void JollyEdge<FP>::apply_mobility(size_t event, FP batch_size, SimulationNode<FP, Sim>& node_from,
                                   SimulationNode<FP, Sim>& node_to)
{
    if (batch_size > node_from.get_result().get_last_value()[event]) {
        batch_size = node_from.get_result().get_last_value()[event];
    }
    node_from.get_result().get_last_value()[event] -= batch_size;
    node_to.get_result().get_last_value()[event] += batch_size;
}

/**
 * edge functor for mobility-based simulation.
 * @see JollyEdge::apply_mobility
 */
template <typename FP, class Sim, class StochasticEdge>
void apply_mobility(size_t event, FP batch_size, StochasticEdge& jollyEdge, SimulationNode<FP, Sim>& node_from,
                    SimulationNode<FP, Sim>& node_to)
{
    jollyEdge.apply_mobility(event, batch_size, node_from, node_to);
}

/**
 * create a mobility-based simulation.
 * After every second time step, for each edge a portion of the population corresponding to the coefficients of the edge
 * changes from one node to the other. In the next timestep, the mobile population returns to their "home" node.
 * Returns are adjusted based on the development in the target node.
 * @param t0 start time of the simulation
 * @param dt time step between mobility
 * @param graph set up for mobility-based simulation
 * @{
 */
template <typename FP, class Sim>
JollyGraphSimulation<FP, Graph<SimulationNode<FP, Sim>, JollyEdge<FP>>>
make_jolly_sim(FP t0, FP dt, const Graph<SimulationNode<FP, Sim>, JollyEdge<FP>>& graph)
{
    return make_graph_sim_jolly<FP>(
        t0, dt, graph, &advance_model<FP, Sim>,
        static_cast<void (*)(size_t, FP, JollyEdge<FP>&, SimulationNode<FP, Sim>&, SimulationNode<FP, Sim>&)>(
            &apply_mobility<FP, Sim, JollyEdge<FP>>));
}

template <typename FP, class Sim>
JollyGraphSimulation<FP, Graph<SimulationNode<FP, Sim>, JollyEdge<FP>>>
make_jolly_sim(FP t0, FP dt, Graph<SimulationNode<FP, Sim>, JollyEdge<FP>>&& graph)
{
    return make_graph_sim_jolly<FP>(
        t0, dt, std::move(graph), &advance_model<FP, Sim>,
        static_cast<void (*)(size_t, FP, JollyEdge<FP>&, SimulationNode<FP, Sim>&, SimulationNode<FP, Sim>&)>(
            &apply_mobility<FP, Sim, JollyEdge<FP>>));
}

/** @} */

} // namespace mio

#endif //METAPOPULATION_MOBILITY_JOLLY_H
