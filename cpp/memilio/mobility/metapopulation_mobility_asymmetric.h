/* 
* Copyright (C) 2020-2025 MEmilio
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
#ifndef METAPOPULATION_MOBILITY_STOCHASTIC_H
#define METAPOPULATION_MOBILITY_STOCHASTIC_H

#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/geography/locations.h"

#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"

#include <boost/numeric/ublas/vector_expression.hpp>
#include <cassert>

namespace mio
{

template <class Sim>
class LocationNode : public SimulationNode<Sim>
{
    using Base = SimulationNode<Sim>;

public:
    template <class... Args, typename = std::enable_if_t<std::is_constructible<Sim, Args...>::value, void>>
    LocationNode(Args&&... args)
        : Base(std::forward<Args>(args)...)
        , m_location(0.000, 0.000)
    {
    }

    auto get_location() const
    {
        return m_location;
    }

private:
    mio::geo::GeographicalLocation m_location; // location of the node
};

/**
 * node functor for mobility-based simulation.
 * @see SimulationNode::advance
 */
template <class Sim>
void advance_model(double t, double dt, LocationNode<Sim>& node)
{
    node.advance(t, dt);
}

/**
 * represents the mobility between two nodes.
 */
class MobilityEdgeDirected
{
public:
    /**
     * create edge with timed movement parameters.
     * @param params mobility rate for each group and compartment
     */
    // MobilityEdgeDirected(const MobilityParametersTimed& params)
    //     : m_parameters(params)
    // {
    // }

    // auto next_event_time() const
    // {
    //     return m_parameters.next_event_time();
    // }

    MobilityEdgeDirected() = default;

    /**
         * compute mobility from node_from to node_to for a given event
         * @param[in] event index specifying which compartment and age group change nodes
         * @param node_from node that people changed from
         * @param node_to node that people changed to
         */
    template <class Sim>
    void apply_mobility(double num_moving, LocationNode<Sim>& node_from, LocationNode<Sim>& node_to);

    // private:
    // MobilityParametersTimed m_parameters;
};

template <class Sim>
void MobilityEdgeDirected::apply_mobility(double num_moving, LocationNode<Sim>& node_from, LocationNode<Sim>& node_to)
{
    // auto next_event = m_parameters.process_next_event();
    // auto num_moving = next_event.number;
    // auto num_available = boost::numeric::ublas::sum(node_from.get_result().get_last_value());
    auto rng          = mio::RandomNumberGenerator();
    auto distribution = DiscreteDistributionInPlace<int>();

    for (int i = 0; i < num_moving; ++i) {
        auto group = distribution(rng, {node_from.get_result().get_last_value()});
        node_from.get_result().get_last_value()[group] -= 1;
        node_to.get_result().get_last_value()[group] += 1;
    }
}

template <class Sim>
void apply_timed_mobility(double t, double num_moving, MobilityEdgeDirected& edge, LocationNode<Sim>& node_from,
                          LocationNode<Sim>& node_to)
{
    // if (edge.next_event_time() >= t + dt) {
    //     return;
    // }
    mio::unused(t);
    edge.apply_mobility(num_moving, node_from, node_to);
}
// /**get_last_value
//      * edge functor for mobility-based simulation.
//      * @see MobilityEdgeDirected::apply_mobility
//      */
// template <class Sim, class StochasticEdge>
// void apply_mobility(StochasticEdge& mobilityEdge, size_t event, SimulationNode<Sim>& node_from,
//                     SimulationNode<Sim>& node_to)
// {
//     mobilityEdge.apply_mobility(event, node_from, node_to);
// }

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
AsymmetricGraphSimulation<Graph<LocationNode<Sim>, MobilityEdgeDirected>>
make_mobility_sim(FP t0, FP dt, const Graph<LocationNode<Sim>, MobilityEdgeDirected>& graph)
{
    return make_asymmetric_graph_sim(
        t0, dt, graph, &advance_model<Sim>,
        static_cast<void (*)(FP, FP, LocationNode<Sim>&, LocationNode<Sim>&)>(apply_timed_mobility<Sim>));
}

template <typename FP, class Sim>
AsymmetricGraphSimulation<Graph<LocationNode<Sim>, MobilityEdgeDirected>>
make_mobility_sim(FP t0, FP dt, Graph<LocationNode<Sim>, MobilityEdgeDirected>&& graph)
{
    using GraphSim = AsymmetricGraphSimulation<Graph<LocationNode<Sim>, MobilityEdgeDirected>, FP, FP,
                                               void (*)(FP, FP, mio::MobilityEdgeDirected&, mio::LocationNode<Sim>&,
                                                        mio::LocationNode<Sim>&),
                                               void (*)(FP, FP, mio::LocationNode<Sim>&)>;
    return GraphSim(t0, dt, std::move(graph), &advance_model<Sim>, &apply_timed_mobility<Sim>);
}

/** @} */

} // namespace mio

#endif //METAPOPULATION_MOBILITY_STOCHASTIC_H
