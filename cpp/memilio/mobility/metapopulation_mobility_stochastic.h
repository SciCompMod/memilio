/*
* Copyright (C) 2020-2026 MEmilio
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
#ifndef METAPOPULATION_MOBILITY_STOCHASTIC_H
#define METAPOPULATION_MOBILITY_STOCHASTIC_H

#include "memilio/compartments/simulation.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"

#include "boost/filesystem.hpp"

#include <cassert>

namespace mio
{

/**
 * status and age dependent mobility coefficients.
 */
template <typename FP>
using MobilityCoefficients = DampingMatrixExpression<FP, VectorDampings<FP>>;

/**
 * parameters that influence mobility.
 */
template <typename FP>
class MobilityParametersStochastic
{
public:
    /**
     * constructor from mobility coefficients.
     * @param coeffs mobility coefficients
     */
    MobilityParametersStochastic(const MobilityCoefficients<FP>& coeffs)
        : m_coefficients(coeffs)
    {
    }

    /**
     * constructor from mobility coefficients.
     * @param coeffs mobility coefficients
     */
    MobilityParametersStochastic(const Eigen::VectorX<FP>& coeffs)
        : m_coefficients(MobilityCoefficients<FP>(coeffs))
    {
    }

    /**
     * equality comparison operators
     */
    //@{
    bool operator==(const MobilityParametersStochastic& other) const
    {
        return m_coefficients == other.m_coefficients;
    }
    bool operator!=(const MobilityParametersStochastic& other) const
    {
        return m_coefficients != other.m_coefficients;
    }
    //@}

    /**
     * Get/Set the mobility coefficients.
     * The coefficients represent the rates for moving
     * from one node to another by age and infection compartment.
     * @{
     */
    /**
     * @return the mobility coefficients.
     */
    const MobilityCoefficients<FP>& get_coefficients() const
    {
        return m_coefficients;
    }
    MobilityCoefficients<FP>& get_coefficients()
    {
        return m_coefficients;
    }
    /**
     * @param coeffs the mobility coefficients.
     */
    void set_coefficients(const MobilityCoefficients<FP>& coeffs)
    {
        m_coefficients = coeffs;
    }

    /**
     * serialize this.
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("MobilityParameters");
        obj.add_element("Coefficients", m_coefficients);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<MobilityParametersStochastic> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("MobilityParameters");
        auto c   = obj.expect_element("Coefficients", Tag<MobilityCoefficients<FP>>{});
        return apply(
            io,
            [](auto&& c_) {
                MobilityParametersStochastic params(c_);
                return params;
            },
            c);
    }

private:
    MobilityCoefficients<FP> m_coefficients; //one per group and compartment
};

/**
 * represents the mobility between two nodes.
 */
template <typename FP>
class MobilityEdgeStochastic
{
public:
    /**
     * create edge with coefficients.
     * @param coeffs mobility rate for each group and compartment
     */
    MobilityEdgeStochastic(const MobilityParametersStochastic<FP>& params)
        : m_parameters(params)
    {
    }

    /**
     * create edge with coefficients.
     * @param coeffs mobility rate for each group and compartment
     */
    MobilityEdgeStochastic(const Eigen::VectorX<FP>& coeffs)
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
    void apply_mobility(size_t event, SimulationNode<FP, Sim>& node_from, SimulationNode<FP, Sim>& node_to);

private:
    MobilityParametersStochastic<FP> m_parameters;
};

template <typename FP>
template <class Sim>
void MobilityEdgeStochastic<FP>::apply_mobility(size_t event, SimulationNode<FP, Sim>& node_from,
                                                SimulationNode<FP, Sim>& node_to)
{
    node_from.get_result().get_last_value()[event] -= 1;
    node_to.get_result().get_last_value()[event] += 1;
}

/**
 * edge functor for mobility-based simulation.
 * @see MobilityEdgeStochastic::apply_mobility
 */
template <typename FP, class Sim, class StochasticEdge>
void apply_mobility(StochasticEdge& mobilityEdge, size_t event, SimulationNode<FP, Sim>& node_from,
                    SimulationNode<FP, Sim>& node_to)
{
    mobilityEdge.apply_mobility(event, node_from, node_to);
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
GraphSimulationStochastic<FP, Graph<SimulationNode<FP, Sim>, MobilityEdgeStochastic<FP>>>
make_mobility_sim(FP t0, FP dt, const Graph<SimulationNode<FP, Sim>, MobilityEdgeStochastic<FP>>& graph)
{
    return make_graph_sim_stochastic<FP>(
        t0, dt, graph, &advance_model<FP, Sim>,
        static_cast<void (*)(MobilityEdgeStochastic<FP>&, size_t, SimulationNode<FP, Sim>&, SimulationNode<FP, Sim>&)>(
            &apply_mobility<FP, Sim, MobilityEdgeStochastic<FP>>));
}

template <typename FP, class Sim>
GraphSimulationStochastic<FP, Graph<SimulationNode<FP, Sim>, MobilityEdgeStochastic<FP>>>
make_mobility_sim(FP t0, FP dt, Graph<SimulationNode<FP, Sim>, MobilityEdgeStochastic<FP>>&& graph)
{
    return make_graph_sim_stochastic<FP>(
        t0, dt, std::move(graph), &advance_model<FP, Sim>,
        static_cast<void (*)(MobilityEdgeStochastic<FP>&, size_t, SimulationNode<FP, Sim>&, SimulationNode<FP, Sim>&)>(
            &apply_mobility<FP, Sim, MobilityEdgeStochastic<FP>>));
}

/** @} */

} // namespace mio

#endif //METAPOPULATION_MOBILITY_STOCHASTIC_H
