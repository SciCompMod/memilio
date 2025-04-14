/* 
* Copyright (C) 2020-2025 MEmilio
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
using MobilityCoefficients = DampingMatrixExpression<VectorDampings>;

/**
 * parameters that influence mobility.
 */
class MobilityParametersStochastic
{
public:
    /**
     * constructor from mobility coefficients.
     * @param coeffs mobility coefficients
     */
    MobilityParametersStochastic(const MobilityCoefficients& coeffs)
        : m_coefficients(coeffs)
    {
    }

    /**
     * constructor from mobility coefficients.
     * @param coeffs mobility coefficients
     */
    MobilityParametersStochastic(const Eigen::VectorXd& coeffs)
        : m_coefficients(MobilityCoefficients(coeffs))
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
    const MobilityCoefficients& get_coefficients() const
    {
        return m_coefficients;
    }
    MobilityCoefficients& get_coefficients()
    {
        return m_coefficients;
    }
    /**
     * @param coeffs the mobility coefficients.
     */
    void set_coefficients(const MobilityCoefficients& coeffs)
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
        auto c   = obj.expect_element("Coefficients", Tag<MobilityCoefficients>{});
        return apply(
            io,
            [](auto&& c_) {
                MobilityParametersStochastic params(c_);
                return params;
            },
            c);
    }

private:
    MobilityCoefficients m_coefficients; //one per group and compartment
};

/** 
 * represents the mobility between two nodes.
 */
class MobilityEdgeStochastic
{
public:
    /**
     * create edge with coefficients.
     * @param coeffs mobility rate for each group and compartment
     */
    MobilityEdgeStochastic(const MobilityParametersStochastic& params)
        : m_parameters(params)
    {
    }

    /**
     * create edge with coefficients.
     * @param coeffs mobility rate for each group and compartment
     */
    MobilityEdgeStochastic(const Eigen::VectorXd& coeffs)
        : m_parameters(coeffs)
    {
    }

    /**
     * get the mobility parameters.
     */
    const MobilityParametersStochastic& get_parameters() const
    {
        return m_parameters;
    }

    /**
     * get the cumulative transition rate of the edge.
    */
    template <class Sim>
    Eigen::VectorXd get_transition_rates(SimulationNode<Sim>& node_from)
    {
        Eigen::VectorXd transitionRates(node_from.get_last_state().size());
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
    void apply_mobility(size_t event, SimulationNode<Sim>& node_from, SimulationNode<Sim>& node_to);

private:
    MobilityParametersStochastic m_parameters;
};

template <class Sim>
void MobilityEdgeStochastic::apply_mobility(size_t event, SimulationNode<Sim>& node_from, SimulationNode<Sim>& node_to)
{
    node_from.get_result().get_last_value()[event] -= 1;
    node_to.get_result().get_last_value()[event] += 1;
}

/**
 * edge functor for mobility-based simulation.
 * @see MobilityEdgeStochastic::apply_mobility
 */
template <class Sim, class StochasticEdge>
void apply_mobility(StochasticEdge& mobilityEdge, size_t event, SimulationNode<Sim>& node_from,
                    SimulationNode<Sim>& node_to)
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
template <class Sim>
GraphSimulationStochastic<Graph<SimulationNode<Sim>, MobilityEdgeStochastic>>
make_mobility_sim(double t0, double dt, const Graph<SimulationNode<Sim>, MobilityEdgeStochastic>& graph)
{
    return make_graph_sim_stochastic(
        t0, dt, graph, &advance_model<Sim>,
        static_cast<void (*)(MobilityEdgeStochastic&, size_t, SimulationNode<Sim>&, SimulationNode<Sim>&)>(
            &apply_mobility<Sim, MobilityEdgeStochastic>));
}

template <class Sim>
GraphSimulationStochastic<Graph<SimulationNode<Sim>, MobilityEdgeStochastic>>
make_mobility_sim(double t0, double dt, Graph<SimulationNode<Sim>, MobilityEdgeStochastic>&& graph)
{
    return make_graph_sim_stochastic(
        t0, dt, std::move(graph), &advance_model<Sim>,
        static_cast<void (*)(MobilityEdgeStochastic&, size_t, SimulationNode<Sim>&, SimulationNode<Sim>&)>(
            &apply_mobility<Sim, MobilityEdgeStochastic>));
}

/** @} */

} // namespace mio

#endif //METAPOPULATION_MOBILITY_STOCHASTIC_H
