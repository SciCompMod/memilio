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
#include "memilio/geography/geolocation.h"
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

#include <cassert>
#include <map>

namespace mio
{

template <typename FP, class Sim>
class LocationNode : public SimulationNode<FP, Sim>
{
    using Base = SimulationNode<FP, Sim>;

public:
    template <class... Args, typename = std::enable_if_t<std::is_constructible<Sim, Args...>::value, void>>
    LocationNode(FP x, FP y, Args&&... args)
        : Base(std::forward<Args>(args)...)
        , m_location(x, y)
        , regional_neighbor_indices{}
    {
    }

    auto get_location() const
    {
        return m_location;
    }

    void set_location(FP x, FP y)
    {
        m_location = mio::geo::GeographicalLocation(x, y);
    }

    FP get_y() const
    {
        return m_location.get_y();
    }

    FP get_x() const
    {
        return m_location.get_x();
    }

    void set_regional_neighbors(const std::vector<std::vector<size_t>>& neighbors)
    {
        regional_neighbor_indices = neighbors;
    }

    auto get_regional_neighbors() const
    {
        return regional_neighbor_indices;
    }

    auto is_quarantined() const
    {
        return m_is_quarantined;
    }

    void set_quarantined(bool quarantine)
    {
        m_is_quarantined = quarantine;
    }

private:
    mio::geo::GeographicalLocation m_location; // location of the node
    std::vector<std::vector<size_t>> regional_neighbor_indices;
    bool m_is_quarantined{false};
};

/**
 * parameters that influence mobility.
 */
template <typename FP>
class BirdFlightParameters
{
public:
    /**
     * constructor from mobility coefficients.
     * @param coeffs mobility coefficients
     */
    BirdFlightParameters()
        : m_coefficients()
    {
    }

    /**
     * constructor from mobility coefficients.
     * @param coeffs mobility coefficients
     */
    BirdFlightParameters(const Eigen::VectorX<FP>& coeffs)
        : m_coefficients()
    {
        for (Eigen::Index i = 0; i < coeffs.size(); ++i) {
            m_coefficients[(size_t)i] = coeffs[i];
        }
    }

    /**
     * equality comparison operators
     */
    //@{
    bool operator==(const BirdFlightParameters& other) const
    {
        return m_coefficients == other.m_coefficients;
    }
    bool operator!=(const BirdFlightParameters& other) const
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
    const std::map<size_t, FP>& get_coefficients() const
    {
        return m_coefficients;
    }
    std::map<size_t, FP>& get_coefficients()
    {
        return m_coefficients;
    }
    /**
     * @param coeffs the mobility coefficients.
     */
    void set_coefficients(const BirdFlightParameters<FP>& coeffs)
    {
        m_coefficients = coeffs;
    }

    void change_coefficient(size_t index, FP value)
    {
        m_coefficients[index] = value;
    }

    bool has_coefficient(size_t index)
    {
        return m_coefficients.contains(index);
    }

    size_t size() const
    {
        return m_coefficients.size();
    }

    bool empty() const
    {
        return m_coefficients.empty();
    }

    /**
     * serialize this.
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("BirdFlightParameters");
        obj.add_element("Coefficients", m_coefficients);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<BirdFlightParameters> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("BirdFlightParameters");
        auto c   = obj.expect_element("Coefficients", Tag<BirdFlightParameters<FP>>{});
        return apply(
            io,
            [](auto&& c_) {
                BirdFlightParameters params(c_);
                return params;
            },
            c);
    }

private:
    std::map<size_t, FP> m_coefficients;
};

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
    JollyEdge(const BirdFlightParameters<FP>& params)
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
    const BirdFlightParameters<FP>& get_parameters() const
    {
        return m_parameters;
    }

    // /**
    //  * get the cumulative transition rate of the edge.
    // */
    // template <class Sim>
    // Eigen::VectorX<FP> get_transition_rates(SimulationNode<FP, Sim>& node_from)
    // {
    //     Eigen::VectorX<FP> transitionRates(node_from.get_last_state().size());
    //     for (Eigen::Index i = 0; i < node_from.get_last_state().size(); ++i) {
    //         transitionRates[i] =
    //             node_from.get_last_state()(i) * m_parameters.get_coefficients().get_baseline()[(size_t)i];
    //     }
    //     return transitionRates;
    // }

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
    BirdFlightParameters<FP> m_parameters;
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
    mio::log_debug("JollyEdge::apply_mobility: moved {} individuals of compartment {} from node a to node b.",
                   batch_size, event);
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
JollyGraphSimulation<FP, Graph<LocationNode<FP, Sim>, JollyEdge<FP>>>
make_jolly_sim(FP t0, FP dt, const Graph<LocationNode<FP, Sim>, JollyEdge<FP>>& graph)
{
    return make_graph_sim_jolly<FP>(
        t0, dt, graph, &advance_model<FP, Sim>,
        static_cast<void (*)(size_t, FP, JollyEdge<FP>&, SimulationNode<FP, Sim>&, SimulationNode<FP, Sim>&)>(
            &apply_mobility<FP, Sim, JollyEdge<FP>>));
}

template <typename FP, class Sim>
JollyGraphSimulation<FP, Graph<LocationNode<FP, Sim>, JollyEdge<FP>>>
make_jolly_sim(FP t0, FP dt, Graph<LocationNode<FP, Sim>, JollyEdge<FP>>&& graph)
{
    return make_graph_sim_jolly<FP>(
        t0, dt, std::move(graph), &advance_model<FP, Sim>,
        static_cast<void (*)(size_t, FP, JollyEdge<FP>&, SimulationNode<FP, Sim>&, SimulationNode<FP, Sim>&)>(
            &apply_mobility<FP, Sim, JollyEdge<FP>>));
}

/** @} */

} // namespace mio

#endif //METAPOPULATION_MOBILITY_JOLLY_H
