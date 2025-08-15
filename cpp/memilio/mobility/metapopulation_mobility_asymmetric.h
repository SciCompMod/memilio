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

#include <algorithm>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <cassert>
#include <numeric>
#include <vector>

namespace mio
{

template <class Sim>
class LocationNode : public SimulationNode<Sim>
{
    using Base = SimulationNode<Sim>;

public:
    template <class... Args, typename = std::enable_if_t<std::is_constructible<Sim, Args...>::value, void>>
    LocationNode(double latitude, double longitude, Args&&... args)
        : Base(std::forward<Args>(args)...)
        , m_location(latitude, longitude)
        , regional_neighbor_indices{}
    {
    }

    auto get_location() const
    {
        return m_location;
    }

    void set_location(double latitude, double longitude)
    {
        m_location = mio::geo::GeographicalLocation(latitude, longitude);
    }

    double get_longitude() const
    {
        return m_location.get_longitude();
    }

    double get_latitude() const
    {
        return m_location.get_latitude();
    }

    void set_regional_neighbors(const std::vector<std::vector<size_t>>& neighbors)
    {
        regional_neighbor_indices = neighbors;
    }

private:
    mio::geo::GeographicalLocation m_location; // location of the node
    std::vector<std::vector<size_t>> regional_neighbor_indices;
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

    MobilityEdgeDirected(size_t size)
        : m_mobility_results(size)
    {
    }

    MobilityEdgeDirected(size_t size, const std::vector<std::vector<size_t>>& saved_compartment_indices)
        : m_mobility_results(size)
        , m_saved_compartment_indices(saved_compartment_indices)
    {
    }

    MobilityEdgeDirected(const std::vector<std::vector<size_t>>& saved_compartment_indices)
        : m_mobility_results(saved_compartment_indices.size() + 1)
        , m_saved_compartment_indices(saved_compartment_indices)
    {
    }

    /**
     * @brief Get the count of exchanges in selected compartments, along with the total number of exchanges.
     *
     * @return A reference to the TimeSeries object representing the mobility results.
     */
    TimeSeries<ScalarType>& get_mobility_results()
    {
        return m_mobility_results;
    }
    const TimeSeries<ScalarType>& get_mobility_results() const
    {
        return m_mobility_results;
    }

    /**
         * compute mobility from node_from to node_to for a given event
         * @param[in] event index specifying which compartment and age group change nodes
         * @param node_from node that people changed from
         * @param node_to node that people changed to
         */
    template <class Sim>
    void apply_mobility(double& t, double& num_moving, LocationNode<Sim>& node_from, LocationNode<Sim>& node_to);

private:
    // MobilityParametersTimed m_parameters;
    TimeSeries<double> m_mobility_results;
    std::vector<std::vector<size_t>> m_saved_compartment_indices;

    void add_mobility_result_time_point(const double t, std::vector<size_t>& travellers)
    {
        const size_t save_indices_size = this->m_saved_compartment_indices.size();
        if (save_indices_size > 0) {

            Eigen::VectorXd condensed_values = Eigen::VectorXd::Zero(save_indices_size + 1);

            // sum up the values of m_saved_compartment_indices for each group (e.g. Age groups)
            std::transform(this->m_saved_compartment_indices.begin(), this->m_saved_compartment_indices.end(),
                           condensed_values.data(), [&travellers](const auto& indices) {
                               return std::accumulate(indices.begin(), indices.end(), 0.0,
                                                      [&travellers](double sum, auto i) {
                                                          return sum + travellers[i];
                                                      });
                           });

            // the last value is the sum of commuters
            condensed_values[save_indices_size] = std::accumulate(travellers.begin(), travellers.end(), 0);

            // Move the condensed values to the m_mobility_results time series
            m_mobility_results.add_time_point(t, std::move(condensed_values));
        }
    }
};

template <class Sim>
void MobilityEdgeDirected::apply_mobility(double& t, double& num_moving, LocationNode<Sim>& node_from,
                                          LocationNode<Sim>& node_to)
{
    // auto next_event = m_parameters.process_next_event();
    // auto num_moving = next_event.number;
    // auto num_available = boost::numeric::ublas::sum(node_from.get_result().get_last_value());
    auto rng          = mio::RandomNumberGenerator();
    auto distribution = DiscreteDistributionInPlace<int>();
    std::vector<size_t> travellers(node_from.get_result().get_last_value().size(), 0);
    for (int i = 0; i < num_moving; ++i) {
        auto group = distribution(rng, {node_from.get_result().get_last_value()});
        node_from.get_result().get_last_value()[group] -= 1;
        travellers[group] += 1;
        node_to.get_result().get_last_value()[group] += 1;
    }
    add_mobility_result_time_point(t, travellers);
}

template <class Sim>
void apply_timed_mobility(double t, double num_moving, MobilityEdgeDirected& edge, LocationNode<Sim>& node_from,
                          LocationNode<Sim>& node_to)
{
    edge.apply_mobility(t, num_moving, node_from, node_to);
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
AsymmetricGraphSimulation<Graph<LocationNode<Sim>, MobilityEdgeDirected>>
make_mobility_sim(FP t0, FP dt, const Graph<LocationNode<Sim>, MobilityEdgeDirected>& graph)
{
    using GraphSim = AsymmetricGraphSimulation<Graph<LocationNode<Sim>, MobilityEdgeDirected>, FP, FP,
                                               void (*)(FP, FP, mio::MobilityEdgeDirected&, mio::LocationNode<Sim>&,
                                                        mio::LocationNode<Sim>&),
                                               void (*)(FP, FP, mio::LocationNode<Sim>&)>;
    return GraphSim(t0, dt, graph, &advance_model<Sim>, &apply_timed_mobility<Sim>);
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
