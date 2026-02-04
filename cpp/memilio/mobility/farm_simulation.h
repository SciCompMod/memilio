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
#ifndef FARM_SIMULATION_H
#define FARM_SIMULATION_H

#include "memilio/utils/random_number_generator.h"
#include "memilio/geography/geolocation.h"
// #include "memilio/mobility/farm_graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"

#include <algorithm>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <cassert>
#include <numeric>
#include <vector>

namespace mio
{

template <typename FP, class Sim>
class FarmNode : public SimulationNode<FP, Sim>
{
    using Base = SimulationNode<FP, Sim>;

public:
    template <class... Args, typename = std::enable_if_t<std::is_constructible<Sim, Args...>::value, void>>
    FarmNode(FP latitude, FP longitude, size_t farm_type, size_t farm_size, int population, int slaughter,
             Args&&... args)
        : Base(std::forward<Args>(args)...)
        , m_location(latitude, longitude)
        , regional_neighbor_indices{}
        , m_farm_type(farm_type)
        , m_farm_size(farm_size)
        , m_population_date(population)
        , m_slaughter_date(slaughter)
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

    FP get_x() const
    {
        return m_location.get_x();
    }
    FP get_longitude() const
    {
        return m_location.get_x();
    }

    FP get_latitude() const
    {
        return m_location.get_y();
    }

    FP get_y() const
    {
        return m_location.get_y();
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

    void set_population_date(FP date)
    {
        m_population_date = date;
    }
    FP get_population_date() const
    {
        return m_population_date;
    }
    void set_slaughter_date(FP date)
    {
        m_slaughter_date = date;
    }
    FP get_slaughter_date() const
    {
        return m_slaughter_date;
    }
    size_t get_farm_type() const
    {
        return m_farm_type;
    }
    size_t get_farm_size() const
    {
        return m_farm_size;
    }
    void slaughter()
    {
        // ...
    }

    bool is_suspicious() const
    {
        return m_date_suspicion > -1;
    }

    int get_date_suspicion() const
    {
        return m_date_suspicion;
    }

    int get_date_confirmation() const
    {
        return m_date_confirmation;
    }
    void set_date_suspicion(int date)
    {
        m_date_suspicion = date;
    }
    void set_date_confirmation(int date)
    {
        m_date_confirmation = date;
    }
    void set_capacity(size_t size)
    {
        m_farm_size = size;
    }
    size_t get_capacity() const
    {
        return m_farm_size;
    }

    size_t get_type() const
    {
        return m_farm_type;
    }

    void set_infection_status(bool infected)
    {
        m_is_infected = infected;
    }
    bool get_infection_status() const
    {
        return m_is_infected;
    }

private:
    mio::geo::GeographicalLocation m_location; // location of the node
    std::vector<std::vector<size_t>> regional_neighbor_indices;
    bool m_is_quarantined{false};
    size_t m_farm_type{0};
    size_t m_farm_size{0};
    FP m_population_date{0};
    FP m_slaughter_date{0};
    FP m_date_suspicion{-1};
    FP m_date_confirmation{-1};
    bool m_is_infected{false};
};

/**
 * node functor for mobility-based simulation.
 * @see SimulationNode::advance
 */
template <typename FP, class Sim>
void advance_model(FP t, FP dt, FarmNode<FP, Sim>& node)
{
    node.advance(t, dt);
}

/**
 * represents the mobility between two nodes.
 */
template <typename FP>
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
    void apply_mobility(const FP t, const FP num_moving, FarmNode<FP, Sim>& node_from, FarmNode<FP, Sim>& node_to,
                        mio::RandomNumberGenerator& rng);

private:
    // MobilityParametersTimed m_parameters;
    TimeSeries<FP> m_mobility_results;
    std::vector<std::vector<size_t>> m_saved_compartment_indices;

    void add_mobility_result_time_point(const FP t, std::vector<size_t>& travellers)
    {
        const size_t save_indices_size = this->m_saved_compartment_indices.size();
        if (save_indices_size > 0) {

            Eigen::VectorXd condensed_values = Eigen::VectorXd::Zero(save_indices_size + 1);

            // sum up the values of m_saved_compartment_indices for each group (e.g. Age groups)
            std::transform(this->m_saved_compartment_indices.begin(), this->m_saved_compartment_indices.end(),
                           condensed_values.data(), [&travellers](const auto& indices) {
                               return std::accumulate(indices.begin(), indices.end(), 0.0,
                                                      [&travellers](FP sum, auto i) {
                                                          return sum + travellers[i];
                                                      });
                           });

            // the last value is the sum of commuters
            condensed_values[save_indices_size] = std::accumulate(travellers.begin(), travellers.end(), 0.0);

            // Move the condensed values to the m_mobility_results time series
            m_mobility_results.add_time_point(t, std::move(condensed_values));
        }
    }
};

template <typename FP>
template <class Sim>
void MobilityEdgeDirected<FP>::apply_mobility(const FP t, const FP num_moving, FarmNode<FP, Sim>& node_from,
                                              FarmNode<FP, Sim>& node_to, mio::RandomNumberGenerator& rng)
{
    // auto next_event = m_parameters.process_next_event();
    // auto num_moving = next_event.number;
    // auto num_available = boost::numeric::ublas::sum(node_from.get_result().get_last_value());
    auto distribution = DiscreteDistributionInPlace<int>();
    std::vector<size_t> travellers(node_from.get_result().get_last_value().size(), 0);
    if (num_moving > std::accumulate(node_from.get_result().get_last_value().begin(),
                                     node_from.get_result().get_last_value().end(), 0.0)) {
        mio::log_warning("Trying to move more individuals than available ({}) at time {}.", num_moving, t);
    }
    else {
        for (int i = 0; i < num_moving; ++i) {
            auto group = distribution(rng, {node_from.get_result().get_last_value()});
            node_from.get_result().get_last_value()[group] -= 1;
            travellers[group] += 1;
            node_to.get_result().get_last_value()[group] += 1;
        }
    }
    add_mobility_result_time_point(t, travellers);
}

template <typename FP, class Sim>
void apply_timed_mobility(const FP t, const FP num_moving, MobilityEdgeDirected<FP>& edge, FarmNode<FP, Sim>& node_from,
                          FarmNode<FP, Sim>& node_to, mio::RandomNumberGenerator& rng)
{
    edge.apply_mobility(t, num_moving, node_from, node_to, rng);
}

/** @} */

} // namespace mio

#endif //METAPOPULATION_MOBILITY_STOCHASTIC_H
