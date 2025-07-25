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

#include "memilio/compartments/simulation.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"

#include "boost/filesystem.hpp"

#include <cassert>
#include <functional>
#include <queue>
#include <streambuf>
#include <vector>

namespace mio
{

class MobilityParametersTimed
{

public:
    /**
     * @brief Construct a new Mobility Parameters Timed object
     * 
     */
    MobilityParametersTimed()
        : _exchanges{} {};
    MobilityParametersTimed(std::filebuf& input_data)
        : _exchanges{}
    {
        insert_input_data(input_data);
    };
    MobilityParametersTimed(int time, int number, int to)
        : _exchanges{}
    {
        _exchanges.push(ExchangeData(time, number, to));
    };

    /**
    * @brief Return the number of exchanged items in the next exchange event.
    * 
    * @return auto 
    */
    auto next_event_number()
    {
        return _exchanges.top().number;
    };
    /**
     * @brief Return the timepoint of the next exchangeclass EdgePropertyT event.
     * 
     * @return auto 
     */
    auto next_event_time() const
    {
        return _exchanges.top().time;
    };
    /**
     * @brief Return the destination node id of the next exchange
     * 
     * @return auto 
     */
    auto next_event_node_to()
    {
        return _exchanges.top().node_to;
    }
    /**
     * @brief Return a const reference to the next event
     * 
     * @return auto 
     */
    auto next_event()
    {
        return _exchanges.top();
    }
    /**
     * @brief Delete the next event from the heap
     * 
     * @return auto 
     */
    auto pop_next_event()
    {
        return _exchanges.pop();
    }
    /**
     * @brief Return the ExchangeData for the next exchange event and delete it from the list.
     * 
     * @return auto 
     */
    auto process_next_event()
    {
        auto next_event = _exchanges.top();
        _exchanges.pop();
        return next_event;
    };

private:
    void insert_input_data(std::filebuf& input_data) {
        //...
    };

    /**
     * @brief Stores Timepoint and number of exchanged items for an exchange process.
     * 
     * @param time Timepoint of the exchange process
     * @param number Number of exchanged items
     */
    struct ExchangeData {
        double time;
        int number;
        int node_to;
    };

    struct CompareExchangeData {
        bool operator()(const ExchangeData& left, const ExchangeData& right)
        {
            return left.time > right.time;
        };
    };

private:
    std::priority_queue<ExchangeData, std::vector<ExchangeData>, CompareExchangeData> _exchanges;
};

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
    MobilityEdgeDirected(const MobilityParametersTimed& params)
        : m_parameters(params)
    {
    }

    /**
         * get the mobility parameters.
         */
    const MobilityParametersTimed& get_parameters() const
    {
        return m_parameters;
    }

    auto next_event_time() const
    {
        return m_parameters.next_event_time();
    }

    /**
         * compute mobility from node_from to node_to for a given event
         * @param[in] event index specifying which compartment and age group change nodes
         * @param node_from node that people changed from
         * @param node_to node that people changed to
         */
    template <class Sim>
    void apply_mobility(SimulationNode<Sim>& node_from, SimulationNode<Sim>& node_to);

private:
    MobilityParametersTimed m_parameters;
};

template <class Sim>
void MobilityEdgeDirected::apply_mobility(SimulationNode<Sim>& node_from, SimulationNode<Sim>& node_to)
{
    auto next_event = m_parameters.process_next_event();
    auto num_moving = next_event.number;
    node_from.get_result().get_last_value()[0] -= num_moving;
    node_to.get_result().get_last_value()[0] += num_moving;
}

template <class Sim>
void apply_timed_mobility(double t, double dt, MobilityEdgeDirected& edge, SimulationNode<Sim>& node_from,
                          SimulationNode<Sim>& node_to)
{
    if (edge.next_event_time() >= t + dt) {
        return;
    }
    edge.apply_mobility(node_from, node_to);
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
AsymmetricGraphSimulation<Graph<SimulationNode<Sim>, MobilityEdgeDirected>>
make_mobility_sim(FP t0, FP dt, const Graph<SimulationNode<Sim>, MobilityEdgeDirected>& graph)
{
    return make_asymmetric_graph_sim(
        t0, dt, graph, &advance_model<Sim>,
        static_cast<void (*)(FP, FP, SimulationNode<Sim>&, SimulationNode<Sim>&)>(apply_mobility<FP, Sim>));
}
// static_cast<void (*)(MobilityEdgeStochastic&, size_t, SimulationNode<Sim>&, SimulationNode<Sim>&)>(&apply_mobility<Sim, MobilityEdgeStochastic>)

template <typename FP, class Sim>
AsymmetricGraphSimulation<Graph<SimulationNode<Sim>, MobilityEdgeDirected>>
make_mobility_sim(FP t0, FP dt, Graph<SimulationNode<Sim>, MobilityEdgeDirected>&& graph)
{
    return make_asymmetric_graph_sim(
        t0, dt, std::move(graph), &advance_model<Sim>,
        static_cast<void (*)(FP, FP, MobilityEdgeDirected&, SimulationNode<Sim>&, SimulationNode<Sim>&)>(
            apply_timed_mobility<Sim>));
}

/** @} */

} // namespace mio

#endif //METAPOPULATION_MOBILITY_STOCHASTIC_H
