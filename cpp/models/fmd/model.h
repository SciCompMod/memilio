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
#ifndef FMD_MODEL_H
#define FMD_MODEL_H

#include "memilio/config.h"
#include "memilio/geography/rtree.h"
#include "memilio/geography/distance.h"
#include "memilio/mobility/metapopulation_mobility_asymmetric.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/graph_builder.h"
#include "smm/simulation.h"
#include "smm/parameters.h"
#include "fmd/infection_state.h"
#include "fmd/adoption_rates.h"

namespace mio
{
namespace fmd
{

using Model = mio::smm::Model<ScalarType, InfectionState, Status, mio::regions::Region>;

using Builder = mio::GraphBuilder<
    mio::LocationNode<ScalarType, mio::smm::Simulation<ScalarType, InfectionState, Status, mio::regions::Region>>,
    mio::MobilityEdgeDirected<ScalarType>>;
using Graph = mio::Graph<
    mio::LocationNode<ScalarType, mio::smm::Simulation<ScalarType, InfectionState, Status, mio::regions::Region>>,
    mio::MobilityEdgeDirected<ScalarType>>;

/**
 * @brief Create a model object with given number of susceptible individuals and adoption rates.
 * 
 * @param num_sus Number of susceptible individuals.
 * @param adoption_rates Vector of adoption rates.
 * @return Model The created model.
 */
Model create_model(size_t num_sus, const std::vector<AR>& adoption_rates = generic_adoption_rates())
{
    auto home = mio::regions::Region(0);
    Model curr_model(Status{InfectionState::Count}, mio::regions::Region(1));
    curr_model.populations[{home, InfectionState::S}]                                              = num_sus;
    curr_model.populations[{home, InfectionState::E}]                                              = 0;
    curr_model.populations[{home, InfectionState::I}]                                              = 0;
    curr_model.populations[{home, InfectionState::INS}]                                            = 0;
    curr_model.populations[{home, InfectionState::ICS}]                                            = 0;
    curr_model.populations[{home, InfectionState::R}]                                              = 0;
    curr_model.populations[{home, InfectionState::V}]                                              = 0;
    curr_model.populations[{home, InfectionState::D}]                                              = 0;
    curr_model.parameters.get<mio::smm::AdoptionRates<ScalarType, Status, mio::regions::Region>>() = adoption_rates;
    return curr_model;
}

/**
 * @brief Generate and set regional neighbors for all nodes in the graph based on given distances.
 * 
 * The function constructs an RTree from the node locations and queries it for neighbors within the specified distances.
 * 
 * @param graph The graph for which to generate neighbors.
 * @param distances The distances to consider for neighbor generation.
 */
void generate_neighbours(Graph& graph, std::vector<mio::geo::Distance> distances)
{
    auto nodes = graph.nodes() | std::views::transform([](const auto& node) {
                     return &node.property;
                 });
    auto tree  = mio::geo::RTree(nodes.begin(), nodes.end());

    for (auto& node : graph.nodes()) {
        node.property.set_regional_neighbors(tree.in_range_indices_query(node.property.get_location(), distances));
    }
}

} // namespace fmd
} // namespace mio

#endif // FMD_MODEL_H