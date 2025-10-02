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
#include "memilio/config.h"
#include "memilio/geography/locations.h"
#include "memilio/geography/tree.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_asymmetric.h"
#include "memilio/mobility/graph.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/random_number_generator.h"
#include "smm/simulation.h"
#include "smm/parameters.h"
#include <ranges>

enum class InfectionState
{
    S,
    E,
    I,
    R,
    Count
};

int main(int /*argc*/, char** /*argv*/)
{
    const auto t0   = 0.;
    const auto tmax = 100.;
    const auto dt   = 1.; //initial time step

    //total compartment sizes
    double num_total = 10000, num_exp = 200, num_ins = 50, num_rec = 0;

    using Model = mio::smm::Model<ScalarType, 1, InfectionState>;
    Model model;

    auto home = mio::regions::Region(0);

    model.populations[{home, InfectionState::E}] = num_exp;
    model.populations[{home, InfectionState::I}] = num_ins;
    model.populations[{home, InfectionState::R}] = num_rec;
    model.populations[{home, InfectionState::S}] = num_total - num_exp - num_ins - num_rec;

    std::vector<mio::AdoptionRate<ScalarType, InfectionState>> adoption_rates;
    adoption_rates.push_back({InfectionState::E, InfectionState::I, home, 0.2, {}});
    adoption_rates.push_back({InfectionState::I, InfectionState::R, home, 0.333, {}});
    adoption_rates.push_back({InfectionState::S, InfectionState::E, home, 0.2, {{InfectionState::I, 0.5}}});
    model.parameters.get<mio::smm::AdoptionRates<ScalarType, InfectionState>>() = adoption_rates;

    auto model2                                   = model;
    model2.populations[{home, InfectionState::I}] = 9000;
    model2.populations[{home, InfectionState::S}] = 800;

    auto rng = mio::RandomNumberGenerator();

    mio::Graph<mio::LocationNode<ScalarType, mio::smm::Simulation<ScalarType, 1, InfectionState>>,
               mio::MobilityEdgeDirected<ScalarType>>
        graph;
    graph.add_node(0, mio::UniformDistribution<ScalarType>::get_instance()(rng, 6.0, 15.0),
                   mio::UniformDistribution<ScalarType>::get_instance()(rng, 48.0, 54.0), model2, t0);
    graph.add_node(1, mio::UniformDistribution<ScalarType>::get_instance()(rng, 6.0, 15.0),
                   mio::UniformDistribution<ScalarType>::get_instance()(rng, 48.0, 54.0), model2, t0);
    size_t num_nodes = 2000;
    for (size_t i = 2; i < num_nodes; i++) {
        auto local_model = model;
        graph.add_node(i, mio::UniformDistribution<ScalarType>::get_instance()(rng, 6.0, 15.0),
                       mio::UniformDistribution<ScalarType>::get_instance()(rng, 48.0, 54.0), local_model, t0);
    }

    std::vector<std::vector<size_t>> interesting_indices;
    interesting_indices.push_back({model.populations.get_flat_index({home, InfectionState::I})});

    auto param = mio::MobilityParametersTimed(2.0, 10, 1);
    graph.add_edge(0, 1, interesting_indices);

    auto distribution = mio::DiscreteDistributionInPlace<size_t>();

    std::vector<ScalarType> uniform_vector(num_nodes, 1.0);
    mio::log_info("Nodes generated");
    for (size_t i = 0; i < 3 * num_nodes; ++i) {
        auto to   = distribution(rng, {uniform_vector});
        auto from = distribution(rng, {uniform_vector});
        graph.add_edge(from, to, interesting_indices);
    }

    mio::log_info("Graph generated");

    auto nodes = graph.nodes() | std::views::transform([](const auto& node) {
                     return &node.property;
                 });
    // mio::log_info("Value: {}", nodes.begin()->get_latitude());
    auto tree = mio::geo::RTree(nodes.begin(), nodes.end());
    mio::log_info("RTree generated");

    for (auto& node : graph.nodes()) {
        node.property.set_regional_neighbors(
            tree.inrange_indices_query(node.property.get_location(), {3.0, 5.0, 10.0}));
    }

    mio::log_info("Neighbors set");

    auto sim = mio::make_mobility_sim(t0, dt, std::move(graph));

    for (size_t i = 0; i < 3 * num_nodes; i++) {
        sim.add_exchange(distribution(rng, {std::vector<double>(100, 1)}) + 1, 10, i);
    }
    for (size_t i = 0; i < 3 * num_nodes; i++) {
        sim.add_exchange(distribution(rng, {std::vector<double>(100, 1)}) + 1, 10, i);
    }
    for (size_t i = 0; i < 3 * num_nodes; i++) {
        sim.add_exchange(distribution(rng, {std::vector<double>(100, 1)}) + 1, 10, i);
    }

    mio::log_info("Number of exchanges: {}", sim.get_parameters().size());

    mio::log_info("Exchanges added");

    sim.advance(tmax);
    mio::log_info("Simulation finished");

    // std::cout << "First table" << std::endl;
    // sim.get_graph().nodes()[0].property.get_result().print_table({"S", "E", "I", "R"});
    // std::cout << "Second Table" << std::endl;
    // sim.get_graph().nodes()[1].property.get_result().print_table({"S", "E", "I", "R"});

    auto& edge_1_0 = sim.get_graph().edges()[1];
    auto& results  = edge_1_0.property.get_mobility_results();
    results.print_table({"Commuter Sick", "Commuter Total"});

    return 0;
}
