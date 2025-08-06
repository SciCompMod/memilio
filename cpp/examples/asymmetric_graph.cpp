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
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_asymmetric.h"
#include "memilio/mobility/graph.h"
#include "memilio/utils/logging.h"
#include "smm/simulation.h"
#include "smm/parameters.h"

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

    using Model = mio::smm::Model<1, InfectionState>;
    Model model;

    auto home = mio::regions::Region(0);

    model.populations[{home, InfectionState::E}] = num_exp;
    model.populations[{home, InfectionState::I}] = num_ins;
    model.populations[{home, InfectionState::R}] = num_rec;
    model.populations[{home, InfectionState::S}] = num_total - num_exp - num_ins - num_rec;

    std::vector<mio::AdoptionRate<InfectionState>> adoption_rates;
    adoption_rates.push_back({InfectionState::E, InfectionState::I, home, 0.2, {}});
    adoption_rates.push_back({InfectionState::I, InfectionState::R, home, 0.333, {}});
    adoption_rates.push_back({InfectionState::S, InfectionState::E, home, 0.2, {{InfectionState::I, 0.5}}});
    model.parameters.get<mio::smm::AdoptionRates<InfectionState>>() = adoption_rates;

    auto model2 = model;

    mio::Graph<mio::LocationNode<mio::smm::Simulation<1, InfectionState>>, mio::MobilityEdgeDirected> graph;
    graph.add_node(0, 12.0, 21.0, model, t0);
    graph.add_node(1, 12.0, 21.0, model2, t0);
    size_t num_nodes = 999;
    for (size_t i = 2; i < num_nodes; i++) {
        auto local_model = model;
        graph.add_node(i, 12.0, 21.0, local_model, t0);
    }

    auto param = mio::MobilityParametersTimed(2.0, 10, 1);
    graph.add_edge(0, 1);

    auto rng          = mio::RandomNumberGenerator();
    auto distribution = mio::DiscreteDistributionInPlace<size_t>();
    std::vector<double> uniform_vector(num_nodes, 1.0);
    mio::log_info("Nodes generated");
    for (size_t i = 0; i < 3 * num_nodes; ++i) {
        auto to   = distribution(rng, {uniform_vector});
        auto from = distribution(rng, {uniform_vector});
        graph.add_edge(from, to);
    }

    mio::log_info("Graph generated");

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

    return 0;
}
