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
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/timer/auto_timer.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/random_number_generator.h"
#include "smm/simulation.h"
#include "smm/parameters.h"
#include "thirdparty/csv.h"
#include <ranges>
#include <omp.h>

enum class InfectionState
{
    S,
    E,
    I,
    INS,
    ICS,
    R,
    D,
    Count
};

int main(int /*argc*/, char** /*argv*/)
{
    const auto t0   = 0.;
    const auto tmax = 100.;
    const auto dt   = 1.; //initial time step

    //total compartment sizes

    using Model = mio::smm::Model<ScalarType, 1, InfectionState>;
    auto home   = mio::regions::Region(0);
    auto S      = InfectionState::S;
    auto E      = InfectionState::E;
    auto I      = InfectionState::I;
    auto INS    = InfectionState::INS;
    auto ICS    = InfectionState::ICS;
    auto R      = InfectionState::R;
    auto D      = InfectionState::D;

    std::vector<mio::AdoptionRate<ScalarType, InfectionState>> adoption_rates;
    // Adoption rates corresponding to our model, paramters are arbitrary
    adoption_rates.push_back({S, E, home, 0.2, {{I, 0.8}, {INS, 0.1}, {ICS, 0.5}}});
    adoption_rates.push_back({E, I, home, 0.2, {}});
    adoption_rates.push_back({I, INS, home, 0.1, {}});
    adoption_rates.push_back({I, ICS, home, 0.1, {}});
    adoption_rates.push_back({ICS, D, home, 0.6, {}});
    adoption_rates.push_back({ICS, R, home, 0.4, {}});
    adoption_rates.push_back({INS, R, home, 0.5, {}});

    mio::Graph<mio::LocationNode<ScalarType, mio::smm::Simulation<ScalarType, 1, InfectionState>>,
               mio::MobilityEdgeDirected<ScalarType>>
        graph;
    mio::log_info("Starting Graph generation");
    {
        mio::timing::AutoTimer<"Graph Nodes Generation"> timer;
<<<<<<< HEAD
        io::CSVReader<4> farms("../../farms.csv");
=======
        io::CSVReader<4> farms("/home/kilian/Documents/projects/memilio-asymmetric-graph/farms.csv");
>>>>>>> dab1c992aa3c028686e0440a7f20cf4de6eb8871
        farms.read_header(io::ignore_extra_column, "farms", "num_cows", "latitude", "longitude");
        int farm_id, num_cows;
        double latitude, longitude;
        while (farms.read_row(farm_id, num_cows, latitude, longitude)) {
            Model curr_model;
            curr_model.populations[{home, InfectionState::S}]                                = num_cows;
            curr_model.populations[{home, InfectionState::E}]                                = 1;
            curr_model.populations[{home, InfectionState::I}]                                = 0;
            curr_model.populations[{home, InfectionState::INS}]                              = 0;
            curr_model.populations[{home, InfectionState::ICS}]                              = 0;
            curr_model.populations[{home, InfectionState::R}]                                = 0;
            curr_model.populations[{home, InfectionState::D}]                                = 0;
            curr_model.parameters.get<mio::smm::AdoptionRates<ScalarType, InfectionState>>() = adoption_rates;
            graph.add_node(farm_id, longitude, latitude, curr_model, t0);
        }
    }
    mio::log_info("Nodes added to Graph");
    auto rng = mio::RandomNumberGenerator();

    std::vector<std::vector<size_t>> interesting_indices;
    interesting_indices.push_back({Model().populations.get_flat_index({home, InfectionState::I})});
    {
        mio::timing::AutoTimer<"Graph Edges Generation"> timer;
<<<<<<< HEAD
        io::CSVReader<2> edges("../../edges.csv");
=======
        io::CSVReader<2> edges("/home/kilian/Documents/projects/memilio-asymmetric-graph/edges.csv");
>>>>>>> dab1c992aa3c028686e0440a7f20cf4de6eb8871
        edges.read_header(io::ignore_extra_column, "from", "to");
        size_t from, to;
        while (edges.read_row(from, to)) {
            graph.add_edge(from, to, interesting_indices);
            // graph.lazy_add_edge(from, to, interesting_indices);
        }
        // graph.sort_edges();
    }
    mio::log_info("Graph generated");

    auto nodes = graph.nodes() | std::views::transform([](const auto& node) {
                     return &node.property;
                 });
    auto tree  = mio::geo::RTree({nodes.begin(), nodes.end()});
    mio::log_info("RTree generated");

    for (auto& node : graph.nodes()) {
        node.property.set_regional_neighbors(tree.inrange_indices_query(node.property.get_location(), {2.0}));
    }

    mio::log_info("Neighbors set");
    auto sim = mio::make_mobility_sim(t0, dt, std::move(graph));

    io::CSVReader<5> exchanges("/home/kilian/Documents/projects/memilio-asymmetric-graph/trade.csv");
    exchanges.read_header(io::ignore_extra_column, "date", "num_animals", "from", "to", "edge");

    int date, num_animals, edge;
    size_t from, to;
    while (exchanges.read_row(date, num_animals, from, to, edge)) {
        sim.add_exchange(date, num_animals, from, to);
    }
    mio::log_info("Exchanges added");

    // #ifdef MEMILIO_ENABLE_OPENMP
    // #pragma omp parallel for
    //     for (size_t i = 0; i < 100; ++i) {
    // #endif

    auto sim2(sim);
    mio::log_info("new Simulation created");

    sim2.advance(tmax);
    mio::log_info("Simulation finished");

    // #ifdef MEMILIO_ENABLE_OPENMP
    //     }
    // #endif

    // sim2.get_graph().nodes()[28].property.get_result().print_table({"S", "E", "I", "R", "D"});
    // std::cout << "Second Table" << std::endl;
    // sim2.get_graph().nodes()[1].property.get_result().print_table({"S", "E", "I", "R", "D"});

    // auto& edge_1_0 = sim2.get_graph().edges()[1];
    // auto& results  = edge_1_0.property.get_mobility_results();
    // results.print_table({"Commuter Sick", "Commuter Total"});

    // auto exchange_results = sim2.sum_exchanges();
    // mio::log_info("Sum of exchanged sick animals: {}", exchange_results[0]);
    // mio::log_info("Sum of exchanged animals: {}", exchange_results[1]);

    auto sth = sim2.exchanges_per_timestep().export_csv("Exchange_statistics.csv", {"Commuter Sick", "Commuter Total"});

    // for (auto node : sim2.get_graph().nodes()) {
    //     if (node.property.get_result().get_num_time_points() < num_time_points) {
    //         mio::log_error("Node with inconsistent number of time points in results.");
    //     }
    // }

    // exchange_results = sim2.sum_nodes();
    // mio::log_info("{}, {}, {}, {}", exchange_results[0], exchange_results[1], exchange_results[2], exchange_results[3]);

    sth = sim2.statistics_per_timestep().export_csv("Simulation_statistics.csv");
    // // auto combined_results = sim2.combine_node_results();
    // // combined_results.print_table({"S", "E", "I", "R", "D"});
    // // auto ioresult = combined_results.export_csv("Simulation_results.csv");

    // sim2.statistics_per_timestep({0, 1, 2, 3, 4}).print_table({"S", "E", "I", "R", "D"});
    mio::log_info("Finished postprocessing");

    return 0;
}
