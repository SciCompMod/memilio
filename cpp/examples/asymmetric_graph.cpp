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
#include "memilio/config.h"
#include "memilio/geography/geolocation.h"
#include "memilio/geography/rtree.h"
#include "memilio/geography/distance.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_asymmetric.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/graph_builder.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/index.h"
#include "memilio/utils/logging.h"
#include "memilio/timer/auto_timer.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/random_number_generator.h"
#include "smm/simulation.h"
#include "smm/parameters.h"
#include "fmd/infection_state.h"
#include "fmd/model.h"
#include "fmd/adoption_rates.h"
#include "thirdparty/csv.h"
#include <ranges>
#include <omp.h>

int main(int /*argc*/, char** /*argv*/)
{
    const auto t0   = 0.;
    const auto tmax = 100.;
    const auto dt   = 1.; //initial time step

    using mio::fmd::InfectionState;
    using Status = mio::Index<InfectionState>;
    using mio::regions::Region;

    //total compartment sizes

    using Model = mio::smm::Model<ScalarType, InfectionState, Status, Region>;
    auto home   = Region(0);

    auto adoption_rates = mio::fmd::generic_adoption_rates();

    mio::fmd::Builder builder;
    mio::log_info("Starting Graph generation");
    {
        mio::timing::AutoTimer<"Graph Nodes Generation"> timer;
        io::CSVReader<4> farms("../../farms200000.csv");
        farms.read_header(io::ignore_extra_column, "farms", "num_cows", "latitude", "longitude");
        int farm_id, num_cows;
        double latitude, longitude;
        while (farms.read_row(farm_id, num_cows, latitude, longitude)) {
            builder.add_node(farm_id, longitude, latitude, mio::fmd::create_model(num_cows, adoption_rates), t0);
        }
    }
    mio::log_info("Nodes added to Graph");
    auto rng = mio::RandomNumberGenerator();

    std::vector<std::vector<size_t>> interesting_indices;
    interesting_indices.push_back(
        {Model(Status{InfectionState::Count}, Region(1)).populations.get_flat_index({home, InfectionState::I})});
    io::CSVReader<2> edges("../../edges200000.csv");
    edges.read_header(io::ignore_extra_column, "from", "to");
    size_t from, to;
    while (edges.read_row(from, to)) {
        builder.add_edge(from, to, interesting_indices);
    }
    auto graph = std::move(builder).build();
    // // mio::log_info("Graph generated");

    auto nodes = graph.nodes() | std::views::transform([](const auto& node) {
                     return &node.property;
                 });
    auto tree  = mio::geo::RTree(nodes.begin(), nodes.end());
    // mio::log_info("RTree generated");

    for (auto& node : graph.nodes()) {
        mio::timing::AutoTimer<"neighbourhood search"> timer;
        node.property.set_regional_neighbors(
            tree.in_range_indices_query(node.property.get_location(), {mio::geo::kilometers(5.0)}));
    }

    mio::log_info("Neighbors set");
    auto sim = mio::make_mobility_sim(t0, dt, std::move(graph));

    io::CSVReader<4> exchanges("../../trade200000.csv");
    exchanges.read_header(io::ignore_extra_column, "date", "num_animals", "from", "to");

    int date, num_animals;
    while (exchanges.read_row(date, num_animals, from, to)) {
        sim.add_exchange(date, num_animals, from, to);
    }
    // mio::log_info("Exchanges added");

    // // #ifdef MEMILIO_ENABLE_OPENMP
    // // #pragma omp parallel for
    // //     for (size_t i = 0; i < 100; ++i) {
    // // #endif

    // // mio::log_info("new Simulation created");
    // sim2.infectionrisk = 0.1;

    sim.advance(tmax);

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
    // // mio::log_info("Sum of exchanged sick animals: {}", exchange_results[0]);
    // // mio::log_info("Sum of exchanged animals: {}", exchange_results[1]);

    // auto sth = sim2.exchanges_per_timestep().export_csv("Exchange_statistics.csv", {"Commuter Sick", "Commuter Total"});

    // // for (auto node : sim2.get_graph().nodes()) {
    // //     if (node.property.get_result().get_num_time_points() < num_time_points) {
    // //         mio::log_error("Node with inconsistent number of time points in results.");
    // //     }
    // // }

    // // exchange_results = sim2.sum_nodes();
    // // // mio::log_info("{}, {}, {}, {}", exchange_results[0], exchange_results[1], exchange_results[2], exchange_results[3]);

    // sth = sim2.statistics_per_timestep().export_csv("Simulation_statistics.csv");
    // // // auto combined_results = sim2.combine_node_results();
    // // // combined_results.print_table({"S", "E", "I", "R", "D"});
    // // // auto ioresult = combined_results.export_csv("Simulation_results.csv");

    // sim2.statistics_per_timestep({0, 1, 2, 3, 4}).print_table({"S", "E", "I", "R", "D"});
    // mio::log_info("Finished postprocessing");

    mio::unused(sim.statistics_per_timestep().export_csv(mio::path_join("AsymmetricParams_single_run.csv")));

    return 0;
}
