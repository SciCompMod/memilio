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
#include "memilio/io/io.h"
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

const std::string farm_file     = "/home/kilian/Documents/data/read_data/farms_for_simulations.csv";
const std::string edge_file     = "/home/kilian/Documents/data/read_data/simulation_edges.csv";
const std::string exchange_file = "/home/kilian/Documents/data/read_data/exchanges.csv";

int main(int /*argc*/, char** /*argv*/)
{
    mio::log_debug("Folder: {}", mio::get_current_dir_name());
    const auto t0   = 0.;
    const auto tmax = 100.;
    const auto dt   = 1.; //initial time step

    using mio::fmd::InfectionState;
    using Status = mio::Index<InfectionState>;
    using mio::regions::Region;

    //total compartment sizes

    using Model = mio::smm::Model<ScalarType, InfectionState, Status, Region>;

    auto home = Region(0);

    auto adoption_rates = mio::fmd::generic_adoption_rates();

    mio::fmd::Builder builder;
    mio::log_info("Starting Graph generation");
    {
        mio::timing::AutoTimer<"Graph Nodes Generation"> timer;
        io::CSVReader<4> farms(farm_file);
        farms.read_header(io::ignore_extra_column, "id_numeric", "x", "y", "Tiere");
        int farm_id, num_cows;
        double latitude, longitude;
        while (farms.read_row(farm_id, longitude, latitude, num_cows)) {
            builder.add_node(farm_id, longitude, latitude, mio::fmd::create_model(num_cows, adoption_rates), t0);
        }
    }
    mio::log_info("Nodes added to Graph");
    auto rng = mio::RandomNumberGenerator();

    std::vector<std::vector<size_t>> interesting_indices;
    interesting_indices.push_back(
        {Model(Status{InfectionState::Count}, Region(1)).populations.get_flat_index({home, InfectionState::I})});
    io::CSVReader<2> edges(edge_file);
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
            tree.in_range_indices_query(node.property.get_location(), {mio::geo::kilometers(2)}));
    }

    mio::log_info("Neighbors set");
    auto sim = mio::make_mobility_sim(t0, dt, std::move(graph));

    io::CSVReader<4> exchanges(exchange_file);
    exchanges.read_header(io::ignore_extra_column, "from_dec", "to_dec", "day", "LOM_length");

    int date, num_animals;
    while (exchanges.read_row(from, to, date, num_animals)) {
        sim.add_exchange(date, num_animals, from, to);
    }
    int index_farm = 0;
    // mio::UniformIntDistribution<int>::get_instance()(sim.get_rng(), 0, int(sim.get_graph().nodes().size() - 1));
    auto E_index = sim.get_graph().nodes()[0].property.get_simulation().get_model().populations.get_flat_index(
        {Region(0), InfectionState::E});
    auto S_index = sim.get_graph().nodes()[0].property.get_simulation().get_model().populations.get_flat_index(
        {Region(0), InfectionState::S});
    auto num_sus     = sim.get_graph().nodes()[index_farm].property.get_result().get_last_value()[S_index];
    auto num_exposed = std::ceil(0.01 * num_sus);
    sim.get_graph().nodes()[index_farm].property.get_result().get_last_value()[E_index] = num_exposed;
    sim.get_graph().nodes()[index_farm].property.get_result().get_last_value()[S_index] = num_sus - num_exposed;
    mio::log_info("Infecting {} individuals at node {}.", num_exposed, index_farm);

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
