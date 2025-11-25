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
#include <chrono>

enum class InfectionState
{
    S,
    E,
    I,
    INS,
    ICS,
    R,
    V,
    D,
    Count
};

int main(int /*argc*/, char** /*argv*/)
{
    const auto t0 = 0.;
    // const auto tmax = 100.;
    const auto dt = 1.; //initial time step

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
    adoption_rates.push_back({S, E, home, 0.0, {}});
    adoption_rates.push_back({E, I, home, 0.2, {}});
    adoption_rates.push_back({I, INS, home, 0.1, {}});
    adoption_rates.push_back({I, ICS, home, 0.1, {}});
    adoption_rates.push_back({ICS, D, home, 0.6, {}});
    adoption_rates.push_back({ICS, R, home, 0.4, {}});
    adoption_rates.push_back({INS, R, home, 0.5, {}});

    mio::log_info("Reading CSVs");
    std::vector<ScalarType> latitudes, longitudes;
    std::vector<int> farm_ids, num_cows_vec, dates, num_animals_exchanges;
    std::vector<int> froms, tos, from_exchanges, to_exchanges;
    io::CSVReader<4> farms("../../../farms16mio.csv");
    farms.read_header(io::ignore_extra_column, "farms", "num_cows", "latitude", "longitude");
    int farm_id, num_cows;
    double latitude, longitude;
    while (farms.read_row(farm_id, num_cows, latitude, longitude)) {
        farm_ids.push_back(farm_id);
        num_cows_vec.push_back(num_cows);
        latitudes.push_back(latitude);
        longitudes.push_back(longitude);
    }
    io::CSVReader<2> edges("../../../edges16mio.csv");
    edges.read_header(io::ignore_extra_column, "from", "to");
    size_t from, to;
    while (edges.read_row(from, to)) {
        froms.push_back(from);
        tos.push_back(to);
    }
    io::CSVReader<4> exchanges("../../../trade16mio.csv");
    exchanges.read_header(io::ignore_extra_column, "date", "num_animals", "from", "to");
    int date, num_animals;
    while (exchanges.read_row(date, num_animals, from, to)) {
        dates.push_back(date);
        num_animals_exchanges.push_back(num_animals);
        from_exchanges.push_back(from);
        to_exchanges.push_back(to);
    }
    std::vector<std::vector<size_t>> interesting_indices;
    interesting_indices.push_back({Model().populations.get_flat_index({home, InfectionState::I})});

    for (size_t num_edges_benchmark = 1024; num_edges_benchmark < 16000000; num_edges_benchmark *= 2) {
        mio::log_info("Running with {} edges", num_edges_benchmark);

        mio::Graph<mio::LocationNode<ScalarType, mio::smm::Simulation<ScalarType, 1, InfectionState>>,
                   mio::MobilityEdgeDirected<ScalarType>>
            graph;
        mio::GraphBuilder<mio::LocationNode<ScalarType, mio::smm::Simulation<ScalarType, 1, InfectionState>>,
                          mio::MobilityEdgeDirected<ScalarType>>
            builder;

        auto start = std::chrono::high_resolution_clock::now();

        for (size_t i = 0; i < farm_ids.size(); ++i) {
            Model curr_model;
            curr_model.populations[{home, InfectionState::S}]                                = num_cows_vec[i];
            curr_model.populations[{home, InfectionState::E}]                                = 0;
            curr_model.populations[{home, InfectionState::I}]                                = 0;
            curr_model.populations[{home, InfectionState::INS}]                              = 0;
            curr_model.populations[{home, InfectionState::ICS}]                              = 0;
            curr_model.populations[{home, InfectionState::R}]                                = 0;
            curr_model.populations[{home, InfectionState::D}]                                = 0;
            curr_model.parameters.get<mio::smm::AdoptionRates<ScalarType, InfectionState>>() = adoption_rates;
            graph.add_node(farm_ids[i], longitudes[i], latitudes[i], curr_model, t0);
        }
        auto rng = mio::RandomNumberGenerator();

        for (size_t i = 0; i < num_edges_benchmark; ++i) {
            graph.add_edge(froms[i], tos[i], interesting_indices);
        }
        auto end      = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        mio::log_info("Classic graph generation finished for {} edges in {} milliseconds", num_edges_benchmark,
                      duration.count());
        start = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < farm_ids.size(); ++i) {
            Model curr_model;
            curr_model.populations[{home, InfectionState::S}]                                = num_cows_vec[i];
            curr_model.populations[{home, InfectionState::E}]                                = 0;
            curr_model.populations[{home, InfectionState::I}]                                = 0;
            curr_model.populations[{home, InfectionState::INS}]                              = 0;
            curr_model.populations[{home, InfectionState::ICS}]                              = 0;
            curr_model.populations[{home, InfectionState::R}]                                = 0;
            curr_model.populations[{home, InfectionState::D}]                                = 0;
            curr_model.parameters.get<mio::smm::AdoptionRates<ScalarType, InfectionState>>() = adoption_rates;
            builder.add_node(farm_ids[i], longitudes[i], latitudes[i], curr_model, t0);
        }
        for (size_t i = 0; i < num_edges_benchmark; ++i) {
            builder.add_edge(from_exchanges[i], to_exchanges[i], interesting_indices);
        }
        auto graph2 = builder.build(true);
        end         = std::chrono::high_resolution_clock::now();
        duration    = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        mio::log_info("New graph construction finished using {} edges in {} milliseconds", num_edges_benchmark,
                      duration.count());

        // auto nodes = graph.nodes() | std::views::transform([](const auto& node) {
        //                  return &node.property;
        //              });
        // auto tree  = mio::geo::RTree(nodes.begin(), nodes.end());

        // for (auto& node : graph.nodes()) {
        //     node.property.set_regional_neighbors(
        //         tree.in_range_indices_query(node.property.get_location(), {mio::geo::kilometers(2.0)}));
        // }

        auto sim = mio::make_mobility_sim(t0, dt, std::move(graph));

        for (size_t i = 0; i < dates.size(); ++i) {
            sim.add_exchange(dates[i], num_animals_exchanges[i], from_exchanges[i], to_exchanges[i]);
        }

        mio::log_info("Starting and immediately finishing simulation");
    }
    return 0;
}
