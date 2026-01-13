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
#include "thirdparty/csv.h"
#include <ranges>
#include <omp.h>

int main(int /*argc*/, char** /*argv*/)
{
    const auto t0   = 0.;
    const auto tmax = 1000.;
    const auto dt   = 1.; //initial time step
    using mio::fmd::InfectionState;
    using Status = mio::Index<InfectionState>;
    using mio::regions::Region;

    //total compartment sizes

    using Model = mio::smm::Model<ScalarType, InfectionState, Status, Region>;
    auto home   = Region(0);
    auto S      = InfectionState::S;
    auto E      = InfectionState::E;
    auto I      = InfectionState::I;
    auto INS    = InfectionState::INS;
    auto ICS    = InfectionState::ICS;
    auto R      = InfectionState::R;
    auto D      = InfectionState::D;

    std::vector<mio::AdoptionRate<ScalarType, Status, Region>> adoption_rates;
    // Adoption rates corresponding to our model, paramters are arbitrary
    adoption_rates.push_back({S, E, home, 0.2, {{I, 0.8}, {INS, 0.1}, {ICS, 0.5}}});
    adoption_rates.push_back({E, I, home, 0.2, {}});
    adoption_rates.push_back({I, INS, home, 0.1, {}});
    adoption_rates.push_back({I, ICS, home, 0.1, {}});
    adoption_rates.push_back({ICS, D, home, 0.6, {}});
    adoption_rates.push_back({ICS, R, home, 0.4, {}});
    adoption_rates.push_back({INS, R, home, 0.5, {}});

    mio::fmd::Builder builder;
    mio::log_info("Starting Graph generation");
    {
        mio::timing::AutoTimer<"Graph Nodes Generation"> timer;
        io::CSVReader<4> farms("/home/kilian/Documents/data/read_data/farms.csv");
        farms.read_header(io::ignore_extra_column, "id_dec", "x", "y", "farm_size");
        size_t farm_id, num_cows;
        double latitude, longitude;
        while (farms.read_row(farm_id, latitude, longitude, num_cows)) {
            auto curr_model = mio::fmd::create_model(num_cows, adoption_rates);
            builder.add_node(farm_id, longitude, latitude, curr_model, t0);
        }
    }
    mio::log_info("Nodes added to Graph");
    auto rng = mio::RandomNumberGenerator();

    std::vector<std::vector<size_t>> interesting_indices;
    interesting_indices.push_back(
        {Model(Status{InfectionState::Count}, Region(1)).populations.get_flat_index({home, InfectionState::I})});
    // graph.reserve_edges(262144);
    {
        mio::timing::AutoTimer<"Graph Edges Generation"> timer;
        io::CSVReader<2> edges("/home/kilian/Documents/data/read_data/edges.csv");
        edges.read_header(io::ignore_extra_column, "from", "to");
        size_t from, to;
        while (edges.read_row(from, to)) {
            builder.add_edge(from, to, interesting_indices);
        }
    }
    auto graph = std::move(builder).build();
    mio::log_info("Graph generated");

    mio::fmd::generate_neighbours(graph, {mio::geo::kilometers(2.0)});

    mio::log_info("Neighbors set");
    auto sim = mio::make_mobility_sim(t0, dt, std::move(graph));

    io::CSVReader<4> exchanges("/home/kilian/Documents/data/read_data/exchanges.csv");
    exchanges.read_header(io::ignore_extra_column, "from_dec", "to_dec", "day", "weight");

    size_t date, num_animals;
    size_t from, to;
    while (exchanges.read_row(from, to, date, num_animals)) {
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

    auto sth = sim2.exchanges_per_timestep().export_csv("Exchange_statistics.csv", {"Commuter Sick", "Commuter Total"});

    sth = sim2.statistics_per_timestep().export_csv("Simulation_statistics.csv");

    mio::log_info("Finished postprocessing");

    return 0;
}
