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
#include "memilio/mobility/farm_graph_simulation.h"
#include "memilio/mobility/farm_simulation.h"
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
#include "thirdparty/csv.h"
#include <ranges>
#include <omp.h>

enum class InfectionState
{
    S,
    E,
    I,
    D,
    Count
};

int main(int /*argc*/, char** /*argv*/)
{
    const auto t0   = 0.;
    const auto tmax = 34.;
    const auto dt   = 1.; //initial time step

    using Status = mio::Index<InfectionState>;
    using mio::regions::Region;
    auto rng = mio::RandomNumberGenerator();
    rng.seed({42, 123});
    //total compartment sizes

    using Model   = mio::smm::Model<ScalarType, InfectionState, Status, Region>;
    using Builder = mio::GraphBuilder<
        mio::FarmNode<ScalarType, mio::smm::Simulation<ScalarType, InfectionState, Status, mio::regions::Region>>,
        mio::MobilityEdgeDirected<ScalarType>>;

    auto home = Region(0);

    using AR = mio::AdoptionRate<ScalarType, Status, Region>;
    std::vector<AR> adoption_rates_0;
    adoption_rates_0.push_back({InfectionState::S, InfectionState::E, home, 1.0, {{InfectionState::I, 0.5}}});
    adoption_rates_0.push_back({InfectionState::E, InfectionState::I, home, 0.3, {}});
    adoption_rates_0.push_back({InfectionState::I, InfectionState::D, home, 0.1, {}});
    adoption_rates_0.push_back({InfectionState::S, InfectionState::E, home, 0.00, {}});

    std::vector<AR> adoption_rates_1;
    adoption_rates_1.push_back({InfectionState::S, InfectionState::E, home, 1.0, {{InfectionState::I, 0.5}}});
    adoption_rates_1.push_back({InfectionState::E, InfectionState::I, home, 0.3, {}});
    adoption_rates_1.push_back({InfectionState::I, InfectionState::D, home, 0.1, {}});
    adoption_rates_1.push_back({InfectionState::S, InfectionState::E, home, 0.00, {}});

    std::vector<AR> adoption_rates_2;
    adoption_rates_2.push_back({InfectionState::S, InfectionState::E, home, 1.0, {{InfectionState::I, 0.5}}});
    adoption_rates_2.push_back({InfectionState::E, InfectionState::I, home, 0.3, {}});
    adoption_rates_2.push_back({InfectionState::I, InfectionState::D, home, 0.1, {}});
    adoption_rates_2.push_back({InfectionState::S, InfectionState::E, home, 0.00, {}});

    std::vector<AR> adoption_rates_3;
    adoption_rates_3.push_back({InfectionState::S, InfectionState::E, home, 1.0, {{InfectionState::I, 0.5}}});
    adoption_rates_3.push_back({InfectionState::E, InfectionState::I, home, 0.3, {}});
    adoption_rates_3.push_back({InfectionState::I, InfectionState::D, home, 0.1, {}});
    adoption_rates_3.push_back({InfectionState::S, InfectionState::E, home, 0.00, {}});

    std::vector<AR> adoption_rates_4;
    adoption_rates_4.push_back({InfectionState::S, InfectionState::E, home, 1.0, {{InfectionState::I, 0.5}}});
    adoption_rates_4.push_back({InfectionState::E, InfectionState::I, home, 0.3, {}});
    adoption_rates_4.push_back({InfectionState::I, InfectionState::D, home, 0.1, {}});
    adoption_rates_4.push_back({InfectionState::S, InfectionState::E, home, 0.00, {}});

    Model curr_model(Status{InfectionState::Count}, mio::regions::Region(1), rng);

    Builder builder;

    io::CSVReader<9> farms("/home/kilian/Documents/projects/jolly/common_data/farms.csv");
    farms.read_header(io::ignore_extra_column, "farm_id", "production_id", "capacity", "population_day",
                      "slaughter_day", "volume", "in_hrz", "x", "y");
    size_t farm_id, type, size, volume;
    double x, y;
    int hrz;
    int population, slaughter;
    while (farms.read_row(farm_id, type, size, population, slaughter, volume, hrz, x, y)) {
        auto model = curr_model;
        if (type == 0)
            model.parameters.get<mio::smm::AdoptionRates<ScalarType, Status, mio::regions::Region>>() =
                adoption_rates_0;
        else if (type == 1)
            model.parameters.get<mio::smm::AdoptionRates<ScalarType, Status, mio::regions::Region>>() =
                adoption_rates_1;
        else if (type == 2)
            model.parameters.get<mio::smm::AdoptionRates<ScalarType, Status, mio::regions::Region>>() =
                adoption_rates_2;
        else if (type == 3)
            model.parameters.get<mio::smm::AdoptionRates<ScalarType, Status, mio::regions::Region>>() =
                adoption_rates_3;
        else if (type == 4)
            model.parameters.get<mio::smm::AdoptionRates<ScalarType, Status, mio::regions::Region>>() =
                adoption_rates_4;

        model.populations[{home, InfectionState::S}] = volume;
        // if (farm_id % 2 == 0)
        //     curr_model.populations[{home, InfectionState::E}] = 1;
        builder.add_node(farm_id, x, y, type, size, population, slaughter, hrz, model, t0);
    }

    std::vector<std::vector<size_t>> interesting_indices;
    interesting_indices.push_back({curr_model.populations.get_flat_index({home, InfectionState::E})});
    interesting_indices.push_back({curr_model.populations.get_flat_index({home, InfectionState::I})});

    auto graph = std::move(builder).build();

    auto nodes = graph.nodes() | std::views::transform([](const auto& node) {
                     return &node.property;
                 });
    auto tree  = mio::geo::RTree(nodes.begin(), nodes.end());

    for (auto& node : graph.nodes()) {
        mio::timing::AutoTimer<"neighbourhood search"> timer;
        node.property.set_regional_neighbors(
            tree.in_range_indices_query(node.property.get_location(),
                                        {mio::geo::kilometers(1), mio::geo::kilometers(10), mio::geo::kilometers(20)}));
    }

    auto sim = mio::make_farm_sim(t0, dt, std::move(graph), rng);

    sim.advance(tmax);

    // sim.return_all_time_series()[7658].print_table({"S", "E", "I", "D"});
    auto result  = sim.get_confirmation_dates();
    auto counter = 0;
    for (auto r : result) {
        if (r >= 0)
            counter++;
    }
    mio::log_info("Number of confirmed farms: {}", counter);

    return 0;
}
