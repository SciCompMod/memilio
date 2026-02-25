/* 
* Copyright (C) 2020-2026 MEmilio
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
#include <string>
#include <omp.h>
#include <vector>

enum class InfectionState
{
    S,
    E,
    I,
    D,
    Count
};

std::vector<int>
simulate(std::string farm_file, ScalarType tmax = 40, ScalarType dt = 1.0, ScalarType suspicion_threshold = 0.05,
         ScalarType sensitivity = 0.95, ScalarType h0 = 0.001, ScalarType r0 = 1000, ScalarType alpha = 2,
         ScalarType A0_SEI = 0.2, ScalarType A0_EI = 0.2, ScalarType A0_ID = 0.2, ScalarType A0_DeathRate = 0.000128,
         ScalarType A1_SEI = 0.2, ScalarType A1_EI = 0.2, ScalarType A1_ID = 0.2, ScalarType A1_DeathRate = 0.000128,
         ScalarType A2_SEI = 0.2, ScalarType A2_EI = 0.2, ScalarType A2_ID = 0.2, ScalarType A2_DeathRate = 0.000588,
         ScalarType A3_SEI = 0.2, ScalarType A3_EI = 0.2, ScalarType A3_ID = 0.2, ScalarType A3_DeathRate = 0.000588,
         ScalarType A4_SEI = 0.2, ScalarType A4_EI = 0.2, ScalarType A4_ID = 0.2, ScalarType A4_DeathRate = 0.000588,
         ScalarType foi_inner_factor0 = 1.0, ScalarType foi_outer_factor0 = 1.0, ScalarType foi_inner_factor1 = 1.0,
         ScalarType foi_outer_factor1 = 1.0, ScalarType foi_inner_factor2 = 1.0, ScalarType foi_outer_factor2 = 1.0,
         ScalarType foi_inner_factor3 = 1.0, ScalarType foi_outer_factor3 = 1.0, ScalarType foi_inner_factor4 = 1.0,
         ScalarType foi_outer_factor4 = 1.0, ScalarType first_infection_day = 0, ScalarType second_infection_day = 2,
         ScalarType third_infection_day = 2, u_int seed = 42)
{
    const auto t0 = 0.;

    using Status = mio::Index<InfectionState>;
    using mio::regions::Region;
    auto rng = mio::RandomNumberGenerator();
    rng.seed({seed});

    //total compartment sizes
    using Model   = mio::smm::Model<ScalarType, InfectionState, Status, Region>;
    using Builder = mio::GraphBuilder<
        mio::FarmNode<ScalarType, mio::smm::Simulation<ScalarType, InfectionState, Status, mio::regions::Region>>,
        mio::MobilityEdgeDirected<ScalarType>>;

    auto home = Region(0);

    using AR = mio::AdoptionRate<ScalarType, Status, Region>;
    //Organic Ducks
    std::vector<AR> adoption_rates_0;
    adoption_rates_0.push_back({InfectionState::S, InfectionState::E, home, 1.0, {{InfectionState::I, A0_SEI}}});
    adoption_rates_0.push_back({InfectionState::E, InfectionState::I, home, A0_EI, {}});
    adoption_rates_0.push_back({InfectionState::I, InfectionState::D, home, A0_ID, {}});
    adoption_rates_0.push_back({InfectionState::S, InfectionState::E, home, 0.00, {}});
    adoption_rates_0.push_back({InfectionState::S, InfectionState::D, home, A0_DeathRate, {}});
    adoption_rates_0.push_back({InfectionState::E, InfectionState::D, home, A0_DeathRate, {}});
    adoption_rates_0.push_back({InfectionState::I, InfectionState::D, home, A0_DeathRate, {}});

    // Conventional Ducks
    std::vector<AR> adoption_rates_1;
    adoption_rates_1.push_back({InfectionState::S, InfectionState::E, home, 1.0, {{InfectionState::I, A1_SEI}}});
    adoption_rates_1.push_back({InfectionState::E, InfectionState::I, home, A1_EI, {}});
    adoption_rates_1.push_back({InfectionState::I, InfectionState::D, home, A1_ID, {}});
    adoption_rates_1.push_back({InfectionState::S, InfectionState::E, home, 0.00, {}});
    adoption_rates_1.push_back({InfectionState::S, InfectionState::D, home, A1_DeathRate, {}});
    adoption_rates_1.push_back({InfectionState::E, InfectionState::D, home, A1_DeathRate, {}});
    adoption_rates_1.push_back({InfectionState::I, InfectionState::D, home, A1_DeathRate, {}});

    // Layer hens
    std::vector<AR> adoption_rates_2;
    adoption_rates_2.push_back({InfectionState::S, InfectionState::E, home, 1.0, {{InfectionState::I, A2_SEI}}});
    adoption_rates_2.push_back({InfectionState::E, InfectionState::I, home, A2_EI, {}});
    adoption_rates_2.push_back({InfectionState::I, InfectionState::D, home, A2_ID, {}});
    adoption_rates_2.push_back({InfectionState::S, InfectionState::E, home, 0.00, {}});
    adoption_rates_2.push_back({InfectionState::S, InfectionState::D, home, A2_DeathRate, {}});
    adoption_rates_2.push_back({InfectionState::E, InfectionState::D, home, A2_DeathRate, {}});
    adoption_rates_2.push_back({InfectionState::I, InfectionState::D, home, A2_DeathRate, {}});

    // Broiler 1
    std::vector<AR> adoption_rates_3;
    adoption_rates_3.push_back({InfectionState::S, InfectionState::E, home, 1.0, {{InfectionState::I, A3_SEI}}});
    adoption_rates_3.push_back({InfectionState::E, InfectionState::I, home, A3_EI, {}});
    adoption_rates_3.push_back({InfectionState::I, InfectionState::D, home, A3_ID, {}});
    adoption_rates_3.push_back({InfectionState::S, InfectionState::E, home, 0.00, {}});
    adoption_rates_3.push_back({InfectionState::S, InfectionState::D, home, A3_DeathRate, {}});
    adoption_rates_3.push_back({InfectionState::E, InfectionState::D, home, A3_DeathRate, {}});
    adoption_rates_3.push_back({InfectionState::I, InfectionState::D, home, A3_DeathRate, {}});

    // Broiler 2
    std::vector<AR> adoption_rates_4;
    adoption_rates_4.push_back({InfectionState::S, InfectionState::E, home, 1.0, {{InfectionState::I, A4_SEI}}});
    adoption_rates_4.push_back({InfectionState::E, InfectionState::I, home, A4_EI, {}});
    adoption_rates_4.push_back({InfectionState::I, InfectionState::D, home, A4_ID, {}});
    adoption_rates_4.push_back({InfectionState::S, InfectionState::E, home, 0.00, {}});
    adoption_rates_4.push_back({InfectionState::S, InfectionState::D, home, A4_DeathRate, {}});
    adoption_rates_4.push_back({InfectionState::E, InfectionState::D, home, A4_DeathRate, {}});
    adoption_rates_4.push_back({InfectionState::I, InfectionState::D, home, A4_DeathRate, {}});

    Model curr_model(Status{InfectionState::Count}, mio::regions::Region(1), rng);

    Builder builder;

    io::CSVReader<9> farms(farm_file);
    farms.read_header(io::ignore_extra_column, "farm_id", "production_id", "capacity", "population_day",
                      "slaughter_day", "volume", "in_hrz", "x", "y");
    size_t farm_id, type, size, volume;
    double x, y;
    int population, slaughter, hrz;
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

    auto graph                                                    = std::move(builder).build();
    graph.nodes()[7658].property.get_result().get_last_value()[1] = 1;

    auto nodes = graph.nodes() | std::views::transform([](const auto& node) {
                     return &node.property;
                 });
    auto tree  = mio::geo::RTree(nodes.begin(), nodes.end());

    for (auto& node : graph.nodes()) {
        node.property.set_regional_neighbors(
            tree.in_range_indices_query(node.property.get_location(),
                                        {mio::geo::kilometers(1), mio::geo::kilometers(10), mio::geo::kilometers(25)}));
    }

    auto sim = mio::make_farm_sim(t0, dt, std::move(graph), rng);
    sim.set_parameters(suspicion_threshold, sensitivity, h0, r0, alpha,
                       {first_infection_day, second_infection_day, third_infection_day},
                       {foi_inner_factor0, foi_inner_factor1, foi_inner_factor2, foi_inner_factor3, foi_inner_factor4},
                       {foi_outer_factor0, foi_outer_factor1, foi_outer_factor2, foi_outer_factor3, foi_outer_factor4});

    sim.advance(tmax);

    return sim.get_confirmation_dates(tmax);
}

#ifndef JOLLY_BINDINGS_SKIP_MAIN
int main()
{
    const std::string farm_file = "/home/kilian/Documents/projects/jolly/common_data/farms.csv";
    auto result                 = simulate(farm_file);
    mio::log_info("Number of infected farms: {}", std::count_if(result.begin(), result.end(), [](int r) {
                      return r >= 0;
                  }));
    return 0;
}
#endif