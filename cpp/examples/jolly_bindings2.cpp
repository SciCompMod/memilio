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
#include <utility>
#include <vector>

enum class InfState
{
    S,
    E,
    I,
    D,
    Count
};

mio::FarmSimulation<mio::Graph<
    mio::FarmNode<ScalarType, mio::smm::Simulation<ScalarType, InfState, mio::Index<InfState>, mio::regions::Region>>,
    mio::MobilityEdgeDirected<ScalarType>>>
simulate_with_init(std::string farm_file, ScalarType tmax, ScalarType dt, ScalarType suspicion_threshold,
                   ScalarType sensitivity, ScalarType h0, ScalarType r0, ScalarType alpha,
                   ScalarType infection_baseline, ScalarType culling_factor, ScalarType A0_SEI, ScalarType A0_EI,
                   ScalarType A0_ID, ScalarType A0_DeathRate, ScalarType A1_SEI, ScalarType A1_EI, ScalarType A1_ID,
                   ScalarType A1_DeathRate, ScalarType A2_SEI, ScalarType A2_EI, ScalarType A2_ID,
                   ScalarType A2_DeathRate, ScalarType A3_SEI, ScalarType A3_EI, ScalarType A3_ID,
                   ScalarType A3_DeathRate, ScalarType A4_SEI, ScalarType A4_EI, ScalarType A4_ID,
                   ScalarType A4_DeathRate, ScalarType foi_inner_factor0, ScalarType foi_outer_factor0,
                   ScalarType foi_inner_factor1, ScalarType foi_outer_factor1, ScalarType foi_inner_factor2,
                   ScalarType foi_outer_factor2, ScalarType foi_inner_factor3, ScalarType foi_outer_factor3,
                   ScalarType foi_inner_factor4, ScalarType foi_outer_factor4, ScalarType damping0, ScalarType damping1,
                   ScalarType damping2, ScalarType damping3, ScalarType damping4, ScalarType first_infection_day,
                   ScalarType second_infection_day, ScalarType third_infection_day, u_int seed)
{
    const auto t0 = 0.;

    using Status = mio::Index<InfState>;
    using mio::regions::Region;
    auto rng = mio::RandomNumberGenerator();
    rng.seed({seed});

    //total compartment sizes
    using Model   = mio::smm::Model<ScalarType, InfState, Status, Region>;
    using Builder = mio::GraphBuilder<
        mio::FarmNode<ScalarType, mio::smm::Simulation<ScalarType, InfState, Status, mio::regions::Region>>,
        mio::MobilityEdgeDirected<ScalarType>>;

    auto home = Region(0);

    using AR = mio::AdoptionRate<ScalarType, Status, Region>;
    //Organic Ducks
    std::vector<AR> adoption_rates_0;
    adoption_rates_0.push_back({InfState::S, InfState::E, home, 1.0, {{InfState::I, A0_SEI}}});
    adoption_rates_0.push_back({InfState::E, InfState::I, home, A0_EI, {}});
    adoption_rates_0.push_back({InfState::I, InfState::D, home, A0_ID, {}});
    adoption_rates_0.push_back({InfState::S, InfState::E, home, 0.00, {}});
    adoption_rates_0.push_back({InfState::S, InfState::D, home, A0_DeathRate, {}});
    adoption_rates_0.push_back({InfState::E, InfState::D, home, A0_DeathRate, {}});
    adoption_rates_0.push_back({InfState::I, InfState::D, home, A0_DeathRate, {}});

    // Conventional Ducks
    std::vector<AR> adoption_rates_1;
    adoption_rates_1.push_back({InfState::S, InfState::E, home, 1.0, {{InfState::I, A1_SEI}}});
    adoption_rates_1.push_back({InfState::E, InfState::I, home, A1_EI, {}});
    adoption_rates_1.push_back({InfState::I, InfState::D, home, A1_ID, {}});
    adoption_rates_1.push_back({InfState::S, InfState::E, home, 0.00, {}});
    adoption_rates_1.push_back({InfState::S, InfState::D, home, A1_DeathRate, {}});
    adoption_rates_1.push_back({InfState::E, InfState::D, home, A1_DeathRate, {}});
    adoption_rates_1.push_back({InfState::I, InfState::D, home, A1_DeathRate, {}});

    // Layer hens
    std::vector<AR> adoption_rates_2;
    adoption_rates_2.push_back({InfState::S, InfState::E, home, 1.0, {{InfState::I, A2_SEI}}});
    adoption_rates_2.push_back({InfState::E, InfState::I, home, A2_EI, {}});
    adoption_rates_2.push_back({InfState::I, InfState::D, home, A2_ID, {}});
    adoption_rates_2.push_back({InfState::S, InfState::E, home, 0.00, {}});
    adoption_rates_2.push_back({InfState::S, InfState::D, home, A2_DeathRate, {}});
    adoption_rates_2.push_back({InfState::E, InfState::D, home, A2_DeathRate, {}});
    adoption_rates_2.push_back({InfState::I, InfState::D, home, A2_DeathRate, {}});

    // Broiler 1
    std::vector<AR> adoption_rates_3;
    adoption_rates_3.push_back({InfState::S, InfState::E, home, 1.0, {{InfState::I, A3_SEI}}});
    adoption_rates_3.push_back({InfState::E, InfState::I, home, A3_EI, {}});
    adoption_rates_3.push_back({InfState::I, InfState::D, home, A3_ID, {}});
    adoption_rates_3.push_back({InfState::S, InfState::E, home, 0.00, {}});
    adoption_rates_3.push_back({InfState::S, InfState::D, home, A3_DeathRate, {}});
    adoption_rates_3.push_back({InfState::E, InfState::D, home, A3_DeathRate, {}});
    adoption_rates_3.push_back({InfState::I, InfState::D, home, A3_DeathRate, {}});

    // Broiler 2
    std::vector<AR> adoption_rates_4;
    adoption_rates_4.push_back({InfState::S, InfState::E, home, 1.0, {{InfState::I, A4_SEI}}});
    adoption_rates_4.push_back({InfState::E, InfState::I, home, A4_EI, {}});
    adoption_rates_4.push_back({InfState::I, InfState::D, home, A4_ID, {}});
    adoption_rates_4.push_back({InfState::S, InfState::E, home, 0.00, {}});
    adoption_rates_4.push_back({InfState::S, InfState::D, home, A4_DeathRate, {}});
    adoption_rates_4.push_back({InfState::E, InfState::D, home, A4_DeathRate, {}});
    adoption_rates_4.push_back({InfState::I, InfState::D, home, A4_DeathRate, {}});

    Model curr_model(Status{InfState::Count}, mio::regions::Region(1), rng);

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

        model.populations[{home, InfState::S}] = volume;
        // if (farm_id % 2 == 0)
        //     curr_model.populations[{home, InfState::E}] = 1;
        builder.add_node(farm_id, x, y, type, size, population, slaughter, hrz, model, t0);
    }

    std::vector<std::vector<size_t>> interesting_indices;
    interesting_indices.push_back({curr_model.populations.get_flat_index({home, InfState::E})});
    interesting_indices.push_back({curr_model.populations.get_flat_index({home, InfState::I})});

    auto graph = std::move(builder).build();

    auto nodes =
        std::ranges::subrange(graph.nodes().begin(), graph.nodes().end()) | std::views::transform([](const auto& node) {
            return &node.property;
        });
    auto tree = mio::geo::RTree(nodes.begin(), nodes.end());

    for (auto& node : graph.nodes()) {
        node.property.set_regional_neighbors(tree.in_range_indices_query(
            node.property.get_location(),
            {mio::geo::kilometers(1), mio::geo::kilometers(3), mio::geo::kilometers(10), mio::geo::kilometers(25)}));
    }

    auto sim = mio::make_farm_sim(t0, dt, std::move(graph), rng);
    sim.set_parameters(suspicion_threshold, sensitivity, h0, r0, alpha, infection_baseline, culling_factor,
                       {first_infection_day, second_infection_day, third_infection_day},
                       {foi_inner_factor0, foi_inner_factor1, foi_inner_factor2, foi_inner_factor3, foi_inner_factor4},
                       {foi_outer_factor0, foi_outer_factor1, foi_outer_factor2, foi_outer_factor3, foi_outer_factor4},
                       {damping0, damping1, damping2, damping3, damping4});
    sim.advance(tmax);

    return sim;
}

mio::FarmSimulation<mio::Graph<
    mio::FarmNode<ScalarType, mio::smm::Simulation<ScalarType, InfState, mio::Index<InfState>, mio::regions::Region>>,
    mio::MobilityEdgeDirected<ScalarType>>>
simulate_continued(
    const mio::FarmSimulation<
        mio::Graph<mio::FarmNode<ScalarType, mio::smm::Simulation<ScalarType, InfState, mio::Index<InfState>,
                                                                  mio::regions::Region>>,
                   mio::MobilityEdgeDirected<ScalarType>>>& sim_to_copy,
    const ScalarType t, ScalarType tmax, ScalarType dt, ScalarType suspicion_threshold, ScalarType sensitivity,
    ScalarType h0, ScalarType r0, ScalarType alpha, ScalarType infection_baseline, ScalarType culling_factor,
    ScalarType A0_SEI, ScalarType A0_EI, ScalarType A0_ID, ScalarType A0_DeathRate, ScalarType A1_SEI, ScalarType A1_EI,
    ScalarType A1_ID, ScalarType A1_DeathRate, ScalarType A2_SEI, ScalarType A2_EI, ScalarType A2_ID,
    ScalarType A2_DeathRate, ScalarType A3_SEI, ScalarType A3_EI, ScalarType A3_ID, ScalarType A3_DeathRate,
    ScalarType A4_SEI, ScalarType A4_EI, ScalarType A4_ID, ScalarType A4_DeathRate, ScalarType foi_inner_factor0,
    ScalarType foi_outer_factor0, ScalarType foi_inner_factor1, ScalarType foi_outer_factor1,
    ScalarType foi_inner_factor2, ScalarType foi_outer_factor2, ScalarType foi_inner_factor3,
    ScalarType foi_outer_factor3, ScalarType foi_inner_factor4, ScalarType foi_outer_factor4, ScalarType damping0,
    ScalarType damping1, ScalarType damping2, ScalarType damping3, ScalarType damping4, ScalarType first_infection_day,
    ScalarType second_infection_day, ScalarType third_infection_day, u_int seed)
{
    using Status = mio::Index<InfState>;
    using mio::regions::Region;
    auto rng = mio::RandomNumberGenerator();
    rng.seed({seed});

    //total compartment sizes
    using Model   = mio::smm::Model<ScalarType, InfState, Status, Region>;
    using Builder = mio::GraphBuilder<
        mio::FarmNode<ScalarType, mio::smm::Simulation<ScalarType, InfState, Status, mio::regions::Region>>,
        mio::MobilityEdgeDirected<ScalarType>>;

    auto home = Region(0);

    using AR = mio::AdoptionRate<ScalarType, Status, Region>;
    //Organic Ducks
    std::vector<AR> adoption_rates_0;
    adoption_rates_0.push_back({InfState::S, InfState::E, home, 1.0, {{InfState::I, A0_SEI}}});
    adoption_rates_0.push_back({InfState::E, InfState::I, home, A0_EI, {}});
    adoption_rates_0.push_back({InfState::I, InfState::D, home, A0_ID, {}});
    adoption_rates_0.push_back({InfState::S, InfState::E, home, 0.00, {}});
    adoption_rates_0.push_back({InfState::S, InfState::D, home, A0_DeathRate, {}});
    adoption_rates_0.push_back({InfState::E, InfState::D, home, A0_DeathRate, {}});
    adoption_rates_0.push_back({InfState::I, InfState::D, home, A0_DeathRate, {}});

    // Conventional Ducks
    std::vector<AR> adoption_rates_1;
    adoption_rates_1.push_back({InfState::S, InfState::E, home, 1.0, {{InfState::I, A1_SEI}}});
    adoption_rates_1.push_back({InfState::E, InfState::I, home, A1_EI, {}});
    adoption_rates_1.push_back({InfState::I, InfState::D, home, A1_ID, {}});
    adoption_rates_1.push_back({InfState::S, InfState::E, home, 0.00, {}});
    adoption_rates_1.push_back({InfState::S, InfState::D, home, A1_DeathRate, {}});
    adoption_rates_1.push_back({InfState::E, InfState::D, home, A1_DeathRate, {}});
    adoption_rates_1.push_back({InfState::I, InfState::D, home, A1_DeathRate, {}});

    // Layer hens
    std::vector<AR> adoption_rates_2;
    adoption_rates_2.push_back({InfState::S, InfState::E, home, 1.0, {{InfState::I, A2_SEI}}});
    adoption_rates_2.push_back({InfState::E, InfState::I, home, A2_EI, {}});
    adoption_rates_2.push_back({InfState::I, InfState::D, home, A2_ID, {}});
    adoption_rates_2.push_back({InfState::S, InfState::E, home, 0.00, {}});
    adoption_rates_2.push_back({InfState::S, InfState::D, home, A2_DeathRate, {}});
    adoption_rates_2.push_back({InfState::E, InfState::D, home, A2_DeathRate, {}});
    adoption_rates_2.push_back({InfState::I, InfState::D, home, A2_DeathRate, {}});

    // Broiler 1
    std::vector<AR> adoption_rates_3;
    adoption_rates_3.push_back({InfState::S, InfState::E, home, 1.0, {{InfState::I, A3_SEI}}});
    adoption_rates_3.push_back({InfState::E, InfState::I, home, A3_EI, {}});
    adoption_rates_3.push_back({InfState::I, InfState::D, home, A3_ID, {}});
    adoption_rates_3.push_back({InfState::S, InfState::E, home, 0.00, {}});
    adoption_rates_3.push_back({InfState::S, InfState::D, home, A3_DeathRate, {}});
    adoption_rates_3.push_back({InfState::E, InfState::D, home, A3_DeathRate, {}});
    adoption_rates_3.push_back({InfState::I, InfState::D, home, A3_DeathRate, {}});

    // Broiler 2
    std::vector<AR> adoption_rates_4;
    adoption_rates_4.push_back({InfState::S, InfState::E, home, 1.0, {{InfState::I, A4_SEI}}});
    adoption_rates_4.push_back({InfState::E, InfState::I, home, A4_EI, {}});
    adoption_rates_4.push_back({InfState::I, InfState::D, home, A4_ID, {}});
    adoption_rates_4.push_back({InfState::S, InfState::E, home, 0.00, {}});
    adoption_rates_4.push_back({InfState::S, InfState::D, home, A4_DeathRate, {}});
    adoption_rates_4.push_back({InfState::E, InfState::D, home, A4_DeathRate, {}});
    adoption_rates_4.push_back({InfState::I, InfState::D, home, A4_DeathRate, {}});

    Model curr_model(Status{InfState::Count}, mio::regions::Region(1), rng);

    Builder builder;

    for (auto& node : sim_to_copy.get_graph().nodes()) {
        auto model = curr_model;
        auto type  = node.property.get_type();
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

        auto last_value                        = node.property.get_result().get_last_value();
        model.populations[{home, InfState::S}] = last_value[model.populations.get_flat_index({home, InfState::S})];
        model.populations[{home, InfState::E}] = last_value[model.populations.get_flat_index({home, InfState::E})];
        model.populations[{home, InfState::I}] = last_value[model.populations.get_flat_index({home, InfState::I})];
        model.populations[{home, InfState::D}] = last_value[model.populations.get_flat_index({home, InfState::D})];
        builder.add_node(node.id, node.property.get_x(), node.property.get_y(), type, node.property.get_farm_size(),
                         node.property.get_population_date(), node.property.get_slaughter_date(),
                         node.property.get_in_hrz(), model, t);

        auto& last_node = builder.get_last_node();
        // Set node properties that were not set via the constructor
        last_node.property.set_date_suspicion(node.property.get_date_suspicion());
        last_node.property.set_date_confirmation(node.property.get_date_confirmation());
        last_node.property.set_quarantined(node.property.is_quarantined());
        last_node.property.set_infection_status(node.property.get_infection_status());
        last_node.property.set_reg_zone_day(node.property.get_reg_zone_day());
        // regional spread rate
        last_node.property.get_simulation()
            .get_model()
            .parameters.get<mio::smm::AdoptionRates<ScalarType, Status, mio::regions::Region>>()[3] =
            node.property.get_simulation()
                .get_model()
                .parameters.get<mio::smm::AdoptionRates<ScalarType, Status, mio::regions::Region>>()[3];
    }

    std::vector<std::vector<size_t>> interesting_indices;
    interesting_indices.push_back({curr_model.populations.get_flat_index({home, InfState::E})});
    interesting_indices.push_back({curr_model.populations.get_flat_index({home, InfState::I})});

    auto graph = std::move(builder).build();

    auto nodes =
        std::ranges::subrange(graph.nodes().begin(), graph.nodes().end()) | std::views::transform([](const auto& node) {
            return &node.property;
        });
    auto tree = mio::geo::RTree(nodes.begin(), nodes.end());

    for (auto& node : graph.nodes()) {
        node.property.set_regional_neighbors(tree.in_range_indices_query(
            node.property.get_location(),
            {mio::geo::kilometers(1), mio::geo::kilometers(3), mio::geo::kilometers(10), mio::geo::kilometers(25)}));
    }

    auto sim = mio::make_farm_sim(t, dt, std::move(graph), rng);
    // Copy graph parameters
    sim.copy_parameters_from_simulation(sim_to_copy);
    sim.set_parameters(suspicion_threshold, sensitivity, h0, r0, alpha, infection_baseline, culling_factor,
                       {first_infection_day, second_infection_day, third_infection_day},
                       {foi_inner_factor0, foi_inner_factor1, foi_inner_factor2, foi_inner_factor3, foi_inner_factor4},
                       {foi_outer_factor0, foi_outer_factor1, foi_outer_factor2, foi_outer_factor3, foi_outer_factor4},
                       {damping0, damping1, damping2, damping3, damping4});

    sim.advance(tmax);

    return sim;
}

#ifndef JOLLY_BINDINGS_SKIP_MAIN
int main()
{
    // Parameters
    ScalarType dt = 1.0, suspicion_threshold = 0.2, sensitivity = 0.95, h0 = 0.0002, r0 = 4000, alpha = 10,
               infection_baseline = 0.0001, culling_factor = 1.0;
    ScalarType A0_SEI = 0.5, A0_EI = 0.3, A0_ID = 0.1, A0_DeathRate = 0.0001;
    ScalarType A1_SEI = 0.5, A1_EI = 0.3, A1_ID = 0.1, A1_DeathRate = 0.0001;
    ScalarType A2_SEI = 0.5, A2_EI = 0.3, A2_ID = 0.1, A2_DeathRate = 0.0001;
    ScalarType A3_SEI = 0.5, A3_EI = 0.3, A3_ID = 0.1, A3_DeathRate = 0.0001;
    ScalarType A4_SEI = 0.5, A4_EI = 0.3, A4_ID = 0.1, A4_DeathRate = 0.0001;
    ScalarType foi_inner_factor0 = 1.0, foi_inner_factor1 = 1.0, foi_inner_factor2 = 1.0, foi_inner_factor3 = 1.0,
               foi_inner_factor4 = 1.0;
    ScalarType foi_outer_factor0 = 1.0, foi_outer_factor1 = 1.0, foi_outer_factor2 = 1.0, foi_outer_factor3 = 1.0,
               foi_outer_factor4 = 1.0;
    ScalarType damping0 = 0.5, damping1 = 1.0, damping2 = 1.0, damping3 = 1.0, damping4 = 1.0;
    ScalarType first_infection_day = 0, second_infection_day = 2, third_infection_day = 2;
    u_int seed = 434;

    const std::string farm_file = "/hpc_data/bick_ju/jolly/farms.csv";
    std::cout << "Running one continous simulation...\n";
    auto sim = simulate_with_init(
        farm_file, 60, dt, suspicion_threshold, sensitivity, h0, r0, alpha, infection_baseline, culling_factor, A0_SEI,
        A0_EI, A0_ID, A0_DeathRate, A1_SEI, A1_EI, A1_ID, A1_DeathRate, A2_SEI, A2_EI, A2_ID, A2_DeathRate, A3_SEI,
        A3_EI, A3_ID, A3_DeathRate, A4_SEI, A4_EI, A4_ID, A4_DeathRate, foi_inner_factor0, foi_outer_factor0,
        foi_inner_factor1, foi_outer_factor1, foi_inner_factor2, foi_outer_factor2, foi_inner_factor3,
        foi_outer_factor3, foi_inner_factor4, foi_outer_factor4, damping0, damping1, damping2, damping3, damping4,
        first_infection_day, second_infection_day, third_infection_day, seed);
    auto result = sim.get_confirmation_dates(60);
    mio::log_info("Number of infected farms: {}", std::count_if(result.begin(), result.end(), [](int r) {
                      return r >= 0;
                  }));
    std::cout << "Running simulation with break...\n";
    double tmax = 10;
    auto sim1   = simulate_with_init(
        farm_file, tmax, dt, suspicion_threshold, sensitivity, h0, r0, alpha, infection_baseline, culling_factor,
        A0_SEI, A0_EI, A0_ID, A0_DeathRate, A1_SEI, A1_EI, A1_ID, A1_DeathRate, A2_SEI, A2_EI, A2_ID, A2_DeathRate,
        A3_SEI, A3_EI, A3_ID, A3_DeathRate, A4_SEI, A4_EI, A4_ID, A4_DeathRate, foi_inner_factor0, foi_outer_factor0,
        foi_inner_factor1, foi_outer_factor1, foi_inner_factor2, foi_outer_factor2, foi_inner_factor3,
        foi_outer_factor3, foi_inner_factor4, foi_outer_factor4, damping0, damping1, damping2, damping3, damping4,
        first_infection_day, second_infection_day, third_infection_day, seed);
    result = sim1.get_confirmation_dates(tmax);
    mio::log_info("Number of infected farms: {}", std::count_if(result.begin(), result.end(), [](int r) {
                      return r >= 0;
                  }));
    auto sim_continue = simulate_continued(
        sim1, tmax, tmax + 50, dt, suspicion_threshold, sensitivity, h0, r0, alpha, infection_baseline, culling_factor,
        A0_SEI, A0_EI, A0_ID, A0_DeathRate, A1_SEI, A1_EI, A1_ID, A1_DeathRate, A2_SEI, A2_EI, A2_ID, A2_DeathRate,
        A3_SEI, A3_EI, A3_ID, A3_DeathRate, A4_SEI, A4_EI, A4_ID, A4_DeathRate, foi_inner_factor0, foi_outer_factor0,
        foi_inner_factor1, foi_outer_factor1, foi_inner_factor2, foi_outer_factor2, foi_inner_factor3,
        foi_outer_factor3, foi_inner_factor4, foi_outer_factor4, damping0, damping1, damping2, damping3, damping4,
        first_infection_day, second_infection_day, third_infection_day, seed);
    result = sim_continue.get_confirmation_dates(tmax + 50);
    mio::log_info("Number of infected farms: {}", std::count_if(result.begin(), result.end(), [](int r) {
                      return r >= 0;
                  }));
    return 0;
}
#endif
