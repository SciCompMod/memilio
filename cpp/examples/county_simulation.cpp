/*
* Copyright (C) 2020-2026 MEmilio
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
#include "memilio/epidemiology/age_group.h"
#include "ode_secir/model.h"
#include "ode_secir/infection_state.h"
#include "ode_secir/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_stochastic.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "models/smm/model.h"
#include "thirdparty/csv.h"

#include <iostream>

enum class InfectionState
{
    S,
    E,
    I,
    C, // Culling
    Empty, // Empty
    D, // Dead (only wildlife)
    Count
};
struct Species : public mio::Index<Species> {
    Species(size_t val)
        : Index<Species>(val)
    {
    }
};

int main(int /*argc*/, char** /*argv*/)
{

    using Status = mio::Index<InfectionState, Species>;
    using mio::regions::Region;
    using enum InfectionState;

    /* Example how to run the stochastic metapopulation models with four regions. Within each region we differentiate by
       age groups, species and infection states. The infection states are S, E, C, I, R, D. For the number of age groups
       and species we choose: */
    const size_t num_regions = 1;
    const size_t num_species = 6;
    using Model              = mio::smm::Model<ScalarType, InfectionState, Status, Region>;

    Model model(Status{Count, Species(num_species)}, Region(num_regions));

    for (auto species = 0; species < 6; ++species) {
        model.populations[{Region(0), S, Species(species)}]     = 0.0;
        model.populations[{Region(0), E, Species(species)}]     = 0.0;
        model.populations[{Region(0), I, Species(species)}]     = 0.0;
        model.populations[{Region(0), C, Species(species)}]     = 0.0;
        model.populations[{Region(0), Empty, Species(species)}] = 0.0;
        model.populations[{Region(0), D, Species(species)}]     = 0.0;
    }

    using AR = mio::smm::AdoptionRates<ScalarType, Status, Region>;

    //Set infection state adoption rates. Adoptions only happen within a region.
    AR::Type adoption_rates;
    for (size_t species = 0; species < 5; ++species) {
        adoption_rates.push_back({{S, Species(species)},
                                  {E, Species(species)},
                                  Region(0),
                                  0.1,
                                  {{{E, Species(5)}, 0.2}, {{I, Species(5)}, 0.5}}});
        adoption_rates.push_back({{E, Species(species)}, {I, Species(species)}, Region(0), 1.0 / 5., {}});
        adoption_rates.push_back({{I, Species(species)}, {C, Species(species)}, Region(0), 0.8 / 3., {}});
        adoption_rates.push_back({{C, Species(species)}, {Empty, Species(species)}, Region(0), 0.99 / 5., {}});
        adoption_rates.push_back({{S, Species(species)}, {Empty, Species(species)}, Region(0), 1. / 21, {}});
        if (species != 3) {
            adoption_rates.push_back({{E, Species(species)}, {S, Species(species)}, Region(0), 1. / 22, {}});
        }
    }
    adoption_rates.push_back(
        {{S, Species(5)}, {E, Species(5)}, Region(0), 0.1, {{{E, Species(5)}, 0.2}, {{I, Species(5)}, 0.5}}});
    adoption_rates.push_back({{E, Species(5)}, {I, Species(5)}, Region(0), 1.0 / 5., {}});
    adoption_rates.push_back({{I, Species(5)}, {D, Species(5)}, Region(0), 1. / 5., {}});

    model.parameters.get<AR>() = adoption_rates;
    const auto t0              = 0.;
    const auto tmax            = 90.;
    const auto dt              = 0.1; //initial time step
    const auto num_nodes       = 556;

    mio::Graph<mio::SimulationNode<ScalarType, mio::Simulation<ScalarType, Model>>,
               mio::MobilityEdgeStochastic<ScalarType>>
        graph;

    io::CSVReader<6> population_data("/home/kilian/Documents/projects/jolly/data/farms_by_county_type.csv");
    size_t county_id, conventional, organic, broiler_1, broiler_2, layer;
    population_data.read_header(io::ignore_extra_column, "county_id", "conventional", "broiler_2", "organic",
                                "broiler_1", "layer");
    size_t node_index = 0;
    while (population_data.read_row(county_id, conventional, broiler_2, organic, broiler_1, layer)) {
        while (node_index < county_id) {
            auto model2                                    = model;
            model2.populations[{Region(0), S, Species(6)}] = 10;
            graph.add_node(node_index, model2, t0);
            ++node_index;
        }
        auto model2                                    = model;
        model2.populations[{Region(0), S, Species(0)}] = conventional;
        model2.populations[{Region(0), S, Species(1)}] = organic;
        model2.populations[{Region(0), S, Species(2)}] = broiler_1;
        model2.populations[{Region(0), S, Species(3)}] = broiler_2;
        model2.populations[{Region(0), S, Species(4)}] = layer;
        model2.populations[{Region(0), S, Species(5)}] = 10;
        graph.add_node(county_id, model2, t0);
        ++node_index;
    }
    while (node_index < num_nodes) {
        auto model2                                    = model;
        model2.populations[{Region(0), S, Species(6)}] = 10;
        graph.add_node(node_index, model2, t0);
        ++node_index;
    }

    auto graph_transition_rates = mio::MobilityCoefficients<ScalarType>(model.populations.numel());
    ScalarType kappa            = 0.01;

    auto coeff_idx                                   = model.populations.get_flat_index({Region(0), S, Species(5)});
    graph_transition_rates.get_baseline()[coeff_idx] = 0.02;
    coeff_idx                                        = model.populations.get_flat_index({Region(0), E, Species(5)});
    graph_transition_rates.get_baseline()[coeff_idx] = 0.01;
    coeff_idx                                        = model.populations.get_flat_index({Region(0), I, Species(5)});
    graph_transition_rates.get_baseline()[coeff_idx] = 0.01;

    graph_transition_rates.get_baseline() *= kappa;

    io::CSVReader<2> adjacency("/home/kilian/Documents/projects/jolly/data/county_adjacency.csv");
    size_t county_one, county_two;
    adjacency.read_header(io::ignore_extra_column, "county_id_1", "county_id_2");
    while (adjacency.read_row(county_one, county_two)) {
        graph.add_edge(county_one, county_two, std::move(graph_transition_rates));
        graph.add_edge(county_two, county_one, std::move(graph_transition_rates));
    }

    auto sim = mio::make_mobility_sim<ScalarType>(t0, dt, std::move(graph));
    {
        mio::timing::AutoTimer<"Simulation"> timer;
        sim.advance(tmax);
    }

    return 0;
}
