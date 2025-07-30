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
#include "smm/simulation.h"
#include "smm/parameters.h"

#include <iostream>

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
    const auto tmax = 10.;
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

    // std::vector<mio::smm::TransitionRate<InfectionState>> transition_rates;
    // for (size_t s = 0; s < static_cast<size_t>(InfectionState::Count); ++s) {
    //     transition_rates.push_back({InfectionState(s), home, home, 0});
    // }
    // model.parameters.get<mio::smm::TransitionRates<InfectionState>>() = transition_rates;
    //modify model for second node
    auto model2 = model;

    mio::Graph<mio::SimulationNode<mio::smm::Simulation<1, InfectionState>>, mio::MobilityEdgeDirected> graph;
    graph.add_node(0, model, t0);

    auto param = mio::MobilityParametersTimed(2, 10, 1);

    graph.add_edge(0, 1, param);

    auto sim = mio::make_mobility_sim(t0, dt, std::move(graph));

    sim.advance(tmax);

    return 0;
}
