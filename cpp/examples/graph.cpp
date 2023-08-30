/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"

#include <iostream>

int main(int argc, char** argv)
{
    const auto t0   = 0.;
    const auto tmax = 10.;
    const auto dt   = 0.5; //time step of migration, daily migration every second step

    mio::oseir::Model model;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Susceptible)}] = 10000;
    model.parameters.set<mio::oseir::TimeExposed>(1);
    model.parameters.get<mio::oseir::ContactPatterns>().get_baseline()(0, 0) = 2.7;
    model.parameters.set<mio::oseir::TimeInfected>(1);

    //two mostly identical groups
    auto model_group1 = model;
    auto model_group2 = model;
    //some contact restrictions in group 1
    model_group1.parameters.get<mio::oseir::ContactPatterns>().add_damping(0.5, mio::SimulationTime(5));
    //infection starts in group 1
    model_group1.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Susceptible)}] = 9990;
    model_group1.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}]     = 10;

    mio::Graph<mio::SimulationNode<mio::Simulation<mio::oseir::Model>>, mio::MigrationEdge> g;
    g.add_node(1001, model_group1, t0);
    g.add_node(1002, model_group2, t0);
    g.add_edge(0, 1, Eigen::VectorXd::Constant((size_t)mio::oseir::InfectionState::Count, 0.01));
    g.add_edge(1, 0, Eigen::VectorXd::Constant((size_t)mio::oseir::InfectionState::Count, 0.01));

    auto sim = mio::make_migration_sim(t0, dt, std::move(g));

    sim.advance(tmax);

    return 0;
}
