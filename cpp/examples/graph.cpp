/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#include "seir/seir.h"
#include "memilio/mobility/migration.h"
#include "memilio/epidemiology/simulation.h"

#include <iostream>

int main(int argc, char** argv)
{
    const auto t0   = 0.;
    const auto tmax = 10.;
    const auto dt   = 0.5; //time step of migration, daily migration every second step

    epi::SeirModel model;
    model.populations[{epi::Index<epi::SeirInfType>(epi::SeirInfType::S)}] = 10000;
    model.parameters.set<epi::StageTimeIncubationInv>(1);
    model.parameters.get<epi::ContactFrequency>().get_baseline()(0, 0) = 2.7;
    model.parameters.set<epi::StageTimeInfectiousInv>(1);

    //two mostly identical groups
    auto model_group1 = model;
    auto model_group2 = model;
    //some contact restrictions in group 1
    model_group1.parameters.get<epi::ContactFrequency>().add_damping(0.5, epi::SimulationTime(5));
    //infection starts in group 1
    model_group1.populations[{epi::Index<epi::SeirInfType>(epi::SeirInfType::S)}] = 9990;
    model_group1.populations[{epi::Index<epi::SeirInfType>(epi::SeirInfType::E)}] = 10;

    epi::Graph<epi::SimulationNode<epi::Simulation<epi::SeirModel>>, epi::MigrationEdge> g;
    g.add_node(1001, model_group1, t0);
    g.add_node(1002, model_group2, t0);
    g.add_edge(0, 1, Eigen::VectorXd::Constant((size_t)epi::SeirInfType::Count, 0.01));
    g.add_edge(1, 0, Eigen::VectorXd::Constant((size_t)epi::SeirInfType::Count, 0.01));

    auto sim = epi::make_migration_sim(t0, dt, std::move(g));

    sim.advance(tmax);

    return 0;
}
