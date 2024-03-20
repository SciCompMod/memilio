/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Henrik Zunker
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
#include "ode_secir/model.h"
#include "ode_secir/infection_state.h"
#include "ode_secir/parameters.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"

#include <iostream>

int main()
{
    const auto t0   = 0.;
    const auto tmax = 30.;
    const auto dt   = 0.5; //time step of migration, daily migration every second step

    mio::osecir::Model model(1);
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}] = 10000;
    model.parameters.set<mio::osecir::StartDay>(0);
    model.parameters.set<mio::osecir::Seasonality>(0.2);

    model.parameters.get<mio::osecir::TimeExposed>()            = 3.2;
    model.parameters.get<mio::osecir::TimeInfectedNoSymptoms>() = 2.0;
    model.parameters.get<mio::osecir::TimeInfectedSymptoms>()   = 5.8;
    model.parameters.get<mio::osecir::TimeInfectedSevere>()     = 9.5;
    model.parameters.get<mio::osecir::TimeInfectedCritical>()   = 7.1;

    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()  = 0.1;
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()    = 0.7;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()    = 0.09;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()    = 0.25;
    model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>() = 0.45;
    model.parameters.get<mio::osecir::TestAndTraceCapacity>()              = 35;
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()         = 0.2;
    model.parameters.get<mio::osecir::CriticalPerSevere>()                 = 0.25;
    model.parameters.get<mio::osecir::DeathsPerCritical>()                 = 0.3;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));

    //two mostly identical groups
    auto model_group1 = model;
    auto model_group2 = model;
    //some contact restrictions in model_group1
    mio::ContactMatrixGroup& contact_matrix_m1 = model_group1.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix_m1[0].add_damping(0.7, mio::SimulationTime(15.));

    //infection starts in group 1
    model_group1.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}] = 9990;
    model_group1.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]     = 100;

    mio::Graph<mio::SimulationNode<mio::osecir::Simulation<>>, mio::MigrationEdge> g;
    g.add_node(1001, model_group1, t0);
    g.add_node(1002, model_group2, t0);
    g.add_edge(0, 1, Eigen::VectorXd::Constant((size_t)mio::osecir::InfectionState::Count, 0.1));
    g.add_edge(1, 0, Eigen::VectorXd::Constant((size_t)mio::osecir::InfectionState::Count, 0.1));

    auto sim = mio::make_migration_sim(t0, dt, std::move(g));

    sim.advance(tmax);

    auto& edge_1_0 = sim.get_graph().edges()[1];
    auto& results  = edge_1_0.property.get_migrated();
    results.print_table({"Commuter INS", "Commuter ISy", "Commuter Total"});

    return 0;
}
