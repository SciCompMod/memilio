/*
* Copyright (C) 2020-2026 MEmilio
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
#include "memilio/config.h"
#include "ode_secir/model.h"
#include "ode_secir/infection_state.h"
#include "ode_secir/parameters.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"

#include <iostream>

int main()
{
    mio::set_log_level(mio::LogLevel::warn);
    const auto t0   = 0.;
    const auto tmax = 30.;
    const auto dt   = 0.5; //time step of Mobility, daily Mobility every second step

    const size_t num_groups = 1;
    mio::osecir::Model<ScalarType> model(num_groups);
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}] = 10000;
    model.parameters.set<mio::osecir::StartDay<ScalarType>>(0);
    model.parameters.set<mio::osecir::Seasonality<ScalarType>>(0.2);

    model.parameters.get<mio::osecir::TimeExposed<ScalarType>>()            = 3.2;
    model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>() = 2.0;
    model.parameters.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()   = 5.8;
    model.parameters.get<mio::osecir::TimeInfectedSevere<ScalarType>>()     = 9.5;
    model.parameters.get<mio::osecir::TimeInfectedCritical<ScalarType>>()   = 7.1;

    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()  = 0.1;
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()    = 0.7;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()    = 0.09;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()    = 0.25;
    model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<ScalarType>>() = 0.45;
    model.parameters.get<mio::osecir::TestAndTraceCapacity<ScalarType>>()              = 35;
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()         = 0.2;
    model.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()                 = 0.25;
    model.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()                 = 0.3;

    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        model.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(1, 1, 10));

    //two mostly identical groups
    auto model_group1 = model;
    auto model_group2 = model;
    //some contact restrictions in model_group1
    mio::ContactMatrixGroup<ScalarType>& contact_matrix_m1 =
        model_group1.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix_m1[0].add_damping(0.7, mio::SimulationTime<ScalarType>(15.0));

    //infection starts in group 1
    model_group1.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}] = 9990;
    model_group1.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]     = 100;

    // get indices of INS and ISy compartments.
    std::vector<std::vector<size_t>> indices_save_edges(2);

    // Reserve Space. The multiplication by 2 is necessary because we have the
    // base and the confirmed compartments for each age group.
    for (auto& vec : indices_save_edges) {
        vec.reserve(2 * num_groups);
    }

    // get indices and write them to the vector
    for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(num_groups); ++i) {
        indices_save_edges[0].emplace_back(
            model.populations.get_flat_index({i, mio::osecir::InfectionState::InfectedNoSymptoms}));
        indices_save_edges[0].emplace_back(
            model.populations.get_flat_index({i, mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}));
        indices_save_edges[1].emplace_back(
            model.populations.get_flat_index({i, mio::osecir::InfectionState::InfectedSymptoms}));
        indices_save_edges[1].emplace_back(
            model.populations.get_flat_index({i, mio::osecir::InfectionState::InfectedSymptomsConfirmed}));
    }

    mio::Graph<mio::SimulationNode<ScalarType, mio::osecir::Simulation<ScalarType>>, mio::MobilityEdge<ScalarType>> g;
    g.add_node(1001, model_group1, t0);
    g.add_node(1002, model_group2, t0);
    g.add_edge(0, 1, Eigen::VectorX<ScalarType>::Constant((size_t)mio::osecir::InfectionState::Count, 0.1),
               indices_save_edges);
    g.add_edge(1, 0, Eigen::VectorX<ScalarType>::Constant((size_t)mio::osecir::InfectionState::Count, 0.1),
               indices_save_edges);

    auto sim = mio::make_mobility_sim<ScalarType>(t0, dt, std::move(g));

    sim.advance(tmax);

    auto& edge_1_0 = sim.get_graph().edges()[1];
    auto& results  = edge_1_0.property.get_mobility_results();
    results.print_table({"Commuter INS", "Commuter ISy", "Commuter Total"});

    return 0;
}
