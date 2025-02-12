/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Anna Wendler
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
#include "ide_secir/infection_state.h"
#include "ide_secir/model.h"
#include "ide_secir/simulation.h"
#include "memilio/config.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"

int main()
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    const ScalarType t0   = 0.;
    const ScalarType tmax = 10.;

    // Time step of graph simulation, after every time step the simulation is restarted in every node.
    const ScalarType dt_graph = 5.;
    // Time step of IDE solver.
    const ScalarType dt_ide_solver = 1.;

    size_t num_agegroups = 1;

    mio::CustomIndexArray<ScalarType, mio::AgeGroup> N_init =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 10000.);
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> deaths_init =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 13.10462213);

    // Create TimeSeries with num_transitions * num_agegroups elements where initial transitions needed for simulation
    // will be stored.
    size_t num_transitions = (size_t)mio::isecir::InfectionTransition::Count;
    mio::TimeSeries<ScalarType> transitions_init(num_transitions * num_agegroups);

    // Add time points for initialization of transitions.
    Vec vec_init(num_transitions * num_agegroups);
    vec_init[(size_t)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 25.;
    vec_init[(size_t)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 15.;
    vec_init[(size_t)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.;
    vec_init[(size_t)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 4.;
    vec_init[(size_t)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]     = 1.;
    vec_init[(size_t)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 4.;
    vec_init[(size_t)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical]     = 1.;
    vec_init[(size_t)mio::isecir::InfectionTransition::InfectedSevereToRecovered]            = 1.;
    vec_init[(size_t)mio::isecir::InfectionTransition::InfectedCriticalToDead]               = 1.;
    vec_init[(size_t)mio::isecir::InfectionTransition::InfectedCriticalToRecovered]          = 1.;

    vec_init = vec_init * dt_ide_solver;
    // Add initial time point to time series.
    transitions_init.add_time_point(-10, vec_init);
    // Add further time points until time t0.
    while (transitions_init.get_last_time() < t0 - dt_ide_solver / 2.) {
        transitions_init.add_time_point(transitions_init.get_last_time() + dt_ide_solver, vec_init);
    }

    // Initialize model.
    mio::isecir::Model model(std::move(transitions_init), N_init, deaths_init, num_agegroups);
    model.check_constraints(dt_ide_solver);

    // Two identical models.
    mio::isecir::Model model1 = model;
    mio::isecir::Model model2 = model;

    // Set up graph with two nodes and no edges.
    mio::Graph<mio::SimulationNode<mio::isecir::Simulation>, mio::MobilityEdge<ScalarType>> g;
    g.add_node(1001, model1, dt_ide_solver);
    g.add_node(1002, model2, dt_ide_solver);

    // Simulate without any mobility between the nodes.
    auto sim = mio::make_no_mobility_sim(t0, dt_graph, std::move(g));

    sim.advance(tmax);

    // Print results of first node.
    auto& node_0  = sim.get_graph().nodes()[0];
    auto result_0 = node_0.property.get_result();
    result_0.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 16, 8);

    return 0;
}
