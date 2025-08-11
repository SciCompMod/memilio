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
    // This is an example for running the IDE-SECIR model in a graph. It demonstrates a simple setup for a graph without
    // mobility that consists of two nodes and no edges.

    using Vec = Eigen::VectorX<ScalarType>;

    const ScalarType t0   = 0.;
    const ScalarType tmax = 10.;

    // Set time step size of IDE solver. Note that the time step size will be constant throughout the simulation.
    const ScalarType dt_ide_solver = 1.;

    size_t num_agegroups = 1;

    mio::CustomIndexArray<ScalarType, mio::AgeGroup> N_init =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 10000.);
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> deaths_init =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 13.10462213);

    // Create TimeSeries with num_transitions * num_agegroups elements where initial transitions needed for simulation
    // will be stored. We require values for the transitions for a sufficient number of time points before the start of
    // the simulation to initialize our model.
    size_t num_transitions = (size_t)mio::isecir::InfectionTransition::Count;
    mio::TimeSeries<ScalarType> transitions_init(num_transitions * num_agegroups);

    // Define vector of transitions that will be added to the time points of the TimeSeries transitions_init.
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
    // Multiply vec_init with dt_ide_solver so that within a time interval of length 1, always the above number of
    // individuals are transitioning from one compartment to another, irrespective of the chosen time step size.
    vec_init = vec_init * dt_ide_solver;

    // In this example, we use the default transition distributions. For these distributions, setting the initial time
    // point of the TimeSeries transitions_init at time -10 will give us a sufficient number of time points before t0=0.
    // For more information on this, we refer to the documentation of TransitionDistributions in
    // models/ide_secir/parameters.h.
    transitions_init.add_time_point(-10, vec_init);
    // Add further time points with distance dt_ide_solver until time t0.
    while (transitions_init.get_last_time() < t0 - dt_ide_solver / 2.) {
        transitions_init.add_time_point(transitions_init.get_last_time() + dt_ide_solver, vec_init);
    }

    // Initialize IDE model that will be used in in each node in the graph below.
    mio::isecir::Model model(std::move(transitions_init), N_init, deaths_init, num_agegroups);
    model.check_constraints(dt_ide_solver);

    // Set up graph with two nodes and no edges. To each node, we pass an id, the above constructed IDE model as well
    // as the time step size that will be used by the numerical solver. For simplicity, we use the same model in
    // both nodes.
    mio::Graph<mio::SimulationNode<ScalarType, mio::isecir::Simulation>, mio::MobilityEdge<ScalarType>> g;
    g.add_node(1001, model, dt_ide_solver);
    g.add_node(1002, model, dt_ide_solver);

    // We use make_no_mobilty_sim to create a GraphSimulation for the model graph we defined above. This function allows
    // us to create a simulation without having to define mobility between nodes. Any edges would be effectively
    // ignored by the simulation.
    auto sim = mio::make_no_mobility_sim(t0, std::move(g));
    // Run the simulation. This advances all graph nodes independently.
    sim.advance(tmax);

    // Print results of first node.
    auto& node_0  = sim.get_graph().nodes()[0];
    auto result_0 = node_0.property.get_result();
    result_0.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 16, 8);

    return 0;
}
