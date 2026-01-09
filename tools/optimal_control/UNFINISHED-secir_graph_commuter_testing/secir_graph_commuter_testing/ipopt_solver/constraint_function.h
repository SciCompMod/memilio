#pragma once

#include <vector>

#include "../optimization_settings/secir_optimization.h"
#include "models/ode_secir/model.h"

#include "../control_parameters/damping_controls.h"
#include "../constraints/update_constraints.h"

#include "../helpers/integrator_selector.h"
#include "../helpers/make_time_grid.h"

#include "../constraints/infection_state_utils.h"

template <typename FP>
void constraint_function(
    mio::Graph<mio::SimulationNode<FP, mio::osecir::Simulation<FP>>, mio::MobilityEdge<FP>> graph_model,
    const SecirOptimization& settings, const FP* ptr_parameters, size_t n, FP* ptr_constraints, size_t m)
{
    // ----------------------------------------------------------------- //
    // Evaluate the constraints on the model.                            //
    // Step 1. Define dampings based on 'FP* ptr_parameters'.            //
    // Step 2. Gather information in 'path_values' in 'terminal_values'. //
    // Step 3. Fill information into 'FP* ptr_constraints'.              //
    // ----------------------------------------------------------------- //
    const size_t num_states      = 400;
    const size_t num_graph_nodes = graph_model.nodes().size();
    // assert(n == 2 * settings.num_dynamic_NPI_parameters() + graph_model.nodes().size());
    assert(n == 2 * settings.num_dynamic_NPI_parameters() + num_states);
    assert(m == settings.num_path_constraints() + settings.num_terminal_constraints());

    std::vector<FP> path_constraint_values(settings.num_path_constraints(), 0.0);
    std::vector<FP> terminal_constraint_values(settings.num_terminal_constraints(), 0.0);

    std::vector<std::pair<FP, FP>> dynamic_NPI_values(settings.num_dynamic_NPI_parameters());
    std::vector<FP> commuter_testing_values(num_states);

    for (const auto& dynamic_NPI : settings.dynamic_NPI_parameters()) {
        size_t index                     = static_cast<size_t>(string_to_control(dynamic_NPI.name()));
        dynamic_NPI_values[index].first  = ptr_parameters[2 * index + 0]; // threshold
        dynamic_NPI_values[index].second = ptr_parameters[2 * index + 1]; // strength
    }
    // for (size_t graph_node = 0; graph_node < num_graph_nodes; graph_node++) {
    //     commuter_testing_values[graph_node] = ptr_parameters[2 * settings.num_dynamic_NPI_parameters() + graph_node];
    // }
    for (size_t state_index = 0; state_index < num_states; state_index++) {
        commuter_testing_values[state_index] = ptr_parameters[2 * settings.num_dynamic_NPI_parameters() + state_index];
    }

    // set_dynamic_NPIs<FP>(settings, graph_model, dynamic_NPI_values);
    set_commuter_testing<FP>(settings, graph_model, commuter_testing_values);

    for (auto& node : graph_model.nodes()) {
        node.property.get_simulation().set_integrator_core(
            std::move(make_integrator<FP>(settings.integrator_type(), settings.dt())));
    }

    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.simulation_days());
    auto graph_sim_mobility    = mio::make_mobility_sim<FP>(settings.t0(), settings.dt(), std::move(graph_model));

    update_path_constraint(settings, graph_sim_mobility, path_constraint_values);
    for (size_t day = 0; day < settings.simulation_days(); day++) {
        graph_sim_mobility.advance(time_steps[day + 1]);
        update_path_constraint(settings, graph_sim_mobility, path_constraint_values);
    }
    update_terminal_constraint(settings, graph_sim_mobility, terminal_constraint_values);

    size_t idx = 0;
    for (size_t path_constraint_index = 0; path_constraint_index < settings.num_path_constraints();
         path_constraint_index++) {
        ptr_constraints[idx] = path_constraint_values[path_constraint_index];
        ptr_constraints[idx] = path_constraint_values[path_constraint_index];
        idx++;
    }
    for (size_t terminal_constraint_index = 0; terminal_constraint_index < settings.num_terminal_constraints();
         terminal_constraint_index++) {
        ptr_constraints[idx] = terminal_constraint_values[terminal_constraint_index];
        ptr_constraints[idx] = terminal_constraint_values[terminal_constraint_index];
        idx++;
    }
    ptr_constraints[idx] = get_tests_used(settings, graph_sim_mobility);
}
