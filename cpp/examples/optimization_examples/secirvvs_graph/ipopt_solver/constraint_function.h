#pragma once

#include <vector>

#include "../optimization_settings/secirvvs_optimization.h"
#include "models/ode_secirvvs/model.h"

#include "../control_parameters/damping_controls.h"
#include "../constraints/update_constraints.h"

#include "../helpers/integrator_selector.h"
#include "../helpers/make_time_grid.h"

#include "../constraints/infection_state_utils.h"

template <typename FP>
void constraint_function(
    mio::Graph<mio::SimulationNode<FP, mio::Simulation<FP, mio::osecirvvs::Model<FP>>>, mio::MobilityEdge<FP>>
        graph_model,
    const SecirvvsOptimization& settings, const FP* ptr_parameters, size_t n, FP* ptr_constraints, size_t m)
{
    // ----------------------------------------------------------------- //
    // Evaluate the constraints on the model.                            //
    // Step 1. Define dampings based on 'FP* ptr_parameters'.            //
    // Step 2. Gather information in 'path_values' in 'terminal_values'. //
    // Step 3. Fill information into 'FP* ptr_constraints'.              //
    // ----------------------------------------------------------------- //
    const size_t num_graph_nodes = graph_model.nodes().size();
    assert(n == settings.num_control_parameters() * settings.num_control_intervals() * num_graph_nodes);
    assert(m == settings.num_path_constraints() + settings.num_terminal_constraints());

    std::vector<FP> path_constraint_values(settings.num_path_constraints(), 0.0);
    std::vector<FP> terminal_constraint_values(settings.num_terminal_constraints(), 0.0);

    std::vector<std::vector<FP>> parameters(
        num_graph_nodes, std::vector<FP>(settings.num_control_parameters() * settings.num_control_intervals()));

    for (size_t control_interval = 0; control_interval < settings.num_control_intervals(); control_interval++) {
        for (size_t control_index = 0; control_index < settings.num_control_parameters(); control_index++) {
            size_t idx = control_index + control_interval * settings.num_control_parameters();
            for (size_t node_index = 0; node_index < num_graph_nodes; node_index++) {
                size_t global_idx =
                    idx + node_index * settings.num_control_intervals() * settings.num_control_parameters();
                parameters[node_index][idx] = settings.activation_function()(ptr_parameters[global_idx]);
            }
        }
    }

    set_control_dampings<FP>(settings, graph_model, parameters);

    for (auto& node : graph_model.nodes()) {
        node.property.get_simulation().set_integrator(make_integrator<FP>(settings.integrator_type(), settings.dt()));
    }

    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());
    auto graph_sim_mobility    = mio::make_mobility_sim<FP>(settings.t0(), settings.dt(), graph_model);

    update_path_constraint<FP>(settings, graph_sim_mobility, path_constraint_values);
    for (size_t interval = 0; interval < settings.num_intervals(); interval++) {
        graph_sim_mobility.advance(time_steps[interval + 1]);
        update_path_constraint<FP>(settings, graph_sim_mobility, path_constraint_values);
    }
    update_terminal_constraint<FP>(settings, graph_sim_mobility, terminal_constraint_values);

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
}
