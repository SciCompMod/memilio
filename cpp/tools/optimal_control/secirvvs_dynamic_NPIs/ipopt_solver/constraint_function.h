#pragma once

#include "../optimization_settings/optimization_settings.h"
#include "models/ode_secirvvs/model.h"

#include "../control_parameters/damping_controls.h"
#include "../constraints/update_constraints.h"

#include "../helpers/integrator_selector.h"
#include "../helpers/make_time_grid.h"

#include "../constraints/infection_state_utils.h"

#include <vector>

template <typename FP>
void constraint_function(
    mio::Graph<mio::SimulationNode<FP, mio::osecirvvs::Simulation<FP>>, mio::MobilityEdge<FP>> graph_model,
    const SecirvvsOptimization& settings, const FP* ptr_parameters, size_t n, FP* ptr_constraints, size_t m)
{
    // ----------------------------------------------------------------- //
    // Evaluate the constraints on the model.                            //
    // Step 1. Define dampings based on 'FP* ptr_parameters'.            //
    // Step 2. Gather information in 'path_values' in 'terminal_values'. //
    // Step 3. Fill information into 'FP* ptr_constraints'.              //
    // ----------------------------------------------------------------- //
    // assert(n == 2 * settings.num_control_parameters());
    // assert(m == settings.num_path_constraints() + settings.num_terminal_constraints());

    std::vector<FP> path_constraint_values(settings.num_path_constraints(), 0.0);
    std::vector<FP> terminal_constraint_values(settings.num_terminal_constraints(), 0.0);

    std::vector<FP> dynamic_NPI_strengths(settings.num_control_parameters());
    for (size_t control_index = 0; control_index < settings.num_control_parameters(); control_index++) {
        dynamic_NPI_strengths[control_index] = ptr_parameters[control_index];
    }

    for (auto& node : graph_model.nodes()) {
        set_dynamic_NPIs<FP>(settings, node.property.get_simulation().get_model(), dynamic_NPI_strengths);
    }

    for (auto& node : graph_model.nodes()) {
        node.property.get_simulation().set_integrator_core(
            std::move(make_integrator<FP>(settings.integrator_type(), settings.dt())));
    }

    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());
    auto graph_simulation      = mio::make_mobility_sim<FP>(settings.t0(), settings.dt(), graph_model);

    update_path_constraint<FP>(settings, graph_simulation.get_graph(), path_constraint_values);
    for (size_t interval = 0; interval < settings.num_intervals(); interval++) {
        graph_simulation.advance(time_steps[interval + 1]);
        update_path_constraint<FP>(settings, graph_simulation.get_graph(), path_constraint_values);
    }
    update_terminal_constraint<FP>(settings, graph_simulation.get_graph(), terminal_constraint_values);

    size_t idx = 0;
    for (size_t path_constraint_index = 0; path_constraint_index < settings.num_path_constraints();
         path_constraint_index++) {
        ptr_constraints[idx] = path_constraint_values[path_constraint_index];
        idx++;
    }
    for (size_t terminal_constraint_index = 0; terminal_constraint_index < settings.num_terminal_constraints();
         terminal_constraint_index++) {
        ptr_constraints[idx] = terminal_constraint_values[terminal_constraint_index];
        idx++;
    }
}
