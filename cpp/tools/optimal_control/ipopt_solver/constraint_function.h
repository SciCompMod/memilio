#pragma once

#include <vector>

#include "tools/optimal_control/control_parameters/damping_controls.h"
#include "tools/optimal_control/constraints/update_constraints.h"

#include "tools/optimal_control/helpers/integrator_selector.h"
#include "tools/optimal_control/helpers/make_time_grid.h"


template <typename FP, class OptimizationSettings>
void constraint_function(const OptimizationSettings& settings, const typename OptimizationSettings::template Graph<FP>& graph,
                         const FP* ptr_parameters, size_t n, FP* ptr_constraints, size_t /*m*/)
{
    // ----------------------------------------------------------------- //
    // Evaluate the constraints on the model.                            //
    // Step 1. Define dampings based on 'FP* ptr_parameters'.            //
    // Step 2. Gather information in 'path_values' in 'terminal_values'. //
    // Step 3. Fill information into 'FP* ptr_constraints'.              //
    // ----------------------------------------------------------------- //
    assert(n == settings.num_control_parameters() * settings.num_control_intervals());
    assert(m == settings.num_path_constraints() + settings.num_terminal_constraints());

    std::vector<FP> parameters(n);
    for (size_t i = 0; i < n; i++) {
        parameters[i] = settings.activation_function(ptr_parameters[i]);
    }

    std::vector<FP> path_constraint_values(settings.num_path_constraints(), 0.0);
    std::vector<FP> terminal_constraint_values(settings.num_terminal_constraints(), 42.0);

    // auto integrator = make_integrator<FP>(settings.integrator_type(), settings.dt());

    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());

    auto graph_copy = graph;
    auto sim = mio::make_mobility_sim<FP>(settings.t0(), settings.dt(), std::move(graph_copy));

    // sim.set_integrator(std::move(integrator));
    
    for (auto& node : sim.get_graph().nodes()) {
        set_control_dampings<FP, OptimizationSettings>(settings, node.property.get_simulation().get_model(), parameters);
    }

    // const auto& final_state = sim.get_result().get_last_value();
    // update_path_constraint<FP, OptimizationSettings>(settings, sim.get_model(), final_state, path_constraint_values);
    for (size_t interval = 0; interval < settings.num_intervals(); interval++) {

        // sim.get_dt() = settings.dt();
        sim.advance(time_steps[interval + 1]);

        for (auto& node : sim.get_graph().nodes()) {
            const auto& final_state = node.property.get_simulation().get_result().get_last_value();

            update_path_constraint<FP, OptimizationSettings>(settings, final_state, path_constraint_values);
        }
    }
    for (auto& node : sim.get_graph().nodes()) {
        const auto& final_state = node.property.get_simulation().get_result().get_last_value();

        update_terminal_constraint<FP, OptimizationSettings>(settings, final_state, terminal_constraint_values);
    }

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
