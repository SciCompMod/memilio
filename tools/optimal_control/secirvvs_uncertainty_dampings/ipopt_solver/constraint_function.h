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
void constraint_function(const SecirvvsOptimization& settings, const FP* ptr_parameters, size_t n, FP* ptr_constraints,
                         size_t m)
{
    std::vector<FP> parameters(settings.num_control_parameters() * settings.num_control_intervals());
    for (size_t control_interval = 0; control_interval < settings.num_control_intervals(); control_interval++) {
        for (size_t control_index = 0; control_index < settings.num_control_parameters(); control_index++) {
            size_t idx      = control_interval * settings.num_control_parameters() + control_index;
            parameters[idx] = settings.activation_function()(ptr_parameters[idx]);
        }
    }

    // Store the maximum constraint values over all simulation runs.
    std::vector<FP> max_path_constraints(settings.num_path_constraints(), 0.0);
    std::vector<FP> max_terminal_constraints(settings.num_terminal_constraints(), 0.0);

    for (std::size_t run = 0; run < settings.num_simulation_runs(); run++) {
        // Set seed for this run
        size_t seed      = settings.base_seed() + run;
        auto graph_model = settings.optimization_model().get_graph_model<FP>(seed);

        // Simulate with these controls and evaluate the constraints.
        std::vector<FP> path_constraint_values(settings.num_path_constraints(), 0.0);
        std::vector<FP> terminal_constraint_values(settings.num_terminal_constraints(), 0.0);

        for (auto& node : graph_model.nodes()) {
            set_control_dampings<FP>(settings, node.property.get_simulation().get_model(), parameters);
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

        // Fill max constraint values over all runs
        for (size_t i = 0; i < settings.num_path_constraints(); ++i)
            max_path_constraints[i] = std::max(max_path_constraints[i], path_constraint_values[i]);
        for (size_t i = 0; i < settings.num_terminal_constraints(); ++i)
            max_terminal_constraints[i] = std::max(max_terminal_constraints[i], terminal_constraint_values[i]);
    }

    size_t idx = 0;
    for (size_t path_constraint_index = 0; path_constraint_index < settings.num_path_constraints();
         path_constraint_index++) {
        ptr_constraints[idx] = max_path_constraints[path_constraint_index];
        idx++;
    }
    for (size_t terminal_constraint_index = 0; terminal_constraint_index < settings.num_terminal_constraints();
         terminal_constraint_index++) {
        ptr_constraints[idx] = max_terminal_constraints[terminal_constraint_index];
        idx++;
    }
}
