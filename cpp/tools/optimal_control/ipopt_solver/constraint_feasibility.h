#pragma once

#include <vector>

#include "tools/optimal_control/control_parameters/damping_controls.h"
#include "tools/optimal_control/constraints/update_constraints.h"

#include "tools/optimal_control/helpers/integrator_selector.h"
#include "tools/optimal_control/helpers/make_time_grid.h"


template <class OptimizationSettings>
void check_constraint_feasibility(const OptimizationSettings& settings, const typename OptimizationSettings::template Graph<double>& graph)
{
    std::vector<double> path_constraint_values(settings.num_path_constraints(), 0.0);
    std::vector<double> terminal_constraint_values(settings.num_terminal_constraints(), 0.0);


    auto integrator = make_integrator<double>(settings.integrator_type(), settings.dt());

    std::vector<double> time_steps = make_time_grid<double>(settings.t0(), settings.tmax(), settings.num_intervals());

    std::vector<double> parameters(settings.num_control_parameters() * settings.num_control_intervals());
    for (size_t control_interval = 0; control_interval < settings.num_control_intervals(); control_interval++) {

        mio::SimulationTime<double> time(time_steps[control_interval * settings.pc_resolution()]);
        
        for(size_t idx_control_parameter = 0; idx_control_parameter < settings.num_control_parameters(); idx_control_parameter++)
        {
            parameters[idx_control_parameter + control_interval * settings.num_control_parameters()] =
                settings.control_parameters()[idx_control_parameter].max();
        }
    }

    auto graph_copy = graph;
    auto sim = mio::make_mobility_sim<double>(settings.t0(), settings.dt(), std::move(graph));

    // sim.set_integrator(std::move(integrator));
    
    for (auto& node : sim.get_graph().nodes()) {
        set_control_dampings<double, OptimizationSettings>(settings, node.property.get_simulation().get_model(), parameters);
    }

    // const auto& final_state = sim.get_result().get_last_value();
    // update_path_constraint<double, OptimizationSettings>(settings, sim.get_model(), final_state, path_constraint_values);
    for (size_t interval = 0; interval < settings.num_intervals(); interval++) {

        // sim.get_dt() = settings.dt();
        sim.advance(time_steps[interval + 1]);
        for (auto& node : sim.get_graph().nodes()) {
            const auto& final_state = node.property.get_simulation().get_result().get_last_value();
            update_path_constraint<double, OptimizationSettings>(settings, final_state, path_constraint_values);
        }
    }
    for (auto& node : sim.get_graph().nodes()) {
        const auto& final_state = node.property.get_simulation().get_result().get_last_value();
        update_terminal_constraint<double, OptimizationSettings>(settings, final_state, terminal_constraint_values);
    }

    double epsilon = 1e-3;

    // Validate and clamp path constraint values
    for (size_t i = 0; i < settings.num_path_constraints(); ++i) {
        double path_constraint_value = path_constraint_values[i];

        Constraint& constraint = settings.path_constraints()[i];
        double min_val         = constraint.min();
        double max_val         = constraint.max();

        if (max_val < path_constraint_value) {
            std::cout << "Path Constraint [" << i << "] \"" << constraint.name() << "\": Increasing from " << max_val
                      << " to " << path_constraint_value * (1 + epsilon) << ".\n";

            assert(0.0 <= path_constraint_value);
            constraint.set_range({min_val, path_constraint_value * (1 + epsilon)});
        }
    }

    // Validate and clamp terminal constraint values
    for (size_t i = 0; i < settings.num_terminal_constraints(); ++i) {
        double terminal_constraint_value = terminal_constraint_values[i];

        Constraint& constraint = settings.terminal_constraints()[i];
        double min_val         = constraint.min();
        double max_val         = constraint.max();

        if (max_val < terminal_constraint_value) {
            std::cout << "Terminal Constraint [" << i << "] \"" << constraint.name() << "\": Increasing from "
                      << max_val << " to " << terminal_constraint_value * (1 + epsilon) << ".\n";

            assert(0.0 <= terminal_constraint_value);
            constraint.set_range({min_val, terminal_constraint_value * (1 + epsilon)});
        }
    }
}