#pragma once

#include "optimization_settings/optimization_settings.h"
#include "control_parameters/damping_controls.h"
#include "constraints/update_constraints.h"
#include "helpers/integrator_selector.h"
#include "helpers/make_time_grid.h"

#include <cstddef>
#include <filesystem>
#include <memory>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <iomanip>

void check_constraint_feasibility(const SecirvvsOptimization& settings)
{
    // Set most restrictive controls to validate if the constraints allow for a feasible solution.
    std::vector<double> parameters(settings.num_control_parameters() * settings.num_control_intervals());
    for (size_t interval = 0; interval < settings.num_control_intervals(); ++interval) {
        for (size_t i = 0; i < settings.num_control_parameters(); ++i) {
            size_t index      = interval * settings.num_control_parameters() + i;
            parameters[index] = settings.control_parameters()[i].max();
        }
    }

    // Simulate with these controls and evaluate the constraints.
    std::vector<double> path_constraint_values(settings.num_path_constraints(), 0.0);
    std::vector<double> terminal_constraint_values(settings.num_terminal_constraints(), 0.0);

    auto graph_model = settings.optimization_model().get_graph_model<double>();

    // Set the control parameters for each node simulation
    for (auto& node : graph_model.nodes()) {
        set_control_dampings<double>(settings, node.property.get_simulation().get_model(), parameters);
    }
    // Set the integrator for each node simulation
    for (auto& node : graph_model.nodes()) {
        node.property.get_simulation().set_integrator_core(
            std::move(make_integrator<double>(settings.integrator_type(), settings.dt())));
    }

    // Create graph simulation
    auto graph_simulation = mio::make_mobility_sim<double>(settings.t0(), settings.dt(), std::move(graph_model));

    // Step through the time grid and evaluate the path constraints
    std::vector<double> time_steps = make_time_grid<double>(settings.t0(), settings.tmax(), settings.num_intervals());

    // Simulate and update constraints
    update_path_constraint<double>(settings, graph_simulation.get_graph(), path_constraint_values);
    for (size_t interval = 0; interval < settings.num_intervals(); interval++) {
        graph_simulation.advance(time_steps[interval + 1]);
        update_path_constraint<double>(settings, graph_simulation.get_graph(), path_constraint_values);
    }
    update_terminal_constraint<double>(settings, graph_simulation.get_graph(), terminal_constraint_values);

    // Print out the path constraint values
    std::cout << "Lower bound for path constraints:\n";
    for (size_t i = 0; i < settings.num_path_constraints(); ++i) {
        const Constraint& constraint = settings.path_constraints()[i];
        std::cout << "  [" << i << "] \"" << constraint.name() << "\": " << path_constraint_values[i] << "\n";
    }
    // Print out the terminal constraint values
    std::cout << "Lower bound on terminal constraints:\n";
    for (size_t i = 0; i < settings.num_terminal_constraints(); ++i) {
        const Constraint& constraint = settings.terminal_constraints()[i];
        std::cout << "  [" << i << "] \"" << constraint.name() << "\": " << terminal_constraint_values[i] << "\n";
    }

    // ---------------------------- //
    // Check constraint feasibility //
    bool violated = false;
    std::ostringstream report;

    report << std::fixed << std::setprecision(6);
    report << "Constraint feasibility report\n";
    report << "=====================================\n\n";

    // --- Path constraints (check only upper bounds) ---
    for (size_t i = 0; i < settings.num_path_constraints(); ++i) {
        const auto& constraint = settings.path_constraints()[i];
        double value           = path_constraint_values[i];
        double max_v           = constraint.max();

        if (value > max_v) {
            violated = true;
            report << "Path constraint violated: [" << i << "] \"" << constraint.name() << "\" value = " << value
                   << " exceeds upper bound " << max_v << "\n";
        }
    }

    // --- Terminal constraints (check only upper bounds) ---
    for (size_t i = 0; i < settings.num_terminal_constraints(); ++i) {
        const auto& constraint = settings.terminal_constraints()[i];
        double value           = terminal_constraint_values[i];
        double max_v           = constraint.max();

        if (value > max_v) {
            violated = true;
            report << "Terminal constraint violated: [" << i << "] \"" << constraint.name() << "\" value = " << value
                   << " exceeds upper bound " << max_v << "\n";
        }
    }

    // ---------------------------------- //
    // Report and throw error if violated //
    if (violated) {
        std::ofstream log_file("constraint_violation.log");
        log_file << report.str();
        log_file.close();

        std::cerr << "Constraint feasibility check failed. See 'constraint_violation.log' for details.\n";
        throw std::runtime_error("Constraint feasibility check failed. See 'constraint_violation.log'.");
    }
    else {
        std::cout << "All constraints are feasible.\n";
    }
}
