
#pragma once

#include <vector>
#include <iostream>

#include "tools/optimal_control/control_parameters/damping_controls.h"

#include "tools/optimal_control/helpers/integrator_selector.h"
#include "tools/optimal_control/helpers/make_time_grid.h"

#include "tools/optimal_control/constraints/infection_state_utils.h"

template <typename FP, class OptimizationSettings>
void save_solution(const OptimizationSettings& settings, const typename OptimizationSettings::template ModelTemplate<FP>& model, size_t n,
                   const FP* ptr_parameters, const FP* z_L, const FP* z_U, size_t m, const FP* ptr_constraints,
                   const FP* lambda, FP obj_value)
{
    std::vector<FP> parameters(n);
    for (size_t i = 0; i < n; i++) {
        parameters[i] = settings.activation_function()(ptr_parameters[i]);
    }

    std::cout << "\nPath Constraints:\n";
    for (size_t i = 0; i < settings.num_path_constraints(); ++i) {
        std::cout << settings.path_constraints()[i].name() << ": " << ptr_constraints[i] << std::endl;
    }
    std::cout << "\nTerminal Constraints:\n";
    for (size_t i = 0; i < settings.num_terminal_constraints(); ++i) {
        std::cout << settings.terminal_constraints()[i].name() << ": "
                  << ptr_constraints[i + settings.num_path_constraints()] << std::endl;
    }

    // Write parameters to CSV file
    std::ofstream file("control_parameters.csv");
    if (!file.is_open()) {
        std::cerr << "Failed to open control_parameters.csv for writing.\n";
        return;
    }

    // Write header
    file << "Time";
    for (const auto& cp : settings.control_parameters()) {
        file << "," << cp.name();
    }
    file << "\n";

    // Write data rows
    size_t num_intervals = settings.num_control_intervals();
    size_t num_controls  = settings.num_control_parameters();
    double dt            = (settings.tmax() - settings.t0()) / num_intervals;

    for (size_t i = 0; i < num_intervals; ++i) {
        double time = settings.t0() + i * dt;
        file << time;

        for (size_t j = 0; j < num_controls; ++j) {
            size_t index = j + i * num_controls;
            file << "," << parameters[index];
        }

        file << "\n";
    }

    file.close();
}
