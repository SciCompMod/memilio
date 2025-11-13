
#pragma once

#include "../optimization_settings/optimization_settings.h"
#include "models/ode_secirvvs/model.h"

#include "../control_parameters/damping_controls.h"
#include "../constraints/update_constraints.h"

#include "../helpers/integrator_selector.h"
#include "../helpers/make_time_grid.h"

#include "../constraints/infection_state_utils.h"

#include <vector>
#include <iostream>

template <typename FP>
void save_solution(
    mio::Graph<mio::SimulationNode<FP, mio::osecirvvs::Simulation<FP>>, mio::MobilityEdge<FP>> graph_model,
    const SecirvvsOptimization& settings, size_t n, const FP* ptr_parameters, const FP* z_L, const FP* z_U, size_t m,
    const FP* ptr_constraints, const FP* lambda, FP obj_value)
{
    std::vector<FP> parameters(settings.num_control_parameters() * settings.num_control_intervals());

    for (size_t control_interval = 0; control_interval < settings.num_control_intervals(); control_interval++) {
        for (size_t control_index = 0; control_index < settings.num_control_parameters(); control_index++) {
            size_t idx = control_interval * settings.num_control_parameters() + control_index;

            parameters[idx] = settings.activation_function()(ptr_parameters[idx]);
        }
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
        file << ", " << cp.name();
    }
    file << "\n";

    size_t num_intervals = settings.num_control_intervals();
    size_t num_controls  = settings.num_control_parameters();
    double dt            = (settings.tmax() - settings.t0()) / num_intervals;

    // Write data rows
    for (size_t i = 0; i < num_intervals; ++i) {
        double time = settings.t0() + i * dt;
        file << time;

        size_t base_index = i * num_controls;
        for (size_t j = 0; j < num_controls; ++j) {
            file << ", " << parameters[base_index + j];
        }
        file << "\n";
    }

    // Repeat last set of parameters at final time
    double final_time = settings.t0() + num_intervals * dt;
    file << final_time;

    size_t last_base_index = (num_intervals - 1) * num_controls;
    for (size_t j = 0; j < num_controls; ++j) {
        file << ", " << parameters[last_base_index + j];
    }
    file << "\n";

    file.close();
}
