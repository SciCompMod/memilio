
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
    std::vector<FP> dynamic_NPI_strengths(settings.num_control_parameters());
    for (size_t control_index = 0; control_index < settings.num_control_parameters(); control_index++) {
        dynamic_NPI_strengths[control_index] = ptr_parameters[control_index];
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

    // Print optimal control parameters
    std::cout << "\nOptimal Control Parameters:\n";
    size_t index = 0; // index into dynamic_NPI_strengths
    for (auto& cp : settings.control_parameters()) {
        std::cout << cp.name() << ", " << "Threshold: " << cp.threshold() << "\n";

        // Print each damping with its corresponding name
        const auto& damping_names = cp.damping_names();
        const auto& dampings      = cp.dampings();

        for (size_t d = 0; d < dampings.size(); ++d) {
            FP strength = std::max(dynamic_NPI_strengths[index++], static_cast<FP>(0));
            double tol  = 1e-7;
            if (strength < tol) {
                strength = 0.0;
            }

            std::cout << "    Damping: " << damping_names[d] << ", Strength: " << strength << "\n";
        }
    }
    // // Write parameters to CSV file
    // std::ofstream file("control_parameters.csv");
    // if (!file.is_open()) {
    //     std::cerr << "Failed to open control_parameters.csv for writing.\n";
    //     return;
    // }
    // // Write header
    // file << "Name, Threshold, Strength\n";
    // for (size_t i = 0; i < settings.num_control_parameters(); ++i) {
    //     const auto& cp = settings.control_parameters()[i];
    //     file << cp.name() << ", " << dynamic_NPI_strengths[i] << "\n";
    // }

    // file.close();

    // // Write parameters to CSV file
    // std::ofstream file("control_parameters.csv");
    // if (!file.is_open()) {
    //     std::cerr << "Failed to open control_parameters.csv for writing.\n";
    //     return;
    // }

    // // Write header
    // file << "Time";
    // for (const auto& cp : settings.control_parameters()) {
    //     file << ", " << cp.name();
    // }
    // file << "\n";

    // size_t num_intervals = settings.num_control_intervals();
    // size_t num_controls  = settings.num_control_parameters();
    // double dt            = (settings.tmax() - settings.t0()) / num_intervals;

    // // Write data rows
    // for (size_t i = 0; i < num_intervals; ++i) {
    //     double time = settings.t0() + i * dt;
    //     file << time;

    //     size_t base_index = i * num_controls;
    //     for (size_t j = 0; j < num_controls; ++j) {
    //         file << ", " << parameters[base_index + j];
    //     }
    //     file << "\n";
    // }

    // // Repeat last set of parameters at final time
    // double final_time = settings.t0() + num_intervals * dt;
    // file << final_time;

    // size_t last_base_index = (num_intervals - 1) * num_controls;
    // for (size_t j = 0; j < num_controls; ++j) {
    //     file << ", " << parameters[last_base_index + j];
    // }
    // file << "\n";

    // file.close();
}
