
#pragma once

#include <vector>
#include <iostream>

#include "../optimization_settings/secirvvs_optimization.h"
#include "models/ode_secirvvs/model.h"

#include "../control_parameters/damping_controls.h"
#include "../constraints/update_constraints.h"

#include "../helpers/integrator_selector.h"
#include "../helpers/make_time_grid.h"

#include "../constraints/infection_state_utils.h"

template <typename FP>
void save_solution(
    mio::Graph<mio::SimulationNode<FP, mio::Simulation<FP, mio::osecirvvs::Model<FP>>>, mio::MobilityEdge<FP>>
        graph_model,
    const SecirvvsOptimization& settings, size_t n, const FP* ptr_parameters, const FP* z_L, const FP* z_U, size_t m,
    const FP* ptr_constraints, const FP* lambda, FP obj_value)
{
    const size_t num_graph_nodes = graph_model.nodes().size();

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

    std::cout << "\nPath Constraints:\n";
    for (size_t i = 0; i < settings.num_path_constraints(); ++i) {
        std::cout << settings.path_constraints()[i].name() << ": " << ptr_constraints[i] << std::endl;
    }
    std::cout << "\nTerminal Constraints:\n";
    for (size_t i = 0; i < settings.num_terminal_constraints(); ++i) {
        std::cout << settings.terminal_constraints()[i].name() << ": "
                  << ptr_constraints[i + settings.num_path_constraints()] << std::endl;
    }

    size_t num_intervals = settings.num_control_intervals();
    size_t num_controls  = settings.num_control_parameters();
    FP dt                = (settings.tmax() - settings.t0()) / num_intervals;

    for (size_t node_index = 0; node_index < num_graph_nodes; ++node_index) {
        std::ofstream file("control_parameters_node_" + std::to_string(node_index) + ".csv");
        if (!file.is_open()) {
            std::cerr << "Failed to open file for node " << node_index << "\n";
            continue;
        }

        // Write header
        file << "Time";
        for (const auto& cp : settings.control_parameters()) {
            file << "," << cp.name();
        }
        file << "\n";

        // Write data rows
        for (size_t i = 0; i < num_intervals; ++i) {
            FP time = settings.t0() + i * dt;
            file << time;

            for (size_t j = 0; j < num_controls; ++j) {
                size_t index = j + i * num_controls;
                file << "," << parameters[node_index][index];
            }

            file << "\n";
        }

        file.close();
    }
}
