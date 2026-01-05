
#pragma once

#include <vector>
#include <iostream>

#include "../optimization_settings/secir_optimization.h"
#include "models/ode_secir/model.h"

#include "../control_parameters/damping_controls.h"
#include "../constraints/update_constraints.h"

#include "../helpers/integrator_selector.h"
#include "../helpers/make_time_grid.h"

#include "../constraints/infection_state_utils.h"

#include <fstream>

template <typename FP>
void save_solution(mio::Graph<mio::SimulationNode<FP, mio::osecir::Simulation<FP>>, mio::MobilityEdge<FP>> graph_model,
                   const SecirOptimization& settings, size_t n, const FP* ptr_parameters, const FP* z_L, const FP* z_U,
                   size_t m, const FP* ptr_constraints, const FP* lambda, FP obj_value)
{
    const size_t num_graph_nodes = graph_model.nodes().size();

    std::cout << "\nPath Constraints:\n";
    for (size_t i = 0; i < settings.num_path_constraints(); ++i) {
        std::cout << settings.path_constraints()[i].name() << ": " << ptr_constraints[i] << std::endl;
    }
    std::cout << "\nTerminal Constraints:\n";
    for (size_t i = 0; i < settings.num_terminal_constraints(); ++i) {
        std::cout << settings.terminal_constraints()[i].name() << ": "
                  << ptr_constraints[i + settings.num_path_constraints()] << std::endl;
    }

    // Reconstruct the simulation to calculate tests used
    const size_t num_states = 400;
    std::vector<std::pair<FP, FP>> dynamic_NPI_values(settings.num_dynamic_NPI_parameters());
    std::vector<FP> commuter_testing_values(num_states);

    for (const auto& dynamic_NPI : settings.dynamic_NPI_parameters()) {
        size_t index                     = static_cast<size_t>(string_to_control(dynamic_NPI.name()));
        dynamic_NPI_values[index].first  = ptr_parameters[2 * index + 0];
        dynamic_NPI_values[index].second = ptr_parameters[2 * index + 1];
    }
    for (size_t state_index = 0; state_index < num_states; state_index++) {
        commuter_testing_values[state_index] = ptr_parameters[2 * settings.num_dynamic_NPI_parameters() + state_index];
    }

    set_commuter_testing<FP>(settings, graph_model, commuter_testing_values);

    for (auto& node : graph_model.nodes()) {
        node.property.get_simulation().set_integrator_core(
            std::move(make_integrator<FP>(settings.integrator_type(), settings.dt())));
    }

    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.simulation_days());
    auto graph_sim_mobility    = mio::make_mobility_sim<FP>(settings.t0(), settings.dt(), std::move(graph_model));

    for (size_t day = 0; day < settings.simulation_days(); day++) {
        graph_sim_mobility.advance(time_steps[day + 1]);
    }

    // Calculate total tests used globally
    FP total_tests_global = 0.0;
    std::vector<FP> tests_per_node(num_graph_nodes, 0.0);

    const auto& graph = graph_sim_mobility.get_graph();
    size_t node_idx   = 0;

    for (const auto& node_entry : graph.nodes()) {
        const auto& node_simulation = node_entry.property.get_simulation();
        const auto& node_model      = node_simulation.get_model();
        const auto& result          = node_simulation.get_result();

        const auto num_age_groups = node_model.parameters.get_num_groups();
        const FP start_detection  = node_model.parameters.get_start_commuter_detection();
        const FP end_detection    = node_model.parameters.get_end_commuter_detection();

        const size_t num_time_points = result.get_num_time_points();

        FP node_tests = 0.0;

        for (size_t t = 1; t < num_time_points; ++t) {
            FP time = result.get_time(t);
            if (time < start_detection || time >= end_detection) {
                continue;
            }

            Eigen::Ref<const Eigen::VectorX<FP>> current_state  = result.get_value(t);
            Eigen::Ref<const Eigen::VectorX<FP>> previous_state = result.get_value(t - 1);

            for (auto age = mio::AgeGroup(0); age < num_age_groups; ++age) {
                auto INSCi =
                    node_model.populations.get_flat_index(mio::Index<mio::AgeGroup, mio::osecir::InfectionState>(
                        age, mio::osecir::InfectionState::InfectedNoSymptomsConfirmed));
                auto ISyCi =
                    node_model.populations.get_flat_index(mio::Index<mio::AgeGroup, mio::osecir::InfectionState>(
                        age, mio::osecir::InfectionState::InfectedSymptomsConfirmed));

                FP new_tests =
                    (current_state[INSCi] - previous_state[INSCi]) + (current_state[ISyCi] - previous_state[ISyCi]);

                if (new_tests > 0) {
                    node_tests += new_tests;
                }
            }
        }

        tests_per_node[node_idx] = node_tests;
        total_tests_global += node_tests;
        node_idx++;
    }

    // Write to output file
    std::ofstream output_file("solution_tests_used.txt");
    if (!output_file.is_open()) {
        std::cerr << "Error: Could not open output file for writing.\n";
        return;
    }

    output_file << "Total tests used (global): " << total_tests_global << "\n\n";

    for (size_t i = 0; i < num_graph_nodes; ++i) {
        output_file << "Node " << i << ": " << tests_per_node[i] << " tests " << commuter_testing_values[i] << "\n";
    }

    output_file.close();
    std::cout << "\nSolution saved to solution_tests_used.txt\n";
    std::cout << "Total tests used: " << total_tests_global << "\n";
}
