
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

    {
        for (auto& node : graph_model.nodes()) {
            set_control_dampings<FP>(settings, node.property.get_simulation().get_model(), parameters);
        }

        for (auto& node : graph_model.nodes()) {
            node.property.get_simulation().set_integrator_core(
                std::move(make_integrator<FP>(settings.integrator_type(), settings.dt())));
        }

        std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());
        auto graph_simulation      = mio::make_mobility_sim<FP>(settings.t0(), settings.dt(), graph_model);

        std::string filename = "population_time_series.csv";
        std::ofstream csv_file(filename);

        csv_file << "Time";
        for (size_t i = 0; i < num_infection_states(); ++i) {
            csv_file << "," << infection_state_to_string(static_cast<mio::osecirvvs::InfectionState>(i));
        }
        csv_file << "\n";

        {
            // Initial totals at t0
            std::vector<FP> totals(num_infection_states(), FP(0.0));

            // Sum over nodes and age groups
            for (auto& node : graph_simulation.get_graph().nodes()) {

                const auto& model = node.property.get_simulation().get_model();
                const auto& ts    = node.property.get_simulation().get_result();
                auto local_state  = ts.get_last_value();

                size_t groups = model.parameters.get_num_groups();

                for (size_t g = 0; g < groups; g++) {
                    for (size_t s = 0; s < num_infection_states(); s++) {
                        size_t idx = g * num_infection_states() + s;
                        totals[s] += local_state[idx];
                    }
                }
            }

            csv_file << std::fixed << std::setprecision(6) << settings.t0();
            for (auto v : totals)
                csv_file << "," << v;
            csv_file << "\n";
        }

        for (size_t k = 0; k < settings.num_intervals(); k++) {
            graph_simulation.advance(time_steps[k + 1]);
            std::vector<FP> totals(num_infection_states(), FP(0.0));

            // Sum over node models and age groups
            for (auto& node : graph_simulation.get_graph().nodes()) {
                const auto& model = node.property.get_simulation().get_model();
                const auto& ts    = node.property.get_simulation().get_result();
                auto local_state  = ts.get_last_value();

                size_t groups = model.parameters.get_num_groups();

                for (size_t g = 0; g < groups; g++) {
                    for (size_t s = 0; s < num_infection_states(); s++) {
                        size_t idx = g * num_infection_states() + s;
                        totals[s] += local_state[idx];
                    }
                }
            }

            csv_file << std::fixed << std::setprecision(6) << time_steps[k + 1];
            for (auto v : totals)
                csv_file << "," << v;
            csv_file << "\n";
        }
        file.close();
    }
}
