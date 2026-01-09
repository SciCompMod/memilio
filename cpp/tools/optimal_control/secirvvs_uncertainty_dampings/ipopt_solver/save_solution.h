
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
void save_solution(const SecirvvsOptimization& settings, size_t n, const FP* ptr_parameters, const FP* z_L,
                   const FP* z_U, size_t m, const FP* ptr_constraints, const FP* lambda, FP obj_value)
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
    FP dt                = (settings.tmax() - settings.t0()) / num_intervals;

    // Write data rows
    for (size_t i = 0; i < num_intervals; ++i) {
        FP time = settings.t0() + i * dt;
        file << time;

        size_t base_index = i * num_controls;
        for (size_t j = 0; j < num_controls; ++j) {
            file << ", " << parameters[base_index + j];
        }
        file << "\n";
    }

    // Repeat last set of parameters at final time
    FP final_time = settings.t0() + num_intervals * dt;
    file << final_time;

    size_t last_base_index = (num_intervals - 1) * num_controls;
    for (size_t j = 0; j < num_controls; ++j) {
        file << ", " << parameters[last_base_index + j];
    }
    file << "\n";

    file.close();
    // Precompute time grid once
    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());

    // --- 2. Loop over all simulation runs (uncertainty) ---
    for (size_t run = 0; run < settings.num_simulation_runs(); run++) {

        size_t seed = settings.base_seed() + run;

        // --- 2a. Construct graph model for this run ---
        auto graph_model = settings.optimization_model().get_graph_model<FP>(seed);

        // --- 2b. Apply control dampings on each node ---
        for (auto& node : graph_model.nodes()) {
            set_control_dampings<FP>(settings, node.property.get_simulation().get_model(), parameters);
        }

        // --- 2c. Set integrators on each node ---
        for (auto& node : graph_model.nodes()) {
            node.property.get_simulation().set_integrator_core(
                std::move(make_integrator<FP>(settings.integrator_type(), settings.dt())));
        }

        // --- 2d. Instantiate the mobility simulator ---
        auto graph_simulation = mio::make_mobility_sim<FP>(settings.t0(), settings.dt(), graph_model);

        // === 3. Prepare CSV output file for this run ===
        std::string filename = "population_time_series_run_" + std::to_string(run) + ".csv";
        std::ofstream csv_file(filename);

        csv_file << "Time";
        for (size_t i = 0; i < num_infection_states(); ++i) {
            csv_file << "," << infection_state_to_string(static_cast<mio::osecirvvs::InfectionState>(i));
        }
        csv_file << "\n";

        // --- Write initial state ---
        {
            std::vector<FP> totals(num_infection_states(), FP(0.0));

            // Sum over nodes and age groups
            for (auto& node : graph_simulation.get_graph().nodes()) {

                const auto& model = node.property.get_simulation().get_model();
                const auto& ts    = node.property.get_simulation().get_result();
                auto local_state  = ts.get_last_value();

                size_t groups = model.parameters.get_num_groups().get();

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

        // --- 4. Advance the mobility simulation over all intervals ---
        for (size_t k = 0; k < settings.num_intervals(); k++) {

            graph_simulation.advance(time_steps[k + 1]);

            std::vector<FP> totals(num_infection_states(), FP(0.0));

            // Sum over node models and age groups
            for (auto& node : graph_simulation.get_graph().nodes()) {

                const auto& model = node.property.get_simulation().get_model();
                const auto& ts    = node.property.get_simulation().get_result();
                auto local_state  = ts.get_last_value();

                size_t groups = model.parameters.get_num_groups().get();

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
    }

    //     // Here is how the simulation runs.
    //     // I have an old script to safe a single solution run which is outdated.
    //     // Please adapt it to the new structure if needed.
    //     // Also now store each different run in a different file.

    // // template <typename FP>
    // // FP objective_function(const SecirvvsOptimization& settings, const FP* ptr_parameters, size_t n)
    // // {
    // //     assert(n == settings.num_control_parameters() * settings.num_control_intervals());

    // //     FP objective_sum = 0.0;

    // //     // Loop over all simulation runs to account for uncertainty
    // //     for (size_t run = 0; run < settings.num_simulation_runs(); run++) {

    // //         // Set seed for this run
    // //         size_t seed = settings.base_seed() + run;

    // //         // Define graph model for this run
    // //         auto graph_model = settings.optimization_model().get_graph_model<FP>(seed);

    // //         // Extract control parameters
    // //         std::vector<FP> parameters(settings.num_control_parameters() * settings.num_control_intervals());
    // //         for (size_t control_interval = 0; control_interval < settings.num_control_intervals(); control_interval++) {
    // //             for (size_t control_index = 0; control_index < settings.num_control_parameters(); control_index++) {
    // //                 size_t idx      = control_interval * settings.num_control_parameters() + control_index;
    // //                 parameters[idx] = settings.activation_function()(ptr_parameters[idx]);
    // //             }
    // //         }

    // //         // Set control dampings for each node
    // //         for (auto& node : graph_model.nodes()) {
    // //             set_control_dampings<FP>(settings, node.property.get_simulation().get_model(), parameters);
    // //         }

    // //         // Set integrators for each node
    // //         for (auto& node : graph_model.nodes()) {
    // //             node.property.get_simulation().set_integrator_core(
    // //                 std::move(make_integrator<FP>(settings.integrator_type(), settings.dt())));
    // //         }

    // //         // Prepare time grid and simulation
    // //         std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());
    // //         auto graph_simulation      = mio::make_mobility_sim<FP>(settings.t0(), settings.dt(), graph_model);

    // //         // Run simulation and compute objective for this run
    // //         FP objective = 0.0;
    // //         for (size_t interval = 0; interval < settings.num_intervals(); interval++) {

    // //             graph_simulation.advance(time_steps[interval + 1]);

    // //             size_t control_interval = interval / settings.pc_resolution();

    // //             FP interval_cost = 0.0;
    // //             for (size_t control_index = 0; control_index < settings.num_control_parameters(); control_index++) {
    // //                 interval_cost += settings.control_parameters()[control_index].cost() *
    // //                                  parameters[control_interval * settings.num_control_parameters() + control_index];
    // //             }
    // //             objective += interval_cost / settings.num_intervals();
    // //         }

    // //         // Accumulate objective
    // //         objective_sum += objective;
    // //     }

    // //     // Return average over all simulation runs
    // //     return objective_sum / static_cast<FP>(settings.num_simulation_runs());
    // // }

    //     mio::osecirvvs::Model<FP> new_model = settings.optimization_model().get_graph_model<FP>(settings.base_seed());
    //     set_control_dampings<FP>(settings, new_model, parameters);

    //     std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());

    //     mio::osecirvvs::Simulation<FP> simulation(new_model, settings.t0(), settings.dt());
    //     simulation.set_integrator_core(std::move(make_integrator<FP>(settings.integrator_type(), settings.dt())));

    //     {
    //         // Open CSV file and write header
    //         std::ofstream csv_file("population_time_series.csv");
    //         csv_file << "Time";
    //         for (size_t i = 0; i < num_infection_states(); ++i) {
    //             csv_file << "," << infection_state_to_string(static_cast<mio::osecirvvs::InfectionState>(i));
    //         }
    //         csv_file << "\n";

    //         // --- Write initial state (time = t0) ---
    //         {
    //             const mio::TimeSeries<FP>& simulation_result       = simulation.get_result();
    //             Eigen::Ref<const Eigen::VectorX<FP>> initial_state = simulation_result.get_last_value();

    //             std::vector<FP> total_by_state(num_infection_states(), FP(0.0));
    //             for (mio::AgeGroup age_group(0); age_group < new_model.parameters.get_num_groups(); ++age_group) {
    //                 for (size_t state_index = 0; state_index < num_infection_states(); ++state_index) {
    //                     size_t idx = age_group.get() * num_infection_states() + state_index;
    //                     total_by_state[state_index] += initial_state[idx];
    //                 }
    //             }

    //             csv_file << std::fixed << std::setprecision(6) << settings.t0();
    //             for (const auto& val : total_by_state) {
    //                 csv_file << "," << val;
    //             }
    //             csv_file << "\n";
    //         }
    //         // ----------------------------------------

    //         for (size_t controlIndex = 0; controlIndex < settings.num_control_intervals(); controlIndex++) {
    //             for (size_t substep = 0; substep < settings.pc_resolution(); substep++) {
    //                 size_t timeStepIndex = controlIndex * settings.pc_resolution() + substep;

    //                 simulation.advance(time_steps[timeStepIndex + 1]);
    //                 simulation.get_dt() = settings.dt();

    //                 const mio::TimeSeries<FP>& simulation_result     = simulation.get_result();
    //                 Eigen::Ref<const Eigen::VectorX<FP>> final_state = simulation_result.get_last_value();

    //                 // Sum over all age groups per infection state
    //                 std::vector<FP> total_by_state(num_infection_states(), FP(0.0));
    //                 for (mio::AgeGroup age_group = 0; age_group < new_model.parameters.get_num_groups(); age_group++) {
    //                     for (size_t state_index = 0; state_index < num_infection_states(); state_index++) {
    //                         size_t idx = age_group.get() * num_infection_states() + state_index;
    //                         total_by_state[state_index] += final_state[idx];
    //                     }
    //                 }

    //                 // Write to CSV
    //                 csv_file << std::fixed << std::setprecision(6) << time_steps[timeStepIndex + 1];
    //                 for (const auto& val : total_by_state) {
    //                     csv_file << "," << val;
    //                 }
    //                 csv_file << "\n";
    //             }
    //         }

    //         // Close CSV file
    //         csv_file.close();
    //     }
}
