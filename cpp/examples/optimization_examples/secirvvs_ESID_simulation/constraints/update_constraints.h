#pragma once

#include <iostream>
#include <vector>
#include <cstddef>
#include <cassert>

#include "constraints.h"
#include "infection_state_utils.h"

#include "../optimization_settings/secirvvs_optimization.h"

#include <vector>

template <typename FP>
void update_path_constraint(const SecirvvsOptimization& optimization_settings,
                            const mio::osecirvvs::Simulation<FP>& simulation, std::vector<FP>& path_constraints)
{
    // ------------------------------------------------------------------------ //
    // This function is called during an intermediate simulation step.         //
    // It queries the current state of the model and updates the values in     //
    // 'path_constraints' based on the maximum encountered so far.             //
    // ------------------------------------------------------------------------ //
    assert(path_constraints.size() == optimization_settings.num_path_constraints());

    const mio::TimeSeries<FP>& simulation_result      = simulation.get_result();
    Eigen::Ref<const Eigen::VectorX<FP>> latest_state = simulation_result.get_last_value();

    const size_t num_age_groups = simulation.get_model().parameters.get_num_groups();
    const size_t num_states     = num_infection_states();

    using std::max;

    for (size_t constraint_idx = 0; constraint_idx < optimization_settings.num_path_constraints(); constraint_idx++) {
        const Constraint& constraint = optimization_settings.path_constraints()[constraint_idx];
        const auto relevant_states   = query_infection_states(constraint.name());

        FP total_value = 0.0;
        for (size_t age_group = 0; age_group < num_age_groups; age_group++) {
            for (const auto& state : relevant_states) {
                size_t index = age_group * num_states + static_cast<size_t>(state);
                total_value += latest_state[index];
            }
        }

        path_constraints[constraint_idx] = max<FP>(path_constraints[constraint_idx], total_value);
    }
}

template <typename FP>
void update_terminal_constraint(const SecirvvsOptimization& optimization_settings,
                                const mio::osecirvvs::Simulation<FP>& simulation, std::vector<FP>& terminal_constraints)
{
    // ------------------------------------------------------------------------ //
    // This function is called at the final simulation step.                    //
    // It queries the final state vector and stores the constraint values in    //
    // 'terminal_constraints'.                                                  //
    // ------------------------------------------------------------------------ //
    assert(terminal_constraints.size() == optimization_settings.num_terminal_constraints());

    const mio::TimeSeries<FP>& simulation_result     = simulation.get_result();
    Eigen::Ref<const Eigen::VectorX<FP>> final_state = simulation_result.get_last_value();

    const size_t num_age_groups = simulation.get_model().parameters.get_num_groups();
    const size_t num_states     = num_infection_states();

    for (size_t constraint_idx = 0; constraint_idx < optimization_settings.num_terminal_constraints();
         constraint_idx++) {
        const Constraint& constraint = optimization_settings.terminal_constraints()[constraint_idx];
        const auto relevant_states   = query_infection_states(constraint.name());

        FP total_value = 0.0;
        for (size_t age_group = 0; age_group < num_age_groups; ++age_group) {
            for (const auto& state : relevant_states) {
                size_t index = age_group * num_states + static_cast<size_t>(state);
                total_value += final_state[index];
            }
        }

        terminal_constraints[constraint_idx] = total_value;
    }
}