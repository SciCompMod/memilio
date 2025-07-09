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
void update_path_constraint(const SecirvvsOptimization& settings, const mio::osecirvvs::Model<FP>& model,
                            std::vector<FP>& path_constraints)
{
    // ------------------------------------------------------------------- //
    // Consider the model to be at an intermediate simulation step.        //
    // Query and fill the constraints of the model into 'path_constraints' //
    // ------------------------------------------------------------------- //
    assert(path_constraints.size() == settings.num_path_constraints());

    using std::max;

    for (size_t constraint_index = 0; constraint_index < settings.num_path_constraints(); constraint_index++) {
        Constraint constraint = settings.path_constraints()[constraint_index];
        auto states           = query_infection_states(constraint.name());

        FP value = 0.0;
        for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
            for (const auto& state : states) {
                value += model.populations[{age_group, state}];
            }
        }
        path_constraints[constraint_index] = max<FP>(path_constraints[constraint_index], value);
    }
}

template <typename FP>
void update_terminal_constraint(const SecirvvsOptimization& settings, const mio::osecirvvs::Model<FP>& model,
                                std::vector<FP>& terminal_constraints)
{
    // ----------------------------------------------------------------------- //
    // Consider the model to be at the final simulation step.                  //
    // Query and fill the constraints of the model into 'terminal_constraints' //
    // ----------------------------------------------------------------------- //
    assert(terminal_constraints.size() == settings.num_terminal_constraints());

    for (size_t constraint_index = 0; constraint_index < settings.num_terminal_constraints(); constraint_index++) {
        Constraint constraint = settings.terminal_constraints()[constraint_index];
        auto states           = query_infection_states(constraint.name());

        FP value = 0.0;
        for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
            for (const auto& state : states) {
                value += model.populations[{age_group, state}];
            }
        }
        terminal_constraints[constraint_index] = value;
    }
}