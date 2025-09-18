#pragma once

#include <iostream>
#include <vector>
#include <cstddef>
#include <cassert>

#include "tools/optimal_control/constraints/constraints.h"
#include "tools/optimal_control/constraints/infection_state_utils.h"

#include <vector>

template <typename FP, class OptimizationSettings>
void update_path_constraint(const OptimizationSettings& settings, const typename OptimizationSettings::template ModelTemplate<FP>& model,
                            const auto& final_state, std::vector<FP>& path_constraints)
{
    // ------------------------------------------------------------------- //
    // Consider the model to be at an intermediate simulation step.        //
    // Query and fill the constraints of the model into 'path_constraints' //
    // ------------------------------------------------------------------- //
    assert(path_constraints.size() == settings.num_path_constraints());

    using std::max;

    for (size_t constraint_index = 0; constraint_index < settings.num_path_constraints(); constraint_index++) {
        Constraint constraint = settings.path_constraints()[constraint_index];
        auto states           = query_infection_states<OptimizationSettings::InfectionState>(constraint.name(), settings.states_strings);

        FP value = 0.0;
        for (mio::AgeGroup agegroup = 0; agegroup < model.parameters.get_num_groups(); agegroup++)
        {
            auto age_group_offset = (int)agegroup * (int)num_infection_states<OptimizationSettings::InfectionState>();

            for (const auto& state : states) {
                value += final_state[(int)state + age_group_offset];
            }
        }

        path_constraints[constraint_index] = max<FP>(path_constraints[constraint_index], value);
    }
}

template <typename FP, class OptimizationSettings>
void update_terminal_constraint(const OptimizationSettings& settings, const typename OptimizationSettings::template ModelTemplate<FP>& model,
                                const auto& final_state, std::vector<FP>& terminal_constraints)
{
    // ----------------------------------------------------------------------- //
    // Consider the model to be at the final simulation step.                  //
    // Query and fill the constraints of the model into 'terminal_constraints' //
    // ----------------------------------------------------------------------- //
    assert(terminal_constraints.size() == settings.num_terminal_constraints());

    for (size_t constraint_index = 0; constraint_index < settings.num_terminal_constraints(); constraint_index++) {
        Constraint constraint = settings.terminal_constraints()[constraint_index];
        auto states           = query_infection_states<OptimizationSettings::InfectionState>(constraint.name(), settings.states_strings);

        FP value = 0.0;
        for (mio::AgeGroup agegroup = 0; agegroup < model.parameters.get_num_groups(); agegroup++)
        {
            auto age_group_offset = (int)agegroup * (int)num_infection_states<OptimizationSettings::InfectionState>();

            for (const auto& state : states) {
                value += final_state[(int)state + age_group_offset];
            }
        }
        terminal_constraints[constraint_index] = value;
    }
}