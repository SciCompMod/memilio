#pragma once

#include <iostream>
#include <vector>
#include <cstddef>
#include <cassert>
#include <string>
#include <utility>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include "tools/optimal_control/constraints/constraints.h"



template <class InfectionState>
inline size_t num_infection_states()
{
    return static_cast<size_t>(InfectionState::Count);
}

template <class InfectionState>
inline std::vector<InfectionState> query_infection_states(const std::string& state_query, const std::vector<std::pair<std::string, InfectionState>>& all_states)
{
    // Split query into separate subqueries
    std::vector<std::string> and_group_tokens;
    std::vector<std::string> or_group_tokens;

    std::istringstream iss(state_query);
    std::string word;
    while (iss >> word) {
        if (word.find('+') != std::string::npos) {
            std::istringstream substream(word);
            std::string sub;
            while (std::getline(substream, sub, '+')) {
                and_group_tokens.push_back(sub);
            }
        }
        else {
            or_group_tokens.push_back(word);
        }
    }

    std::vector<InfectionState> result;

    if (!and_group_tokens.empty()) {
        for (const auto& [name, state] : all_states) {
            bool match_all =
                std::all_of(and_group_tokens.begin(), and_group_tokens.end(), [&](const std::string& token) {
                    return name.find(token) != std::string::npos;
                });
            if (match_all) {
                result.push_back(state);
            }
        }
    }

    if (!or_group_tokens.empty()) {
        for (const auto& [name, state] : all_states) {
            bool match_any = std::any_of(or_group_tokens.begin(), or_group_tokens.end(), [&](const std::string& token) {
                return name.find(token) != std::string::npos;
            });
            if (match_any && std::find(result.begin(), result.end(), state) == result.end()) {
                result.push_back(state);
            }
        }
    }

    if (result.empty()) {
        throw std::invalid_argument("No matching InfectionState for query: " + state_query);
    }

    return result;
}

template <typename FP, class OptimizationSettings>
void update_path_constraint(const OptimizationSettings& settings, const auto& final_state, std::vector<FP>& path_constraints)
{
    // ------------------------------------------------------------------- //
    // Consider the model to be at an intermediate simulation step.        //
    // Query and fill the constraints of the model into 'path_constraints' //
    // ------------------------------------------------------------------- //
    assert(path_constraints.size() == settings.num_path_constraints());

    using std::max;

    for (size_t constraint_index = 0; constraint_index < settings.num_path_constraints(); constraint_index++) {
        Constraint constraint = settings.path_constraints()[constraint_index];
        auto states           = query_infection_states<typename OptimizationSettings::InfectionState>(constraint.name(), settings.states_strings());

        FP value = 0.0;
        for (int agegroup = 0; agegroup < settings.optimization_model().num_age_groups(); agegroup++)
        {
            auto age_group_offset = agegroup * (int)num_infection_states<typename OptimizationSettings::InfectionState>();

            for (const auto& state : states) {
                value += final_state[(int)state + age_group_offset];
            }
        }

        path_constraints[constraint_index] = max<FP>(path_constraints[constraint_index], value);
    }
}

template <typename FP, class OptimizationSettings>
void update_terminal_constraint(const OptimizationSettings& settings, const auto& final_state, std::vector<FP>& terminal_constraints)
{
    // ----------------------------------------------------------------------- //
    // Consider the model to be at the final simulation step.                  //
    // Query and fill the constraints of the model into 'terminal_constraints' //
    // ----------------------------------------------------------------------- //
    assert(terminal_constraints.size() == settings.num_terminal_constraints());

    for (size_t constraint_index = 0; constraint_index < settings.num_terminal_constraints(); constraint_index++) {
        Constraint constraint = settings.terminal_constraints()[constraint_index];
        auto states           = query_infection_states<typename OptimizationSettings::InfectionState>(constraint.name(), settings.states_strings());

        FP value = 0.0;
        for (int agegroup = 0; agegroup < settings.optimization_model().num_age_groups(); agegroup++)
        {
            auto age_group_offset = agegroup * (int)num_infection_states<typename OptimizationSettings::InfectionState>();

            for (const auto& state : states) {
                value += final_state[(int)state + age_group_offset];
            }
        }
        terminal_constraints[constraint_index] = value;
    }
}
