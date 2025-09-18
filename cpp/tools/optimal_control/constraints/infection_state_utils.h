#pragma once

#include <string>
#include <vector>
#include <utility>
#include <sstream>
#include <algorithm>
#include <stdexcept>

template <class InfectionState>
inline size_t num_infection_states()
{
    return static_cast<size_t>(InfectionState::Count);
}

template <class InfectionState>
inline std::vector<InfectionState> query_infection_states(const std::string& state_query, std::vector<std::pair<std::string, InfectionState>>& all_states)
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