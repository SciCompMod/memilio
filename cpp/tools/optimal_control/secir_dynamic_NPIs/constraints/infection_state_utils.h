#pragma once

#include "models/ode_secir/model.h"

#include <string>
#include <vector>
#include <utility>
#include <sstream>
#include <algorithm>
#include <stdexcept>

inline size_t num_infection_states()
{
    return static_cast<size_t>(mio::osecir::InfectionState::Count);
}

inline std::string infection_state_to_string(mio::osecir::InfectionState state)
{
    using InfectionState = mio::osecir::InfectionState;

    switch (state) {
    case InfectionState::Susceptible:
        return "Susceptible";
    case InfectionState::Exposed:
        return "Exposed";
    case InfectionState::InfectedNoSymptoms:
        return "InfectedNoSymptoms";
    case InfectionState::InfectedNoSymptomsConfirmed:
        return "InfectedNoSymptomsConfirmed";
    case InfectionState::InfectedSymptoms:
        return "InfectedSymptoms";
    case InfectionState::InfectedSymptomsConfirmed:
        return "InfectedSymptomsConfirmed";
    case InfectionState::InfectedSevere:
        return "InfectedSevere";
    case InfectionState::InfectedCritical:
        return "InfectedCritical";
    case InfectionState::Recovered:
        return "Recovered";
    case InfectionState::Dead:
        return "Dead";
    case InfectionState::Count:
        return "Count";
    default:
        return "InfectionState Unknown";
    }
}

inline std::vector<mio::osecir::InfectionState> query_infection_states(const std::string& state_query)
{
    using InfectionState = mio::osecir::InfectionState;

    static const std::vector<std::pair<std::string, InfectionState>> all_states = {
        {"Susceptible", InfectionState::Susceptible},
        {"Exposed", InfectionState::Exposed},
        {"InfectedNoSymptoms", InfectionState::InfectedNoSymptoms},
        {"InfectedNoSymptomsConfirmed", InfectionState::InfectedNoSymptomsConfirmed},
        {"InfectedSymptoms", InfectionState::InfectedSymptoms},
        {"InfectedSymptomsConfirmed", InfectionState::InfectedSymptomsConfirmed},
        {"InfectedSevere", InfectionState::InfectedSevere},
        {"InfectedCritical", InfectionState::InfectedCritical},
        {"Recovered", InfectionState::Recovered},
        {"Dead", InfectionState::Dead}};

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
