#pragma once

#include "models/ode_secirvvs/model.h"

#include <string>
#include <vector>
#include <utility>
#include <sstream>
#include <algorithm>
#include <stdexcept>

inline size_t num_infection_states()
{
    return static_cast<size_t>(mio::osecirvvs::InfectionState::Count);
}

inline std::string infection_state_to_string(mio::osecirvvs::InfectionState state)
{
    using InfectionState = mio::osecirvvs::InfectionState;

    switch (state) {
    case InfectionState::SusceptibleNaive:
        return "SusceptibleNaive";
    case InfectionState::SusceptiblePartialImmunity:
        return "SusceptiblePartialImmunity";
    case InfectionState::SusceptibleImprovedImmunity:
        return "SusceptibleImprovedImmunity";
    case InfectionState::ExposedNaive:
        return "ExposedNaive";
    case InfectionState::ExposedPartialImmunity:
        return "ExposedPartialImmunity";
    case InfectionState::ExposedImprovedImmunity:
        return "ExposedImprovedImmunity";
    case InfectionState::InfectedNoSymptomsNaive:
        return "InfectedNoSymptomsNaive";
    case InfectionState::InfectedNoSymptomsPartialImmunity:
        return "InfectedNoSymptomsPartialImmunity";
    case InfectionState::InfectedNoSymptomsImprovedImmunity:
        return "InfectedNoSymptomsImprovedImmunity";
    case InfectionState::InfectedNoSymptomsNaiveConfirmed:
        return "InfectedNoSymptomsNaiveConfirmed";
    case InfectionState::InfectedNoSymptomsPartialImmunityConfirmed:
        return "InfectedNoSymptomsPartialImmunityConfirmed";
    case InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed:
        return "InfectedNoSymptomsImprovedImmunityConfirmed";
    case InfectionState::InfectedSymptomsNaive:
        return "InfectedSymptomsNaive";
    case InfectionState::InfectedSymptomsPartialImmunity:
        return "InfectedSymptomsPartialImmunity";
    case InfectionState::InfectedSymptomsImprovedImmunity:
        return "InfectedSymptomsImprovedImmunity";
    case InfectionState::InfectedSymptomsNaiveConfirmed:
        return "InfectedSymptomsNaiveConfirmed";
    case InfectionState::InfectedSymptomsPartialImmunityConfirmed:
        return "InfectedSymptomsPartialImmunityConfirmed";
    case InfectionState::InfectedSymptomsImprovedImmunityConfirmed:
        return "InfectedSymptomsImprovedImmunityConfirmed";
    case InfectionState::InfectedSevereNaive:
        return "InfectedSevereNaive";
    case InfectionState::InfectedSeverePartialImmunity:
        return "InfectedSeverePartialImmunity";
    case InfectionState::InfectedSevereImprovedImmunity:
        return "InfectedSevereImprovedImmunity";
    case InfectionState::InfectedCriticalNaive:
        return "InfectedCriticalNaive";
    case InfectionState::InfectedCriticalPartialImmunity:
        return "InfectedCriticalPartialImmunity";
    case InfectionState::InfectedCriticalImprovedImmunity:
        return "InfectedCriticalImprovedImmunity";
    case InfectionState::DeadNaive:
        return "DeadNaive";
    case InfectionState::DeadPartialImmunity:
        return "DeadPartialImmunity";
    case InfectionState::DeadImprovedImmunity:
        return "DeadImprovedImmunity";
    case InfectionState::Count:
        return "Count";
    default:
        return "InfectionState Unknown";
    }
}

inline std::vector<mio::osecirvvs::InfectionState> query_infection_states(const std::string& state_query)
{
    using InfectionState = mio::osecirvvs::InfectionState;

    static const std::vector<std::pair<std::string, InfectionState>> all_states = {
        {"SusceptibleNaive", InfectionState::SusceptibleNaive},
        {"SusceptiblePartialImmunity", InfectionState::SusceptiblePartialImmunity},
        {"SusceptibleImprovedImmunity", InfectionState::SusceptibleImprovedImmunity},
        {"ExposedNaive", InfectionState::ExposedNaive},
        {"ExposedPartialImmunity", InfectionState::ExposedPartialImmunity},
        {"ExposedImprovedImmunity", InfectionState::ExposedImprovedImmunity},
        {"InfectedNoSymptomsNaive", InfectionState::InfectedNoSymptomsNaive},
        {"InfectedNoSymptomsPartialImmunity", InfectionState::InfectedNoSymptomsPartialImmunity},
        {"InfectedNoSymptomsImprovedImmunity", InfectionState::InfectedNoSymptomsImprovedImmunity},
        {"InfectedNoSymptomsNaiveConfirmed", InfectionState::InfectedNoSymptomsNaiveConfirmed},
        {"InfectedNoSymptomsPartialImmunityConfirmed", InfectionState::InfectedNoSymptomsPartialImmunityConfirmed},
        {"InfectedNoSymptomsImprovedImmunityConfirmed", InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed},
        {"InfectedSymptomsNaive", InfectionState::InfectedSymptomsNaive},
        {"InfectedSymptomsPartialImmunity", InfectionState::InfectedSymptomsPartialImmunity},
        {"InfectedSymptomsImprovedImmunity", InfectionState::InfectedSymptomsImprovedImmunity},
        {"InfectedSymptomsNaiveConfirmed", InfectionState::InfectedSymptomsNaiveConfirmed},
        {"InfectedSymptomsPartialImmunityConfirmed", InfectionState::InfectedSymptomsPartialImmunityConfirmed},
        {"InfectedSymptomsImprovedImmunityConfirmed", InfectionState::InfectedSymptomsImprovedImmunityConfirmed},
        {"InfectedSevereNaive", InfectionState::InfectedSevereNaive},
        {"InfectedSeverePartialImmunity", InfectionState::InfectedSeverePartialImmunity},
        {"InfectedSevereImprovedImmunity", InfectionState::InfectedSevereImprovedImmunity},
        {"InfectedCriticalNaive", InfectionState::InfectedCriticalNaive},
        {"InfectedCriticalPartialImmunity", InfectionState::InfectedCriticalPartialImmunity},
        {"InfectedCriticalImprovedImmunity", InfectionState::InfectedCriticalImprovedImmunity},
        {"DeadNaive", InfectionState::DeadNaive},
        {"DeadPartialImmunity", InfectionState::DeadPartialImmunity},
        {"DeadImprovedImmunity", InfectionState::DeadImprovedImmunity}};

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
        for (const auto& [n, state] : all_states) {
            const auto& name = n; // create a named reference or copy
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
        for (const auto& [n, state] : all_states) {
            const auto& name = n;
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
