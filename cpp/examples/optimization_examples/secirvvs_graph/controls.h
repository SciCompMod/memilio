// clang-format off
#pragma once


#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>

// ------------------------- //
// Handle Control Parameters //
// ------------------------- //

enum class ControlParameter {
    SchoolClosure,
    HomeOffice,
    PhysicalDistancingSchool,
    PhysicalDistancingWork,
    PhysicalDistancingOther,
    Count
};

inline ControlParameter string_to_control(const std::string& control_name) {
    static const std::vector<std::pair<std::string, ControlParameter>> control_pairs = {
        {"SchoolClosure", ControlParameter::SchoolClosure},
        {"HomeOffice", ControlParameter::HomeOffice},
        {"PhysicalDistancingSchool", ControlParameter::PhysicalDistancingSchool},
        {"PhysicalDistancingWork", ControlParameter::PhysicalDistancingWork},
        {"PhysicalDistancingOther", ControlParameter::PhysicalDistancingOther},
    };
    auto it = std::find_if(
        control_pairs.begin(), control_pairs.end(),
        [&](const auto& pair) {
            return pair.first == control_name;
        }
    );
    if (it != control_pairs.end()) {
        return it->second;
    }
    throw std::runtime_error("Invalid control name: " + control_name);
}

template <typename FP>
struct ControlPolicy {
    std::string name;
    std::pair<FP, FP> range;
    FP effectiveness;
    FP cost;

    constexpr const std::string& get_name() const noexcept {
        return name;
    }

    constexpr std::pair<FP, FP> get_range() const noexcept {
        return range;
    }

    constexpr FP get_min() const noexcept {
        return range.first;
    }

    constexpr FP get_max() const noexcept {
        return range.second;
    }

    constexpr FP get_effectiveness() const noexcept {
        return effectiveness;
    }

    constexpr FP get_cost() const noexcept {
        return cost;
    }
};

