#pragma once
#include "memilio/io/result_io.h"
#include "abm/simulation.h"
#include "abm/time.h"
#include "abm/world.h"
#include <Eigen/Dense>

// DEBUGGING
struct LogInfectionStatePerAgeGroup : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    static Type log(const mio::abm::Simulation& sim);
};

struct LogInfectionPerLocationTypePerAgeGroup : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    static Type log(const mio::abm::Simulation& sim);
};

struct LogLocationTypeAndId : mio::LogAlways {
    using Type = std::vector<std::tuple<uint32_t, mio::abm::LocationType>>;
    static Type log(const mio::abm::Simulation& sim);
};

struct LogLocationIdAndPersonId : mio::LogAlways {
    using Type = std::vector<std::tuple<uint32_t, uint32_t>>;
    static Type log(const mio::abm::Simulation& sim);
};

struct LogWhoInfected : mio::LogAlways {
    using Type = std::vector<uint32_t>;
    static Type log(const mio::abm::Simulation& sim);
};

struct LogHouseholdId : mio::LogOnce {
    using Type = std::vector<std::tuple<uint32_t, uint32_t>>;
    static Type log(const mio::abm::Simulation& sim)
    {
        Type household_information{};
        household_information.reserve(sim.get_world().get_persons().size());
        for (auto&& person : sim.get_world().get_persons()) {
            household_information.push_back(std::make_tuple(
                person.get_person_id(), person.get_assigned_location_index(mio::abm::LocationType::Home)));
        }
        return household_information;
    }
};

struct LogInfectionDetailed : mio::LogAlways {
    using Type = std::vector<std::tuple<uint32_t, uint32_t, mio::abm::LocationType>>;
    static Type log(const mio::abm::Simulation& sim);
};

struct LogAmountOfInfections : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    static Type log(const mio::abm::Simulation& sim);
};