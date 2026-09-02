#ifndef LOCATION_SPLIT_CUSTOM_LOGGERS_H
#define LOCATION_SPLIT_CUSTOM_LOGGERS_H

#include "abm/infection_state.h"
#include "abm/location_id.h"
#include "abm/location_type.h"
#include "abm/person_id.h"
#include "abm/simulation.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/io/history.h"

#include <Eigen/Dense>
#include <cstdint>
#include <tuple>
#include <vector>

/**
 * @file custom_loggers.h
 * @brief Loggers used by the location split simulation.
 *
 * The panvXabmSim branch logged the site of every new Infection with Location::get_infected_persons(),
 * which only exists on that branch. Here LogPersonStates records the Location and the InfectionState of
 * every Person at every time step instead, and the infection events are reconstructed from that log in
 * multi_run_simulator.cpp. This keeps the simulation free of changes to the ABM itself.
 */

/// @brief One record of LogPersonStates: who was where and in which InfectionState.
struct PersonState {
    mio::abm::PersonId person_id;
    mio::abm::LocationId location_id;
    mio::abm::LocationType location_type;
    mio::abm::InfectionState infection_state;
    mio::AgeGroup age;
};

/// @brief Number of Person%s per InfectionState and AgeGroup, flattened as age * Count + state.
struct LogInfectionStatePerAgeGroup : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorX<ScalarType>>;
    static Type log(const mio::abm::Simulation<>& sim);
};

/// @brief Number of Person%s that are not Susceptible, i.e. the cumulative number of infections.
struct LogTotalInfections : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorX<ScalarType>>;
    static Type log(const mio::abm::Simulation<>& sim);
};

/// @brief Location and InfectionState of every Person at the current time step.
struct LogPersonStates : mio::LogAlways {
    using Type = std::vector<PersonState>;
    static Type log(const mio::abm::Simulation<>& sim);
};

/// @brief Person%s that are infectious, i.e. infected and currently shedding virus.
struct LogInfectiousPersons : mio::LogAlways {
    using Type = std::vector<mio::abm::PersonId>;
    static Type log(const mio::abm::Simulation<>& sim);
};

/// @brief The Home assigned to each Person. Logged once.
struct LogHouseholdId : mio::LogOnce {
    using Type = std::vector<std::tuple<mio::abm::PersonId, mio::abm::LocationId>>;
    static Type log(const mio::abm::Simulation<>& sim);
};

/// @brief The Work assigned to each Person, or an invalid LocationId. Logged once.
struct LogWorkId : mio::LogOnce {
    using Type = std::vector<std::tuple<mio::abm::PersonId, mio::abm::LocationId>>;
    static Type log(const mio::abm::Simulation<>& sim);
};

/// @brief The LocationType of every Location. Logged once.
struct LogLocationIdAndType : mio::LogOnce {
    using Type = std::vector<std::tuple<mio::abm::LocationId, mio::abm::LocationType>>;
    static Type log(const mio::abm::Simulation<>& sim);
};

#endif // LOCATION_SPLIT_CUSTOM_LOGGERS_H
