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
 * Only aggregated counts are logged here. The site of every new Infection is detected hour by hour in
 * MultiRunSimulator::run_single_simulation by comparing every Person with its state of the previous
 * hour, which needs memory for one hour instead of the whole simulation.
 */

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

#endif // LOCATION_SPLIT_CUSTOM_LOGGERS_H
