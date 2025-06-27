#pragma once
#include "memilio/io/result_io.h"
#include "abm/simulation.h"
#include "abm/time.h"
#include "abm/world.h"
#include <Eigen/Dense>

struct LogInfectionStatePerAgeGroup : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    static Type log(const mio::abm::Simulation& sim);
};

struct LogInfectionPerLocationTypePerAgeGroup : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    static Type log(const mio::abm::Simulation& sim);
};