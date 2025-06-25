#pragma once
#include "abm/abm.h"
#include <Eigen/Dense>

struct LogInfectionStatePerAgeGroup : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    static Type log(const mio::abm::Simulation& sim);
};

struct LogInfectionPerLocationTypePerAgeGroup : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    static Type log(const mio::abm::Simulation& sim);
};