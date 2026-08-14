#pragma once
#include <Eigen/Core>
#include <map>
#include <memory>
#include <string>
#include "abm/pcr_surveillance.h"
#include "memilio/config.h"

// Pre-built ABM population: German-demographics-based households, schools, workplaces,
// shops, events, and hospital/ICU (via CityBuilder)
struct ABMPopulation {
    struct Impl;
    std::shared_ptr<Impl> impl;
    explicit ABMPopulation(int total_population = 500, ScalarType quarantine_compliance = 1.0);
};

// Runs one forward pass
std::map<std::string, Eigen::MatrixXd> forward_pass(const ABMPopulation& population,
                                                    const std::map<std::string, ScalarType>& params,
                                                    const mio::abm::TestingBudget& design);
