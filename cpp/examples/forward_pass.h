#pragma once
#include <Eigen/Core>
#include <memory>
#include <utility>
#include "memilio/config.h"

// Pre-built ABM population: German-demographics-based households, schools, workplaces,
// shops, events, and hospital/ICU (via CityBuilder), with all agents initially susceptible.
// Construct once with ABMPopulation(total_population), then reuse across forward_pass() calls.
struct ABMPopulation {
    struct Impl;
    std::shared_ptr<Impl> impl;
    explicit ABMPopulation(int total_population = 500);
};

// Runs one forward pass: copies the population, sets beta/kappa, randomises initial
// infections, and simulates. Reusing the same ABMPopulation avoids rebuilding the
// household/location structure on every call.
// Returns (histogram_school, histogram_work):
//   histogram_school: (n_days, 42) — [day, count_ct0..count_ct40] for age 5–14 at School
//   histogram_work:   (n_days, 83) — [day, count_ct0..count_ct40 (age 15–34), count_ct0..count_ct40 (age 35–59)]
// Encoding: 0 = max viral load, 40 = not detected.
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> forward_pass(const ABMPopulation& population,
                                                          ScalarType beta, ScalarType kappa);

// Convenience wrapper: builds a fresh ABMPopulation on every call (original behaviour).
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> forward_pass(ScalarType beta, ScalarType kappa);
