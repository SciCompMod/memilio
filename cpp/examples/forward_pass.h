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
// Returns (histogram, cohort):
//   histogram:    (n_days, 42)              — [day, count_ct0..count_ct40]
//   fixed_cohort: (n_days, cohort_budget+1) — [day, ct_student0..ct_student_{k-1}]
// Encoding: 0 = max viral load, 40 = not detected, 255 = unused cohort slot.
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> forward_pass(const ABMPopulation& population,
                                                          ScalarType beta, ScalarType kappa,
                                                          int cohort_budget = 50);

// Convenience wrapper: builds a fresh ABMPopulation on every call (original behaviour).
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> forward_pass(ScalarType beta, ScalarType kappa,
                                                          int cohort_budget = 50);
