/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Khoa Nguyen, Nils Waßmuth
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
// Minimal standalone driver for forward_pass(): builds a small city, runs one 30-day forward
// pass with budget-constrained PCR surveillance + diagnostic care-seeking, and writes the raw
// resolved tests (the positives table and the aggregated negatives) to files.
#include "abm/pcr_surveillance.h"
#include "forward_pass.h"
#include "memilio/io/directories.h"
#include "memilio/utils/logging.h"

#include <fstream>
#include <iostream>

int main()
{
    mio::set_log_level(mio::LogLevel::warn);

    // Build a small city
    ABMPopulation population(/*total_population=*/500);

    // Random surveillance: sample 5% of the eligible population per surveillance event,
    // every day. Results are reported and isolation started on a positive with  a 1-day delay.
    // Diagnostic care-seeking runs (symptomatic individuals) alongside surveillance with the hazard fixed inside forward_pass().
    mio::abm::TestingBudget design;
    design.event_budget_fraction = 0.05;
    design.test_period_days      = 1;

    const ScalarType beta  = 1.0;
    const ScalarType t_exposed = 5.0;
    auto out               = forward_pass(population, {{"beta", beta}, {"t_exposed", t_exposed}}, design);

    const Eigen::MatrixXd& positives = out.at("positives"); // [day, person_id, age, loc, ct, source]
    const Eigen::MatrixXd& negatives = out.at("negatives"); // [day, source, age, loc, count]
    const auto n_ever_infected       = static_cast<long>(out.at("n_ever_infected")(0, 0));

    auto outdir      = mio::create_directories_or_exit(mio::example_results_dir("abm_minimal"));
    auto pos_path    = outdir / "positives.txt";
    auto neg_path    = outdir / "negatives.txt";
    std::ofstream(pos_path.string()) << positives << "\n";
    std::ofstream(neg_path.string()) << negatives << "\n";

    // Column 5 of the positives table is the TestSource (0 = survey, 1 = diagnostic).
    long survey_positives = 0;
    for (Eigen::Index i = 0; i < positives.rows(); ++i) {
        survey_positives += (static_cast<int>(positives(i, 5)) ==
                             static_cast<int>(mio::abm::TestSource::Survey))
                                ? 1
                                : 0;
    }

    const long n_positives = static_cast<long>(positives.rows());
    std::cout << "positives table (" << positives.rows() << "x" << positives.cols() << ") -> " << pos_path << "\n";
    std::cout << "negatives table (" << negatives.rows() << "x" << negatives.cols() << ") -> " << neg_path << "\n";
    std::cout << "detected positives: " << n_positives << " (survey: " << survey_positives
              << ", diagnostic: " << n_positives - survey_positives << ")\n";
    std::cout << "Persons ever infected: " << n_ever_infected << std::endl;

    return 0;
}
