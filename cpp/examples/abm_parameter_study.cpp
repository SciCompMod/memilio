/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding, Sascha Korf
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
#include "abm/result_simulation.h"
#include "abm/household.h"
#include "abm/lockdown_rules.h"
#include "abm/model.h"
#include "abm/time.h"
#include "memilio/timer/auto_timer.h"

#include "memilio/compartments/parameter_studies.h"
#include "memilio/data/analyze_result.h"
#include "memilio/io/io.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/base_dir.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/stl_util.h"
#include "abm/city_builder.h"

#include <string>

// constexpr size_t num_age_groups = 5;

/// An ABM setup taken from abm_minimal.cpp.
mio::abm::Model make_model(size_t num_persons, mio::RandomNumberGenerator& rng)
{

    auto model = CityBuilder::build_world(CityConfig{static_cast<int>(num_persons)}, rng);

    //infections and masks
    for (auto& person : model.get_persons()) {
        auto prng = mio::abm::PersonalRandomNumberGenerator(person);
        //some % of people are infected, large enough to have some infection activity without everyone dying
        auto pct_infected = 0.0005;
        if (mio::UniformDistribution<ScalarType>::get_instance()(prng, 0.0, 1.0) < pct_infected) {
            auto infection =
                mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype, person.get_age(), model.parameters,
                                    mio::abm::TimePoint(0), mio::abm::InfectionState::Exposed);
            person.add_new_infection(std::move(infection));
        }
    }

    return model;
}

int main()
{
    mio::mpi::init();

    mio::set_log_level(mio::LogLevel::off);

    // Set start and end time for the simulation.
    auto t0   = mio::abm::TimePoint(0);
    auto tmax = t0 + mio::abm::days(14);
    // Set the number of simulations to run in the study
    const size_t num_runs = 128 * 128;

    // Create a parameter study.
    // Note that the study for the ABM currently does not make use of the arguments "parameters" or "dt", as we create
    // a new model for each simulation. Hence we set both arguments to 0.
    // This is mostly due to https://github.com/SciCompMod/memilio/issues/1400
    mio::ParameterStudy study(0, t0, tmax, mio::abm::TimeSpan(0), num_runs);

    mio::timing::AutoTimer<"abm_timer"> timer;

    // Optional: set seeds to get reproducable results
    study.get_rng().seed({12341234, 53456, 63451, 5232576, 84586, 52345});

    const std::string result_dir = mio::path_join(mio::base_dir(), "example_results");
    std::cout << "Writing results to " << result_dir << std::endl;
    if (!mio::create_directory(result_dir)) {
        mio::log_error("Could not create result directory \"{}\".", result_dir);
        return 1;
    }

    auto ensemble_results = study.run(
        [](auto, auto t0_, auto, size_t) {
            return mio::abm::ResultSimulation(make_model(2'000'000, mio::thread_local_rng()), t0_);
        },
        [result_dir](auto&& sim, auto&& run_idx) {
            auto interpolated_result = mio::interpolate_simulation_result(sim.get_result());
            std::string outpath = mio::path_join(result_dir, "abm_minimal_run_" + std::to_string(run_idx) + ".txt");
            std::ofstream outfile_run(outpath);
            sim.get_result().print_table(outfile_run, {"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4);
            std::cout << "Results written to " << outpath << std::endl;
            auto params = std::vector<mio::abm::Model>{};
            return std::vector{interpolated_result};
        });

    if (ensemble_results.size() > 0) {
        auto ensemble_results_p05 = ensemble_percentile(ensemble_results, 0.05);
        auto ensemble_results_p25 = ensemble_percentile(ensemble_results, 0.25);
        auto ensemble_results_p50 = ensemble_percentile(ensemble_results, 0.50);
        auto ensemble_results_p75 = ensemble_percentile(ensemble_results, 0.75);
        auto ensemble_results_p95 = ensemble_percentile(ensemble_results, 0.95);

        mio::unused(save_result(ensemble_results_p05, {0}, 8,
                                mio::path_join(result_dir, "Results_" + std::string("p05") + ".h5")));
        mio::unused(save_result(ensemble_results_p25, {0}, 8,
                                mio::path_join(result_dir, "Results_" + std::string("p25") + ".h5")));
        mio::unused(save_result(ensemble_results_p50, {0}, 8,
                                mio::path_join(result_dir, "Results_" + std::string("p50") + ".h5")));
        mio::unused(save_result(ensemble_results_p75, {0}, 8,
                                mio::path_join(result_dir, "Results_" + std::string("p75") + ".h5")));
        mio::unused(save_result(ensemble_results_p95, {0}, 8,
                                mio::path_join(result_dir, "Results_" + std::string("p95") + ".h5")));
    }

    mio::mpi::finalize();

    return 0;
}
