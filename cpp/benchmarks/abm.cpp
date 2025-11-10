/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele
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
#include "abm/simulation.h"
#include "abm/common_abm_loggers.h"
#include "benchmark/benchmark.h"

mio::abm::Simulation<> make_simulation(size_t num_persons, std::initializer_list<uint32_t> seeds)
{
    auto rng = mio::RandomNumberGenerator();
    rng.seed(seeds);
    auto model      = mio::abm::Model(5);
    model.get_rng() = rng;

    //create persons at home
    // const auto mean_home_size    = 5.0;
    // const auto min_home_size     = 1;
    // auto& home_size_distribution = mio::PoissonDistribution<int>::get_instance();
    // auto home                    = model.add_location(mio::abm::LocationType::Home);
    // auto planned_home_size       = home_size_distribution(model.get_rng(), mean_home_size);
    // auto home_size               = 0;
    // for (size_t i = 0; i < num_persons; ++i) {
    //     if (home_size >= std::max(min_home_size, planned_home_size)) {
    //         home              = model.add_location(mio::abm::LocationType::Home);
    //         planned_home_size = home_size_distribution(model.get_rng(), mean_home_size);
    //         home_size         = 0;
    //     }

    //     auto age    = mio::AgeGroup(mio::UniformIntDistribution<size_t>::get_instance()(
    //         model.get_rng(), size_t(0), model.parameters.get_num_groups() - 1));
    //     auto person = model.add_person(home, age);
    //     model.assign_location(uint32_t(i), home);
    //     home_size++;
    // }

    //create other locations
    // for (auto loc_type :
    //      {mio::abm::LocationType::School, mio::abm::LocationType::Work, mio::abm::LocationType::SocialEvent,
    //       mio::abm::LocationType::BasicsShop, mio::abm::LocationType::Hospital, mio::abm::LocationType::ICU}) {

    //     const auto num_locs = std::max(size_t(1), num_persons / 100);
    //     std::vector<mio::abm::LocationId> locs(num_locs);
    //     std::generate(locs.begin(), locs.end(), [&] {
    //         return model.add_location(loc_type);
    //     });
    //     for (size_t p = 0; p < num_persons; ++p) {
    //         auto loc_idx =
    //             mio::UniformIntDistribution<size_t>::get_instance()(model.get_rng(), size_t(0), num_locs - 1);
    //         model.assign_location(uint32_t(p), locs[loc_idx]);
    //     }
    // }

    //infections and masks
    for (auto& person : model.get_persons()) {
        auto prng = mio::abm::PersonalRandomNumberGenerator(person);
        //some % of people are infected, large enough to have some infection activity without everyone dying
        auto pct_infected = 0.05;
        if (mio::UniformDistribution<ScalarType>::get_instance()(prng, 0.0, 1.0) < pct_infected) {
            auto infection =
                mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype, person.get_age(), model.parameters,
                                    mio::abm::TimePoint(0), mio::abm::InfectionState::Exposed);
            person.add_new_infection(std::move(infection));
        }
    }

    return mio::abm::Simulation(mio::abm::TimePoint(0), std::move(model));
}

/**
 * Benchmark for the ABM simulation.
 * @param num_persons Number of persons in the simulation.
 * @param seeds Seeds for the random number generator.
 */
void abm_benchmark(benchmark::State& state, size_t num_persons, std::initializer_list<uint32_t> seeds)
{
    mio::set_log_level(mio::LogLevel::warn);

    for (auto&& _ : state) {
        state.PauseTiming(); //exclude the setup from the benchmark
        auto sim = make_simulation(num_persons, seeds);
        state.ResumeTiming();

        //simulated time should be long enough to have full infection runs and mobility to every location
        auto final_time = sim.get_time() + mio::abm::days(5);
        mio::History<mio::DataWriterToMemory, mio::abm::LogDataForMobility> history;
        sim.advance(final_time);

        //debug output can be enabled to check for unexpected results (e.g. infections dieing out)
        //normally should have no significant effect on runtime
        // const bool monitor_infection_activity = false;
        // if constexpr (monitor_infection_activity) {
        //     std::cout << "num_persons = " << num_persons << "\n";
        //     for (auto inf_state = 0; inf_state < (int)mio::abm::InfectionState::Count; inf_state++) {
        //         std::cout << "inf_state = " << inf_state << ", sum = "
        //                   << sim.get_model().get_subpopulation_combined(sim.get_time(),
        //                                                                 mio::abm::InfectionState(inf_state))
        //                   << "\n";
        //     }
        // }
    }
}

//Measure ABM simulation run time with different sizes and different seeds.
//Fixed RNG seeds to make runs comparable. When there are code changes, the simulation will still
//run differently due to different sequence of random numbers being drawn. But for large enough sizes
//RNG should average out, so runs should be comparable even with code changes.
//We run a few different benchmarks to hopefully catch abnormal cases. Then seeds may
//have to be adjusted to get the benchmark back to normal.
//For small sizes (e.g. 10k) extreme cases are too likely, i.e. infections die out
//or overwhelm everything, so we don't benchmark these. Results should be mostly transferrable.

int main(int argc, char** argv)
{
    // Default problem size
    size_t num_persons = 1000000;

    //print omp_threads
#ifdef MEMILIO_ENABLE_OPENMP
    int omp_threads = 1;
#pragma omp parallel
    {
#pragma omp single
        omp_threads = omp_get_num_threads();
    }
    std::cout << "Running ABM benchmark with " << omp_threads << " OpenMP threads.\n";
#else
    std::cout << "Running ABM benchmark without OpenMP.\n";
#endif

    // Parse custom arguments for problem size BEFORE benchmark::Initialize
    // Remove custom args from argv to prevent benchmark from seeing them
    std::vector<char*> filtered_argv;
    filtered_argv.push_back(argv[0]); // Keep program name

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.find("--num_persons=") == 0) {
            num_persons = std::stoul(arg.substr(14));
            // Don't add this to filtered_argv
        }
        else {
            // Keep other arguments for benchmark
            filtered_argv.push_back(argv[i]);
        }
    }

    // Update argc to reflect filtered arguments
    int filtered_argc = static_cast<int>(filtered_argv.size());

    // Register the benchmark with the specified problem size
    std::string benchmark_name = "abm_benchmark_" + std::to_string(num_persons);
    benchmark::RegisterBenchmark(benchmark_name.c_str(), [num_persons](benchmark::State& state) {
        abm_benchmark(state, num_persons, {1415921265u, 35897932u});
    })->Unit(benchmark::kMillisecond);

    // Initialize and run benchmarks with filtered arguments
    benchmark::Initialize(&filtered_argc, filtered_argv.data());
    if (benchmark::ReportUnrecognizedArguments(filtered_argc, filtered_argv.data())) {
        return 1;
    }
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}