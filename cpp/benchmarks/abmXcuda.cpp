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
#include "memilio/io/history.h"
#include "benchmark/benchmark.h"

template <typename T>
void write_log_to_file(const T& history)
{
    auto logg = history.get_log();
    // Write the results to a file.
    auto time_at_location = std::get<0>(logg);
    std::string input;
    std::ofstream myfile("test_output2.txt");
    myfile << "time_at_location:\n";
    for (auto i = 0; i < time_at_location.size(); ++i) {
        myfile << "Next Timepoint " << i << "\n";
        for (auto j = 0; j< time_at_location[i].size(); ++j) {
            auto id = time_at_location[i][j];
            myfile << id << "\n";
        }
        
    }
        
    myfile << input << "\n";
    std::cout << "Looger logged to file\n";

    myfile.close();
}



mio::abm::Simulation<> make_simulation(size_t num_persons, size_t num_p_p_loc, std::initializer_list<uint32_t> seeds, double pct_infected)
{
    auto rng = mio::RandomNumberGenerator();
    rng.seed(seeds);
    

    //create persons at home
    const size_t mean_home_size    = 5;
    const size_t amount_of_home_size = num_persons / mean_home_size;

    auto model      = mio::abm::Model(1);
    model.get_rng() = rng;


    auto home = model.add_location(mio::abm::LocationType::Home);
    auto const age = mio::AgeGroup( model.parameters.get_num_groups() - 1);
    for (size_t i = 0; i < num_persons; ++i) {    
        if (i % amount_of_home_size == 0 && i < amount_of_home_size) {
            home = model.add_location(mio::abm::LocationType::Home);
        }
        auto person = model.add_person(home, age);
        model.assign_location(uint32_t(i), home);   
    }

    // create other locations
    for (auto loc_type :
         {mio::abm::LocationType::SocialEvent, mio::abm::LocationType::Hospital, mio::abm::LocationType::ICU, mio::abm::LocationType::Work, mio::abm::LocationType::School, mio::abm::LocationType::BasicsShop}) {
        const auto num_locs = std::max(size_t(1), num_persons / num_p_p_loc);
        std::vector<mio::abm::LocationId> locs(num_locs);
        std::generate(locs.begin(), locs.end(), [&] {
            return model.add_location(loc_type);
        });
        for (size_t p = 0; p < num_persons; ++p) {
            auto loc_idx =
                mio::UniformIntDistribution<size_t>::get_instance()(model.get_rng(), size_t(0), num_locs - 1);
            model.assign_location(uint32_t(p), locs[loc_idx]);
        }
    }

    // infections and masks
    for (auto& person : model.get_persons()) {
        auto prng = mio::abm::PersonalRandomNumberGenerator(person);
        if (mio::UniformDistribution<double>::get_instance()(prng, 0.0, 1.0) < pct_infected) {
            auto state = mio::abm::InfectionState(
                mio::UniformIntDistribution<int>::get_instance()(prng, 1, int(mio::abm::InfectionState::Count) - 1));
            auto infection = mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                 model.parameters, mio::abm::TimePoint(0), state);
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
void abm_benchmark(benchmark::State& state, size_t num_persons, size_t num_p_p_loc, std::initializer_list<uint32_t> seeds)
{
    mio::set_log_level(mio::LogLevel::warn);

    for (auto&& _ : state) {
        state.PauseTiming(); //exclude the setup from the benchmark
        auto sim = make_simulation(num_persons,num_p_p_loc, seeds, 0.5);
        state.ResumeTiming();

        mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
            Eigen::Index(mio::abm::InfectionState::Count)};

        //simulated time should be long enough to have full infection runs and mobility to every location
        auto final_time = sim.get_time() + mio::abm::hours(24);
        sim.advance(final_time);
        // std::ofstream outfile("abm_minimal.txt");
        // std::get<0>(historyTimeSeries.get_log())
        //     .print_table({"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4, outfile);
        // std::cout << "Results written to abm_minimal.txt" << std::endl;
    }
}

// void abm_benchmark(size_t num_persons, size_t num_p_p_loc, std::initializer_list<uint32_t> seeds)
// {
//     mio::set_log_level(mio::LogLevel::warn);
//     auto sim = make_simulation(num_persons,num_p_p_loc, seeds);

//     //simulated time should be long enough to have full infection runs and mobility to every location
//     auto final_time = sim.get_time() + mio::abm::hours(24);
//     sim.advance(final_time);
// }

// int main(){
//     abm_benchmark(100,100, {28841971u, 69399375u});
//     return 0;
// }

//Measure ABM simulation run time with different sizes and different seeds.
//Fixed RNG seeds to make runs comparable. When there are code changes, the simulation will still
//run differently due to different sequence of random numbers being drawn. But for large enough sizes
//RNG should average out, so runs should be comparable even with code changes.
//We run a few different benchmarks to hopefully catch abnormal cases. Then seeds may
//have to be adjusted to get the benchmark back to normal.
//For small sizes (e.g. 10k) extreme cases are too likely, i.e. infections die out
//or overwhelm everything, so we don't benchmark these. Results should be mostly transferrable.
// BENCHMARK_CAPTURE(abm_benchmark, abm_benchmark_50k, 50000,100, {14159265u, 35897932u})->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(abm_benchmark, abm_benchmark_100k, 100000,100, {38462643u, 38327950u})->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(abm_benchmark, abm_benchmark_200k, 1000000,100, {28841971u, 69399375u})->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
