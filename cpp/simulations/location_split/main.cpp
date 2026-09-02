/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Sascha Korf
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

/**
 * @file main.cpp
 * @brief Multi run ABM simulation on a synthetic city.
 *
 * The city is generated from the demographic parameters of a region and the epidemic is started from
 * a number of randomly drawn infected Person%s. This is the port of cpp/examples/panvXabmSim from the
 * abmXpanvadere branch onto the ABM of main, without the Panvadere coupling: the location that the
 * branch split off and handed to an external microsimulation is not modelled separately here yet, the
 * initialization is random.
 */

#include "city_builder.h"
#include "defaults.h"
#include "file_utils.h"
#include "multi_run_simulator.h"

#include "memilio/utils/logging.h"
#include "memilio/utils/random_number_generator.h"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>

namespace
{

void print_help(const char* program_name)
{
    std::cout << "Usage: " << program_name << " [options]\n";
    std::cout << "Options:\n";
    std::cout << "  --region <name>            germany (default), france or usa\n";
    std::cout << "  --runs <number>            Number of runs (default: " << Config::DEFAULT_RUNS
              << ", range: " << Config::MIN_RUNS << "-" << Config::MAX_RUNS << ")\n";
    std::cout << "  --days <number>            Simulated days (default: " << Config::DEFAULT_DAYS
              << ", range: " << Config::MIN_DAYS << "-" << Config::MAX_DAYS << ")\n";
    std::cout << "  --n_persons <number>       Total population (default: " << Config::DEFAULT_POPULATION
              << ", range: " << Config::MIN_POPULATION << "-" << Config::MAX_POPULATION << ")\n";
    std::cout << "  --initial_infections <n>   Randomly drawn initially infected persons (default: "
              << Config::DEFAULT_INITIAL_INFECTIONS << ")\n";
    std::cout << "  --infection_k <number>     Scale of the viral shed (default: " << Config::DEFAULT_INFECTION_K
              << ", range: " << Config::MIN_INFECTION_K << "-" << Config::MAX_INFECTION_K << ")\n";
    std::cout << "  --output_dir <path>        Output directory (default: " << Config::DEFAULT_OUTPUT_DIR
              << "/run_<timestamp>)\n";
    std::cout << "  --seed <number>            Seed for reproducibility (default: built in seed)\n";
    std::cout << "  --no_infection_events      Do not write the per run infection event csv\n";
    std::cout << "  --write_contacts           Write the pairwise contact hours per run (expensive)\n";
    std::cout << "  --print_city               Print the derived city infrastructure and exit\n";
    std::cout << "  --help                     Show this message\n";
    std::cout << "\nExamples:\n";
    std::cout << "  " << program_name << " --runs 50 --days 14\n";
    std::cout << "  " << program_name << " --region france --n_persons 10000 --initial_infections 25\n";
}

/// @brief Parse an argument that must lie in [min, max]. Exits on error.
template <class T>
T parse_number(const std::string& option, const char* text, T min, T max)
{
    T value;
    try {
        if constexpr (std::is_floating_point_v<T>) {
            value = static_cast<T>(std::stod(text));
        }
        else {
            value = static_cast<T>(std::stoll(text));
        }
    }
    catch (const std::exception&) {
        std::cerr << "Error: invalid number for " << option << ": " << text << "\n";
        std::exit(1);
    }
    if (value < min || value > max) {
        std::cerr << "Error: " << option << " must be between " << min << " and " << max << "\n";
        std::exit(1);
    }
    return value;
}

struct CommandLine {
    MultiRunConfig config{};
    bool print_city_only = false;
};

CommandLine parse_command_line(int argc, char* argv[])
{
    CommandLine cli;
    auto& config = cli.config;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        const bool has_value  = i + 1 < argc;

        if (arg == "--help") {
            print_help(argv[0]);
            std::exit(0);
        }
        else if (arg == "--region" && has_value) {
            if (!CityParameters::region_from_string(argv[++i], config.city_config.region)) {
                std::cerr << "Error: unknown region '" << argv[i] << "'. Valid: germany, france, usa\n";
                std::exit(1);
            }
        }
        else if (arg == "--runs" && has_value) {
            config.num_runs = parse_number<int>(arg, argv[++i], Config::MIN_RUNS, Config::MAX_RUNS);
        }
        else if (arg == "--days" && has_value) {
            config.simulation_days = parse_number<int>(arg, argv[++i], Config::MIN_DAYS, Config::MAX_DAYS);
        }
        else if (arg == "--n_persons" && has_value) {
            config.city_config.total_population =
                parse_number<int>(arg, argv[++i], Config::MIN_POPULATION, Config::MAX_POPULATION);
        }
        else if (arg == "--initial_infections" && has_value) {
            config.num_initial_infections = parse_number<int>(arg, argv[++i], 0, Config::MAX_POPULATION);
        }
        else if (arg == "--infection_k" && has_value) {
            config.infection_parameter_k =
                parse_number<double>(arg, argv[++i], Config::MIN_INFECTION_K, Config::MAX_INFECTION_K);
        }
        else if (arg == "--output_dir" && has_value) {
            config.output_base_dir = argv[++i];
        }
        else if (arg == "--seed" && has_value) {
            config.custom_seed = static_cast<uint32_t>(
                parse_number<long long>(arg, argv[++i], 1, std::numeric_limits<uint32_t>::max()));
        }
        else if (arg == "--no_infection_events") {
            config.write_infection_events = false;
        }
        else if (arg == "--write_contacts") {
            config.write_contacts = true;
        }
        else if (arg == "--print_city") {
            cli.print_city_only = true;
        }
        else {
            std::cerr << "Error: unknown argument '" << arg << "'\n";
            print_help(argv[0]);
            std::exit(1);
        }
    }

    if (config.num_initial_infections > config.city_config.total_population) {
        std::cerr << "Error: --initial_infections must not exceed --n_persons\n";
        std::exit(1);
    }

    return cli;
}

void print_config_summary(const MultiRunConfig& config)
{
    std::cout << "\n=== Multi-Run Simulation Setup ===" << std::endl;
    std::cout << "Region: " << CityParameters::region_to_string(config.city_config.region) << std::endl;
    std::cout << "Population: " << config.city_config.total_population << std::endl;
    std::cout << "Number of runs: " << config.num_runs << std::endl;
    std::cout << "Simulation days: " << config.simulation_days << std::endl;
    std::cout << "Initial infections: " << config.num_initial_infections << std::endl;
    std::cout << "Infection parameter K: " << config.infection_parameter_k << std::endl;
    std::cout << "Output directory: " << config.output_base_dir << std::endl;
    std::cout << "==================================" << std::endl;
}

mio::IOResult<void> main_flow(int argc, char* argv[])
{
    mio::set_log_level(mio::LogLevel::critical);

    auto cli     = parse_command_line(argc, argv);
    auto& config = cli.config;

    if (cli.print_city_only) {
        CityBuilder::print_city_summary(config.city_config);
        return mio::success();
    }

    auto rng           = mio::RandomNumberGenerator();
    const uint32_t s   = (config.custom_seed != 0) ? config.custom_seed : 12345678u;
    rng.seed({s, s + 1, s + 2, s + 3, s + 4, s + 5});
    rng.synchronize();

    if (config.output_base_dir == Config::DEFAULT_OUTPUT_DIR) {
        config.output_base_dir = config.output_base_dir + "/run_" + current_date_time();
    }

    print_config_summary(config);
    const auto start_time = std::chrono::high_resolution_clock::now();

    BOOST_OUTCOME_TRY(auto results, MultiRunSimulator::run_multi_simulation(config, rng));
    BOOST_OUTCOME_TRY(MultiRunSimulator::save_multi_run_results(results, config.output_base_dir, config));

    const auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(
        std::chrono::high_resolution_clock::now() - start_time);

    std::cout << "\n=== Simulation Summary ===" << std::endl;
    std::cout << "Total execution time: " << duration.count() << " seconds" << std::endl;
    std::cout << "Successful runs: " << results.successful_runs << "/" << config.num_runs << std::endl;
    std::cout << "Results saved to: " << config.output_base_dir << std::endl;
    std::cout << "==========================" << std::endl;

    return mio::success();
}

} // namespace

int main(int argc, char* argv[])
{
    auto result = main_flow(argc, argv);
    if (result.has_error()) {
        std::cerr << "Simulation failed: " << result.error().message() << std::endl;
        return 1;
    }
    return 0;
}
