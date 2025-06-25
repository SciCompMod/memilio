#include <iostream>
#include <string>
#include <chrono>
#include "boost/filesystem.hpp"
#include "memilio/io/result_io.h"

#include "include/multi_run_simulator.h"
#include "include/event_simulator.h"
#include "include/file_utils.h"

namespace fs = boost::filesystem;

MultiRunConfig parse_multi_run_config(int argc, char* argv[])
{
    MultiRunConfig config;

    // Set defaults
    config.city_config                     = CityConfig{};
    config.event_config.type               = EventType::Restaurant;
    config.event_config.use_panvadere_init = true;
    config.event_config.panvadere_file     = "./data/restaurant/infections.txt";
    config.num_runs                        = 100;
    config.simulation_days                 = 30;

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--event" && i + 1 < argc) {
            std::string event_type = argv[++i];
            if (event_type == "restaurant")
                config.event_config.type = EventType::Restaurant;
            else if (event_type == "work")
                config.event_config.type = EventType::WorkMeeting;
            else if (event_type == "choir")
                config.event_config.type = EventType::Choir;
        }
        else if (arg == "--init" && i + 1 < argc) {
            std::string init_type = argv[++i];
            if (init_type == "panvadere") {
                config.event_config.use_panvadere_init = true;
                config.event_config.panvadere_file     = "./data/restaurant/infections.txt";
            }
            else if (init_type == "uniform") {
                config.event_config.use_panvadere_init = false;
            }
        }
        else if (arg == "--runs" && i + 1 < argc) {
            config.num_runs = std::stoi(argv[++i]);
        }
        else if (arg == "--days" && i + 1 < argc) {
            config.simulation_days = std::stoi(argv[++i]);
        }
        else if (arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [options]\n";
            std::cout << "Options:\n";
            std::cout << "  --event <type>        Set the event type (restaurant, work, choir)\n";
            std::cout << "  --init <type>         Set the initialization type (panvadere, uniform)\n";
            std::cout << "  --runs <number>       Set the number of runs\n";
            std::cout << "  --days <number>       Set the number of simulation days\n";
            std::cout << "  --help                Show this help message\n";
            return config;
        }
    }

    return config;
}

void print_summary(const MultiRunResults& results)
{
    std::cout << "\n=== Simulation Summary ===" << std::endl;
    std::cout << "Event Type: " << EventSimulator::event_type_to_string(results.event_type) << std::endl;
    std::cout << "Initialization: " << (results.used_panvadere_init ? "Panvadere" : "Event Simulation") << std::endl;
    std::cout << "Infection Parameter K: " << results.infection_parameter_k << std::endl;
    std::cout << "Successful Runs: " << results.successful_runs << std::endl;
    std::cout << "=========================" << std::endl;
}

mio::IOResult<void> main_flow(int argc, char* argv[])
{
    auto start_time = std::chrono::high_resolution_clock::now();

    // Parse configuration
    auto config = parse_multi_run_config(argc, argv);

    std::cout << "=== Multi-Run Simulation Setup ===" << std::endl;
    std::cout << "Event Type: " << EventSimulator::event_type_to_string(config.event_config.type) << std::endl;
    std::cout << "Initialization: " << (config.event_config.use_panvadere_init ? "Panvadere" : "Event Simulation")
              << std::endl;
    std::cout << "Number of runs: " << config.num_runs << std::endl;
    std::cout << "Population: " << config.city_config.total_population << std::endl;
    std::cout << "Simulation days: " << config.simulation_days << std::endl;

    // Run multi-simulation
    BOOST_OUTCOME_TRY(auto results, MultiRunSimulator::run_multi_simulation(config));

    // Save results
    std::string result_dir = config.output_base_dir + "/results_" + currentDateTime();
    BOOST_OUTCOME_TRY(MultiRunSimulator::save_multi_run_results(results, result_dir));

    auto end_time       = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

    print_summary(results);
    std::cout << "\n✓ Multi-run simulation completed!" << std::endl;
    std::cout << "⏱️  Total time: " << total_duration.count() << "s" << std::endl;
    std::cout << "📁 Results: " << result_dir << std::endl;

    return mio::success();
}

int main(int argc, char* argv[])
{
    auto result = main_flow(argc, argv);
    if (result.has_error()) {
        std::cerr << "❌ Multi-run simulation failed: " << result.error().message() << std::endl;
        return 1;
    }
    return 0;
}