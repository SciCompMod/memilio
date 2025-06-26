#include <iostream>
#include <string>
#include <chrono>
#include <unordered_map>
#include "boost/filesystem.hpp"
#include "memilio/io/result_io.h"

#include "include/multi_run_simulator.h"
#include "include/event_simulator.h"
#include "include/file_utils.h"

namespace fs = boost::filesystem;

namespace
{
// Event type mapping for cleaner parsing
const std::unordered_map<std::string, EventType> EVENT_TYPE_MAP = {
    {"restaurant_table_equals_household", EventType::Restaurant_Table_Equals_Household},
    {"restaurant_table_equals_half_household", EventType::Restaurant_Table_Equals_Half_Household},
    {"work_meeting_many", EventType::WorkMeeting_Many_Meetings},
    {"work_meeting_few", EventType::WorkMeeting_Few_Meetings}};

// Simulation type mapping
const std::unordered_map<std::string, SimType> SIM_TYPE_MAP = {
    {"panvadere", SimType::Panvadere}, {"memilio", SimType::Memilio}, {"both", SimType::Both}};
} // namespace

void print_help(const char* program_name)
{
    std::cout << "Usage: " << program_name << " [options]\n";
    std::cout << "Options:\n";
    std::cout << "  --event <type>        Set the event type\n";
    std::cout << "                        Valid types: restaurant_table_equals_household,\n";
    std::cout << "                                    restaurant_table_equals_half_household,\n";
    std::cout << "                                    work_meeting_many, work_meeting_few\n";
    std::cout << "  --sim <type>          Set the simulation type (panvadere, memilio, both)\n";
    std::cout << "  --runs <number>       Set the number of runs (default: 10)\n";
    std::cout << "  --days <number>       Set the number of simulation days (default: 30)\n";
    std::cout << "  --n_persons <number>  Set the total population\n";
    std::cout << "  --help                Show this help message\n";
}

bool parse_event_type(const std::string& event_type, MultiRunConfig& config)
{
    auto it = EVENT_TYPE_MAP.find(event_type);
    if (it != EVENT_TYPE_MAP.end()) {
        config.event_config.type = it->second;
        return true;
    }

    std::cerr << "Error: Unknown event type '" << event_type << "'\n";
    std::cerr << "Valid types: ";
    for (const auto& [key, value] : EVENT_TYPE_MAP) {
        std::cerr << key << " ";
    }
    std::cerr << std::endl;
    return false;
}

bool parse_sim_type(const std::string& sim_type, MultiRunConfig& config)
{
    auto it = SIM_TYPE_MAP.find(sim_type);
    if (it != SIM_TYPE_MAP.end()) {
        config.simulation_type = it->second;
        // Note: panvadere_file is now automatically determined based on event type
        return true;
    }

    std::cerr << "Error: Unknown simulation type '" << sim_type << "'\n";
    std::cerr << "Valid types: panvadere, memilio, both\n";
    return false;
}

MultiRunConfig parse_multi_run_config(int argc, char* argv[])
{
    MultiRunConfig config;

    // Set defaults
    config.city_config       = CityConfig{};
    config.event_config.type = EventType::Restaurant_Table_Equals_Household;
    config.simulation_type   = SimType::Both;
    config.num_runs          = 10;
    config.simulation_days   = 30;

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--event" && i + 1 < argc) {
            if (!parse_event_type(argv[++i], config)) {
                exit(1);
            }
        }
        else if (arg == "--sim" && i + 1 < argc) {
            if (!parse_sim_type(argv[++i], config)) {
                exit(1);
            }
        }
        else if (arg == "--runs" && i + 1 < argc) {
            try {
                config.num_runs = std::stoi(argv[++i]);
                if (config.num_runs <= 0) {
                    std::cerr << "Error: Number of runs must be positive\n";
                    exit(1);
                }
            }
            catch (const std::exception& e) {
                std::cerr << "Error: Invalid number for --runs\n";
                exit(1);
            }
        }
        else if (arg == "--days" && i + 1 < argc) {
            try {
                config.simulation_days = std::stoi(argv[++i]);
                if (config.simulation_days <= 0) {
                    std::cerr << "Error: Number of days must be positive\n";
                    exit(1);
                }
            }
            catch (const std::exception& e) {
                std::cerr << "Error: Invalid number for --days\n";
                exit(1);
            }
        }
        else if (arg == "--n_persons" && i + 1 < argc) {
            try {
                config.city_config.total_population = std::stoi(argv[++i]);
                if (config.city_config.total_population <= 0) {
                    std::cerr << "Error: Population must be positive\n";
                    exit(1);
                }
            }
            catch (const std::exception& e) {
                std::cerr << "Error: Invalid number for --n_persons\n";
                exit(1);
            }
        }
        else if (arg == "--help") {
            print_help(argv[0]);
            exit(0);
        }
        else {
            std::cerr << "Error: Unknown argument '" << arg << "'\n";
            print_help(argv[0]);
            exit(1);
        }
    }

    return config;
}

void print_config_summary(const MultiRunConfig& config)
{
    std::cout << "\n=== Multi-Run Simulation Setup ===" << std::endl;
    std::cout << "Event Type: " << EventSimulator::event_type_to_string(config.event_config.type) << std::endl;
    std::cout << "Simulation Type: " << EventSimulator::simulation_type_to_string(config.simulation_type) << std::endl;
    std::cout << "Number of runs: " << config.num_runs << std::endl;
    std::cout << "Population: " << config.city_config.total_population << std::endl;
    std::cout << "Simulation days: " << config.simulation_days << std::endl;
    std::cout << "===================================" << std::endl;
}

void print_summary(const MultiRunResults& results)
{
    std::cout << "\n=== Simulation Summary ===" << std::endl;
    std::cout << "Event Type: " << EventSimulator::event_type_to_string(results.event_type) << std::endl;
    std::cout << "Simulation Type: " << EventSimulator::simulation_type_to_string(results.simulation_type) << std::endl;
    std::cout << "Infection Parameter K: " << results.infection_parameter_k << std::endl;
    std::cout << "=========================" << std::endl;
}

mio::IOResult<void> main_flow(int argc, char* argv[])
{
    auto start_time = std::chrono::high_resolution_clock::now();

    // Parse configuration
    auto config = parse_multi_run_config(argc, argv);
    print_config_summary(config);

    // Run multi-simulation
    BOOST_OUTCOME_TRY(auto results, MultiRunSimulator::run_multi_simulation(config));

    // Save results
    std::string result_dir = config.output_base_dir + "/results_" + currentDateTime();
    BOOST_OUTCOME_TRY(MultiRunSimulator::save_multi_run_results(results, result_dir));

    // Print final results
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