#include <iostream>
#include <string>
#include <chrono>
#include "boost/filesystem.hpp"
#include "memilio/io/result_io.h"

#include "include/config_validation.h"
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
    {"restaurant_table_equals_random_household", EventType::Restaurant_Table_Equals_Random},
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
    std::cout << "  --runs <number>       Set the number of runs (default: " << Config::DEFAULT_RUNS
              << ", range: " << Config::MIN_RUNS << "-" << Config::MAX_RUNS << ")\n";
    std::cout << "  --days <number>       Set the number of simulation days (default: " << Config::DEFAULT_DAYS
              << ", range: " << Config::MIN_DAYS << "-" << Config::MAX_DAYS << ")\n";
    std::cout << "  --n_persons <number>  Set the total population (default: " << Config::DEFAULT_POPULATION
              << ", range: " << Config::MIN_POPULATION << "-" << Config::MAX_POPULATION << ")\n";
    std::cout << "  --output_dir <path>   Set the output directory (default: " << Config::DEFAULT_OUTPUT_DIR << ")\n";
    std::cout << "                        The directory will be created if it does not exist.\n";
    std::cout << "                        If not set, results will be saved in the current directory.\n";
    std::cout << "  --help                Show this help message\n";
    std::cout << "\nExamples:\n";
    std::cout << "  " << program_name << " --event restaurant_table_equals_household --runs 50\n";
    std::cout << "  " << program_name << " --sim panvadere --days 14 --n_persons 10000\n";
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
    config.city_config                       = CityConfig{};
    config.city_config.total_population      = Config::DEFAULT_POPULATION;
    config.event_config.type                 = EventType::Restaurant_Table_Equals_Household;
    config.event_config.event_duration_hours = Config::DEFAULT_EVENT_HOURS;
    config.simulation_type                   = SimType::Memilio; // Default to Panvadere simulatio
    config.num_runs                          = Config::DEFAULT_RUNS;
    config.simulation_days                   = Config::DEFAULT_DAYS;
    config.output_base_dir                   = Config::DEFAULT_OUTPUT_DIR;

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
                if (config.num_runs < Config::MIN_RUNS || config.num_runs > Config::MAX_RUNS) {
                    std::cerr << "Error: Number of runs must be between " << Config::MIN_RUNS << " and "
                              << Config::MAX_RUNS << "\n";
                    exit(1);
                }
            }
            catch (const std::exception& e) {
                std::cerr << "Error: Invalid number for --runs: " << argv[i] << "\n";
                exit(1);
            }
        }
        else if (arg == "--days" && i + 1 < argc) {
            try {
                config.simulation_days = std::stoi(argv[++i]);
                if (config.simulation_days < Config::MIN_DAYS || config.simulation_days > Config::MAX_DAYS) {
                    std::cerr << "Error: Number of days must be between " << Config::MIN_DAYS << " and "
                              << Config::MAX_DAYS << "\n";
                    exit(1);
                }
            }
            catch (const std::exception& e) {
                std::cerr << "Error: Invalid number for --days: " << argv[i] << "\n";
                exit(1);
            }
        }
        else if (arg == "--n_persons" && i + 1 < argc) {
            try {
                config.city_config.total_population = std::stoi(argv[++i]);
                if (config.city_config.total_population < Config::MIN_POPULATION ||
                    config.city_config.total_population > Config::MAX_POPULATION) {
                    std::cerr << "Error: Population must be between " << Config::MIN_POPULATION << " and "
                              << Config::MAX_POPULATION << "\n";
                    exit(1);
                }
            }
            catch (const std::exception& e) {
                std::cerr << "Error: Invalid number for --n_persons: " << argv[i] << "\n";
                exit(1);
            }
        }
        else if (arg == "--output_dir" && i + 1 < argc) {
            config.output_base_dir = argv[++i];
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

bool validate_config(const MultiRunConfig& config)
{
    // Validate event configuration
    // if (!config.event_config.is_valid()) {
    //     std::cerr << "Error: Invalid event configuration\n";
    //     return false;
    // }

    // Validate simulation parameters
    if (config.num_runs < Config::MIN_RUNS || config.num_runs > Config::MAX_RUNS) {
        std::cerr << "Error: Invalid number of runs\n";
        return false;
    }

    if (config.simulation_days < Config::MIN_DAYS || config.simulation_days > Config::MAX_DAYS) {
        std::cerr << "Error: Invalid number of simulation days\n";
        return false;
    }

    if (config.city_config.total_population < Config::MIN_POPULATION ||
        config.city_config.total_population > Config::MAX_POPULATION) {
        std::cerr << "Error: Invalid population size\n";
        return false;
    }

    return true;
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
    std::cout << "Infection Parameter K: " << results.infection_parameter_k << std::endl;
    std::cout << "=========================" << std::endl;
}

mio::IOResult<void> main_flow(int argc, char* argv[])
{
    mio::set_log_level(mio::LogLevel::critical);
    auto start_time = std::chrono::high_resolution_clock::now();

    // Parse configuration
    auto config = parse_multi_run_config(argc, argv);

    // Validate complete configuration
    if (!validate_config(config)) {
        return mio::failure(mio::StatusCode::InvalidValue, "Configuration validation failed");
    }

    print_config_summary(config);

    // Run multi-simulation
    BOOST_OUTCOME_TRY(auto results, MultiRunSimulator::run_multi_simulation(config));

    // Save results
    if (config.output_base_dir == Config::DEFAULT_OUTPUT_DIR) {
        config.output_base_dir = config.output_base_dir + "/results_" + currentDateTime();
    }
    BOOST_OUTCOME_TRY(MultiRunSimulator::save_multi_run_results(results, config.output_base_dir));

    // Print final results
    auto end_time       = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

    print_summary(results);

    std::cout << "\n✓ Multi-run simulation completed successfully!" << std::endl;
    std::cout << "Total execution time: " << total_duration.count() << " seconds" << std::endl;
    std::cout << "Results saved to: " << config.output_base_dir << std::endl;
    std::cout << "Successful runs: " << results.successful_runs << "/" << config.num_runs << std::endl;

    return mio::success();
}

int main(int argc, char* argv[])
{
    auto result = main_flow(argc, argv);
    if (result.has_error()) {
        std::cerr << "Multi-run simulation failed: " << result.error().message() << std::endl;
        return 1;
    }
    return 0;
}