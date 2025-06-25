#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include "boost/filesystem.hpp"
#include "memilio/io/result_io.h"

// Your project includes
#include "include/config_manager.h"
#include "include/simulation_runner.h"
#include "include/file_utils.h"

namespace fs     = boost::filesystem;
using ScalarType = double;

mio::IOResult<void> main_flow(int argc, char* argv[])
{
    auto start_time = std::chrono::high_resolution_clock::now();

    // Load configuration
    auto config = ConfigManager::load_config_from_args(argc, argv);
    mio::set_log_level(mio::LogLevel::err);

    // Quick validation
    if (config.n_persons <= 0 || config.simulation_days <= 0) {
        return mio::failure(mio::StatusCode::InvalidValue, "Invalid simulation parameters");
    }

    if (!boost::filesystem::exists(config.infection_data_file)) {
        return mio::failure(mio::StatusCode::InvalidFileFormat,
                            "Infection data file not found: " + config.infection_data_file);
    }

    std::cout << "=== Simulation Setup ===" << std::endl;
    std::cout << "Input file: " << config.infection_data_file << std::endl;
    std::cout << "Persons: " << config.n_persons << std::endl;
    std::cout << "Days: " << config.simulation_days << std::endl;

    std::cout << "\n[1/4] Setting up directories..." << std::endl;
    // Setup result directories
    std::string result_dir = config.input_dir + "/results_" + currentDateTime();
    BOOST_OUTCOME_TRY(create_result_folders(result_dir));

    std::cout << "[2/4] Running simulation..." << std::endl;
    // Run simulation
    BOOST_OUTCOME_TRY(auto results, SimulationRunner::run_simulation(config));

    std::cout << "[3/4] Saving results..." << std::endl;
    // Save results
    BOOST_OUTCOME_TRY(SimulationRunner::save_simulation_results(results, result_dir));

    std::cout << "[4/4] Finalizing..." << std::endl;
    // Copy to last run directory
    std::string last_run_dir = config.input_dir + "/results_last_run";
    BOOST_OUTCOME_TRY(copy_result_folder(result_dir, last_run_dir));

    std::string summary_file = result_dir + "/panvXabm_results.txt";
    BOOST_OUTCOME_TRY(SimulationRunner::write_summary_output(results, summary_file));

    auto end_time       = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

    std::cout << "\n✓ Simulation completed successfully!" << std::endl;
    std::cout << "Total time: " << total_duration.count() << "s" << std::endl;
    std::cout << "Results: " << result_dir << std::endl;

    return mio::success();
}

int main(int argc, char* argv[])
{
    auto result = main_flow(argc, argv);
    if (result.has_error()) {
        std::cerr << "Simulation failed: " << result.error().message() << std::endl;
        std::cerr << "Try --help for usage information" << std::endl;
        return 1;
    }
    return 0;
}