#include "../include/multi_run_simulator.h"
#include "../include/file_utils.h"
#include <algorithm>
#include <iostream>
#include <fstream>

mio::IOResult<MultiRunResults> MultiRunSimulator::run_multi_simulation(const MultiRunConfig& config)
{
    MultiRunResults results;
    results.event_type      = config.event_config.type;
    results.simulation_type = config.event_config.simulation_type;

    // Step 1: Build city (done once)
    std::cout << "Building city..." << std::endl;
    BOOST_OUTCOME_TRY(auto base_world, CityBuilder::build_world(config.city_config));

    // Step 2: Get map from specific event to ids of persons in simulation
    std::cout << "Mapping events to person IDs..." << std::endl;
    BOOST_OUTCOME_TRY(auto event_map, EventSimulator::map_events_to_persons(base_world, config.event_config));
    results.event_person_map = event_map;

    // Step 3: Calculate K parameter
    std::cout << "Calculating infection parameter K..." << std::endl;
    BOOST_OUTCOME_TRY(auto k_value, EventSimulator::calculate_infection_parameter_k(config.event_config));
    results.infection_parameter_k = k_value;

    // Step 3: Get initial infections
    std::map<uint32_t, bool> initial_infections;
    if (config.event_config.use_panvadere_init) {
        std::cout << "Loading infections from Panvadere..." << std::endl;
        BOOST_OUTCOME_TRY(auto infections, EventSimulator::initialize_from_panvadere(config.event_config.panvadere_file,
                                                                                     config.event_config.type));
        initial_infections = infections;
    }
    else {
        std::cout << "Simulating event infections..." << std::endl;
        BOOST_OUTCOME_TRY(auto infections,
                          EventSimulator::initialize_from_event_simulation(config.event_config, base_world));
        initial_infections = infections;
    }

    // Step 4: Run multiple simulations
    std::cout << "Running " << config.num_runs << " simulations..." << std::endl;
    results.all_runs.reserve(config.num_runs);

    for (int run = 0; run < config.num_runs; ++run) {
        if (run % 10 == 0) {
            std::cout << "Run " << run << "/" << config.num_runs << std::endl;
        }

        auto single_result =
            run_single_simulation_with_infections(base_world, initial_infections, k_value, config.simulation_days);

        if (single_result.has_value()) {
            results.all_runs.push_back(single_result.value());
            results.successful_runs++;
        }
        else {
            std::cerr << "Run " << run << " failed: " << single_result.error().message() << std::endl;
        }
    }

    // Step 5: Calculate percentiles
    if (results.successful_runs > 0) {
        auto final_results =
            calculate_percentiles(results.all_runs, results.event_type, results.used_panvadere_init, k_value);
        results.percentile_25 = final_results.percentile_25;
        results.percentile_50 = final_results.percentile_50;
        results.percentile_75 = final_results.percentile_75;
    }

    return mio::success(results);
}

mio::IOResult<void> MultiRunSimulator::save_multi_run_results(const MultiRunResults& results,
                                                              const std::string& base_dir)
{
    // Create directory structure
    std::string event_dir   = base_dir + "/" + EventSimulator::event_type_to_string(results.event_type);
    std::string init_method = results.used_panvadere_init ? "panvadere" : "simulation";
    std::string full_dir    = event_dir + "_" + init_method;

    BOOST_OUTCOME_TRY(create_result_folders(full_dir));

    // Save percentile results
    std::string summary_file = full_dir + "/summary.txt";
    std::ofstream summary(summary_file);
    if (summary.is_open()) {
        summary << "Multi-Run Simulation Summary\n";
        summary << "Event Type: " << EventSimulator::event_type_to_string(results.event_type) << "\n";
        summary << "Initialization: " << init_method << "\n";
        summary << "Infection Parameter K: " << results.infection_parameter_k << "\n";
        summary << "Total Runs: " << results.all_runs.size() << "\n";
        summary.close();
    }

    // TODO: Save detailed percentile data
    // Save percentile 25, 50, 75 time series data

    return mio::success();
}

mio::IOResult<SimulationRunner::SimulationResults>
MultiRunSimulator::run_single_simulation_with_infections(const mio::abm::World& base_world,
                                                         const std::map<uint32_t, bool>& initial_infections,
                                                         double k_parameter, int simulation_days)
{
    // TODO: Implement single simulation run with initial infections
    // 1. Copy base world
    // 2. Apply initial infections to specified people/households
    // 3. Set infection parameter K
    // 4. Run simulation for specified days
    // 5. Return results

    SimulationRunner::SimulationResults results;

    // Placeholder implementation
    return mio::success(results);
}