#include "../include/multi_run_simulator.h"
#include "../include/file_utils.h"
#include "abm/time.h"
#include <algorithm>
#include <iostream>
#include <fstream>

mio::IOResult<MultiRunResults> MultiRunSimulator::run_multi_simulation(const MultiRunConfig& config)
{
    MultiRunResults results;
    results.event_type      = config.event_config.type;
    results.simulation_type = config.simulation_type;

    // Step 1: Build city (done once)
    std::cout << "Building city..." << std::endl;
    BOOST_OUTCOME_TRY(auto base_world, CityBuilder::build_world(config.city_config));

    // Step 2: Get map from specific event to ids of persons in simulation
    std::cout << "Mapping events to person IDs..." << std::endl;
    //TODO: Implement this function to map events to persons
    // BOOST_OUTCOME_TRY(auto event_map, EventSimulator::map_events_to_persons(base_world, config.event_config));
    // results.event_person_map = event_map;

    // Step 3: Calculate K parameter
    std::cout << "Calculating infection parameter K..." << std::endl;
    //TODO: Implement this function to calculate K parameter based on event type
    BOOST_OUTCOME_TRY(auto k_value, EventSimulator::calculate_infection_parameter_k(config.event_config));
    results.infection_parameter_k = k_value;

    // Step 3: Get initial infections
    std::map<uint32_t, bool> initial_infections{};
    initial_infections = {}; // Initialize empty map
    // TODO: Implement this.
    // For now just infect the first 10 persons in the base world
    const auto& persons = base_world.get_persons();
    uint32_t count      = 0;
    for (const auto& person : persons) {
        initial_infections[person.get_person_id()] = true;
        count++;
        if (count >= 10) {
            break; // Limit to first 10 persons
        }
    }
    // if (config.event_config.use_panvadere_init) {
    //     std::cout << "Loading infections from Panvadere..." << std::endl;
    //     BOOST_OUTCOME_TRY(auto infections, EventSimulator::initialize_from_panvadere(config.event_config.panvadere_file,
    //                                                                                  config.event_config.type));
    //     initial_infections = infections;
    // }
    // else {
    //     std::cout << "Simulating event infections..." << std::endl;
    //     BOOST_OUTCOME_TRY(auto infections,
    //                       EventSimulator::initialize_from_event_simulation(config.event_config, base_world));
    //     initial_infections = infections;
    // }

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
    return mio::success(results);
}

mio::IOResult<void> MultiRunSimulator::save_multi_run_results(const MultiRunResults& results,
                                                              const std::string& base_dir)
{
    BOOST_OUTCOME_TRY(create_result_folders(base_dir));

    // Save percentile results
    std::string summary_file = base_dir + "/summary.txt";
    std::ofstream summary(summary_file);
    if (summary.is_open()) {
        summary << "Multi-Run Simulation Summary\n";
        summary << "Event Type: " << EventSimulator::event_type_to_string(results.event_type) << "\n";
        summary << "Initialization: " << EventSimulator::simulation_type_to_string(results.simulation_type) << "\n";
        summary << "Infection Parameter K: " << results.infection_parameter_k << "\n";
        summary << "Total Runs: " << results.all_runs.size() << "\n";
        summary.close();
    }

    //Save percentile results
    //TODO: Add percentile calculation like in paper Simulation

    return mio::success();
}

void assign_infection_state(mio::abm::World& world, const std::map<uint32_t, bool>& initial_infections,
                            mio::abm::TimePoint t)
{
    for (const auto& [person_id, infected] : initial_infections) {
        if (infected) {
            auto rng = mio::abm::Person::RandomNumberGenerator(world.get_rng(), world.get_person(person_id));
            world.get_person(person_id).add_new_infection(
                mio::abm::Infection(rng, mio::abm::VirusVariant::Alpha, world.get_person(person_id).get_age(),
                                    world.parameters, t, mio::abm::InfectionState::Exposed));
        }
    }
}

mio::IOResult<SimulationResults>
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
    SimulationResults results;
    auto t0   = mio::abm::TimePoint(0); // Start time per simulation
    auto tmax = mio::abm::TimePoint(0) + mio::abm::days(simulation_days); // End time per simulation

    auto world = base_world; // Copy the base world
    // Apply initial infections
    assign_infection_state(world, initial_infections, t0);
    // Set infection parameter K
    world.parameters.get<mio::abm::InfectionRateFromViralShed>() = k_parameter;

    auto sim = mio::abm::Simulation(t0, std::move(world));

    mio::History<mio::abm::TimeSeriesWriter, LogInfectionPerLocationTypePerAgeGroup> historyInfectionPerLocationType{
        Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups())};

    sim.advance(tmax, historyInfectionPerLocationType);

    results.infection_per_loc_type =
        std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionPerLocationType.get_log())};

    // Placeholder implementation
    return mio::success(results);
}