#include "../include/multi_run_simulator.h"
#include "../include/file_utils.h"
#include "abm/time.h"
#include "memilio/io/io.h"
#include "abm/analyze_result.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <tuple>

std::string location_type_to_string(mio::abm::LocationType type)
{
    switch (type) {
    case mio::abm::LocationType::Home:
        return "Home";
    case mio::abm::LocationType::Work:
        return "Work";
    case mio::abm::LocationType::School:
        return "School";
    case mio::abm::LocationType::SocialEvent:
        return "SocialEvent";
    case mio::abm::LocationType::BasicsShop:
        return "BasicsShop";
    case mio::abm::LocationType::Hospital:
        return "Hospital";
    case mio::abm::LocationType::ICU:
        return "ICU";
    default:
        return "Unknown";
    }
}

mio::IOResult<MultiRunResults> MultiRunSimulator::run_multi_simulation(const MultiRunConfig& config,
                                                                       mio::RandomNumberGenerator rng)
{
    MultiRunResults results;
    for (int run = 0; run < config.num_runs; ++run) {
        if (run % 10 == 0 || run == config.num_runs - 1) {
            std::cout << "Run " << run + 1 << "/" << config.num_runs << std::endl;
        }

        results.event_type      = config.event_config.type;
        results.simulation_type = config.simulation_type;

        // Step 1: Build city
        // std::cout << "Building city..." << std::endl;
        auto run_rng_counter =
            mio::rng_totalsequence_counter<uint64_t>(static_cast<uint32_t>(0), mio::Counter<uint32_t>(0));
        rng.set_counter(run_rng_counter);
        BOOST_OUTCOME_TRY(auto base_world, CityBuilder::build_world(config.city_config, rng));
        for (auto& person : base_world.get_persons()) {
            // Reset the infection state for each person
            person.get_rng_counter() =
                mio::Counter<uint32_t>(run * 1000000 + run); // Ensure unique counter for each run
        }
        run_rng_counter =
            mio::rng_totalsequence_counter<uint64_t>(static_cast<uint32_t>(run), mio::Counter<uint32_t>(0));
        rng.set_counter(run_rng_counter);
        // We could reset it to zero if we dont want random assignment of infections.

        // CityBuilder::print_city_summary(config.city_config);

        // Step 2: Get map from specific event to ids of persons in simulation
        // std::cout << "Mapping events to person IDs..." << std::endl;
        BOOST_OUTCOME_TRY(auto event_map, EventSimulator::map_events_to_persons(base_world, config.event_config.type));

        // Step 3: Calculate K parameter
        if (results.infection_parameter_k < 0.0) {
            std::cout << "Calculating infection parameter K..." << std::endl;
            // BOOST_OUTCOME_TRY(results.infection_parameter_k, EventSimulator::calculate_infection_parameter_k(
            //                                                      config.event_config, base_world, event_map));
            results.infection_parameter_k = 20; // Placeholder value for K parameter
        }

        // Step 3: Get initial infections
        std::vector<uint32_t> initial_infections;
        if (config.simulation_type == SimType::Panvadere) {
            // std::cout << "Loading infections from Panvadere..." << std::endl;
            BOOST_OUTCOME_TRY(auto infections,
                              EventSimulator::initialize_from_panvadere(config.event_config.type, event_map));
            initial_infections = infections;
        }
        else if (config.simulation_type == SimType::Memilio) {
            // std::cout << "Simulating event infections..." << std::endl;
            BOOST_OUTCOME_TRY(auto infections,
                              EventSimulator::initialize_from_event_simulation(config.event_config.type, event_map));
            initial_infections = infections;
        }

        // Step 4: Run multiple simulations
        // std::cout << "Running " << config.num_runs << " simulations..." << std::endl;
        results.all_runs.reserve(config.num_runs);

        auto single_result = run_single_simulation_with_infections(
            base_world, initial_infections, results.infection_parameter_k, config.simulation_days);

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

    //Save LocationType and ID just take the first run and write out toa txt file
    std::string location_type_and_id_file = base_dir + "/location_type_and_id.txt";
    std::ofstream location_file(location_type_and_id_file);
    if (location_file.is_open()) {
        location_file << "Person ID, Location Type\n";
        for (const auto& timestep : results.all_runs[0].infection_per_location_type_and_id) {
            for (const auto& person_info : timestep) {
                location_file << std::get<0>(person_info) << ", " << location_type_to_string(std::get<1>(person_info))
                              << "\n";
            }
            location_file << "\n"; // Separate time steps with a newline
        }
        location_file.close();
    }

    //Save percentile results
    auto ensembl_inf_per_loc_type  = std::vector<std::vector<mio::TimeSeries<ScalarType>>>{};
    auto ensembl_inf_per_age_group = std::vector<std::vector<mio::TimeSeries<ScalarType>>>{};

    auto ensemble_params = std::vector<std::vector<mio::abm::World>>{};

    for (const auto& run : results.all_runs) {
        ensembl_inf_per_loc_type.push_back(run.infection_per_loc_type);
        ensembl_inf_per_age_group.push_back(run.infection_state_per_age_group);
        ensemble_params.push_back(run.ensemble_params);
    }

    BOOST_OUTCOME_TRY(save_results(ensembl_inf_per_loc_type, ensemble_params, {0},
                                   base_dir + "/infection_per_location_type_per_age_group", false, true));
    BOOST_OUTCOME_TRY(save_results(ensembl_inf_per_age_group, ensemble_params, {0},
                                   base_dir + "/infection_state_per_age_group", false, true));

    return mio::success();
}

void assign_infection_state(mio::abm::World& world, const std::vector<uint32_t>& initial_infections,
                            mio::abm::TimePoint t)
{
    for (const auto& person_id : initial_infections) {
        auto rng = mio::abm::Person::RandomNumberGenerator(world.get_rng(), world.get_person(person_id));
        world.get_person(person_id).add_new_infection(
            mio::abm::Infection(rng, mio::abm::VirusVariant::Alpha, world.get_person(person_id).get_age(),
                                world.parameters, t, mio::abm::InfectionState::InfectedNoSymptoms));
    }
}

mio::IOResult<SimulationResults>
MultiRunSimulator::run_single_simulation_with_infections(mio::abm::World& base_world,
                                                         const std::vector<uint32_t>& initial_infections,
                                                         double k_parameter, int simulation_days)
{
    // TODO: Implement single simulation run with initial infections
    // 1. Apply initial infections to specified people/households
    // 2. Set infection parameter K
    // 3. Run simulation for specified days
    // 4. Return results
    SimulationResults results;
    auto t0   = mio::abm::TimePoint(0); // Start time per simulation
    auto tmax = mio::abm::TimePoint(0) + mio::abm::days(simulation_days); // End time per simulation

    // Apply initial infections
    assign_infection_state(base_world, initial_infections, t0);
    // Set infection parameter K
    base_world.parameters.get<mio::abm::InfectionRateFromViralShed>() = k_parameter;

    auto sim = mio::abm::Simulation(t0, std::move(base_world));

    mio::History<mio::abm::TimeSeriesWriter, LogInfectionPerLocationTypePerAgeGroup> historyInfectionPerLocationType{
        Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups())};
    mio::History<mio::abm::TimeSeriesWriter, LogInfectionStatePerAgeGroup> historyInfectionStatePerAgeGroup{
        Eigen::Index((size_t)mio::abm::InfectionState::Count * sim.get_world().parameters.get_num_groups())};

    mio::History<mio::abm::DataWriterToMemoryDelta, LogLocationTypeAndId> historyLocationTypeAndId;

    sim.advance(tmax, historyInfectionPerLocationType, historyInfectionStatePerAgeGroup, historyLocationTypeAndId);

    results.infection_per_loc_type =
        std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionPerLocationType.get_log())};
    results.infection_state_per_age_group =
        std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionStatePerAgeGroup.get_log())};
    results.ensemble_params                    = std::vector<mio::abm::World>{sim.get_world()};
    results.infection_per_location_type_and_id = std::get<0>(historyLocationTypeAndId.get_log());

    // Placeholder implementation
    return mio::success(results);
}