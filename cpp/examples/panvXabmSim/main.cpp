#include <fstream>
#include <vector>
#include <iostream>
#include <chrono>
#include <string>
#include <map>
#include "abm/abm.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/uncertain_value.h"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "memilio/utils/miompi.h"
#include "memilio/io/binary_serializer.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"
#include "abm/common_abm_loggers.h"
#include "memilio/utils/random_number_generator.h"

#include "include/constants.h"
#include "include/custom_loggers.h"
#include "include/world_creator.h"
#include "include/file_utils.h"

namespace fs     = boost::filesystem;
using ScalarType = double;

mio::IOResult<void> main_flow()
{
    // Default infection data file path
    std::string infection_data_file =
        "/Users/saschakorf/Nosynch/Arbeit/memilio/memilio/cpp/examples/PanVadere/restaurant/simulation_runs/"
        "lu-2020_airflow_inlet_right_outlet_left/infections.txt";
    std::string input_dir = "/Users/saschakorf/Nosynch/Arbeit/memilio/memilio/cpp/examples/results";

    std::string precomputed_dir = input_dir + "/results";
    std::string result_dir      = input_dir + "/results_" + currentDateTime();
    auto created                = create_result_folders(result_dir);

    std::cout << "Creating restaurant simulation from: " << infection_data_file << std::endl;

    int n_persons = 1000;

    // Create the simulation world
    auto world = create_world_from_file(infection_data_file, n_persons);

    // Set up simulation timeframe
    auto t0   = mio::abm::TimePoint(0);
    auto tmax = t0 + mio::abm::days(20);

    auto ensemble_infection_per_loc_type =
        std::vector<std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of infection per location type results
    ensemble_infection_per_loc_type.reserve(size_t(1));

    auto ensemble_infection_state_per_age_group =
        std::vector<std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of infection state per age group results
    ensemble_infection_state_per_age_group.reserve(size_t(1));

    auto ensemble_params = std::vector<std::vector<mio::abm::World>>{}; // Vector of all worlds
    ensemble_params.reserve(size_t(1));

    // Initialize simulation
    auto sim = mio::abm::Simulation(t0, std::move(world));

    mio::History<mio::abm::TimeSeriesWriter, LogInfectionPerLocationTypePerAgeGroup> historyInfectionPerLocationType{
        Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups())};
    mio::History<mio::abm::TimeSeriesWriter, LogInfectionStatePerAgeGroup> historyInfectionStatePerAgeGroup{
        Eigen::Index((size_t)mio::abm::InfectionState::Count * sim.get_world().parameters.get_num_groups())};
    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
        Eigen::Index(mio::abm::InfectionState::Count)};

    // Run simulation and collect data
    sim.advance(tmax, historyInfectionPerLocationType, historyInfectionStatePerAgeGroup, historyTimeSeries);

    auto temp_sim_infection_per_loc_tpye =
        std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionPerLocationType.get_log())};
    auto temp_sim_infection_state_per_age_group =
        std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionStatePerAgeGroup.get_log())};
    // Push result of the simulation back to the result vector
    ensemble_infection_per_loc_type.emplace_back(temp_sim_infection_per_loc_tpye);
    ensemble_infection_state_per_age_group.emplace_back(temp_sim_infection_state_per_age_group);
    ensemble_params.emplace_back(std::vector<mio::abm::World>{sim.get_world()});

    BOOST_OUTCOME_TRY(save_results(ensemble_infection_state_per_age_group, ensemble_params, {0},
                                   (fs::path)result_dir / "infection_state_per_age_group" / "0", true));
    BOOST_OUTCOME_TRY(save_results(ensemble_infection_per_loc_type, ensemble_params, {0},
                                   (fs::path)result_dir / "infection_per_location_type_per_age_group" / "0", true));

    // copy results into a fixed name folder to have easier access
    std::string last_run_dir = input_dir + "/results_last_run";
    auto copied              = copy_result_folder(result_dir, last_run_dir);

    // Write results to file
    std::ofstream outfile("panvXabm_results.txt");
    std::get<0>(historyTimeSeries.get_log())
        .print_table({"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4, outfile);

    std::cout << "Results written to panvXabm_results.txt" << std::endl;

    return mio::success();
}

int main()
{
    auto result = main_flow();
    if (result.has_error()) {
        std::cerr << "Error: " << result.error().message() << std::endl;
        return 1;
    }
    return 0;
}