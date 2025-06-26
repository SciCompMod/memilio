#include "../include/simulation_runner.h"
#include "../include/file_utils.h"
#include "boost/filesystem.hpp"
#include "abm/common_abm_loggers.h"
#include <fstream>

namespace fs = boost::filesystem;

mio::IOResult<SimulationRunner::SimulationResults> SimulationRunner::run_simulation(const SimulationConfig& config)
{
    SimulationResults results;

    // Set up simulation timeframe
    auto t0   = mio::abm::TimePoint(0);
    auto tmax = t0 + mio::abm::days(config.simulation_days);

    // Initialize result containers
    results.infection_per_loc_type.reserve(1);
    results.infection_state_per_age_group.reserve(1);
    results.ensemble_params.reserve(1);

    // Initialize simulation
    auto sim = mio::abm::Simulation(t0, std::move(world));

    // Create history loggers
    mio::History<mio::abm::TimeSeriesWriter, LogInfectionPerLocationTypePerAgeGroup> historyInfectionPerLocationType{
        Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups())};

    mio::History<mio::abm::TimeSeriesWriter, LogInfectionStatePerAgeGroup> historyInfectionStatePerAgeGroup{
        Eigen::Index((size_t)mio::abm::InfectionState::Count * sim.get_world().parameters.get_num_groups())};

    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
        Eigen::Index(mio::abm::InfectionState::Count)};

    // Run simulation
    sim.advance(tmax, historyInfectionPerLocationType, historyInfectionStatePerAgeGroup, historyTimeSeries);

    // Collect results
    auto temp_infection_per_loc_type =
        std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionPerLocationType.get_log())};
    auto temp_infection_state_per_age_group =
        std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionStatePerAgeGroup.get_log())};

    results.infection_per_loc_type.emplace_back(temp_infection_per_loc_type);
    results.infection_state_per_age_group.emplace_back(temp_infection_state_per_age_group);
    results.ensemble_params.emplace_back(std::vector<mio::abm::World>{sim.get_world()});
    results.time_series = std::get<0>(historyTimeSeries.get_log());

    return mio::success(results);
}

mio::IOResult<void> SimulationRunner::save_simulation_results(const SimulationResults& results,
                                                              const std::string& result_dir)
{
    BOOST_OUTCOME_TRY(save_results(results.infection_state_per_age_group, results.ensemble_params, {0},
                                   (fs::path)result_dir / "infection_state_per_age_group" / "0", true));

    BOOST_OUTCOME_TRY(save_results(results.infection_per_loc_type, results.ensemble_params, {0},
                                   (fs::path)result_dir / "infection_per_location_type_per_age_group" / "0", true));

    return mio::success();
}
