#pragma once
#include "city_builder.h"
#include "event_simulator.h"
#include <vector>

struct MultiRunConfig {
    CityConfig city_config;
    EventSimulationConfig event_config;
    SimType simulation_type     = SimType::Both;
    int num_runs                = 100;
    int simulation_days         = 30;
    std::string output_base_dir = "./results";
};

struct SimulationResults {
    std::vector<std::vector<mio::TimeSeries<ScalarType>>> infection_per_loc_type;
    std::vector<std::vector<mio::TimeSeries<ScalarType>>> infection_state_per_age_group;
    std::vector<std::vector<mio::abm::World>> ensemble_params;
    mio::TimeSeries<ScalarType> time_series;

    // Add a constructor to initialize time_series
    SimulationResults()
        : time_series(Eigen::Index(mio::abm::InfectionState::Count))
    {
    }
};

struct MultiRunResults {
    std::vector<SimulationResults> all_runs;
    EventType event_type;
    SimType simulation_type;
    double infection_parameter_k = -1.0; // Default value indicating not set
    int successful_runs          = 0;
};

class MultiRunSimulator
{
public:
    static mio::IOResult<MultiRunResults> run_multi_simulation(const MultiRunConfig& config);
    static mio::IOResult<void> save_multi_run_results(const MultiRunResults& results, const std::string& base_dir);

private:
    static mio::IOResult<SimulationResults>
    run_single_simulation_with_infections(const mio::abm::World& base_world,
                                          const std::map<uint32_t, bool>& initial_infections, double k_parameter,
                                          int simulation_days);
};