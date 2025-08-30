#pragma once
#include "custom_loggers.h"
#include "abm/common_abm_loggers.h"
#include "city_builder.h"
#include "event_simulator.h"
#include <vector>

struct MultiRunConfig {
    CityConfig city_config;
    EventSimulationConfig event_config;
    SimType simulation_type;
    int num_runs;
    int simulation_days;
    double infection_parameter_k = -1.0; // Default value indicating not set
    std::string output_base_dir;
    uint32_t custom_seed = 0; // Default value indicating use predefined seeds
};

struct SimulationResults {
    std::vector<mio::TimeSeries<ScalarType>> infection_per_loc_type;
    std::vector<mio::TimeSeries<ScalarType>> infection_state_per_age_group;
    std::vector<std::vector<std::tuple<uint32_t, mio::abm::LocationType>>> infection_per_location_type_and_id;
    std::vector<mio::abm::World> ensemble_params;
    std::vector<mio::abm::World> ensemble_params_no_agegroups;
    // std::vector<std::vector<std::tuple<uint32_t, mio::abm::LocationType>>> infection_per_location_type_and_id;

    std::vector<mio::TimeSeries<ScalarType>> history_infected_amount; // Figure 1
    std::vector<std::tuple<uint32_t, uint32_t>> history_household_id; // Figure 2
    std::vector<std::tuple<uint32_t, uint32_t>> history_work_id; // Figure 2
    std::vector<std::tuple<uint32_t, mio::abm::LocationType>> loc_id_and_type;
    std::vector<std::vector<std::tuple<uint32_t, uint32_t>>> history_person_and_location_id; // Figure 2 / 3
    std::vector<std::vector<std::tuple<uint32_t, uint32_t, mio::abm::LocationType>>>
        history_detailed_infection; // Figure 2/3
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
    static mio::IOResult<MultiRunResults> run_multi_simulation(const MultiRunConfig& config,
                                                               mio::RandomNumberGenerator rng);
    static mio::IOResult<void> save_multi_run_results(const MultiRunResults& results, const std::string& base_dir,
                                                      const MultiRunConfig& config);

private:
    static mio::IOResult<SimulationResults>
    run_single_simulation_with_infections(mio::abm::World& base_world, const std::vector<uint32_t>& initial_infections,
                                          double k_parameter, int simulation_days);
};