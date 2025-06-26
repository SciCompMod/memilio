#pragma once
#include "simulation_runner.h"
#include "city_builder.h"
#include "event_simulator.h"
#include <vector>

struct MultiRunConfig {
    CityConfig city_config;
    EventSimulationConfig event_config;
    int num_runs                = 100;
    int simulation_days         = 30;
    std::string output_base_dir = "./results";
};

struct MultiRunResults {
    std::vector<SimulationRunner::SimulationResults> all_runs;
    EventType event_type;
    SimType simulation_type;
    double infection_parameter_k;
};

class MultiRunSimulator
{
public:
    static mio::IOResult<MultiRunResults> run_multi_simulation(const MultiRunConfig& config);
    static mio::IOResult<void> save_multi_run_results(const MultiRunResults& results, const std::string& base_dir);

private:
    static mio::IOResult<SimulationRunner::SimulationResults>
    run_single_simulation_with_infections(const mio::abm::World& base_world,
                                          const std::map<uint32_t, bool>& initial_infections, double k_parameter,
                                          int simulation_days);
};