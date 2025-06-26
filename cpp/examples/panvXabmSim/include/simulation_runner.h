#pragma once
#include "abm/abm.h"
#include "memilio/io/result_io.h"
#include "custom_loggers.h"
#include <vector>
#include <string>

using ScalarType = double;

enum class SimType
{
    Panvadere,
    Memilio,
    Both
};

class SimulationRunner
{
public:
    struct SimulationConfig {
        SimType simulation_type = SimType::Both;
        std::string infection_data_file;
        std::string input_dir;
        int n_persons       = 1000;
        int simulation_days = 20;
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

    static mio::IOResult<SimulationResults> run_simulation(const SimulationConfig& config);
    static mio::IOResult<void> save_simulation_results(const SimulationResults& results, const std::string& result_dir);
    static mio::IOResult<void> write_summary_output(const SimulationResults& results, const std::string& filename);
};