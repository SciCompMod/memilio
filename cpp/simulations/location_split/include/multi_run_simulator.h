#ifndef LOCATION_SPLIT_MULTI_RUN_SIMULATOR_H
#define LOCATION_SPLIT_MULTI_RUN_SIMULATOR_H

#include "city_builder.h"
#include "custom_loggers.h"
#include "defaults.h"

#include "abm/model.h"
#include "memilio/io/io.h"
#include "memilio/utils/time_series.h"

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

/**
 * @brief Configuration of a multi run simulation.
 */
struct MultiRunConfig {
    CityConfig city_config{};
    int num_runs                = Config::DEFAULT_RUNS;
    int simulation_days         = Config::DEFAULT_DAYS;
    int num_initial_infections  = Config::DEFAULT_INITIAL_INFECTIONS;
    double infection_parameter_k = Config::DEFAULT_INFECTION_K;
    std::string output_base_dir  = Config::DEFAULT_OUTPUT_DIR;
    uint32_t custom_seed         = 0; ///< 0 means "use the built in seed".
    bool write_infection_events  = true; ///< Write one csv of infection events per run.
    bool write_contacts          = false; ///< Write the pairwise contact hours. Quadratic in the location size.
    bool write_person_locations  = false; ///< Write where every Person was at every hour. Large.
};

/// @brief One new Infection, with the Location the Person was at when it happened.
struct InfectionEvent {
    int time_step; ///< Hours since the start of the simulation.
    mio::abm::PersonId person_id;
    mio::abm::LocationId location_id;
    mio::abm::LocationType location_type;
    mio::AgeGroup age;
};

/**
 * @brief Everything one simulation run produces.
 */
struct SimulationResults {
    std::vector<mio::TimeSeries<ScalarType>> infection_state_per_age_group;
    std::vector<mio::TimeSeries<ScalarType>> infection_per_location_type_per_age_group;
    std::vector<mio::TimeSeries<ScalarType>> total_infections;

    std::vector<std::tuple<mio::abm::PersonId, mio::abm::LocationId>> household_id;
    std::vector<std::tuple<mio::abm::PersonId, mio::abm::LocationId>> work_id;
    std::vector<std::tuple<mio::abm::LocationId, mio::abm::LocationType>> location_id_and_type;

    std::vector<InfectionEvent> infection_events;
    std::vector<mio::abm::PersonId> initial_infections;
    /// (person 1, person 2, location, hours spent together). Only filled if requested.
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> contact_hours;
    /// Per time step, the Location of every Person. Only filled if requested. Needed to reconstruct
    /// who infected whom, because the infector must have shared the Location with the infectee.
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> person_locations;

    double average_household_size_of_initial_infections           = 0.0;
    double average_persons_above_age_group_0_in_initial_households = 0.0;
};

/**
 * @brief Results of all runs of a multi run simulation.
 */
struct MultiRunResults {
    std::vector<SimulationResults> all_runs;
    double infection_parameter_k = -1.0;
    int successful_runs          = 0;
};

/**
 * @brief Runs the city ABM several times and writes the aggregated results.
 */
class MultiRunSimulator
{
public:
    /// @brief Run all runs configured in @p config.
    static mio::IOResult<MultiRunResults> run_multi_simulation(const MultiRunConfig& config,
                                                               mio::RandomNumberGenerator rng);

    /// @brief Write all results of a multi run simulation below @p base_dir.
    static mio::IOResult<void> save_multi_run_results(const MultiRunResults& results, const std::string& base_dir,
                                                      const MultiRunConfig& config);

private:
    static mio::IOResult<SimulationResults> run_single_simulation(mio::abm::Model&& model,
                                                                  const std::vector<mio::abm::PersonId>& initial_infections,
                                                                  const MultiRunConfig& config);
};

/// @brief Human readable name of a LocationType, used in the csv output.
std::string location_type_to_string(mio::abm::LocationType type);

#endif // LOCATION_SPLIT_MULTI_RUN_SIMULATOR_H
