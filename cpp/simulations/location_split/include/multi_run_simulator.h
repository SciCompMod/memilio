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
    int num_threads              = 0; ///< Runs executed in parallel with OpenMP. 0 means the OpenMP default.
};

/// @brief One new Infection, with the Location the Person was at when it happened.
struct InfectionEvent {
    int time_step; ///< Hours since the start of the simulation.
    mio::abm::PersonId person_id;
    mio::abm::LocationId location_id;
    mio::abm::LocationType location_type;
    mio::AgeGroup age;
};

/// @brief The Location%s a Person is assigned to in one run. Invalid ids mean "none".
struct PersonAssignment {
    mio::abm::PersonId person_id;
    mio::AgeGroup age;
    mio::abm::LocationId household;
    mio::abm::LocationId work;
    mio::abm::LocationId school; ///< The school, or the room of the school if schools are split.
};

/**
 * @brief Everything one simulation run produces.
 *
 * Every run builds its own city, so the assignments below belong to this run only.
 */
struct SimulationResults {
    int run = -1; ///< Index of the run. Runs with the same index and seed share the city across arms.

    std::vector<mio::TimeSeries<ScalarType>> infection_state_per_age_group;
    /// New Infection%s per day, per LocationType and AgeGroup.
    std::vector<mio::TimeSeries<ScalarType>> infection_per_location_type_per_age_group;
    std::vector<mio::TimeSeries<ScalarType>> total_infections;

    std::vector<PersonAssignment> persons;
    std::vector<std::tuple<mio::abm::LocationId, mio::abm::LocationType>> location_id_and_type;

    std::vector<InfectionEvent> infection_events;
    std::vector<mio::abm::PersonId> initial_infections;
    /// (person 1, person 2, location, hours spent together). Only filled if requested.
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> contact_hours;

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
    /**
     * @brief Simulate one run hour by hour.
     * Infection events, contact hours and person locations are derived on the fly from the state of
     * the previous hour, so the memory needed does not grow with the simulated time (apart from the
     * contact hours, which grow with the number of distinct contacts).
     */
    static mio::IOResult<SimulationResults> run_single_simulation(mio::abm::Model&& model,
                                                                  const std::vector<mio::abm::PersonId>& initial_infections,
                                                                  const MultiRunConfig& config, int run);
};

/// @brief Human readable name of a LocationType, used in the csv output.
std::string location_type_to_string(mio::abm::LocationType type);

#endif // LOCATION_SPLIT_MULTI_RUN_SIMULATOR_H
