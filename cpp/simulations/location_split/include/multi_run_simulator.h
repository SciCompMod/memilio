#ifndef LOCATION_SPLIT_MULTI_RUN_SIMULATOR_H
#define LOCATION_SPLIT_MULTI_RUN_SIMULATOR_H

#include "city_builder.h"
#include "constants.h"
#include "defaults.h"

#include "abm/infection_state.h"
#include "abm/location_id.h"
#include "abm/location_type.h"
#include "abm/person_id.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/io/io.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/time_series.h"

#include <cstdint>
#include <string>
#include <vector>

/// @brief Configuration of a multi run simulation.
struct MultiRunConfig {
    CityConfig city_config{};
    int num_runs                 = Config::DEFAULT_RUNS;
    int simulation_days          = Config::DEFAULT_DAYS;
    int num_initial_infections   = Config::DEFAULT_INITIAL_INFECTIONS;
    double infection_parameter_k = Config::DEFAULT_INFECTION_K;
    std::string output_base_dir  = Config::DEFAULT_OUTPUT_DIR;
    uint32_t custom_seed         = 0; ///< 0 means "use the built in seed".
    bool write_person_locations  = false; ///< Write where every Person was at every hour. Large.
    bool redraw_school_rooms     = false; ///< Assign the pupils of a school to its rooms anew every day.
    bool write_school_rooms      = false; ///< Write persons, susceptible, infectious and new infections per school room and hour.
    int num_threads              = 0; ///< Threads used; 0 means the OpenMP default.
};

/// @brief One new Infection, with the Location the Person was at when it happened.
struct InfectionEvent {
    int time_step; ///< Hours since the start; the Infection happened in the hour before.
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
    mio::abm::LocationId school; ///< The school, or its room if schools are split.
};

/// @brief School room hours with the same day and the same persons present, aggregated over all school rooms.
struct SchoolRoomHours {
    int day;
    int persons; ///< Persons present in the room.
    int susceptible; ///< Susceptible Person%s present.
    int infectious; ///< Person%s present in InfectedNoSymptoms, InfectedSymptoms, InfectedSevere or InfectedCritical.
    int64_t room_hours; ///< Number of (room, hour) pairs with these counts.
    int64_t new_infections; ///< Of the susceptible present, how many were infected in that hour.
};

/// @brief Everything one run produces. Every run builds its own city.
struct RunResult {
    int run = -1; ///< Runs with the same index and seed have the same city in every arm.
    std::vector<PersonAssignment> persons;
    std::vector<mio::abm::LocationType> location_types; ///< Type of every Location, indexed by LocationId.
    std::vector<mio::abm::PersonId> initial_infections;
    std::vector<InfectionEvent> infection_events;
    std::vector<SchoolRoomHours> school_rooms; ///< Only filled with MultiRunConfig::write_school_rooms.

    /// Persons per InfectionState and AgeGroup every hour, flattened as age * InfectionState::Count + state.
    mio::TimeSeries<ScalarType> infection_states{Eigen::Index((size_t)mio::abm::InfectionState::Count * num_age_groups)};
    /// Persons that are not Susceptible any more, i.e. the cumulative number of infections, every hour.
    mio::TimeSeries<ScalarType> total_infections{Eigen::Index(1)};
    /// New Infection%s per day, LocationType and AgeGroup, flattened as age * LocationType::Count + type.
    mio::TimeSeries<ScalarType> infections_per_location_type{
        Eigen::Index((size_t)mio::abm::LocationType::Count * num_age_groups)};
};

/// @brief Simulate all runs of @p config. Runs are independent and executed in parallel if possible.
mio::IOResult<std::vector<RunResult>> run_simulations(const MultiRunConfig& config, mio::RandomNumberGenerator rng);

/// @brief Write the results of all runs to config.output_base_dir.
mio::IOResult<void> save_results(const std::vector<RunResult>& runs, const MultiRunConfig& config);

/// @brief Name of a LocationType as used in the csv output.
std::string location_type_to_string(mio::abm::LocationType type);

#endif // LOCATION_SPLIT_MULTI_RUN_SIMULATOR_H
