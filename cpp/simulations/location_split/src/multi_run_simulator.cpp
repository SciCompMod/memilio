#include "multi_run_simulator.h"
#include "file_utils.h"

#include "abm/infection.h"
#include "abm/model.h"
#include "abm/personal_rng.h"
#include "abm/simulation.h"
#include "memilio/config.h"
#include "memilio/data/analyze_result.h"
#include "memilio/io/result_io.h"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef MEMILIO_ENABLE_OPENMP
// The ABM library parallelises every Model itself. Its cache updates use `omp single`, which binds to
// the innermost enclosing parallel region, i.e. to the loop over runs below if that were parallel:
// only one of the runs would update its caches and the others would crash. The runs are therefore
// executed one after another and each Model uses the threads.
constexpr bool parallel_runs = false;
#else
constexpr bool parallel_runs = true;
#endif

using mio::abm::InfectionState;
using mio::abm::LocationType;

std::string location_type_to_string(LocationType type)
{
    switch (type) {
    case LocationType::Home:
        return "Home";
    case LocationType::School:
        return "School";
    case LocationType::Work:
        return "Work";
    case LocationType::SocialEvent:
        return "SocialEvent";
    case LocationType::BasicsShop:
        return "BasicsShop";
    case LocationType::Hospital:
        return "Hospital";
    case LocationType::ICU:
        return "ICU";
    case LocationType::Car:
        return "Car";
    case LocationType::PublicTransport:
        return "PublicTransport";
    case LocationType::TransportWithoutContact:
        return "TransportWithoutContact";
    case LocationType::Cemetery:
        return "Cemetery";
    default:
        return "Unknown";
    }
}

namespace
{

// ---------------------------------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------------------------------

/// @brief Draw @p count different Person%s uniformly at random.
std::vector<mio::abm::PersonId> draw_initial_infections(const mio::abm::Model& model, int count,
                                                        mio::RandomNumberGenerator& rng)
{
    const auto population = model.get_persons().size();
    const auto n          = std::min(static_cast<size_t>(std::max(0, count)), population);

    std::unordered_set<uint64_t> drawn;
    std::vector<mio::abm::PersonId> drawn_ids;
    while (drawn_ids.size() < n) {
        const auto index = mio::UniformIntDistribution<size_t>::get_instance()(rng, size_t(0), population - 1);
        const auto id    = model.get_persons()[index].get_id();
        if (drawn.insert(id.get()).second) {
            drawn_ids.push_back(id);
        }
    }
    return drawn_ids;
}

// Loggers. Simulation::advance calls log(sim) once at the start and after every hour.

/// @brief Persons per InfectionState and AgeGroup, and persons that are not Susceptible any more.
struct InfectionStateLogger {
    mio::TimeSeries<ScalarType> states{Eigen::Index((size_t)InfectionState::Count * num_age_groups)};
    mio::TimeSeries<ScalarType> total_infections{Eigen::Index(1)};

    void log(const mio::abm::Simulation<>& sim)
    {
        const auto t          = sim.get_time();
        const auto num_states = (size_t)InfectionState::Count;
        auto counts           = states.add_time_point(t.days(), Eigen::VectorX<ScalarType>::Zero(states.get_num_elements()));
        auto infected         = total_infections.add_time_point(t.days(), Eigen::VectorX<ScalarType>::Zero(1));
        for (const auto& person : sim.get_model().get_persons()) {
            const auto state = person.get_infection_state(t);
            counts[Eigen::Index(num_states * (size_t)person.get_age().get() + (size_t)state)] += 1;
            if (state != InfectionState::Susceptible) {
                infected[0] += 1;
            }
        }
    }
};

/**
 * @brief Every new Infection with the Location it happened at.
 *
 * A Person that was Susceptible in the previous hour and is not any more was infected at the Location
 * it occupied in the previous hour, because the ABM lets Person%s interact before it moves them. Only
 * the previous hour is kept, so the memory does not grow with the simulated time.
 */
struct InfectionEventLogger {
    std::vector<InfectionEvent> events;

    void log(const mio::abm::Simulation<>& sim)
    {
        const auto persons = sim.get_model().get_persons();
        previous.resize(persons.size());
        for (size_t i = 0; i < persons.size(); ++i) {
            const auto& person = persons[i];
            const auto state   = person.get_infection_state(sim.get_time());
            if (step > 0 && previous[i].state == InfectionState::Susceptible && state != InfectionState::Susceptible) {
                events.push_back({step, person.get_id(), previous[i].location, previous[i].type, person.get_age()});
            }
            previous[i] = {state, person.get_location(), person.get_location_type()};
        }
        ++step;
    }

private:
    struct PreviousHour {
        InfectionState state;
        mio::abm::LocationId location;
        LocationType type;
    };
    std::vector<PreviousHour> previous;
    int step = 0;
};

/// @brief Streams the Location of every Person at every hour to a csv file.
struct PersonLocationLogger {
    std::ofstream file;

    void log(const mio::abm::Simulation<>& sim)
    {
        for (const auto& person : sim.get_model().get_persons()) {
            file << step << "," << person.get_id().get() << "," << person.get_location().get() << "\n";
        }
        ++step;
    }

private:
    int step = 0;
};

/// @brief New Infection%s per day, LocationType and AgeGroup. An event in step k happened in hour [k - 1, k).
mio::TimeSeries<ScalarType> infections_per_day(const std::vector<InfectionEvent>& events, int num_days)
{
    const auto num_types = (size_t)LocationType::Count;
    mio::TimeSeries<ScalarType> series(Eigen::Index(num_types * num_age_groups));
    for (int day = 0; day < num_days; ++day) {
        series.add_time_point(day, Eigen::VectorX<ScalarType>::Zero(series.get_num_elements()));
    }
    for (const auto& event : events) {
        const int day = (event.time_step - 1) / 24;
        if (day < num_days) {
            series.get_value(day)[Eigen::Index(num_types * (size_t)event.age.get() + (size_t)event.location_type)] +=
                1;
        }
    }
    return series;
}

/// @brief Infect the initial Person%s, then simulate one run in steps of one hour.
mio::IOResult<RunResult> simulate_run(mio::abm::Model&& model, const std::vector<mio::abm::PersonId>& initial,
                                      const MultiRunConfig& config, int run)
{
    RunResult result;
    result.run                = run;
    result.initial_infections = initial;

    const auto t0 = mio::abm::TimePoint(0);
    for (const auto& id : initial) {
        auto& person = model.get_person(id);
        auto rng     = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
        person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                     model.parameters, t0, InfectionState::InfectedNoSymptoms));
    }
    model.parameters.get<mio::abm::InfectionRateFromViralShed>()[{mio::abm::VirusVariant::Wildtype}] =
        config.infection_parameter_k;

    for (const auto& person : model.get_persons()) {
        result.persons.push_back({person.get_id(), person.get_age(), person.get_assigned_location(LocationType::Home),
                                  person.get_assigned_location(LocationType::Work),
                                  person.get_assigned_location(LocationType::School)});
    }
    for (const auto& location : model.get_locations()) {
        result.location_types.push_back(location.get_type());
    }

    InfectionStateLogger state_logger;
    InfectionEventLogger event_logger;
    auto sim        = mio::abm::Simulation(t0, std::move(model)); // steps of one hour
    const auto tmax = t0 + mio::abm::days(config.simulation_days);

    if (config.write_person_locations) {
        PersonLocationLogger location_logger;
        const auto filename = config.output_base_dir + "/runs/run_" + std::to_string(run) + "_person_locations.csv";
        location_logger.file.open(filename);
        if (!location_logger.file.is_open()) {
            return mio::failure(mio::StatusCode::FileNotFound, "Could not open " + filename);
        }
        location_logger.file << "Timestep,Person_ID,Location_ID\n";
        sim.advance(tmax, state_logger, event_logger, location_logger);
    }
    else {
        sim.advance(tmax, state_logger, event_logger);
    }

    result.infection_states             = std::move(state_logger.states);
    result.total_infections             = std::move(state_logger.total_infections);
    result.infection_events             = std::move(event_logger.events);
    result.infections_per_location_type = infections_per_day(result.infection_events, config.simulation_days);
    return mio::success(std::move(result));
}

} // namespace

mio::IOResult<std::vector<RunResult>> run_simulations(const MultiRunConfig& config, mio::RandomNumberGenerator rng)
{
    if (config.write_person_locations) {
        BOOST_OUTCOME_TRY(mio::create_directory(config.output_base_dir + "/runs", /*create_parents=*/true));
    }
#ifdef _OPENMP
    if (config.num_threads > 0) {
        omp_set_num_threads(config.num_threads);
    }
#endif

    // Every run writes only its own slot, so the results do not depend on the number of threads.
    std::vector<mio::IOResult<RunResult>> results(static_cast<size_t>(config.num_runs),
                                                  mio::failure(mio::StatusCode::UnknownError, "Run not executed."));
    int finished = 0;

#pragma omp parallel for schedule(dynamic, 1) if (parallel_runs)
    for (int run = 0; run < config.num_runs; ++run) {
        // Every run gets its own stream of the counter based generator, so it is reproducible
        // regardless of the order in which the runs are executed.
        auto run_rng = rng;
        run_rng.set_counter(
            mio::rng_totalsequence_counter<uint64_t>(static_cast<uint32_t>(run), mio::Counter<uint32_t>(0)));
        auto model = CityBuilder::build_model(config.city_config, run_rng);

        // The city is built from the Model's generator, so a copy of it taken now is a stream that neither
        // the city construction nor the simulation uses.
        auto draw_rng      = model.get_rng();
        const auto initial = draw_initial_infections(model, config.num_initial_infections, draw_rng);

        // The same counter offset for every Person of this run, as on the panvXabmSim branch; the Person
        // index still separates the individual streams.
        model.reset_rng(mio::Counter<uint32_t>(static_cast<uint32_t>(run)));

        results[static_cast<size_t>(run)] = simulate_run(std::move(model), initial, config, run);

#pragma omp critical(location_split_progress)
        {
            ++finished;
            if (finished % 10 == 0 || finished == config.num_runs) {
                std::cout << "Finished " << finished << "/" << config.num_runs << " runs" << std::endl;
            }
        }
    }

    std::vector<RunResult> runs;
    for (auto& result : results) {
        BOOST_OUTCOME_TRY(auto&& run, std::move(result));
        runs.push_back(std::move(run));
    }
    return mio::success(std::move(runs));
}

// ---------------------------------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------------------------------

namespace
{

mio::IOResult<std::ofstream> open_file(const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        return mio::failure(mio::StatusCode::FileNotFound, "Could not open " + filename);
    }
    return mio::success(std::move(file));
}

/// @brief The location id, or an empty field for "none".
std::string id_or_empty(mio::abm::LocationId id)
{
    return id == mio::abm::LocationId::invalid_id() ? std::string() : std::to_string(id.get());
}

mio::IOResult<void> write_summary(const MultiRunConfig& config, size_t num_runs)
{
    BOOST_OUTCOME_TRY(auto&& file, open_file(config.output_base_dir + "/summary.txt"));
    const auto& city = config.city_config;
    file << "Multi-Run Simulation Summary\n"
         << "Region: " << CityParameters::region_to_string(city.region) << "\n"
         << "Population: " << city.total_population << "\n"
         << "Simulation days: " << config.simulation_days << "\n"
         << "Initial infections: " << config.num_initial_infections << "\n"
         << "Infection Parameter K: " << config.infection_parameter_k << "\n"
         << "Successful runs: " << num_runs << "/" << config.num_runs << "\n"
         << "Seed: " << config.custom_seed << "\n"
         << "Contact weights: " << contact_weights_to_string(city.contact_weights) << "\n"
         << "Room sizes (0 = no split): school " << city.school_room_size << ", work " << city.work_room_size
         << ", social event " << city.event_room_size << "\n";
    return mio::success();
}

/// @brief Which arm of the location split study this is, read by the comparison plots in viz/.
mio::IOResult<void> write_scenario(const MultiRunConfig& config)
{
    BOOST_OUTCOME_TRY(auto&& file, open_file(config.output_base_dir + "/scenario.json"));
    const auto& city = config.city_config;
    const std::vector<std::tuple<std::string, std::string, int>> splits{
        {"School", "school", city.school_room_size},
        {"Work", "work", city.work_room_size},
        {"SocialEvent", "event", city.event_room_size}};

    // The label names every split; the first split location type is reported as the axis of the arm.
    std::string label, split_type = "none";
    int room_size = 0;
    for (const auto& [type, word, size] : splits) {
        if (size > 0) {
            label += (label.empty() ? "" : ", ") + word + " rooms of " + std::to_string(size);
            if (split_type == "none") {
                split_type = type;
                room_size  = size;
            }
        }
    }
    auto json_size = [](int size) {
        return size > 0 ? std::to_string(size) : std::string("null");
    };
    file << "{\"label\": \"" << (label.empty() ? std::string("no split") : label)
         << "\", \"split_location_type\": \"" << split_type << "\", \"room_size\": " << json_size(room_size)
         << ", \"school_room_size\": " << json_size(city.school_room_size)
         << ", \"work_room_size\": " << json_size(city.work_room_size)
         << ", \"event_room_size\": " << json_size(city.event_room_size) << "}\n";
    return mio::success();
}

/// @brief Every run builds its own city, but the Location%s are the same in every run.
mio::IOResult<void> write_location_types(const RunResult& run, const std::string& filename)
{
    BOOST_OUTCOME_TRY(auto&& file, open_file(filename));
    file << "Location_ID,Location_Type\n";
    for (size_t id = 0; id < run.location_types.size(); ++id) {
        file << id << "," << location_type_to_string(run.location_types[id]) << "\n";
    }
    return mio::success();
}

mio::IOResult<void> write_persons(const RunResult& run, const std::string& filename)
{
    BOOST_OUTCOME_TRY(auto&& file, open_file(filename));
    file << "Person_ID,Age_Group,Household_ID,Work_ID,School_ID\n";
    for (const auto& p : run.persons) {
        file << p.person_id.get() << "," << p.age.get() << "," << id_or_empty(p.household) << ","
             << id_or_empty(p.work) << "," << id_or_empty(p.school) << "\n";
    }
    return mio::success();
}

mio::IOResult<void> write_infection_events(const RunResult& run, const std::string& filename)
{
    BOOST_OUTCOME_TRY(auto&& file, open_file(filename));
    file << "Timestep,Person_ID,Location_ID,Location_Type,Age_Group\n";
    // The initially infected are infected before the first step and have no infecting Location.
    for (const auto& id : run.initial_infections) {
        file << "0," << id.get() << ",,Initial,\n";
    }
    for (const auto& e : run.infection_events) {
        file << e.time_step << "," << e.person_id.get() << "," << e.location_id.get() << ","
             << location_type_to_string(e.location_type) << "," << (size_t)e.age.get() << "\n";
    }
    return mio::success();
}

/**
 * @brief Write one time series per run (results_run<i>.h5) and the p05 to p95 percentiles over the runs.
 * mio::save_results is not used because it also serializes the Model parameters as a graph, which
 * abm::Model does not support.
 */
mio::IOResult<void> write_time_series(const std::vector<RunResult>& runs,
                                      mio::TimeSeries<ScalarType> RunResult::*series, int num_groups,
                                      const std::string& dir)
{
    std::vector<std::vector<mio::TimeSeries<ScalarType>>> ensemble;
    for (const auto& run : runs) {
        ensemble.push_back({run.*series});
        BOOST_OUTCOME_TRY(mio::save_result(ensemble.back(), {0}, num_groups,
                                           dir + "/results_run" + std::to_string(run.run) + ".h5"));
    }
    for (const auto& [name, p] : {std::pair<const char*, ScalarType>{"p05", 0.05},
                                  {"p25", 0.25},
                                  {"p50", 0.50},
                                  {"p75", 0.75},
                                  {"p95", 0.95}}) {
        BOOST_OUTCOME_TRY(mio::create_directory(dir + "/" + name));
        BOOST_OUTCOME_TRY(
            mio::save_result(mio::ensemble_percentile(ensemble, p), {0}, num_groups, dir + "/" + name + "/Results.h5"));
    }
    return mio::success();
}

} // namespace

mio::IOResult<void> save_results(const std::vector<RunResult>& runs, const MultiRunConfig& config)
{
    const auto& dir = config.output_base_dir;
    BOOST_OUTCOME_TRY(create_result_folders(dir));

    BOOST_OUTCOME_TRY(write_summary(config, runs.size()));
    BOOST_OUTCOME_TRY(write_scenario(config));
    BOOST_OUTCOME_TRY(write_location_types(runs.front(), dir + "/location_id_and_type.csv"));
    CityBuilder::save_city_to_file(config.city_config, dir + "/city_config.csv");

    // Per run files; the person locations were already written while the runs were simulated.
    BOOST_OUTCOME_TRY(mio::create_directory(dir + "/runs"));
    for (const auto& run : runs) {
        const auto prefix = dir + "/runs/run_" + std::to_string(run.run);
        BOOST_OUTCOME_TRY(write_persons(run, prefix + "_persons.csv"));
        BOOST_OUTCOME_TRY(write_infection_events(run, prefix + "_infection_events.csv"));
    }

    const auto groups = static_cast<int>(num_age_groups);
    BOOST_OUTCOME_TRY(
        write_time_series(runs, &RunResult::infection_states, groups, dir + "/infection_state_per_age_group"));
    BOOST_OUTCOME_TRY(write_time_series(runs, &RunResult::infections_per_location_type, groups,
                                        dir + "/infection_per_location_type_per_age_group"));
    BOOST_OUTCOME_TRY(write_time_series(runs, &RunResult::total_infections, 1, dir + "/total_infections"));
    return mio::success();
}
