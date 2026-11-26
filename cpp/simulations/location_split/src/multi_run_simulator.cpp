#include "multi_run_simulator.h"
#include "constants.h"
#include "file_utils.h"

#include "abm/common_abm_loggers.h"
#include "abm/infection.h"
#include "abm/personal_rng.h"
#include "abm/simulation.h"
#include "abm/time.h"
#include "memilio/io/history.h"
#include "memilio/data/analyze_result.h"
#include "memilio/io/result_io.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "memilio/config.h"

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

std::string location_type_to_string(mio::abm::LocationType type)
{
    switch (type) {
    case mio::abm::LocationType::Home:
        return "Home";
    case mio::abm::LocationType::School:
        return "School";
    case mio::abm::LocationType::Work:
        return "Work";
    case mio::abm::LocationType::SocialEvent:
        return "SocialEvent";
    case mio::abm::LocationType::BasicsShop:
        return "BasicsShop";
    case mio::abm::LocationType::Hospital:
        return "Hospital";
    case mio::abm::LocationType::ICU:
        return "ICU";
    case mio::abm::LocationType::Car:
        return "Car";
    case mio::abm::LocationType::PublicTransport:
        return "PublicTransport";
    case mio::abm::LocationType::TransportWithoutContact:
        return "TransportWithoutContact";
    case mio::abm::LocationType::Cemetery:
        return "Cemetery";
    default:
        return "Unknown";
    }
}

namespace
{

/**
 * @brief Draw the Person%s that start the epidemic.
 * @param[in] model The city.
 * @param[in] count Number of Person%s to infect.
 * @param[in,out] rng Random number generator.
 * @return The drawn PersonIds, without duplicates.
 */
std::vector<mio::abm::PersonId> draw_initial_infections(const mio::abm::Model& model, int count,
                                                        mio::RandomNumberGenerator& rng)
{
    const auto population = model.get_persons().size();
    const auto n          = std::min(static_cast<size_t>(std::max(0, count)), population);

    std::unordered_set<uint64_t> drawn;
    std::vector<mio::abm::PersonId> initial_infections;
    initial_infections.reserve(n);

    while (initial_infections.size() < n) {
        const auto index =
            mio::UniformIntDistribution<size_t>::get_instance()(rng, size_t(0), population - 1);
        const auto person_id = model.get_persons()[index].get_id();
        if (drawn.insert(person_id.get()).second) {
            initial_infections.push_back(person_id);
        }
    }
    return initial_infections;
}

/// @brief Infect the given Person%s at time @p t.
void assign_infection_state(mio::abm::Model& model, const std::vector<mio::abm::PersonId>& initial_infections,
                            mio::abm::TimePoint t)
{
    for (const auto& person_id : initial_infections) {
        auto& person = model.get_person(person_id);
        auto rng     = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
        person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                     model.parameters, t,
                                                     mio::abm::InfectionState::InfectedNoSymptoms));
    }
}

/**
 * @brief Statistics of the households of the initially infected Person%s.
 * @return (average household size, average number of members above age group 0).
 */
std::pair<double, double> initial_household_statistics(const mio::abm::Model& model,
                                                       const std::vector<mio::abm::PersonId>& initial_infections)
{
    if (initial_infections.empty()) {
        return {0.0, 0.0};
    }

    // Members and non-toddler members per Home.
    std::unordered_map<uint32_t, std::pair<int, int>> members_per_home;
    for (const auto& person : model.get_persons()) {
        const auto home = person.get_assigned_location(mio::abm::LocationType::Home);
        if (home == mio::abm::LocationId::invalid_id()) {
            continue;
        }
        auto& entry = members_per_home[home.get()];
        entry.first++;
        if (person.get_age() != age_group_0_to_4) {
            entry.second++;
        }
    }

    double total_size          = 0.0;
    double total_above_group_0 = 0.0;
    for (const auto& person_id : initial_infections) {
        const auto home = model.get_person(person_id).get_assigned_location(mio::abm::LocationType::Home);
        const auto it   = members_per_home.find(home.get());
        if (it != members_per_home.end()) {
            total_size += it->second.first;
            total_above_group_0 += it->second.second;
        }
    }

    const auto n = static_cast<double>(initial_infections.size());
    return {total_size / n, total_above_group_0 / n};
}

/// @brief New Infection%s per day, per LocationType and AgeGroup, flattened as age * Count + type.
mio::TimeSeries<ScalarType> infections_per_location_type(const std::vector<InfectionEvent>& events, int num_days,
                                                         size_t num_groups)
{
    const auto num_types = (size_t)mio::abm::LocationType::Count;
    const auto dim       = Eigen::Index(num_types * num_groups);

    mio::TimeSeries<ScalarType> series(dim);
    for (int day = 0; day < num_days; ++day) {
        series.add_time_point(day, Eigen::VectorX<ScalarType>::Zero(dim));
    }
    for (const auto& event : events) {
        // An event in time step k happened in the hour [k - 1, k).
        const int day = (event.time_step - 1) / 24;
        if (day >= 0 && day < num_days) {
            const auto index = num_types * (size_t)event.age.get() + (size_t)event.location_type;
            series.get_value(day)[Eigen::Index(index)] += 1;
        }
    }
    return series;
}

/// @brief Add one hour to every pair of Person%s that are at the same Location.
void add_contact_hours(const mio::abm::Model& model,
                       std::map<std::tuple<uint32_t, uint32_t, uint32_t>, uint32_t>& contact_hours)
{
    std::unordered_map<uint32_t, std::vector<uint32_t>> location_to_persons;
    for (const auto& person : model.get_persons()) {
        location_to_persons[person.get_location().get()].push_back(static_cast<uint32_t>(person.get_id().get()));
    }
    for (const auto& [location_id, persons] : location_to_persons) {
        for (size_t i = 0; i < persons.size(); ++i) {
            for (size_t j = i + 1; j < persons.size(); ++j) {
                const auto p1 = std::min(persons[i], persons[j]);
                const auto p2 = std::max(persons[i], persons[j]);
                contact_hours[{p1, p2, location_id}]++;
            }
        }
    }
}

/// @brief Append the Location of every Person at time step @p step to the open csv.
void write_person_locations(const mio::abm::Model& model, int step, std::ofstream& file)
{
    for (const auto& person : model.get_persons()) {
        file << step << "," << person.get_id().get() << "," << person.get_location().get() << "\n";
    }
}

/// @brief The Location%s each Person is assigned to.
std::vector<PersonAssignment> person_assignments(const mio::abm::Model& model)
{
    std::vector<PersonAssignment> assignments;
    assignments.reserve(model.get_persons().size());
    for (const auto& person : model.get_persons()) {
        assignments.push_back(PersonAssignment{person.get_id(), person.get_age(),
                                               person.get_assigned_location(mio::abm::LocationType::Home),
                                               person.get_assigned_location(mio::abm::LocationType::Work),
                                               person.get_assigned_location(mio::abm::LocationType::School)});
    }
    return assignments;
}

/// @brief Directory of the per run files.
std::string runs_directory(const MultiRunConfig& config)
{
    return config.output_base_dir + "/runs";
}

} // namespace

mio::IOResult<SimulationResults>
MultiRunSimulator::run_single_simulation(mio::abm::Model&& model,
                                         const std::vector<mio::abm::PersonId>& initial_infections,
                                         const MultiRunConfig& config, int run)
{
    SimulationResults results;
    results.run = run;

    const auto t0   = mio::abm::TimePoint(0);
    const auto tmax = t0 + mio::abm::days(config.simulation_days);
    const auto dt   = mio::abm::hours(1);

    assign_infection_state(model, initial_infections, t0);
    std::tie(results.average_household_size_of_initial_infections,
             results.average_persons_above_age_group_0_in_initial_households) =
        initial_household_statistics(model, initial_infections);
    results.initial_infections = initial_infections;

    model.parameters.get<mio::abm::InfectionRateFromViralShed>()[{mio::abm::VirusVariant::Wildtype}] =
        config.infection_parameter_k;

    const auto num_groups = (size_t)model.parameters.get_num_groups();
    auto sim              = mio::abm::Simulation(t0, std::move(model)); // steps of one hour
    auto persons          = sim.get_model().get_persons();

    results.persons = person_assignments(sim.get_model());
    results.location_id_and_type.reserve(sim.get_model().get_locations().size());
    for (const auto& location : sim.get_model().get_locations()) {
        results.location_id_and_type.emplace_back(location.get_id(), location.get_type());
    }

    mio::History<mio::abm::TimeSeriesWriter, LogInfectionStatePerAgeGroup> history_infection_state{
        Eigen::Index((size_t)mio::abm::InfectionState::Count * num_groups)};
    mio::History<mio::abm::TimeSeriesWriter, LogTotalInfections> history_total_infections{Eigen::Index(1)};

    std::map<std::tuple<uint32_t, uint32_t, uint32_t>, uint32_t> contact_hours;
    std::ofstream location_file;
    if (config.write_person_locations) {
        const auto filename = runs_directory(config) + "/run_" + std::to_string(run) + "_person_locations.csv";
        location_file.open(filename);
        if (!location_file.is_open()) {
            return mio::failure(mio::StatusCode::FileNotFound, "Could not open " + filename);
        }
        location_file << "Timestep,Person_ID,Location_ID\n";
    }

    // State of the previous hour. A Person that was Susceptible then and is not any more now was
    // infected in this hour, at the Location it occupied in the previous hour, because the ABM lets
    // Person%s interact before it moves them.
    std::vector<mio::abm::InfectionState> previous_state(persons.size());
    std::vector<mio::abm::LocationId> previous_location(persons.size());
    std::vector<mio::abm::LocationType> previous_location_type(persons.size());
    auto remember_current_hour = [&]() {
        for (size_t i = 0; i < persons.size(); ++i) {
            previous_state[i]         = persons[i].get_infection_state(sim.get_time());
            previous_location[i]      = persons[i].get_location();
            previous_location_type[i] = persons[i].get_location_type();
        }
    };

    int step = 0;
    history_infection_state.log(sim);
    history_total_infections.log(sim);
    if (config.write_contacts) {
        add_contact_hours(sim.get_model(), contact_hours);
    }
    if (config.write_person_locations) {
        write_person_locations(sim.get_model(), step, location_file);
    }
    remember_current_hour();

    while (sim.get_time() < tmax) {
        sim.advance(sim.get_time() + dt);
        ++step;

        history_infection_state.log(sim);
        history_total_infections.log(sim);

        for (size_t i = 0; i < persons.size(); ++i) {
            if (previous_state[i] == mio::abm::InfectionState::Susceptible &&
                persons[i].get_infection_state(sim.get_time()) != mio::abm::InfectionState::Susceptible) {
                results.infection_events.push_back(InfectionEvent{step, persons[i].get_id(), previous_location[i],
                                                                  previous_location_type[i], persons[i].get_age()});
            }
        }
        if (config.write_contacts) {
            add_contact_hours(sim.get_model(), contact_hours);
        }
        if (config.write_person_locations) {
            write_person_locations(sim.get_model(), step, location_file);
        }
        remember_current_hour();
    }

    results.infection_state_per_age_group =
        std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(history_infection_state.get_log())};
    results.total_infections = std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(history_total_infections.get_log())};
    results.infection_per_location_type_per_age_group = std::vector<mio::TimeSeries<ScalarType>>{
        infections_per_location_type(results.infection_events, config.simulation_days, num_groups)};

    if (config.write_contacts) {
        results.contact_hours.reserve(contact_hours.size());
        for (const auto& [key, hours] : contact_hours) {
            results.contact_hours.emplace_back(std::get<0>(key), std::get<1>(key), std::get<2>(key), hours);
        }
    }

    return mio::success(results);
}

mio::IOResult<MultiRunResults> MultiRunSimulator::run_multi_simulation(const MultiRunConfig& config,
                                                                       mio::RandomNumberGenerator rng)
{
    MultiRunResults results;
    results.infection_parameter_k = config.infection_parameter_k;

    if (config.write_person_locations) {
        BOOST_OUTCOME_TRY(mio::create_directory(runs_directory(config), /*create_parents=*/true));
    }

#ifdef _OPENMP
    if (config.num_threads > 0) {
        omp_set_num_threads(config.num_threads);
    }
#endif

    // Every run writes only its own slot, so the result does not depend on the number of threads or
    // on the order in which the runs finish.
    std::vector<mio::IOResult<SimulationResults>> run_results(
        static_cast<size_t>(config.num_runs),
        mio::failure(mio::StatusCode::UnknownError, "Run was not executed."));
    int finished_runs = 0;

#pragma omp parallel for schedule(dynamic, 1) if (parallel_runs)
    for (int run = 0; run < config.num_runs; ++run) {
        // Every run gets its own stream of the counter based generator so that the runs are
        // independent and reproducible regardless of the order they are executed in.
        auto run_rng = rng;
        run_rng.set_counter(
            mio::rng_totalsequence_counter<uint64_t>(static_cast<uint32_t>(run), mio::Counter<uint32_t>(0)));

        auto model = CityBuilder::build_model(config.city_config, run_rng);

        // The city is built from the Model's generator, so a copy of it taken after the build is a
        // stream that neither the city construction nor the simulation itself uses.
        auto draw_rng                 = model.get_rng();
        const auto initial_infections = draw_initial_infections(model, config.num_initial_infections, draw_rng);

        // Give every Person of this run the same counter offset, as on the panvXabmSim branch. The
        // Person index still separates the individual streams.
        model.reset_rng(mio::Counter<uint32_t>(static_cast<uint32_t>(run)));

        run_results[static_cast<size_t>(run)] = run_single_simulation(std::move(model), initial_infections, config, run);

#pragma omp critical(location_split_progress)
        {
            ++finished_runs;
            if (finished_runs % 10 == 0 || finished_runs == config.num_runs) {
                std::cout << "Finished " << finished_runs << "/" << config.num_runs << " runs" << std::endl;
            }
        }
    }

    results.all_runs.reserve(run_results.size());
    for (size_t run = 0; run < run_results.size(); ++run) {
        if (run_results[run]) {
            results.all_runs.push_back(std::move(run_results[run].value()));
            results.successful_runs++;
        }
        else {
            std::cerr << "Run " << run << " failed: " << run_results[run].error().message() << std::endl;
        }
    }

    if (results.successful_runs == 0) {
        return mio::failure(mio::StatusCode::InvalidValue, "All simulation runs failed.");
    }
    return mio::success(results);
}

namespace
{

/// @brief Write the Location%s every Person of one run is assigned to. An empty field means "none".
void save_person_assignments(const std::vector<PersonAssignment>& persons, const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: could not open " << filename << " for writing.\n";
        return;
    }
    auto write_id = [&file](mio::abm::LocationId id) {
        if (id != mio::abm::LocationId::invalid_id()) {
            file << id.get();
        }
    };
    file << "Person_ID,Age_Group,Household_ID,Work_ID,School_ID\n";
    for (const auto& person : persons) {
        file << person.person_id.get() << "," << person.age.get() << ",";
        write_id(person.household);
        file << ",";
        write_id(person.work);
        file << ",";
        write_id(person.school);
        file << "\n";
    }
}

/// @brief Write the infection events of one run.
void save_infection_events(const std::vector<InfectionEvent>& events, const std::vector<mio::abm::PersonId>& initial,
                           const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: could not open " << filename << " for writing.\n";
        return;
    }
    file << "Timestep,Person_ID,Location_ID,Location_Type,Age_Group\n";
    // The initially infected are seeded before the first time step and have no infecting Location.
    for (const auto& person_id : initial) {
        file << "0," << person_id.get() << ",,Initial,\n";
    }
    for (const auto& event : events) {
        file << event.time_step << "," << event.person_id.get() << "," << event.location_id.get() << ","
             << location_type_to_string(event.location_type) << "," << (size_t)event.age.get() << "\n";
    }
}

/// @brief Write the pairwise contact hours of one run.
void save_contact_hours(const std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>>& contacts,
                        const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: could not open " << filename << " for writing.\n";
        return;
    }
    file << "Person_1,Person_2,Location_ID,Contact_Hours\n";
    for (const auto& [p1, p2, location_id, hours] : contacts) {
        file << p1 << "," << p2 << "," << location_id << "," << hours << "\n";
    }
}

/**
 * @brief Write one h5 per run plus the p05, p25, p50, p75 and p95 percentiles over the runs.
 *
 * mio::save_results is not used here because it also serializes the Model parameters as a graph,
 * which is not supported for abm::Model.
 */
mio::IOResult<void> save_ensemble(const std::vector<std::vector<mio::TimeSeries<ScalarType>>>& ensemble,
                                  const std::vector<int>& run_indices, int num_groups, const std::string& dir)
{
    for (size_t i = 0; i < ensemble.size(); ++i) {
        BOOST_OUTCOME_TRY(mio::save_result(ensemble[i], {0}, num_groups,
                                           dir + "/results_run" + std::to_string(run_indices[i]) + ".h5"));
    }

    for (const auto& [name, p] : {std::pair<const char*, ScalarType>{"p05", 0.05},
                                  {"p25", 0.25},
                                  {"p50", 0.50},
                                  {"p75", 0.75},
                                  {"p95", 0.95}}) {
        const auto percentile_dir = dir + "/" + name;
        BOOST_OUTCOME_TRY(mio::create_directory(percentile_dir));
        BOOST_OUTCOME_TRY(mio::save_result(mio::ensemble_percentile(ensemble, p), {0}, num_groups,
                                           percentile_dir + "/Results.h5"));
    }
    return mio::success();
}

} // namespace

mio::IOResult<void> MultiRunSimulator::save_multi_run_results(const MultiRunResults& results,
                                                              const std::string& base_dir,
                                                              const MultiRunConfig& config)
{
    BOOST_OUTCOME_TRY(create_result_folders(base_dir));

    const auto& first_run = results.all_runs.front();

    // Run configuration.
    {
        std::ofstream summary(base_dir + "/summary.txt");
        if (summary.is_open()) {
            summary << "Multi-Run Simulation Summary\n";
            summary << "Region: " << CityParameters::region_to_string(config.city_config.region) << "\n";
            summary << "Population: " << config.city_config.total_population << "\n";
            summary << "Simulation days: " << config.simulation_days << "\n";
            summary << "Initial infections: " << config.num_initial_infections << "\n";
            summary << "Infection Parameter K: " << results.infection_parameter_k << "\n";
            summary << "Successful runs: " << results.successful_runs << "/" << config.num_runs << "\n";
            summary << "Seed: " << config.custom_seed << "\n";
            summary << "Contact weights: " << contact_weights_to_string(config.city_config.contact_weights) << "\n";
            summary << "Room sizes (0 = no split): school " << config.city_config.school_room_size << ", work "
                    << config.city_config.work_room_size << ", social event " << config.city_config.event_room_size
                    << "\n";
        }
    }

    // Which arm of the location split study this is, read by the comparison plots in viz/.
    {
        std::ofstream scenario(base_dir + "/scenario.json");
        if (scenario.is_open()) {
            const auto& city = config.city_config;
            // One location type is split per arm; with several, the first one is reported as the axis.
            const std::vector<std::tuple<std::string, std::string, int>> splits{
                {"School", "school", city.school_room_size},
                {"Work", "work", city.work_room_size},
                {"SocialEvent", "event", city.event_room_size}};
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
            scenario << "{\"label\": \"" << (label.empty() ? std::string("no split") : label)
                     << "\", \"split_location_type\": \"" << split_type << "\", \"room_size\": " << json_size(room_size)
                     << ", \"school_room_size\": " << json_size(city.school_room_size)
                     << ", \"work_room_size\": " << json_size(city.work_room_size)
                     << ", \"event_room_size\": " << json_size(city.event_room_size) << "}\n";
        }
    }

    // Statistics of the households the epidemic starts in.
    {
        std::ofstream file(base_dir + "/initial_infection_households.csv");
        if (file.is_open()) {
            file << "Run,Average_Household_Size,Average_Persons_Above_Age_Group_0\n";
            for (const auto& run_result : results.all_runs) {
                file << run_result.run << "," << run_result.average_household_size_of_initial_infections << ","
                     << run_result.average_persons_above_age_group_0_in_initial_households << "\n";
            }
        }
    }

    // Locations and their types. Every run builds its own city, but the Location%s are created in the
    // same order and number in every run, so this list holds for all runs. Which Person is assigned to
    // which Location differs between runs and is written per run below.
    {
        std::ofstream file(base_dir + "/location_id_and_type.csv");
        if (file.is_open()) {
            file << "Location_ID,Location_Type\n";
            for (const auto& [location_id, type] : first_run.location_id_and_type) {
                file << location_id.get() << "," << location_type_to_string(type) << "\n";
            }
        }
    }

    CityBuilder::save_city_to_file(config.city_config, base_dir + "/city_config.csv");

    // Per run detail. The person locations were already written during the runs.
    {
        const std::string runs_dir = base_dir + "/runs";
        BOOST_OUTCOME_TRY(mio::create_directory(runs_dir));
        for (const auto& run_result : results.all_runs) {
            const auto prefix = runs_dir + "/run_" + std::to_string(run_result.run);
            save_person_assignments(run_result.persons, prefix + "_persons.csv");
            if (config.write_infection_events) {
                save_infection_events(run_result.infection_events, run_result.initial_infections,
                                      prefix + "_infection_events.csv");
            }
            if (config.write_contacts) {
                save_contact_hours(run_result.contact_hours, prefix + "_contact_hours.csv");
            }
        }
    }

    // Time series, saved as single runs plus percentiles.
    std::vector<std::vector<mio::TimeSeries<ScalarType>>> ensemble_infection_state;
    std::vector<std::vector<mio::TimeSeries<ScalarType>>> ensemble_infection_per_location_type;
    std::vector<std::vector<mio::TimeSeries<ScalarType>>> ensemble_total_infections;

    std::vector<int> run_indices;
    for (const auto& run : results.all_runs) {
        run_indices.push_back(run.run);
        ensemble_infection_state.push_back(run.infection_state_per_age_group);
        ensemble_infection_per_location_type.push_back(run.infection_per_location_type_per_age_group);
        ensemble_total_infections.push_back(run.total_infections);
    }

    BOOST_OUTCOME_TRY(save_ensemble(ensemble_infection_state, run_indices, static_cast<int>(num_age_groups),
                                    base_dir + "/infection_state_per_age_group"));
    BOOST_OUTCOME_TRY(save_ensemble(ensemble_infection_per_location_type, run_indices, static_cast<int>(num_age_groups),
                                    base_dir + "/infection_per_location_type_per_age_group"));
    BOOST_OUTCOME_TRY(save_ensemble(ensemble_total_infections, run_indices, 1, base_dir + "/total_infections"));

    return mio::success();
}
