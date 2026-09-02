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

/**
 * @brief Reconstruct the new Infection%s from the per time step log of Person states.
 *
 * A Person is newly infected in time step k if it was Susceptible in step k-1 and is not any more in
 * step k. Because the ABM lets Person%s interact before it moves them, the Infection happened at the
 * Location the Person occupied in step k-1.
 */
std::vector<InfectionEvent> derive_infection_events(const std::vector<std::vector<PersonState>>& person_states)
{
    std::vector<InfectionEvent> events;
    for (size_t k = 1; k < person_states.size(); ++k) {
        const auto& prev = person_states[k - 1];
        const auto& curr = person_states[k];
        assert(prev.size() == curr.size() && "The population must not change during a run.");

        for (size_t i = 0; i < curr.size(); ++i) {
            if (prev[i].infection_state == mio::abm::InfectionState::Susceptible &&
                curr[i].infection_state != mio::abm::InfectionState::Susceptible) {
                events.push_back(InfectionEvent{static_cast<int>(k), curr[i].person_id, prev[i].location_id,
                                                prev[i].location_type, curr[i].age});
            }
        }
    }
    return events;
}

/// @brief Number of new Infection%s per LocationType and AgeGroup over time, flattened as age * Count + type.
mio::TimeSeries<ScalarType> infections_per_location_type(const std::vector<InfectionEvent>& events, int num_time_points,
                                                         size_t num_groups)
{
    const auto num_types = (size_t)mio::abm::LocationType::Count;
    const auto dim       = Eigen::Index(num_types * num_groups);

    mio::TimeSeries<ScalarType> series(dim);
    for (int k = 0; k < num_time_points; ++k) {
        series.add_time_point(k / 24.0, Eigen::VectorX<ScalarType>::Zero(dim));
    }
    for (const auto& event : events) {
        if (event.time_step < num_time_points) {
            const auto index = num_types * (size_t)event.age.get() + (size_t)event.location_type;
            series.get_value(event.time_step)[Eigen::Index(index)] += 1;
        }
    }
    return series;
}

/// @brief Hours that each pair of Person%s spent at the same Location.
std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>>
calculate_contact_hours(const std::vector<std::vector<PersonState>>& person_states)
{
    std::map<std::tuple<uint32_t, uint32_t, uint32_t>, uint32_t> contact_hours;

    for (const auto& step : person_states) {
        std::unordered_map<uint32_t, std::vector<uint32_t>> location_to_persons;
        for (const auto& state : step) {
            location_to_persons[state.location_id.get()].push_back(static_cast<uint32_t>(state.person_id.get()));
        }

        for (const auto& [location_id, persons] : location_to_persons) {
            for (size_t i = 0; i < persons.size(); ++i) {
                for (size_t j = i + 1; j < persons.size(); ++j) {
                    auto p1 = persons[i];
                    auto p2 = persons[j];
                    if (p1 > p2) {
                        std::swap(p1, p2);
                    }
                    contact_hours[{p1, p2, location_id}]++;
                }
            }
        }
    }

    std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> result;
    result.reserve(contact_hours.size());
    for (const auto& [key, hours] : contact_hours) {
        result.emplace_back(std::get<0>(key), std::get<1>(key), std::get<2>(key), hours);
    }
    return result;
}

} // namespace

mio::IOResult<SimulationResults>
MultiRunSimulator::run_single_simulation(mio::abm::Model&& model,
                                         const std::vector<mio::abm::PersonId>& initial_infections,
                                         const MultiRunConfig& config)
{
    SimulationResults results;

    const auto t0   = mio::abm::TimePoint(0);
    const auto tmax = t0 + mio::abm::days(config.simulation_days);

    assign_infection_state(model, initial_infections, t0);
    std::tie(results.average_household_size_of_initial_infections,
             results.average_persons_above_age_group_0_in_initial_households) =
        initial_household_statistics(model, initial_infections);
    results.initial_infections = initial_infections;

    model.parameters.get<mio::abm::InfectionRateFromViralShed>()[{mio::abm::VirusVariant::Wildtype}] =
        config.infection_parameter_k;

    const auto num_groups = (size_t)model.parameters.get_num_groups();
    auto sim              = mio::abm::Simulation(t0, std::move(model));

    mio::History<mio::abm::TimeSeriesWriter, LogInfectionStatePerAgeGroup> history_infection_state{
        Eigen::Index((size_t)mio::abm::InfectionState::Count * num_groups)};
    mio::History<mio::abm::TimeSeriesWriter, LogTotalInfections> history_total_infections{Eigen::Index(1)};
    mio::History<mio::DataWriterToMemory, LogPersonStates> history_person_states;
    mio::History<mio::DataWriterToMemory, LogHouseholdId, LogWorkId, LogLocationIdAndType> history_topology;

    sim.advance(tmax, history_infection_state, history_total_infections, history_person_states, history_topology);

    results.infection_state_per_age_group =
        std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(history_infection_state.get_log())};
    results.total_infections = std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(history_total_infections.get_log())};

    results.household_id         = std::get<0>(history_topology.get_log()).front();
    results.work_id              = std::get<1>(history_topology.get_log()).front();
    results.location_id_and_type = std::get<2>(history_topology.get_log()).front();

    const auto& person_states = std::get<0>(history_person_states.get_log());
    results.infection_events  = derive_infection_events(person_states);
    results.infection_per_location_type_per_age_group = std::vector<mio::TimeSeries<ScalarType>>{
        infections_per_location_type(results.infection_events, static_cast<int>(person_states.size()), num_groups)};

    if (config.write_contacts) {
        results.contact_hours = calculate_contact_hours(person_states);
    }

    return mio::success(results);
}

mio::IOResult<MultiRunResults> MultiRunSimulator::run_multi_simulation(const MultiRunConfig& config,
                                                                       mio::RandomNumberGenerator rng)
{
    MultiRunResults results;
    results.infection_parameter_k = config.infection_parameter_k;
    results.all_runs.reserve(static_cast<size_t>(config.num_runs));

    for (int run = 0; run < config.num_runs; ++run) {
        if (run % 10 == 0 || run == config.num_runs - 1) {
            std::cout << "Run " << run + 1 << "/" << config.num_runs << std::endl;
        }

        // Every run gets its own stream of the counter based generator so that the runs are
        // independent and reproducible regardless of the order they are executed in.
        rng.set_counter(mio::rng_totalsequence_counter<uint64_t>(static_cast<uint32_t>(run), mio::Counter<uint32_t>(0)));

        auto model = CityBuilder::build_model(config.city_config, rng);

        // The city is built from the Model's generator, so a copy of it taken after the build is a
        // stream that neither the city construction nor the simulation itself uses.
        auto draw_rng                 = model.get_rng();
        const auto initial_infections = draw_initial_infections(model, config.num_initial_infections, draw_rng);

        // Give every Person of this run the same counter offset, as on the panvXabmSim branch. The
        // Person index still separates the individual streams.
        model.reset_rng(mio::Counter<uint32_t>(static_cast<uint32_t>(run)));

        auto single_result = run_single_simulation(std::move(model), initial_infections, config);
        if (single_result) {
            results.all_runs.push_back(std::move(single_result.value()));
            results.successful_runs++;
        }
        else {
            std::cerr << "Run " << run << " failed: " << single_result.error().message() << std::endl;
        }
    }

    if (results.successful_runs == 0) {
        return mio::failure(mio::StatusCode::InvalidValue, "All simulation runs failed.");
    }
    return mio::success(results);
}

namespace
{

/// @brief Write a two column csv of (Person, Location) pairs.
void save_person_location_csv(const std::vector<std::tuple<mio::abm::PersonId, mio::abm::LocationId>>& entries,
                              const std::string& filename, const std::string& location_column)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: could not open " << filename << " for writing.\n";
        return;
    }
    file << "Person_ID," << location_column << "\n";
    for (const auto& [person_id, location_id] : entries) {
        file << person_id.get() << ",";
        if (location_id == mio::abm::LocationId::invalid_id()) {
            file << "\n";
        }
        else {
            file << location_id.get() << "\n";
        }
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
                                  int num_groups, const std::string& dir)
{
    for (size_t run = 0; run < ensemble.size(); ++run) {
        BOOST_OUTCOME_TRY(
            mio::save_result(ensemble[run], {0}, num_groups, dir + "/results_run" + std::to_string(run) + ".h5"));
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
        }
    }

    // Statistics of the households the epidemic starts in.
    {
        std::ofstream file(base_dir + "/initial_infection_households.csv");
        if (file.is_open()) {
            file << "Run,Average_Household_Size,Average_Persons_Above_Age_Group_0\n";
            for (size_t run = 0; run < results.all_runs.size(); ++run) {
                file << run << "," << results.all_runs[run].average_household_size_of_initial_infections << ","
                     << results.all_runs[run].average_persons_above_age_group_0_in_initial_households << "\n";
            }
        }
    }

    // Topology of the city. Identical for all runs up to the random household composition, so the
    // first run is written as a reference.
    save_person_location_csv(first_run.household_id, base_dir + "/household_id.csv", "Household_ID");
    save_person_location_csv(first_run.work_id, base_dir + "/work_id.csv", "Work_ID");
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

    // Per run detail.
    if (config.write_infection_events || config.write_contacts) {
        const std::string runs_dir = base_dir + "/runs";
        BOOST_OUTCOME_TRY(mio::create_directory(runs_dir));
        for (size_t run = 0; run < results.all_runs.size(); ++run) {
            const auto prefix = runs_dir + "/run_" + std::to_string(run);
            if (config.write_infection_events) {
                save_infection_events(results.all_runs[run].infection_events,
                                      results.all_runs[run].initial_infections, prefix + "_infection_events.csv");
            }
            if (config.write_contacts) {
                save_contact_hours(results.all_runs[run].contact_hours, prefix + "_contact_hours.csv");
            }
        }
    }

    // Time series, saved as single runs plus percentiles.
    std::vector<std::vector<mio::TimeSeries<ScalarType>>> ensemble_infection_state;
    std::vector<std::vector<mio::TimeSeries<ScalarType>>> ensemble_infection_per_location_type;
    std::vector<std::vector<mio::TimeSeries<ScalarType>>> ensemble_total_infections;

    for (const auto& run : results.all_runs) {
        ensemble_infection_state.push_back(run.infection_state_per_age_group);
        ensemble_infection_per_location_type.push_back(run.infection_per_location_type_per_age_group);
        ensemble_total_infections.push_back(run.total_infections);
    }

    BOOST_OUTCOME_TRY(save_ensemble(ensemble_infection_state, static_cast<int>(num_age_groups),
                                    base_dir + "/infection_state_per_age_group"));
    BOOST_OUTCOME_TRY(save_ensemble(ensemble_infection_per_location_type, static_cast<int>(num_age_groups),
                                    base_dir + "/infection_per_location_type_per_age_group"));
    BOOST_OUTCOME_TRY(save_ensemble(ensemble_total_infections, 1, base_dir + "/total_infections"));

    return mio::success();
}
