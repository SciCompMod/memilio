/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Henrik Zunker
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#include "abm/household.h"
#include "abm/infection_state.h"
#include "abm/location_type.h"
#include "abm/lockdown_rules.h"
#include "abm/model.h"
#include "abm/person_id.h"
#include "abm/time.h"
#include "graph_abm/graph_abm_mobility.h"
#include "graph_abm/graph_abmodel.h"
#include "memilio/io/history.h"
#include "memilio/mobility/graph.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/random_number_generator.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Constants
constexpr size_t NUM_AGE_GROUPS       = 6;
constexpr size_t NUM_INFECTION_STATES = static_cast<size_t>(mio::abm::InfectionState::Count); // 8

// Age group indices (0-4, 5-14, 15-34, 35-59, 60-79, 80+)
const auto AGE_0_4   = mio::AgeGroup(0);
const auto AGE_5_14  = mio::AgeGroup(1);
const auto AGE_15_34 = mio::AgeGroup(2);
const auto AGE_35_59 = mio::AgeGroup(3);
const auto AGE_60_79 = mio::AgeGroup(4);
const auto AGE_80    = mio::AgeGroup(5);

// Approximate German age distribution (fractions)
const double AGE_DISTRIBUTION[NUM_AGE_GROUPS] = {0.046, 0.094, 0.213, 0.283, 0.223, 0.141};

// Custom Logger: Per-node age-stratified infection state counts
/**
 * @brief Logger that records per-age-group, per-infection-state counts.
 *
 * Output vector layout (length = NUM_AGE_GROUPS * NUM_INFECTION_STATES = 48):
 *   [age0_S, age0_E, age0_INS, age0_ISy, age0_ISev, age0_ICrit, age0_R, age0_D,
 *    age1_S, age1_E, ..., age5_D]
 *
 * Compatible with TimeSeriesWriter for direct TimeSeries output.
 */
struct LogAgeStratifiedState : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorX<ScalarType>>;

    static Type log(const mio::abm::Simulation<mio::GraphABModel>& sim)
    {
        Eigen::VectorX<ScalarType> counts = Eigen::VectorX<ScalarType>::Zero(NUM_AGE_GROUPS * NUM_INFECTION_STATES);

        auto t = sim.get_time();
        for (const auto& person : sim.get_model().get_persons()) {
            size_t age   = person.get_age().get();
            size_t state = static_cast<size_t>(person.get_infection_state(t));
            if (age < NUM_AGE_GROUPS && state < NUM_INFECTION_STATES) {
                counts[static_cast<Eigen::Index>(age * NUM_INFECTION_STATES + state)] += 1.0;
            }
        }
        return {t, counts};
    }
};

// Configuration
struct Config {
    size_t num_nodes               = 5;
    size_t persons_per_node        = 2000;
    size_t num_days                = 30;
    size_t run_id                  = 0;
    uint64_t seed                  = 42;
    size_t max_dampings            = 3;
    double fixed_commuter_fraction = 0.1; // used only when param_sampling is disabled
    bool param_sampling            = true; // sample calibration parameters per run
    std::string output_dir         = "abm_gnn_data";
};

// Calibration Parameters (sampled per run, used as GNN input)
/**
 * @brief The calibration parameter vector  this is what the GNN learns to map
 *        to the per-node time series. At inference time, gradient descent on
 *        these parameters (with the trained surrogate) recovers the values
 *        that best explain observed real-world data.
 *
 * Keep the order in sync with PARAM_NAMES below and with the Python loader.
 */
struct SimParameters {
    double transmission_rate_scale; // multiplier on AerosolTransmissionRates [0.5, 2.0]
    double contact_rate_work; // MaximumContacts at workplaces          [5, 20]
    double contact_rate_school; // MaximumContacts at schools             [10, 30]
    double contact_rate_social; // MaximumContacts at social events       [5, 25]
    double contact_rate_shop; // MaximumContacts at shops               [5, 25]
    double initial_infected_fraction; // symptomatic seed prevalence            [0.001, 0.05]
    double initial_recovered_fraction; // pre-existing immunity fraction         [0.0, 0.3]
    double commuter_fraction; // fraction of workers commuting          [0.05, 0.25]
};

constexpr size_t NUM_CALIB_PARAMS               = 8;
const char* const PARAM_NAMES[NUM_CALIB_PARAMS] = {
    "transmission_rate_scale", "contact_rate_work",         "contact_rate_school",        "contact_rate_social",
    "contact_rate_shop",       "initial_infected_fraction", "initial_recovered_fraction", "commuter_fraction"};

/**
 * @brief Samples a calibration parameter vector uniformly from the configured ranges.
 *
 * Uses the run-local RNG so that each run_id reproducibly maps to one parameter
 * vector.
 */
SimParameters sample_parameters(std::mt19937& rng)
{
    auto u = [&](double lo, double hi) {
        return std::uniform_real_distribution<double>(lo, hi)(rng);
    };
    SimParameters p;
    p.transmission_rate_scale    = u(0.5, 2.0);
    p.contact_rate_work          = u(5.0, 20.0);
    p.contact_rate_school        = u(10.0, 30.0);
    p.contact_rate_social        = u(5.0, 25.0);
    p.contact_rate_shop          = u(5.0, 25.0);
    p.initial_infected_fraction  = u(0.001, 0.05);
    p.initial_recovered_fraction = u(0.0, 0.3);
    p.commuter_fraction          = u(0.05, 0.25);
    return p;
}

/**
 * @brief Returns a deterministic default parameter set (used when sampling is disabled).
 */
SimParameters default_parameters(double commuter_fraction)
{
    SimParameters p;
    p.transmission_rate_scale    = 1.0;
    p.contact_rate_work          = 10.0;
    p.contact_rate_school        = 20.0;
    p.contact_rate_social        = 10.0;
    p.contact_rate_shop          = 20.0;
    p.initial_infected_fraction  = 0.01;
    p.initial_recovered_fraction = 0.05;
    p.commuter_fraction          = commuter_fraction;
    return p;
}

Config parse_args(int argc, char* argv[])
{
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--num_nodes" && i + 1 < argc)
            cfg.num_nodes = std::stoul(argv[++i]);
        else if (arg == "--persons_per_node" && i + 1 < argc)
            cfg.persons_per_node = std::stoul(argv[++i]);
        else if (arg == "--num_days" && i + 1 < argc)
            cfg.num_days = std::stoul(argv[++i]);
        else if (arg == "--run_id" && i + 1 < argc)
            cfg.run_id = std::stoul(argv[++i]);
        else if (arg == "--seed" && i + 1 < argc)
            cfg.seed = std::stoull(argv[++i]);
        else if (arg == "--max_dampings" && i + 1 < argc)
            cfg.max_dampings = std::stoul(argv[++i]);
        else if (arg == "--fixed_commuter_fraction" && i + 1 < argc)
            cfg.fixed_commuter_fraction = std::stod(argv[++i]);
        else if (arg == "--no_param_sampling")
            cfg.param_sampling = false;
        else if (arg == "--output_dir" && i + 1 < argc)
            cfg.output_dir = argv[++i];
        else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: graph_abm_data_generation [options]\n"
                      << "  --num_nodes N                Number of graph nodes/regions (default: 5)\n"
                      << "  --persons_per_node P         Approximate population per node (default: 2000)\n"
                      << "  --num_days D                 Simulation days (default: 30)\n"
                      << "  --run_id R                   Run identifier for output naming (default: 0)\n"
                      << "  --seed S                     Random seed (default: 42)\n"
                      << "  --max_dampings M             Max NPI interventions (0-3, default: 3)\n"
                      << "  --fixed_commuter_fraction F  Commuter fraction when sampling is off (default: 0.1)\n"
                      << "  --no_param_sampling          Disable per-run parameter sampling\n"
                      << "  --output_dir DIR             Output directory (default: abm_gnn_data)\n";
            std::exit(0);
        }
    }
    return cfg;
}

// Model Setup Helpers

/**
 * @brief Sets COVID-19 infection parameters (matching the ODE paper parameters).
 *
 * The transmission_rate_scale from the calibration parameter vector is applied
 * here as a multiplier on AerosolTransmissionRates.
 */
void set_infection_parameters(mio::GraphABModel& model, const SimParameters& params)
{
    // Time distributions for disease progression
    // Based on Kuehn et al. (2021), same as the ODE surrogate paper
    model.parameters.get<mio::abm::TimeExposedToNoSymptoms>() =
        mio::ParameterDistributionLogNormal(std::log(3.335), 0.3);
    model.parameters.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>() =
        mio::ParameterDistributionLogNormal(std::log(2.59), 0.3);
    model.parameters.get<mio::abm::TimeInfectedNoSymptomsToRecovered>() =
        mio::ParameterDistributionLogNormal(std::log(5.0), 0.3);
    model.parameters.get<mio::abm::TimeInfectedSymptomsToRecovered>() =
        mio::ParameterDistributionLogNormal(std::log(6.95), 0.3);
    model.parameters.get<mio::abm::TimeInfectedSymptomsToSevere>() =
        mio::ParameterDistributionLogNormal(std::log(5.0), 0.3);
    model.parameters.get<mio::abm::TimeInfectedSevereToRecovered>() =
        mio::ParameterDistributionLogNormal(std::log(7.28), 0.3);
    model.parameters.get<mio::abm::TimeInfectedSevereToCritical>() =
        mio::ParameterDistributionLogNormal(std::log(7.28), 0.3);
    model.parameters.get<mio::abm::TimeInfectedCriticalToRecovered>() =
        mio::ParameterDistributionLogNormal(std::log(13.07), 0.3);
    model.parameters.get<mio::abm::TimeInfectedCriticalToDead>() =
        mio::ParameterDistributionLogNormal(std::log(13.07), 0.3);

    // Transmission parameters (scaled by calibration parameter)
    model.parameters.get<mio::abm::AerosolTransmissionRates>() = params.transmission_rate_scale;

    // Age groups that go to school / work
    model.parameters.get<mio::abm::AgeGroupGotoSchool>()           = false;
    model.parameters.get<mio::abm::AgeGroupGotoSchool>()[AGE_5_14] = true;

    model.parameters.get<mio::abm::AgeGroupGotoWork>()            = false;
    model.parameters.get<mio::abm::AgeGroupGotoWork>()[AGE_15_34] = true;
    model.parameters.get<mio::abm::AgeGroupGotoWork>()[AGE_35_59] = true;

    model.parameters.check_constraints();
}

/**
 * @brief Creates household members for the 6 age groups.
 */
struct HouseholdMembers {
    mio::abm::HouseholdMember child;
    mio::abm::HouseholdMember young_adult;
    mio::abm::HouseholdMember adult;
    mio::abm::HouseholdMember senior;

    HouseholdMembers()
        : child(NUM_AGE_GROUPS)
        , young_adult(NUM_AGE_GROUPS)
        , adult(NUM_AGE_GROUPS)
        , senior(NUM_AGE_GROUPS)
    {
        // Children: 0-4 or 5-14
        child.set_age_weight(AGE_0_4, 1);
        child.set_age_weight(AGE_5_14, 2); // more school-age children

        // Young adults: 15-34
        young_adult.set_age_weight(AGE_15_34, 1);

        // Adults: 35-59
        adult.set_age_weight(AGE_35_59, 1);

        // Seniors: 60-79 or 80+
        senior.set_age_weight(AGE_60_79, 3);
        senior.set_age_weight(AGE_80, 1);
    }
};

/**
 * @brief Adds households and persons to a model to reach target population.
 *
 * Household composition (approximate German census distribution):
 *   - Single senior:        15%
 *   - Single adult:         25%
 *   - Two adults:           20%
 *   - Family (2+1 child):   15%
 *   - Family (2+2 children):10%
 *   - Single parent + child: 5%
 *   - Senior couple:        10%
 */
void create_population(mio::GraphABModel& model, size_t target_population, const HouseholdMembers& members)
{
    // Average household size ~ 2.0 persons
    size_t total_households = std::max<size_t>(target_population / 2, 1);

    // Distribution of household types (fractions)
    size_t n_single_senior = static_cast<size_t>(total_households * 0.15);
    size_t n_single_adult  = static_cast<size_t>(total_households * 0.25);
    size_t n_two_adults    = static_cast<size_t>(total_households * 0.20);
    size_t n_family_1child = static_cast<size_t>(total_households * 0.15);
    size_t n_family_2child = static_cast<size_t>(total_households * 0.10);
    size_t n_single_parent = static_cast<size_t>(total_households * 0.05);
    size_t n_senior_couple = total_households - n_single_senior - n_single_adult - n_two_adults - n_family_1child -
                             n_family_2child - n_single_parent;

    auto add_group = [&](auto household, size_t count) {
        if (count == 0)
            return;
        auto group = mio::abm::HouseholdGroup();
        group.add_households(household, static_cast<int>(count));
        add_household_group_to_model(model, group);
    };

    // Single senior
    auto single_senior = mio::abm::Household();
    single_senior.add_members(members.senior, 1);
    add_group(single_senior, n_single_senior);

    // Single adult
    auto single_adult = mio::abm::Household();
    single_adult.add_members(members.young_adult, 1);
    add_group(single_adult, n_single_adult);

    // Two adults
    auto two_adults = mio::abm::Household();
    two_adults.add_members(members.young_adult, 1);
    two_adults.add_members(members.adult, 1);
    add_group(two_adults, n_two_adults);

    // Family with 1 child
    auto family_1 = mio::abm::Household();
    family_1.add_members(members.child, 1);
    family_1.add_members(members.young_adult, 1);
    family_1.add_members(members.adult, 1);
    add_group(family_1, n_family_1child);

    // Family with 2 children
    auto family_2 = mio::abm::Household();
    family_2.add_members(members.child, 2);
    family_2.add_members(members.young_adult, 1);
    family_2.add_members(members.adult, 1);
    add_group(family_2, n_family_2child);

    // Single parent + child
    auto single_parent = mio::abm::Household();
    single_parent.add_members(members.child, 1);
    single_parent.add_members(members.adult, 1);
    add_group(single_parent, n_single_parent);

    // Senior couple
    auto senior_couple = mio::abm::Household();
    senior_couple.add_members(members.senior, 2);
    add_group(senior_couple, n_senior_couple);
}

/**
 * @brief Creates shared locations (School, Work, Shop, etc.) for a model node.
 *
 * Per-location MaximumContacts values come from the sampled calibration vector,
 * so each run explores a different point in contact-rate space.
 *
 * @returns Map from LocationType to LocationId for later reference.
 */
std::map<mio::abm::LocationType, mio::abm::LocationId>
create_locations(mio::GraphABModel& model, size_t /*population_size*/, const SimParameters& params)
{
    std::map<mio::abm::LocationType, mio::abm::LocationId> locations;

    // School (only children go)
    auto school = model.add_location(mio::abm::LocationType::School);
    model.get_location(school).get_infection_parameters().set<mio::abm::MaximumContacts>(params.contact_rate_school);
    locations[mio::abm::LocationType::School] = school;

    // Work (adults go)
    auto work = model.add_location(mio::abm::LocationType::Work);
    model.get_location(work).get_infection_parameters().set<mio::abm::MaximumContacts>(params.contact_rate_work);
    locations[mio::abm::LocationType::Work] = work;

    // Social event
    auto event = model.add_location(mio::abm::LocationType::SocialEvent);
    model.get_location(event).get_infection_parameters().set<mio::abm::MaximumContacts>(params.contact_rate_social);
    locations[mio::abm::LocationType::SocialEvent] = event;

    // Shop
    auto shop = model.add_location(mio::abm::LocationType::BasicsShop);
    model.get_location(shop).get_infection_parameters().set<mio::abm::MaximumContacts>(params.contact_rate_shop);
    locations[mio::abm::LocationType::BasicsShop] = shop;

    // Hospital
    auto hospital = model.add_location(mio::abm::LocationType::Hospital);
    model.get_location(hospital).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    locations[mio::abm::LocationType::Hospital] = hospital;

    // ICU
    auto icu = model.add_location(mio::abm::LocationType::ICU);
    model.get_location(icu).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    locations[mio::abm::LocationType::ICU] = icu;

    return locations;
}

/**
 * @brief Builds the initial infection-state distribution from the sampled
 *        calibration parameters.
 *
 * The two free parameters from the calibration vector are:
 *   - initial_infected_fraction  (symptomatic seed prevalence)
 *   - initial_recovered_fraction (pre-existing immunity)
 * The remaining compartments are scaled relative to the symptomatic fraction
 * (mirroring the ratios used in the ODE paper's persistent-threat regime).
 */
std::vector<double> build_infection_distribution(const SimParameters& params)
{
    double p_infected  = params.initial_infected_fraction;
    double p_exposed   = p_infected * 0.5;
    double p_no_symp   = p_infected * 0.3;
    double p_severe    = p_infected * 0.02;
    double p_critical  = p_infected * 0.005;
    double p_recovered = params.initial_recovered_fraction;
    double p_dead      = 0.0; // start with no deaths

    double p_susceptible = 1.0 - (p_exposed + p_no_symp + p_infected + p_severe + p_critical + p_recovered + p_dead);
    if (p_susceptible < 0.3) {
        // rescale to ensure enough susceptibles
        double scale = 0.5 / (1.0 - p_susceptible);
        p_exposed *= scale;
        p_no_symp *= scale;
        p_infected *= scale;
        p_severe *= scale;
        p_critical *= scale;
        p_recovered *= scale;
        p_susceptible = 1.0 - (p_exposed + p_no_symp + p_infected + p_severe + p_critical + p_recovered + p_dead);
    }

    return {p_susceptible, p_exposed, p_no_symp, p_infected, p_severe, p_critical, p_recovered, p_dead};
}

/**
 * @brief Assigns infection states and non-home locations to persons of a model.
 */
void initialize_persons(mio::GraphABModel& model,
                        const std::map<mio::abm::LocationType, mio::abm::LocationId>& locations,
                        const std::vector<double>& infection_dist, mio::abm::TimePoint start_date)
{
    for (auto& person : model.get_persons()) {
        // Assign infection state
        mio::abm::InfectionState infection_state = mio::abm::InfectionState(
            mio::DiscreteDistribution<size_t>::get_instance()(model.get_rng(), infection_dist));

        if (infection_state != mio::abm::InfectionState::Susceptible) {
            auto rng = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         model.parameters, start_date, infection_state));
        }

        // Assign shared locations (within same model)
        auto mid = model.get_id();

        person.set_assigned_location(mio::abm::LocationType::SocialEvent,
                                     locations.at(mio::abm::LocationType::SocialEvent), mid);
        person.set_assigned_location(mio::abm::LocationType::BasicsShop,
                                     locations.at(mio::abm::LocationType::BasicsShop), mid);
        person.set_assigned_location(mio::abm::LocationType::Hospital, locations.at(mio::abm::LocationType::Hospital),
                                     mid);
        person.set_assigned_location(mio::abm::LocationType::ICU, locations.at(mio::abm::LocationType::ICU), mid);

        // School for children (age groups 0 and 1)
        if (person.get_age() == AGE_5_14) {
            person.set_assigned_location(mio::abm::LocationType::School, locations.at(mio::abm::LocationType::School),
                                         mid);
        }

        // Work for adults (age groups 2 and 3) - commuting handled separately
        if (person.get_age() == AGE_15_34 || person.get_age() == AGE_35_59) {
            person.set_assigned_location(mio::abm::LocationType::Work, locations.at(mio::abm::LocationType::Work), mid);
        }
    }
}

/**
 * @brief Sets up cross-node commuting: a fraction of workers are assigned to
 *        a work location in a neighboring node.
 *
 * @param models Vector of all graph models.
 * @param all_locations Vector of location maps per node.
 * @param adjacency The adjacency list (neighbors per node).
 * @param commuter_fraction Fraction of workers that commute.
 * @param rng Random number generator.
 */
void setup_commuting(std::vector<mio::GraphABModel>& models,
                     const std::vector<std::map<mio::abm::LocationType, mio::abm::LocationId>>& all_locations,
                     const std::vector<std::vector<size_t>>& adjacency, double commuter_fraction, std::mt19937& rng)
{
    for (size_t node_idx = 0; node_idx < models.size(); ++node_idx) {
        auto& model     = models[node_idx];
        auto& neighbors = adjacency[node_idx];

        if (neighbors.empty())
            continue;

        std::uniform_int_distribution<size_t> neighbor_dist(0, neighbors.size() - 1);
        std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

        for (auto& person : model.get_persons()) {
            // Only working-age adults can commute
            if (person.get_age() != AGE_15_34 && person.get_age() != AGE_35_59)
                continue;

            if (prob_dist(rng) < commuter_fraction) {
                // Assign this worker to a neighboring node's workplace
                size_t target_node = neighbors[neighbor_dist(rng)];
                person.set_assigned_location(mio::abm::LocationType::Work,
                                             all_locations[target_node].at(mio::abm::LocationType::Work),
                                             models[target_node].get_id());
            }
        }
    }
}

/**
 * @brief Generates a random sparse adjacency list for the graph.
 *
 * Creates a connected graph where each node is connected to ~3-5 neighbors.
 * Connection probability decreases with "distance" (node index difference).
 *
 * @returns Adjacency list and weighted adjacency matrix.
 */
std::pair<std::vector<std::vector<size_t>>, std::vector<std::vector<double>>> generate_graph_topology(size_t num_nodes,
                                                                                                      std::mt19937& rng)
{
    std::vector<std::vector<size_t>> adj_list(num_nodes);
    std::vector<std::vector<double>> adj_matrix(num_nodes, std::vector<double>(num_nodes, 0.0));

    if (num_nodes <= 1) {
        return {adj_list, adj_matrix};
    }

    // First: create a ring to ensure connectivity
    for (size_t i = 0; i < num_nodes; ++i) {
        size_t j = (i + 1) % num_nodes;
        if (i != j) {
            adj_list[i].push_back(j);
            adj_list[j].push_back(i);
            double weight    = std::uniform_real_distribution<double>(0.3, 1.0)(rng);
            adj_matrix[i][j] = weight;
            adj_matrix[j][i] = weight;
        }
    }

    // Add random additional edges for denser connectivity
    size_t extra_edges = std::max<size_t>(num_nodes, 2);
    std::uniform_int_distribution<size_t> node_dist(0, num_nodes - 1);

    for (size_t e = 0; e < extra_edges; ++e) {
        size_t i = node_dist(rng);
        size_t j = node_dist(rng);
        if (i == j || adj_matrix[i][j] > 0)
            continue;

        double weight = std::uniform_real_distribution<double>(0.1, 0.5)(rng);
        adj_list[i].push_back(j);
        adj_list[j].push_back(i);
        adj_matrix[i][j] = weight;
        adj_matrix[j][i] = weight;
    }

    return {adj_list, adj_matrix};
}

/**
 * @brief Generates random damping (NPI intervention) schedule.
 *
 * Returns pairs of (day, severity) where severity is in [0.3, 0.9].
 * Days are drawn uniformly from [2, num_days-2], sorted ascending.
 */
std::vector<std::pair<size_t, double>> generate_dampings(size_t num_days, size_t max_dampings, std::mt19937& rng)
{
    std::uniform_int_distribution<size_t> n_damp_dist(0, max_dampings);
    size_t n_dampings = n_damp_dist(rng);

    if (n_dampings == 0 || num_days < 5) {
        return {};
    }

    std::vector<std::pair<size_t, double>> dampings;
    std::uniform_int_distribution<size_t> day_dist(2, std::max<size_t>(num_days - 2, 3));
    std::uniform_real_distribution<double> severity_dist(0.3, 0.9);

    std::vector<size_t> days;
    for (size_t i = 0; i < n_dampings; ++i) {
        size_t day = day_dist(rng);
        // Ensure minimum distance of 3 days between dampings
        bool too_close = false;
        for (auto d : days) {
            if (std::abs(static_cast<int>(day) - static_cast<int>(d)) < 3) {
                too_close = true;
                break;
            }
        }
        if (!too_close) {
            days.push_back(day);
        }
    }

    std::sort(days.begin(), days.end());

    for (auto day : days) {
        dampings.push_back({day, severity_dist(rng)});
    }

    return dampings;
}

// Output Functions

/**
 * @brief Writes per-node time series to CSV.
 *
 * Format: time, age0_S, age0_E, age0_INS, age0_ISy, age0_ISev, age0_ICrit, age0_R, age0_D, age1_S, ...
 */
void write_timeseries(const mio::TimeSeries<ScalarType>& ts, const std::string& filepath)
{
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open " << filepath << " for writing.\n";
        return;
    }

    // Header
    out << "time";
    const char* state_names[] = {"S", "E", "INS", "ISy", "ISev", "ICrit", "R", "D"};
    for (size_t age = 0; age < NUM_AGE_GROUPS; ++age) {
        for (size_t s = 0; s < NUM_INFECTION_STATES; ++s) {
            out << ",age" << age << "_" << state_names[s];
        }
    }
    out << "\n";

    // Data
    out << std::fixed << std::setprecision(4);
    for (Eigen::Index t = 0; t < ts.get_num_time_points(); ++t) {
        out << ts.get_time(t);
        for (Eigen::Index c = 0; c < ts.get_num_elements(); ++c) {
            out << "," << ts.get_value(t)[c];
        }
        out << "\n";
    }
}

void write_adjacency(const std::vector<std::vector<double>>& adj_matrix, const std::string& filepath)
{
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open " << filepath << " for writing.\n";
        return;
    }

    size_t n = adj_matrix.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (j > 0)
                out << ",";
            out << std::fixed << std::setprecision(4) << adj_matrix[i][j];
        }
        out << "\n";
    }
}

void write_parameters(const SimParameters& params, const std::string& filepath)
{
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open " << filepath << " for writing.\n";
        return;
    }
    // CSV: name,value (one row per parameter  easy to load with pandas)
    out << "name,value\n";
    out << std::fixed << std::setprecision(6);
    const double values[NUM_CALIB_PARAMS] = {params.transmission_rate_scale,    params.contact_rate_work,
                                             params.contact_rate_school,        params.contact_rate_social,
                                             params.contact_rate_shop,          params.initial_infected_fraction,
                                             params.initial_recovered_fraction, params.commuter_fraction};
    for (size_t i = 0; i < NUM_CALIB_PARAMS; ++i) {
        out << PARAM_NAMES[i] << "," << values[i] << "\n";
    }
}

void write_dampings(const std::vector<std::pair<size_t, double>>& dampings, const std::string& filepath)
{
    std::ofstream out(filepath);
    out << "day,severity\n";
    for (auto& [day, severity] : dampings) {
        out << day << "," << std::fixed << std::setprecision(4) << severity << "\n";
    }
}

void write_metadata(const Config& cfg, const SimParameters& params, size_t actual_total_population,
                    double runtime_seconds, const std::string& filepath)
{
    std::ofstream out(filepath);
    out << "num_nodes=" << cfg.num_nodes << "\n"
        << "persons_per_node=" << cfg.persons_per_node << "\n"
        << "actual_total_population=" << actual_total_population << "\n"
        << "num_days=" << cfg.num_days << "\n"
        << "run_id=" << cfg.run_id << "\n"
        << "seed=" << cfg.seed << "\n"
        << "max_dampings=" << cfg.max_dampings << "\n"
        << "param_sampling=" << (cfg.param_sampling ? "true" : "false") << "\n"
        << "commuter_fraction=" << params.commuter_fraction << "\n"
        << "num_age_groups=" << NUM_AGE_GROUPS << "\n"
        << "num_infection_states=" << NUM_INFECTION_STATES << "\n"
        << "num_calib_params=" << NUM_CALIB_PARAMS << "\n"
        << "runtime_seconds=" << std::fixed << std::setprecision(3) << runtime_seconds << "\n";
}

int main(int argc, char* argv[])
{
    mio::set_log_level(mio::LogLevel::warn);

    Config cfg = parse_args(argc, argv);

    std::cout << "=== Graph ABM Data Generation for GNN Surrogate ===\n"
              << "  Nodes: " << cfg.num_nodes << "\n"
              << "  Persons/node: " << cfg.persons_per_node << "\n"
              << "  Days: " << cfg.num_days << "\n"
              << "  Run ID: " << cfg.run_id << "\n"
              << "  Seed: " << cfg.seed << "\n"
              << "  Max dampings: " << cfg.max_dampings << "\n"
              << "  Param sampling: " << (cfg.param_sampling ? "on" : "off") << "\n"
              << "  Output: " << cfg.output_dir << "\n"
              << std::endl;

    // Create output directory
    std::ostringstream run_dir_ss;
    run_dir_ss << cfg.output_dir << "/run_" << std::setw(4) << std::setfill('0') << cfg.run_id;
    std::string run_dir = run_dir_ss.str();
    std::filesystem::create_directories(run_dir);

    // Seed RNG
    std::mt19937 rng(cfg.seed + cfg.run_id);
    mio::thread_local_rng().seed({static_cast<uint32_t>(cfg.seed + cfg.run_id)});

    auto start_date = mio::abm::TimePoint(0);
    auto end_date   = start_date + mio::abm::days(static_cast<int>(cfg.num_days));

    // --- Step 0: Sample (or set) calibration parameters for this run ---
    SimParameters sim_params =
        cfg.param_sampling ? sample_parameters(rng) : default_parameters(cfg.fixed_commuter_fraction);
    write_parameters(sim_params, run_dir + "/parameters.csv");
    std::cout << "  Calibration parameters sampled (transmission_rate_scale=" << std::fixed << std::setprecision(3)
              << sim_params.transmission_rate_scale << ", initial_infected=" << sim_params.initial_infected_fraction
              << ")." << std::endl;

    // --- Step 1: Generate graph topology ---
    auto [adj_list, adj_matrix] = generate_graph_topology(cfg.num_nodes, rng);
    write_adjacency(adj_matrix, run_dir + "/adjacency.csv");
    std::cout << "  Graph topology generated." << std::endl;

    // --- Step 2: Generate random dampings ---
    auto dampings = generate_dampings(cfg.num_days, cfg.max_dampings, rng);
    write_dampings(dampings, run_dir + "/dampings.csv");
    std::cout << "  Dampings: " << dampings.size() << " interventions." << std::endl;

    // --- Step 3: Create models for each node ---
    HouseholdMembers hh_members;
    std::vector<mio::GraphABModel> models;
    std::vector<std::map<mio::abm::LocationType, mio::abm::LocationId>> all_locations;
    models.reserve(cfg.num_nodes);
    all_locations.reserve(cfg.num_nodes);

    // Initial infection-state distribution from sampled parameters
    // (same regime across all nodes, like the ODE paper)
    auto infection_dist = build_infection_distribution(sim_params);

    for (size_t i = 0; i < cfg.num_nodes; ++i) {
        models.emplace_back(NUM_AGE_GROUPS, static_cast<int>(i));
        auto& model = models.back();

        // Set infection parameters (transmission rate from sim_params)
        set_infection_parameters(model, sim_params);

        // Create population via households
        size_t node_pop = cfg.persons_per_node;
        create_population(model, node_pop, hh_members);

        // Create shared locations (contact rates from sim_params)
        auto locs = create_locations(model, node_pop, sim_params);
        all_locations.push_back(locs);

        // Initialize persons (infection states + location assignments)
        initialize_persons(model, locs, infection_dist, start_date);

        // Apply dampings (interventions)
        for (auto& [day, severity] : dampings) {
            auto t_damping = start_date + mio::abm::days(static_cast<int>(day));
            // Apply home office + social event restrictions
            mio::abm::set_home_office(t_damping, severity, model.parameters);
            mio::abm::close_social_events(t_damping, severity, model.parameters);
            if (severity > 0.5) {
                mio::abm::set_school_closure(t_damping, severity * 0.8, model.parameters);
            }
        }

        std::cout << "  Node " << i << ": " << model.get_persons().size() << " persons created." << std::endl;
    }

    // --- Step 4: Set up cross-node commuting (commuter fraction from sim_params) ---
    setup_commuting(models, all_locations, adj_list, sim_params.commuter_fraction, rng);
    std::cout << "  Commuting setup complete." << std::endl;

    // --- Step 5: Build graph simulation ---
    using HistoryType = mio::History<mio::DataWriterToMemory, LogAgeStratifiedState>;

    mio::Graph<mio::ABMSimulationNode<HistoryType>, mio::ABMMobilityEdge<HistoryType>> graph;

    // Add nodes
    for (size_t i = 0; i < cfg.num_nodes; ++i) {
        graph.add_node(models[i].get_id(), HistoryType{}, start_date, std::move(models[i]));
    }

    // Add edges (bidirectional for each adjacency entry)
    for (size_t i = 0; i < cfg.num_nodes; ++i) {
        for (size_t j : adj_list[i]) {
            if (j > i) { // avoid duplicate edges
                graph.add_edge(static_cast<int>(i), static_cast<int>(j));
                graph.add_edge(static_cast<int>(j), static_cast<int>(i));
            }
        }
    }

    auto exchange_time_span = mio::abm::hours(12);
    auto sim                = mio::make_abm_graph_sim<HistoryType>(start_date, exchange_time_span, std::move(graph));

    // --- Step 6: Run simulation ---
    std::cout << "  Running simulation for " << cfg.num_days << " days..." << std::flush;
    auto t_start = std::chrono::high_resolution_clock::now();

    sim.advance(end_date);

    auto t_end          = std::chrono::high_resolution_clock::now();
    double runtime_secs = std::chrono::duration<double>(t_end - t_start).count();
    std::cout << " done (" << std::fixed << std::setprecision(2) << runtime_secs << "s)" << std::endl;

    // --- Step 7: Extract and write results ---
    // DataWriterToMemory stores a vector of (TimePoint, VectorX) pairs per logger.
    // Convert to TimeSeries for output.
    size_t total_population = 0;

    auto graph_nodes = sim.get_graph().nodes();
    for (size_t i = 0; i < graph_nodes.size(); ++i) {
        auto& node = graph_nodes[i];
        // get_history() returns tuple<History...>; get first History, then its log
        auto& history           = std::get<0>(node.property.get_history());
        const auto& log_data    = history.get_log();
        const auto& logged_data = std::get<0>(log_data); // vector of (TimePoint, VectorX)

        // Convert logged data to TimeSeries
        Eigen::Index num_elements = static_cast<Eigen::Index>(NUM_AGE_GROUPS * NUM_INFECTION_STATES);
        mio::TimeSeries<ScalarType> ts(num_elements);

        for (const auto& entry : logged_data) {
            ts.add_time_point(entry.first.days(), entry.second);
        }

        std::ostringstream ts_path;
        ts_path << run_dir << "/timeseries_node_" << std::setw(3) << std::setfill('0') << i << ".csv";
        write_timeseries(ts, ts_path.str());

        // Count final population in this node
        size_t node_pop = 0;
        if (ts.get_num_time_points() > 0) {
            auto last = ts.get_last_value();
            for (Eigen::Index c = 0; c < last.size(); ++c) {
                node_pop += static_cast<size_t>(last[c]);
            }
        }
        total_population += node_pop;
    }

    write_metadata(cfg, sim_params, total_population, runtime_secs, run_dir + "/metadata.txt");

    std::cout << "\n=== Run " << cfg.run_id << " complete ===" << std::endl;
    std::cout << "  Total population: " << total_population << std::endl;
    std::cout << "  Runtime: " << std::fixed << std::setprecision(2) << runtime_secs << "s" << std::endl;
    std::cout << "  Output: " << run_dir << std::endl;

    return 0;
}
