/*
* Copyright (C) 2020-2024 MEmilio
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

#include "memilio/compartments/parameter_studies.h"
#include "memilio/geography/regions.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/result_io.h"
#include "memilio/io/mobility_io.h"
#include "memilio/utils/date.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/random_number_generator.h"
#include "ode_secirts/parameters.h"
#include "ode_secirts/parameters_io.h"
#include "ode_secirts/parameter_space.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/mobility/metapopulation_mobility_detailed.h"
#include "memilio/utils/stl_util.h"
#include "boost/filesystem.hpp"
#include <cstdio>
#include <iomanip>

namespace fs = boost::filesystem;

/**
 * Assigns a uniform distribution to an UncertainValue with a specified range.
 * The value is set to the average of min and max, and the distribution is UNIFORM(min, max).
 * @param[in,out] param UncertainValue to configure.
 * @param[in] min Lower bound of the uniform distribution.
 * @param[in] max Upper bound of the uniform distribution.
 */
void assign_uniform_distribution(mio::UncertainValue<double>& param, double min, double max)
{
    param = mio::UncertainValue<double>(0.5 * (min + max));
    param.set_distribution(mio::ParameterDistributionUniform(min, max));
}

/**
 * Assigns uniform distributions to an array of UncertainValues.
 * Each element i is set to the average of min[i] and max[i] with a UNIFORM(min[i], max[i]) distribution.
 * @param[in,out] array Array of UncertainValues to configure.
 * @param[in] min Array of lower bounds for each element.
 * @param[in] max Array of upper bounds for each element.
 * @tparam N Size of the array, must match the number of elements in min and max.
 */
template <size_t N>
void assign_uniform_distribution_array(mio::CustomIndexArray<mio::UncertainValue<double>, mio::AgeGroup>& array,
                                       const double (&min)[N], const double (&max)[N])
{
    assert(N == array.numel() && "Array size must match the number of elements in min and max.");
    for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(N); ++i) {
        assign_uniform_distribution(array[i], min[static_cast<size_t>(i)], max[static_cast<size_t>(i)]);
    }
}

/**
 * Assigns a uniform distribution to all elements of an array of UncertainValues using a single range.
 * Each element is set to the average of min and max with a UNIFORM(min, max) distribution.
 * @param[in,out] array Array of UncertainValues to configure.
 * @param[in] min Lower bound of the uniform distribution applied to all elements.
 * @param[in] max Upper bound of the uniform distribution applied to all elements.
 */
void assign_uniform_distribution_array(mio::CustomIndexArray<mio::UncertainValue<double>, mio::AgeGroup>& array,
                                       double min, double max)
{
    for (auto i = mio::AgeGroup(0); i < array.size<mio::AgeGroup>(); ++i) {
        assign_uniform_distribution(array[i], min, max);
    }
}

/**
 * Configures epidemiological parameters for COVID-19 model (Omicron variant) based on literature.
 * @param[in,out] params Object that the parameters will be added to.
 * @return IOResult<void> indicating success (currently no failure cases defined).
 */
mio::IOResult<void> set_covid_parameters(mio::osecirts::Parameters<double>& params)
{
    constexpr size_t num_age_groups = 6;

    // --- Transition Times ---
    // Incubation and infectious periods sourced from literature
    // Sources: doi.org/10.1016/j.lanepe.2022.100446, doi.org/10.3201/eid2806.220158
    const double time_exposed_min              = 1.66;
    const double time_exposed_max              = 1.66;
    const double time_infected_no_symptoms_min = 1.44;
    const double time_infected_no_symptoms_max = 1.44;

    // Symptomatic period: doi.org/10.1016/S0140-6736(22)00327-0
    const double time_infected_symptoms_min = 6.58;
    const double time_infected_symptoms_max = 7.16;

    // Severe and critical periods: doi.org/10.1186/s12879-022-07971-6
    const double time_infected_severe_min[num_age_groups]   = {1.8, 1.8, 1.8, 2.5, 3.5, 4.91};
    const double time_infected_severe_max[num_age_groups]   = {2.3, 2.3, 2.3, 3.67, 5, 7.01};
    const double time_infected_critical_min[num_age_groups] = {9.29, 9.29, 9.29, 10.842, 11.15, 11.07};
    const double time_infected_critical_max[num_age_groups] = {10.57, 10.57, 10.57, 12.86, 13.23, 13.25};

    assign_uniform_distribution_array(params.get<mio::osecirts::TimeExposed<double>>(), time_exposed_min,
                                      time_exposed_max);
    assign_uniform_distribution_array(params.get<mio::osecirts::TimeInfectedNoSymptoms<double>>(),
                                      time_infected_no_symptoms_min, time_infected_no_symptoms_max);
    assign_uniform_distribution_array(params.get<mio::osecirts::TimeInfectedSymptoms<double>>(),
                                      time_infected_symptoms_min, time_infected_symptoms_max);
    assign_uniform_distribution_array(params.get<mio::osecirts::TimeInfectedSevere<double>>(), time_infected_severe_min,
                                      time_infected_severe_max);
    assign_uniform_distribution_array(params.get<mio::osecirts::TimeInfectedCritical<double>>(),
                                      time_infected_critical_min, time_infected_critical_max);

    // --- Transmission Probabilities ---
    // Adjusted for Omicron variant: doi.org/10.7554/eLife.78933
    const double variant_factor                                = 1.94;
    const double transmission_prob_contact_min[num_age_groups] = {0.02 * variant_factor, 0.05 * variant_factor,
                                                                  0.05 * variant_factor, 0.05 * variant_factor,
                                                                  0.08 * variant_factor, 0.10 * variant_factor};
    const double transmission_prob_contact_max[num_age_groups] = {0.04 * variant_factor, 0.07 * variant_factor,
                                                                  0.07 * variant_factor, 0.07 * variant_factor,
                                                                  0.10 * variant_factor, 0.15 * variant_factor};

    // Relative transmission from asymptomatic cases (fixed for simplicity, could use DOI: 10.1097/INF.0000000000003791)
    const double rel_transmission_no_symptoms_min = 0.5;
    const double rel_transmission_no_symptoms_max = 0.5;

    // Risk of infection from symptomatic cases (depends on incidence and testing capacity)
    const double risk_infection_symptomatic_min     = 0.0;
    const double risk_infection_symptomatic_max     = 0.2;
    const double max_risk_infection_symptomatic_min = 0.4;
    const double max_risk_infection_symptomatic_max = 0.5;

    assign_uniform_distribution_array(params.get<mio::osecirts::TransmissionProbabilityOnContact<double>>(),
                                      transmission_prob_contact_min, transmission_prob_contact_max);
    assign_uniform_distribution_array(params.get<mio::osecirts::RelativeTransmissionNoSymptoms<double>>(),
                                      rel_transmission_no_symptoms_min, rel_transmission_no_symptoms_max);
    assign_uniform_distribution_array(params.get<mio::osecirts::RiskOfInfectionFromSymptomatic<double>>(),
                                      risk_infection_symptomatic_min, risk_infection_symptomatic_max);
    assign_uniform_distribution_array(params.get<mio::osecirts::MaxRiskOfInfectionFromSymptomatic<double>>(),
                                      max_risk_infection_symptomatic_min, max_risk_infection_symptomatic_max);

    // --- Disease Progression Probabilities ---
    // Recovery from asymptomatic infection: doi.org/10.1101/2022.05.05.22274697
    const double recovered_per_infected_no_symptoms_min[num_age_groups] = {0.2, 0.25, 0.2, 0.2, 0.175, 0.1};
    const double recovered_per_infected_no_symptoms_max[num_age_groups] = {0.4, 0.45, 0.35, 0.3, 0.25, 0.15};

    // Severe cases from symptomatic infection (2021 data with factors): doi.org/10.1016/S0140-6736(22)00462-7
    const double severe_per_infected_symptoms_min[num_age_groups] = {1 * 0.006,   0.8 * 0.006, 0.4 * 0.015,
                                                                     0.3 * 0.049, 0.25 * 0.15, 0.35 * 0.2};
    const double severe_per_infected_symptoms_max[num_age_groups] = {1 * 0.009,   0.8 * 0.009, 0.4 * 0.023,
                                                                     0.3 * 0.074, 0.25 * 0.18, 0.35 * 0.25};

    // Critical cases from severe (Delta-adjusted, factor 0.52): doi.org/10.1177/14034948221108548
    const double critical_factor                         = 0.52;
    const double critical_per_severe_min[num_age_groups] = {0.05 * critical_factor, 0.05 * critical_factor,
                                                            0.05 * critical_factor, 0.10 * critical_factor,
                                                            0.25 * critical_factor, 0.35 * critical_factor};
    const double critical_per_severe_max[num_age_groups] = {0.10 * critical_factor, 0.10 * critical_factor,
                                                            0.10 * critical_factor, 0.20 * critical_factor,
                                                            0.35 * critical_factor, 0.45 * critical_factor};

    // Deaths from critical cases (factor 0.39): doi.org/10.1136/bmjgh-2023-012328
    const double death_factor                            = 0.39;
    const double deaths_per_critical_min[num_age_groups] = {0.00 * death_factor, 0.00 * death_factor,
                                                            0.10 * death_factor, 0.10 * death_factor,
                                                            0.30 * death_factor, 0.50 * death_factor};
    const double deaths_per_critical_max[num_age_groups] = {0.10 * death_factor, 0.10 * death_factor,
                                                            0.18 * death_factor, 0.18 * death_factor,
                                                            0.50 * death_factor, 0.70 * death_factor};

    assign_uniform_distribution_array(params.get<mio::osecirts::RecoveredPerInfectedNoSymptoms<double>>(),
                                      recovered_per_infected_no_symptoms_min, recovered_per_infected_no_symptoms_max);
    assign_uniform_distribution_array(params.get<mio::osecirts::SeverePerInfectedSymptoms<double>>(),
                                      severe_per_infected_symptoms_min, severe_per_infected_symptoms_max);
    assign_uniform_distribution_array(params.get<mio::osecirts::CriticalPerSevere<double>>(), critical_per_severe_min,
                                      critical_per_severe_max);
    assign_uniform_distribution_array(params.get<mio::osecirts::DeathsPerCritical<double>>(), deaths_per_critical_min,
                                      deaths_per_critical_max);

    // --- Immunity Parameters ---
    // Exposure reduction (no reduction assumed here)
    const double reduc_exposed_partial_immunity_min  = 1.0;
    const double reduc_exposed_partial_immunity_max  = 1.0;
    const double reduc_exposed_improved_immunity_min = 1.0;
    const double reduc_exposed_improved_immunity_max = 1.0;

    // Symptom reduction: doi.org/10.1056/NEJMoa2119451
    const double reduc_infected_symptoms_partial_min  = 0.746;
    const double reduc_infected_symptoms_partial_max  = 0.961;
    const double reduc_infected_symptoms_improved_min = 0.295;
    const double reduc_infected_symptoms_improved_max = 0.344;

    // Severe/critical/death reduction
    // Partial immunity: doi.org/10.1056/NEJMoa2119451 (week 4 report)
    const double reduc_severe_critical_dead_partial_min = 0.52;
    const double reduc_severe_critical_dead_partial_max = 0.82;
    // Improved immunity: doi.org/10.1136/bmj-2022-071502
    const double reduc_severe_critical_dead_improved_min = 0.1;
    const double reduc_severe_critical_dead_improved_max = 0.19;

    // Time reduction for mild infections: doi.org/10.1101/2021.09.24.21263978
    const double reduc_time_infected_mild = 0.5;

    assign_uniform_distribution_array(params.get<mio::osecirts::ReducExposedPartialImmunity<double>>(),
                                      reduc_exposed_partial_immunity_min, reduc_exposed_partial_immunity_max);
    assign_uniform_distribution_array(params.get<mio::osecirts::ReducExposedImprovedImmunity<double>>(),
                                      reduc_exposed_improved_immunity_min, reduc_exposed_improved_immunity_max);
    assign_uniform_distribution_array(params.get<mio::osecirts::ReducInfectedSymptomsPartialImmunity<double>>(),
                                      reduc_infected_symptoms_partial_min, reduc_infected_symptoms_partial_max);
    assign_uniform_distribution_array(params.get<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<double>>(),
                                      reduc_infected_symptoms_improved_min, reduc_infected_symptoms_improved_max);
    assign_uniform_distribution_array(
        params.get<mio::osecirts::ReducInfectedSevereCriticalDeadPartialImmunity<double>>(),
        reduc_severe_critical_dead_partial_min, reduc_severe_critical_dead_partial_max);
    assign_uniform_distribution_array(
        params.get<mio::osecirts::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>(),
        reduc_severe_critical_dead_improved_min, reduc_severe_critical_dead_improved_max);
    assign_uniform_distribution_array(params.get<mio::osecirts::ReducTimeInfectedMild<double>>(),
                                      reduc_time_infected_mild, reduc_time_infected_mild);

    // --- Seasonality ---
    // Seasonal variation
    const double seasonality_min = 0.1;
    const double seasonality_max = 0.3;
    assign_uniform_distribution(params.get<mio::osecirts::Seasonality<double>>(), seasonality_min, seasonality_max);

    // --- Variant-Specific Parameters ---
    params.get<mio::osecirts::StartDayNewVariant>() = mio::get_day_in_year(mio::Date(2022, 6, 6));

    // --- Waning Immunity Durations ---
    // Temporary immunity periods: doi.org/10.1016/S1473-3099(22)00801-5
    const double immunity_interval_partial_min  = 60;
    const double immunity_interval_partial_max  = 60;
    const double immunity_interval_improved_min = 60;
    const double immunity_interval_improved_max = 60;

    assign_uniform_distribution_array(params.get<mio::osecirts::TimeTemporaryImmunityPI<double>>(),
                                      immunity_interval_partial_min, immunity_interval_partial_max);
    assign_uniform_distribution_array(params.get<mio::osecirts::TimeTemporaryImmunityII<double>>(),
                                      immunity_interval_improved_min, immunity_interval_improved_max);

    // Waning immunity duration: doi.org/10.1016/S1473-3099(22)00801-5
    params.get<mio::osecirts::TimeWaningPartialImmunity<double>>()  = 365.0;
    params.get<mio::osecirts::TimeWaningImprovedImmunity<double>>() = 365.0;

    return mio::success();
}

enum class ContactLocation
{
    Home = 0,
    School,
    Work,
    Other,
    Count,
};

/**
 * different types of NPI, used as DampingType.
 */
enum class Intervention
{
    Home,
    SchoolClosure,
    HomeOffice,
    GatheringBanFacilitiesClosure,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
};

/**
 * different level of NPI, used as DampingLevel.
 */
enum class InterventionLevel
{
    Main,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
    Holidays,
};

static const std::map<ContactLocation, std::string> contact_locations = {{ContactLocation::Home, "home"},
                                                                         {ContactLocation::School, "school_pf_eig"},
                                                                         {ContactLocation::Work, "work"},
                                                                         {ContactLocation::Other, "other"}};

/**
 * Configures contact matrices for the model by reading data from files and applying scaling adjustments.
 * @param[in] data_dir Directory containing contact data files.
 * @param[in,out] params Parameters object to store the configured contact matrices.
 * @param[in] avg_transport_time Average transport time (default: 0.0), used to scale contact frequencies.
 * @param[in] share_staying_local Proportion of individuals staying local (default: 1.0), affects contact scaling.
 * @return IOResult indicating success or an IO error if file reading fails.
 */
mio::IOResult<void> set_contact_matrices(const fs::path& data_dir, mio::osecirts::Parameters<double>& params,
                                         double avg_transport_time = 0.0, double share_staying_local = 1.0)
{
    // Read transport contact matrix
    BOOST_OUTCOME_TRY(auto&& transport_matrix,
                      mio::read_mobility_plain((data_dir / "contacts" / "contacts_transport.txt").string()));

    // Init contact matrices for all locations
    constexpr size_t num_age_groups = 6;
    mio::ContactMatrixGroup contact_matrices(contact_locations.size(), num_age_groups);

    for (const auto& [location_id, location_name] : contact_locations) {
        // Read baseline contact matrix for this location
        BOOST_OUTCOME_TRY(
            auto&& baseline,
            mio::read_mobility_plain((data_dir / "contacts" / ("baseline_" + location_name + ".txt")).string()));

        // Adjust "other" location by subtracting transport contacts
        if (location_name == "other") {
            baseline = (baseline - transport_matrix).cwiseAbs(); // Ensure non-negative values
        }

        // Scale contacts based on travel time and local staying proportion
        auto scaled_baseline =
            (1.0 - share_staying_local) * baseline / (1.0 - avg_transport_time) + share_staying_local * baseline;

        // Adjust contacts: only scale for commuting age groups (2-5), others remain unchanged
        Eigen::MatrixXd adjusted_baseline = Eigen::MatrixXd::Zero(num_age_groups, num_age_groups);
        for (size_t i = 0; i < num_age_groups; ++i) {
            for (size_t j = 0; j < num_age_groups; ++j) {
                if ((i >= 2 && i <= 5) || (j >= 2 && j <= 5)) {
                    adjusted_baseline(i, j) = scaled_baseline(i, j);
                }
                else {
                    adjusted_baseline(i, j) = baseline(i, j);
                }
            }
        }

        // Assign adjusted baseline and zero minimum to the contact matrix
        contact_matrices[static_cast<size_t>(location_id)].get_baseline() = adjusted_baseline;
        contact_matrices[static_cast<size_t>(location_id)].get_minimum() =
            Eigen::MatrixXd::Zero(num_age_groups, num_age_groups);
    }

    // Store the configured contact matrices in the parameters object
    params.get<mio::osecirts::ContactPatterns<double>>() = mio::UncertainContactMatrix<double>(contact_matrices);

    return mio::success();
}

/**
 * Configures a transport-specific contact matrix for the model by reading data from a file.
 * @param[in] data_dir Directory containing the transport contact data file ("contacts/contacts_transport.txt").
 * @param[in,out] params Parameters object to store the configured transport contact matrix.
 * @return IOResult<void> indicating success or an IO error if file reading fails.
 */
mio::IOResult<void> set_contact_matrices_transport(const fs::path& data_dir, mio::osecirts::Parameters<double>& params)
{
    // Define the path to the transport contacts file
    fs::path transport_file = data_dir / "contacts" / "contacts_transport.txt";

    // number of age groups
    const auto num_age_groups = static_cast<size_t>(params.get_num_groups());

    // Read the transport contact matrix from the file
    BOOST_OUTCOME_TRY(auto&& transport_matrix, mio::read_mobility_plain(transport_file.string()));

    // Initialize contact matrices with a single group for transport
    auto contact_matrices = mio::ContactMatrixGroup(1, num_age_groups);

    // Assign the transport matrix as the baseline and set a zero minimum
    contact_matrices[0].get_baseline() = transport_matrix;
    contact_matrices[0].get_minimum()  = Eigen::MatrixXd::Zero(num_age_groups, num_age_groups);

    // Store the configured contact matrix in the parameters object
    params.get<mio::osecirts::ContactPatterns<double>>() = mio::UncertainContactMatrix<double>(contact_matrices);

    return mio::success();
}

/**
 * Scales contact matrices for all nodes in a graph based on commuting data and travel time.
 * Computes average travel time and commuter share, then applies these to adjust local contact patterns.
 * @param[in,out] params_graph Graph containing model parameters for each node, updated with scaled contact matrices.
 * @param[in] data_dir Directory containing baseline contact data files.
 * @param[in] mobility_data_dir Path to the file with commuter mobility data.
 * @param[in] commuting_weights Weights for each age group to compute total population (size must match age groups).
 * @return IOResult<void> indicating success or an IO error if file reading or contact matrix setup fails.
 */
mio::IOResult<void> scale_contacts_local(mio::ExtendedGraph<mio::osecirts::Model<double>>& params_graph,
                                         const fs::path& data_dir, const std::string& mobility_data_dir,
                                         const std::vector<double>& commuting_weights)
{
    // Read commuter mobility data
    BOOST_OUTCOME_TRY(auto&& mobility_data_commuter, mio::read_mobility_plain(mobility_data_dir));

    // Calculate average travel time across all edges
    const double avg_travel_time = std::accumulate(params_graph.edges().begin(), params_graph.edges().end(), 0.0,
                                                   [](double sum, const auto& edge) {
                                                       return sum + edge.property.travel_time;
                                                   }) /
                                   static_cast<double>(params_graph.edges().size());

    // Calculate total population weighted by age group commuting factors
    const double total_population =
        std::accumulate(params_graph.nodes().begin(), params_graph.nodes().end(), 0.0,
                        [&commuting_weights](double sum, const auto& node) {
                            const auto& populations = node.property.base_sim.populations;
                            return sum + populations.get_group_total(mio::AgeGroup(0)) * commuting_weights[0] +
                                   populations.get_group_total(mio::AgeGroup(1)) * commuting_weights[1] +
                                   populations.get_group_total(mio::AgeGroup(2)) * commuting_weights[2] +
                                   populations.get_group_total(mio::AgeGroup(3)) * commuting_weights[3] +
                                   populations.get_group_total(mio::AgeGroup(4)) * commuting_weights[4] +
                                   populations.get_group_total(mio::AgeGroup(5)) * commuting_weights[5];
                        });

    // Calculate total number of commuters from mobility data
    const double total_commuters =
        std::accumulate(params_graph.edges().begin(), params_graph.edges().end(), 0.0,
                        [&mobility_data_commuter](double sum, const auto& edge) {
                            return sum + mobility_data_commuter(edge.start_node_idx, edge.end_node_idx);
                        });

    // Compute average commuter share relative to total population
    const double avg_commuter_share = total_commuters / total_population;

    // Log results on the root MPI process
    if (mio::mpi::is_root()) {
        std::cout << "Average commuter share: " << avg_commuter_share << "\n";
        std::cout << "Average travel time: " << avg_travel_time << "\n";
    }

    // Scale contact matrices for each node using computed averages
    for (auto& node : params_graph.nodes()) {
        BOOST_OUTCOME_TRY(
            set_contact_matrices(data_dir, node.property.base_sim.parameters, avg_travel_time, avg_commuter_share));
    }

    return mio::success();
}

/**
 * @brief Configures graph nodes for counties with epidemiological data.
 * Reads node IDs and population data from files and initializes models for each node.
 * @param[in] params Model parameters applied to all nodes.
 * @param[in] start_date Start date for reading epidemiological data.
 * @param[in] end_date End date for reading epidemiological data.
 * @param[in] data_dir Directory containing data files (e.g., contact matrices).
 * @param[in] population_data_path Path to JSON file with population data.
 * @param[in] stay_times_data_path Path to TXT file with stay times for local entities.
 * @param[in] is_node_for_county True for county IDs, false for district IDs.
 * @param[in,out] params_graph Graph to populate with nodes.
 * @param[in] read_func Function to read input data and set model compartments.
 * @param[in] node_func Function to get node IDs from population data.
 * @param[in] scaling_factor_inf Factors for confirmed cases to adjust for undetected infections.
 * @param[in] scaling_factor_icu Factor for ICU cases to adjust for underreporting.
 * @param[in] tnt_capacity_factor Factor for test-and-trace capacity relative to population.
 * @param[in] immunity_population Immunity distribution per age group across.
 * @param[in] num_days Number of simulation days (default: 0), used for vaccination data.
 * @param[in] export_time_series If true, exports daily time series data to the input directory.
 * @param[in] rki_age_groups If true, uses RKI age group definitions.
 * @param[in] masks If true, applies mask-related transmission reductions in mobility nodes.
 * @param[in] ffp2 If true, uses FFP2 mask factors; otherwise, uses surgical mask factors.
 * @return IOResult<void> indicating success or an IO error if file reading fails.
 */
template <class ReadFunction, class NodeIdFunction, typename FP = double>
mio::IOResult<void>
set_nodes(const mio::osecirts::Parameters<FP>& params, mio::Date start_date, mio::Date end_date,
          const fs::path& data_dir, const std::string& population_data_path, const std::string& stay_times_data_path,
          bool is_node_for_county, mio::ExtendedGraph<mio::osecirts::Model<FP>>& params_graph, ReadFunction&& read_func,
          NodeIdFunction&& node_func, const std::vector<FP>& scaling_factor_inf, FP scaling_factor_icu,
          FP tnt_capacity_factor, const std::vector<std::vector<double>>& immunity_population, int num_days = 0,
          bool export_time_series = false, bool rki_age_groups = true, bool masks = false, bool ffp2 = false)
{
    // Read stay times for mobility
    BOOST_OUTCOME_TRY(auto&& duration_stay, mio::read_duration_stay(stay_times_data_path));

    // Get node IDs (counties or districts)
    BOOST_OUTCOME_TRY(auto&& node_ids, node_func(population_data_path, is_node_for_county, rki_age_groups));

    // Initialize models for each node
    std::vector<mio::osecirts::Model<FP>> nodes(node_ids.size(),
                                                mio::osecirts::Model<FP>(static_cast<size_t>(params.get_num_groups())));
    for (auto& node : nodes) {
        node.parameters = params; // Assign base parameters to each node
    }

    // Populate nodes with epidemiological data
    BOOST_OUTCOME_TRY(read_func(nodes, start_date, node_ids, scaling_factor_inf, scaling_factor_icu, data_dir.string(),
                                num_days, immunity_population, export_time_series));

    // Configure each node
    for (size_t node_idx = 0; node_idx < nodes.size(); ++node_idx) {
        auto& node = nodes[node_idx];

        // Set test-and-trace capacity with uncertainty
        const FP tnt_capacity = node.populations.get_total() * tnt_capacity_factor;
        auto& tnt_value       = node.parameters.template get<mio::osecirts::TestAndTraceCapacity<FP>>();
        assign_uniform_distribution(tnt_value, 0.8 * tnt_capacity, 1.2 * tnt_capacity);

        // Configure school holidays based on node's state
        const int county_id = static_cast<int>(mio::regions::CountyId(node_ids[node_idx]));
        const auto holiday_periods =
            mio::regions::get_holidays(mio::regions::get_state_id(county_id), start_date, end_date);
        auto& contacts = node.parameters.template get<mio::osecirts::ContactPatterns<FP>>();
        contacts.get_school_holidays().resize(holiday_periods.size());
        std::transform(holiday_periods.begin(), holiday_periods.end(), contacts.get_school_holidays().begin(),
                       [start_date](const auto& period) {
                           return std::make_pair(
                               mio::SimulationTime(mio::get_offset_in_days(period.first, start_date)),
                               mio::SimulationTime(mio::get_offset_in_days(period.second, start_date)));
                       });

        // Add uncertainty to population compartments
        for (auto i = mio::AgeGroup(0); i < params.get_num_groups(); i++) {
            for (auto j = mio::Index<typename mio::osecirts::Model<FP>::Compartments>(0);
                 j < mio::osecirts::Model<FP>::Compartments::Count; ++j) {
                auto& compartment_value = nodes[node_idx].populations[{i, j}];
                compartment_value =
                    mio::UncertainValue(0.5 * (1.1 * double(compartment_value) + 0.9 * double(compartment_value)));
                compartment_value.set_distribution(mio::ParameterDistributionUniform(0.9 * double(compartment_value),
                                                                                     1.1 * double(compartment_value)));
            }
        }

        // Create mobility node with zero population
        auto mobility_model = node;
        mobility_model.populations.set_total(0);

        // Adjust transmission probability in mobility node if masks are enabled
        if (masks) {
            constexpr double surgical_mask_factor = 0.1; // Reduction factor for surgical masks
            constexpr double ffp2_mask_factor     = 0.001; // Reduction factor for FFP2 masks
            double mask_factors[]                 = {1.0,
                                     surgical_mask_factor,
                                     ffp2 ? ffp2_mask_factor : surgical_mask_factor,
                                     ffp2 ? ffp2_mask_factor : surgical_mask_factor,
                                     ffp2 ? ffp2_mask_factor : surgical_mask_factor,
                                     ffp2 ? ffp2_mask_factor : surgical_mask_factor};

            constexpr double variant_factor = 1.94; // Omicron adjustment: doi.org/10.7554/eLife.78933
            double trans_prob_min[]         = {0.02 * variant_factor, 0.05 * variant_factor, 0.05 * variant_factor,
                                       0.05 * variant_factor, 0.08 * variant_factor, 0.10 * variant_factor};
            double trans_prob_max[]         = {0.04 * variant_factor, 0.07 * variant_factor, 0.07 * variant_factor,
                                       0.07 * variant_factor, 0.10 * variant_factor, 0.15 * variant_factor};

            for (int i = 0; i < 6; ++i) {
                trans_prob_min[i] *= mask_factors[i];
                trans_prob_max[i] *= mask_factors[i];
            }
            assign_uniform_distribution_array(
                mobility_model.parameters.template get<mio::osecirts::TransmissionProbabilityOnContact<double>>(),
                trans_prob_min, trans_prob_max);
        }

        // Disable vaccinations in mobility node
        for (int day = 0; day < num_days; ++day) {
            auto t = mio::SimulationDay(static_cast<size_t>(day));
            for (auto age = mio::AgeGroup(0); age < params.get_num_groups(); ++age) {
                mobility_model.parameters.template get<mio::osecirts::DailyPartialVaccinations<double>>()[{age, t}] = 0;
                mobility_model.parameters.template get<mio::osecirts::DailyFullVaccinations<double>>()[{age, t}]    = 0;
                mobility_model.parameters.template get<mio::osecirts::DailyBoosterVaccinations<double>>()[{age, t}] = 0;
            }
        }

        // Add node to graph with its mobility counterpart and stay duration
        params_graph.add_node(node_ids[node_idx], node, mobility_model,
                              duration_stay(static_cast<Eigen::Index>(node_idx)));
    }

    return mio::success();
}

/**
 * @brief Configures graph edges for commuting between counties or districts.
 * Reads commuting matrices, travel times, and paths from files, then creates edges between node pairs with sufficient mobility.
 * @param[in] travel_times_dir Path to TXT file containing travel times between counties.
 * @param[in] mobility_data_dir Path to TXT file containing commuting matrices.
 * @param[in] travel_path_dir Path to TXT file containing paths between counties.
 * @param[in,out] params_graph Graph to add with edges.
 * @param[in] migrating_compartments Infection states that participate in commuting.
 * @param[in] contact_locations_size Number of contact locations (e.g., home, work) for mobility coefficients.
 * @param[in] commuting_weights Weights for each age groups commuting contribution (default: empty, filled with 1.0).
 * @param[in] threshold_edges Minimum commuter coefficient to create an edge (default: 4e-5).
 * @return IOResult<void> indicating success or an IO error if file reading fails.
 */
mio::IOResult<void> set_edges(const std::string& travel_times_dir, const std::string& mobility_data_dir,
                              const std::string& travel_path_dir,
                              mio::ExtendedGraph<mio::osecirts::Model<double>>& params_graph,
                              const std::vector<mio::osecirts::InfectionState>& migrating_compartments,
                              size_t contact_locations_size, std::vector<double> commuting_weights = {},
                              const double threshold_edges = 4e-5)
{
    // Read mobility, travel times, and path data
    BOOST_OUTCOME_TRY(auto&& mobility_data_commuter, mio::read_mobility_plain(mobility_data_dir));
    BOOST_OUTCOME_TRY(auto&& travel_times, mio::read_mobility_plain(travel_times_dir));
    BOOST_OUTCOME_TRY(auto&& path_mobility, mio::read_path_mobility(travel_path_dir));

    // Iterate over all pairs of nodes to set edges
    for (size_t county_idx_i = 0; county_idx_i < params_graph.nodes().size(); ++county_idx_i) {
        for (size_t county_idx_j = 0; county_idx_j < params_graph.nodes().size(); ++county_idx_j) {
            auto& populations = params_graph.nodes()[county_idx_i].property.base_sim.populations;
            const size_t num_age_groups =
                static_cast<size_t>(params_graph.nodes()[county_idx_i].property.base_sim.parameters.get_num_groups());

            // Initialize mobility coefficients aligned with contact locations
            auto mobility_coeffs = mio::MobilityCoefficientGroup(contact_locations_size, populations.numel());

            // Default commuting weights to 1.0 for all age groups if not provided
            if (commuting_weights.empty()) {
                commuting_weights.resize(num_age_groups, 1.0);
            }

            // Calculate working population weighted by commuting factors
            double working_population   = 0.0;
            const auto min_commuter_age = mio::AgeGroup(2); // Youngest commuting age
            const auto max_commuter_age = mio::AgeGroup(4); // Partially retired group
            for (auto age = min_commuter_age; age <= max_commuter_age; ++age) {
                working_population += populations.get_group_total(age) * commuting_weights[static_cast<size_t>(age)];
            }

            // Compute commuter coefficient from i to j
            const double commuter_coeff_ij = mobility_data_commuter(county_idx_i, county_idx_j) / working_population;

            // Assign mobility coefficients for commuting compartments at work location
            for (auto age = min_commuter_age; age <= max_commuter_age; ++age) {
                for (auto compartment : migrating_compartments) {
                    const size_t coeff_index = populations.get_flat_index({age, compartment});
                    mobility_coeffs[static_cast<size_t>(ContactLocation::Work)].get_baseline()[coeff_index] =
                        commuter_coeff_ij * commuting_weights[static_cast<size_t>(age)];
                }
            }

            // Add edge if commuter coefficient exceeds threshold
            if (commuter_coeff_ij > threshold_edges) {
                params_graph.add_edge(county_idx_i, county_idx_j, std::move(mobility_coeffs),
                                      travel_times(county_idx_i, county_idx_j),
                                      path_mobility[county_idx_i][county_idx_j]);
            }
        }
    }

    return mio::success();
}

/**
 * @brief Creates an input graph for a parameter study by reading data from files.
 * @param[in] start_date Start date of the simulation.
 * @param[in] end_date End date of the simulation.
 * @param[in] num_days Number of days to simulate.
 * @param[in] data_dir Directory containing input data files.
 * @param[in] masks If true, applies mask-related transmission reductions.
 * @param[in] ffp2 If true, uses FFP2 mask factors; otherwise, uses surgical masks (requires masks=true).
 * @param[in] edges If true, configures graph edges for mobility.
 * @return IOResult containing the created graph or an IO error if file reading fails.
 */
mio::IOResult<mio::ExtendedGraph<mio::osecirts::Model<double>>>
get_graph(mio::Date start_date, const mio::Date end_date, const int num_days, const std::string& data_dir,
          const bool masks, const bool ffp2, const std::vector<std::vector<double>> immunity_population,
          const bool edges)
{
    // file paths
    const std::string travel_times_dir = mio::path_join(data_dir, "mobility", "travel_times_pathes.txt");
    const std::string durations_dir    = mio::path_join(data_dir, "mobility", "activity_duration_work.txt");
    const std::string mobility_data    = mio::path_join(data_dir, "mobility", "commuter_migration_with_locals.txt");
    const std::string travel_paths     = mio::path_join(data_dir, "mobility", "wegketten_ohne_komma.txt");
    const std::string population_file = mio::path_join(data_dir, "pydata", "Germany", "county_current_population.json");

    // Initialize global parameters
    constexpr int num_age_groups = 6;
    mio::osecirts::Parameters<double> params(num_age_groups);
    params.get<mio::osecirts::StartDay>() = mio::get_day_in_year(start_date);

    // Configure parameters and contact matrices
    BOOST_OUTCOME_TRY(set_covid_parameters(params));
    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, params));

    // Create graph
    mio::ExtendedGraph<mio::osecirts::Model<double>> params_graph;
    const std::vector<double> scaling_factor_infected(num_age_groups, 1.7);
    constexpr double scaling_factor_icu  = 1.0;
    constexpr double tnt_capacity_factor = 1.43 / 100000.0;

    BOOST_OUTCOME_TRY(set_nodes(params, start_date, end_date, data_dir, population_file, durations_dir, true,
                                params_graph, mio::osecirts::read_input_data_county<mio::osecirts::Model<double>>,
                                mio::get_node_ids, scaling_factor_infected, scaling_factor_icu, tnt_capacity_factor,
                                immunity_population, num_days, false, true, masks, ffp2));

    // Set transport contact matrices for mobility models in all nodes
    for (auto& node : params_graph.nodes()) {
        BOOST_OUTCOME_TRY(set_contact_matrices_transport(data_dir, node.property.mobility_sim.parameters));
    }

    // Add edges if enabled
    if (edges) {
        const auto migrating_compartments = {mio::osecirts::InfectionState::SusceptibleNaive,
                                             mio::osecirts::InfectionState::ExposedNaive,
                                             mio::osecirts::InfectionState::InfectedNoSymptomsNaive,
                                             mio::osecirts::InfectionState::InfectedSymptomsNaive,
                                             mio::osecirts::InfectionState::SusceptibleImprovedImmunity,
                                             mio::osecirts::InfectionState::SusceptiblePartialImmunity,
                                             mio::osecirts::InfectionState::ExposedPartialImmunity,
                                             mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunity,
                                             mio::osecirts::InfectionState::InfectedSymptomsPartialImmunity,
                                             mio::osecirts::InfectionState::ExposedImprovedImmunity,
                                             mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunity,
                                             mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunity};

        const std::vector<double> commuting_weights = {0.0, 0.0, 1.0, 1.0, 0.33, 0.0};
        BOOST_OUTCOME_TRY(set_edges(travel_times_dir, mobility_data, travel_paths, params_graph, migrating_compartments,
                                    contact_locations.size(), commuting_weights));
        BOOST_OUTCOME_TRY(scale_contacts_local(params_graph, data_dir, mobility_data, commuting_weights));
    }

    return params_graph;
}

/**
 * @brief Runs a parameter study simulation and saves results.
 * Configures a graph, runs multiple simulations, and exports time series and flow data.
 * @param[in] data_dir Directory containing input data files.
 * @param[in] res_dir Directory to save simulation results.
 * @param[in] num_runs Number of simulation runs for the parameter study.
 * @param[in] num_days Number of days to simulate.
 * @param[in] masks Whether masks are enabled.
 * @param[in] ffp2 Whether FFP2 is enabled.
 * @param[in] edges Whether edges should be used.
 * @return IOResult<void> indicating success or an error during execution.
 */
mio::IOResult<void> run(const std::string& data_dir, std::string res_dir, int num_runs, int num_days, bool masks,
                        bool ffp2, bool edges)
{
    // Uncomment to suppress logs
    // mio::set_log_level(mio::LogLevel::critical);
    const mio::Date start_date(2022, 8, 1);
    const mio::Date end_date(2022, 11, 1);

    if (!masks && ffp2) {
        return mio::failure(mio::StatusCode::InvalidValue, "FFP2 only possible with masks enabled");
    }

    const auto immunity_population =
        std::vector<std::vector<double>>{{0.04, 0.61, 0.35}, {0.04, 0.61, 0.35},   {0.075, 0.62, 0.305},
                                         {0.08, 0.62, 0.3},  {0.035, 0.58, 0.385}, {0.01, 0.41, 0.58}};

    BOOST_OUTCOME_TRY(auto params_graph,
                      get_graph(start_date, end_date, num_days, data_dir, masks, ffp2, immunity_population, edges));

    // Construct result directory name
    res_dir += "/masks_" + std::to_string(masks) + (ffp2 ? "_ffp2" : "") + (!edges ? "_no_edges" : "");
    if (mio::mpi::is_root()) {
        std::cout << "Result directory: " << res_dir << "\n";
    }

    // Ensure result directory exists
    if (!fs::exists(res_dir)) {
        fs::create_directories(res_dir);
    }

    // Extract county IDs from graph nodes
    std::vector<int> county_ids;
    county_ids.reserve(params_graph.nodes().size());
    std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), std::back_inserter(county_ids),
                   [](const auto& node) {
                       return node.id;
                   });

    // Run parameter study
    auto parameter_study = mio::ParameterStudy<
        mio::osecirts::Simulation<double, mio::FlowSimulation<double, mio::osecirts::Model<double>>>,
        mio::ExtendedGraph<mio::osecirts::Model<double>>,
        mio::ExtendedGraph<
            mio::osecirts::Simulation<double, mio::FlowSimulation<double, mio::osecirts::Model<double>>>>>(
        params_graph, 0.0, num_days, 0.01, num_runs);

    if (mio::mpi::is_root()) {
        std::cout << "Seeds: ";
        for (auto seed : parameter_study.get_rng().get_seeds()) {
            std::cout << seed << ", ";
        }
        std::cout << "\n";
    }

    auto ensemble = parameter_study.run(
        [](auto&& graph) {
            return draw_sample(graph);
        },
        [&](auto&& results_graph, size_t run_idx) {
            // get base simulation results
            std::vector<mio::TimeSeries<double>> base_results;
            base_results.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(base_results),
                           [](auto& node) {
                               return mio::interpolate_simulation_result(node.property.base_sim.get_result());
                           });

            std::vector<mio::osecirts::Model<double>> base_params;
            base_params.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(base_params),
                           [](auto& node) {
                               return node.property.base_sim.get_model();
                           });

            std::vector<mio::TimeSeries<double>> base_flows;
            base_flows.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(base_flows),
                           [](auto& node) {
                               return mio::interpolate_simulation_result(node.property.base_sim.get_flows());
                           });

            // get mobility simulation results
            std::vector<mio::osecirts::Model<double>> mobility_params;
            mobility_params.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(),
                           std::back_inserter(mobility_params), [](auto& node) {
                               return node.property.mobility_sim.get_model();
                           });

            std::vector<mio::TimeSeries<double>> mobility_flows;
            mobility_flows.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(),
                           std::back_inserter(mobility_flows), [](auto& node) {
                               return mio::interpolate_simulation_result(node.property.mobility_sim.get_flows());
                           });

            if (mio::mpi::is_root()) {
                std::cout << "Run " << run_idx << " complete.\n";
            }

            return std::make_tuple(base_results, base_params, base_flows, mobility_params, mobility_flows);
        });

    // Save results if ensemble is non-empty
    if (!ensemble.empty()) {
        std::vector<std::vector<mio::TimeSeries<double>>> ensemble_results, ensemble_flows, ensemble_flows_mobility;
        std::vector<std::vector<mio::osecirts::Model<double>>> ensemble_params, ensemble_params_mobility;

        ensemble_results.reserve(ensemble.size());
        ensemble_params.reserve(ensemble.size());
        ensemble_flows.reserve(ensemble.size());
        ensemble_params_mobility.reserve(ensemble.size());
        ensemble_flows_mobility.reserve(ensemble.size());

        for (auto&& run : ensemble) {
            ensemble_results.emplace_back(std::move(std::get<0>(run)));
            ensemble_params.emplace_back(std::move(std::get<1>(run)));
            ensemble_flows.emplace_back(std::move(std::get<2>(run)));
            ensemble_params_mobility.emplace_back(std::move(std::get<3>(run)));
            ensemble_flows_mobility.emplace_back(std::move(std::get<4>(run)));
        }

        // Save base results
        BOOST_OUTCOME_TRY(mio::save_results(ensemble_results, ensemble_params, county_ids, res_dir, false));

        // Save base flows
        const std::string flows_dir = res_dir + "/flows";
        if (mio::mpi::is_root()) {
            fs::create_directories(flows_dir);
            std::cout << "Saving flow results to \"" << flows_dir << "\".\n";
        }
        BOOST_OUTCOME_TRY(mio::save_results(ensemble_flows, ensemble_params, county_ids, flows_dir, false));

        // Save mobility flows
        const std::string mobility_dir       = res_dir + "/mobility";
        const std::string mobility_flows_dir = mobility_dir + "/flows";
        if (mio::mpi::is_root()) {
            fs::create_directories(mobility_flows_dir);
            std::cout << "Saving mobility flow results to \"" << mobility_flows_dir << "\".\n";
        }
        BOOST_OUTCOME_TRY(mio::save_results(ensemble_flows_mobility, ensemble_params_mobility, county_ids,
                                            mobility_flows_dir, false));
    }

    return mio::success();
}

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::warn);
    mio::mpi::init();

    std::string data_dir    = "";
    std::string results_dir = "";
    int num_days            = 10;
    int num_runs            = 2;
    bool masks              = true;
    bool ffp2               = true;
    bool edges              = true;

    if (argc == 8) {
        data_dir    = argv[1];
        results_dir = argv[2];
        num_runs    = std::atoi(argv[3]);
        num_days    = std::atoi(argv[4]);
        masks       = std::atoi(argv[5]) != 0;
        ffp2        = std::atoi(argv[6]) != 0;
        edges       = std::atoi(argv[7]) != 0;

        if (mio::mpi::is_root()) {
            printf("Reading data from \"%s\", saving results to \"%s\".\n", data_dir.c_str(), results_dir.c_str());
            printf("Number of runs: %d, Number of days: %d.\n", num_runs, num_days);
            printf("Masks: %s, FFP2: %s, Edges: %s.\n", masks ? "true" : "false", ffp2 ? "true" : "false",
                   edges ? "true" : "false");
        }
    }
    else if (argc == 5) {
        data_dir    = argv[1];
        results_dir = argv[2];
        num_runs    = std::atoi(argv[3]);
        num_days    = std::atoi(argv[4]);

        if (mio::mpi::is_root()) {
            printf("Reading data from \"%s\", saving results to \"%s\".\n", data_dir.c_str(), results_dir.c_str());
            printf("Using default values - Number of runs: %d, Number of days: %d, Masks: true, FFP2: true, Edges: "
                   "true.\n",
                   num_runs, num_days);
        }
    }
    else {
        if (mio::mpi::is_root()) {
            printf("Usage:\n");
            printf("2022_omicron_late_phase_mobility <data_dir> <results_dir> <num_runs> <num_days> <masks> <ffp2> "
                   "<edges>\n");
            printf("\tRun simulation with data from <data_dir>, saving results to <results_dir>.\n");
            printf("\t<num_runs>: Number of simulation runs.\n");
            printf("\t<num_days>: Number of days to simulate.\n");
            printf("\t<masks>: 0 or 1 to disable or enable masks.\n");
            printf("\t<ffp2>: 0 or 1 to disable or enable FFP2 (requires masks).\n");
            printf("\t<edges>: 0 or 1 to disable or enable edges.\n");
        }
        mio::mpi::finalize();
        return 0;
    }

    auto result = run(data_dir, results_dir, num_runs, num_days, masks, ffp2, edges);
    if (!result) {
        if (mio::mpi::is_root()) {
            printf("%s\n", result.error().formatted_message().c_str());
        }
        mio::mpi::finalize();
        return -1;
    }

    mio::mpi::finalize();
    return 0;
}
