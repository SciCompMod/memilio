/*
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors:
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
#include "ode_secirvvs/parameters_io.h"
#include "ode_secirvvs/parameter_space.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/stl_util.h"
#include "boost/filesystem.hpp"
#include <cstddef>
#include <cstdio>
#include <iomanip>

namespace fs = boost::filesystem;

/**
 * indices of contact matrix corresponding to locations where contacts occur.
 */
enum class ContactLocation
{
    Home = 0,
    School,
    Work,
    Other,
    Count,
};

/**
 * Set a value and distribution of an UncertainValue.
 * Assigns average of min and max as a value and UNIFORM(min, max) as a distribution.
 * @param p uncertain value to set.
 * @param min minimum of distribution.
 * @param max minimum of distribution.
 */
void assign_uniform_distribution(mio::UncertainValue<double>& p, double min, double max)
{
    p = mio::UncertainValue(0.5 * (max + min));
    p.set_distribution(mio::ParameterDistributionUniform(min, max));
}

/**
 * Set a value and distribution of an array of UncertainValues.
 * Assigns average of min[i] and max[i] as a value and UNIFORM(min[i], max[i]) as a distribution for
 * each element i of the array.
 * @param array array of UncertainValues to set.
 * @param min minimum of distribution for each element of array.
 * @param max minimum of distribution for each element of array.
 */
template <size_t N>
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue<double>, mio::AgeGroup>& array,
                                       const double (&min)[N], const double (&max)[N])
{
    assert(N == array.numel());
    for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(N); ++i) {
        assign_uniform_distribution(array[i], min[size_t(i)], max[size_t(i)]);
    }
}

/**
 * Set a value and distribution of an array of UncertainValues.
 * Assigns average of min and max as a value and UNIFORM(min, max) as a distribution to every element of the array.
 * @param array array of UncertainValues to set.
 * @param min minimum of distribution.
 * @param max minimum of distribution.
 */
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue<double>, mio::AgeGroup>& array,
                                       double min, double max)
{
    for (auto i = mio::AgeGroup(0); i < array.size<mio::AgeGroup>(); ++i) {
        assign_uniform_distribution(array[i], min, max);
    }
}

/**
 * Set epidemiological parameters of Covid19.
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
mio::IOResult<void> set_covid_parameters(mio::osecirvvs::Parameters<double>& params, bool long_time)
{
    //times
    const double timeExposedMin            = 2.67;
    const double timeExposedMax            = 4.;
    const double timeInfectedNoSymptomsMin = 1.2;
    const double timeInfectedNoSymptomsMax = 2.53;
    const double timeInfectedSymptomsMin[] = {5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465};
    const double timeInfectedSymptomsMax[] = {8.427, 8.427, 8.4684, 8.3139, 8.169, 8.085};
    const double timeInfectedSevereMin[]   = {3.925, 3.925, 4.85, 6.4, 7.2, 9.};
    const double timeInfectedSevereMax[]   = {6.075, 6.075, 7., 8.7, 9.8, 13.};
    const double timeInfectedCriticalMin[] = {4.95, 4.95, 4.86, 14.14, 14.4, 10.};
    const double timeInfectedCriticalMax[] = {8.95, 8.95, 8.86, 20.58, 19.8, 13.2};

    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeExposed<double>>(), timeExposedMin,
                                      timeExposedMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedNoSymptoms<double>>(),
                                      timeInfectedNoSymptomsMin, timeInfectedNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSymptoms<double>>(),
                                      timeInfectedSymptomsMin, timeInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSevere<double>>(), timeInfectedSevereMin,
                                      timeInfectedSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedCritical<double>>(),
                                      timeInfectedCriticalMin, timeInfectedCriticalMax);

    //probabilities
    double fac_variant                                 = 1.4;
    const double transmissionProbabilityOnContactMin[] = {0.02 * fac_variant, 0.05 * fac_variant, 0.05 * fac_variant,
                                                          0.05 * fac_variant, 0.08 * fac_variant, 0.1 * fac_variant};

    const double transmissionProbabilityOnContactMax[] = {0.04 * fac_variant, 0.07 * fac_variant, 0.07 * fac_variant,
                                                          0.07 * fac_variant, 0.10 * fac_variant, 0.15 * fac_variant};
    const double relativeTransmissionNoSymptomsMin     = 0.5;
    const double relativeTransmissionNoSymptomsMax     = 0.5;
    // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
    // depends on incidence and test and trace capacity
    const double riskOfInfectionFromSymptomaticMin    = 0.0;
    const double riskOfInfectionFromSymptomaticMax    = 0.2;
    const double maxRiskOfInfectionFromSymptomaticMin = 0.4;
    const double maxRiskOfInfectionFromSymptomaticMax = 0.5;
    const double recoveredPerInfectedNoSymptomsMin[]  = {0.2, 0.2, 0.15, 0.15, 0.15, 0.15};
    const double recoveredPerInfectedNoSymptomsMax[]  = {0.3, 0.3, 0.25, 0.25, 0.25, 0.25};
    const double severePerInfectedSymptomsMin[]       = {0.006, 0.006, 0.015, 0.049, 0.15, 0.20};
    const double severePerInfectedSymptomsMax[]       = {0.009, 0.009, 0.023, 0.074, 0.18, 0.25};
    const double criticalPerSevereMin[]               = {0.05, 0.05, 0.05, 0.10, 0.25, 0.35};
    const double criticalPerSevereMax[]               = {0.10, 0.10, 0.10, 0.20, 0.35, 0.45};
    const double deathsPerCriticalMin[]               = {0.00, 0.00, 0.10, 0.10, 0.30, 0.5};
    const double deathsPerCriticalMax[]               = {0.10, 0.10, 0.18, 0.18, 0.50, 0.7};

    const double reducExposedPartialImmunityMin                     = 0.75;
    const double reducExposedPartialImmunityMax                     = 0.85;
    const double reducExposedImprovedImmunityMin                    = 0.281;
    const double reducExposedImprovedImmunityMax                    = 0.381;
    const double reducInfectedSymptomsPartialImmunityMin            = 0.6;
    const double reducInfectedSymptomsPartialImmunityMax            = 0.7;
    const double reducInfectedSymptomsImprovedImmunityMin           = 0.193;
    const double reducInfectedSymptomsImprovedImmunityMax           = 0.293;
    const double reducInfectedSevereCriticalDeadPartialImmunityMin  = 0.05;
    const double reducInfectedSevereCriticalDeadPartialImmunityMax  = 0.15;
    const double reducInfectedSevereCriticalDeadImprovedImmunityMin = 0.041;
    const double reducInfectedSevereCriticalDeadImprovedImmunityMax = 0.141;

    double temp_reducTimeInfectedMild;
    if (long_time) {
        temp_reducTimeInfectedMild = 1.0;
    }
    else {
        temp_reducTimeInfectedMild = 0.5;
    }
    const double reducTimeInfectedMild = temp_reducTimeInfectedMild;

    array_assign_uniform_distribution(params.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>(),
                                      transmissionProbabilityOnContactMin, transmissionProbabilityOnContactMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RelativeTransmissionNoSymptoms<double>>(),
                                      relativeTransmissionNoSymptomsMin, relativeTransmissionNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic<double>>(),
                                      riskOfInfectionFromSymptomaticMin, riskOfInfectionFromSymptomaticMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<double>>(),
                                      maxRiskOfInfectionFromSymptomaticMin, maxRiskOfInfectionFromSymptomaticMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<double>>(),
                                      recoveredPerInfectedNoSymptomsMin, recoveredPerInfectedNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::SeverePerInfectedSymptoms<double>>(),
                                      severePerInfectedSymptomsMin, severePerInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::CriticalPerSevere<double>>(), criticalPerSevereMin,
                                      criticalPerSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::DeathsPerCritical<double>>(), deathsPerCriticalMin,
                                      deathsPerCriticalMax);

    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducExposedPartialImmunity<double>>(),
                                      reducExposedPartialImmunityMin, reducExposedPartialImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducExposedImprovedImmunity<double>>(),
                                      reducExposedImprovedImmunityMin, reducExposedImprovedImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<double>>(),
                                      reducInfectedSymptomsPartialImmunityMin, reducInfectedSymptomsPartialImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<double>>(),
                                      reducInfectedSymptomsImprovedImmunityMin,
                                      reducInfectedSymptomsImprovedImmunityMax);
    array_assign_uniform_distribution(
        params.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<double>>(),
        reducInfectedSevereCriticalDeadPartialImmunityMin, reducInfectedSevereCriticalDeadPartialImmunityMax);
    array_assign_uniform_distribution(
        params.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>(),
        reducInfectedSevereCriticalDeadImprovedImmunityMin, reducInfectedSevereCriticalDeadImprovedImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducTimeInfectedMild<double>>(),
                                      reducTimeInfectedMild, reducTimeInfectedMild);

    //sasonality
    const double seasonality_min = 0.1;
    const double seasonality_max = 0.3;

    assign_uniform_distribution(params.get<mio::osecirvvs::Seasonality<double>>(), seasonality_min, seasonality_max);

    return mio::success();
}

static const std::map<ContactLocation, std::string> contact_locations = {{ContactLocation::Home, "home"},
                                                                         {ContactLocation::School, "school_pf_eig"},
                                                                         {ContactLocation::Work, "work"},
                                                                         {ContactLocation::Other, "other"}};

/**
 * Set contact matrices.
 * Reads contact matrices from files in the data directory.
 * @param data_dir data directory.
 * @param params Object that the contact matrices will be added to.
 * @returns any io errors that happen during reading of the files.
 */
mio::IOResult<void> set_contact_matrices(const fs::path& data_dir, mio::osecirvvs::Parameters<double>& params)
{
    //TODO: io error handling
    auto contact_matrices = mio::ContactMatrixGroup(contact_locations.size(), size_t(params.get_num_groups()));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(auto&& baseline,
                          mio::read_mobility_plain(
                              (data_dir / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));

        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline;
        contact_matrices[size_t(contact_location.first)].get_minimum()  = Eigen::MatrixXd::Zero(6, 6);
    }
    params.get<mio::osecirvvs::ContactPatterns<double>>() = mio::UncertainContactMatrix<double>(contact_matrices);

    return mio::success();
}

/**
 * @brief Get the initial number deaths from the RKI data.
 * 
 * @param data_dir data directory used to get the cases data.
 * @param start_date start date of the simulation.
 * @param vregion vector of regions to get the initial deaths for.
 * @return double initial number of deaths stratified by age group and region.
 */
std::vector<std::vector<double>> get_initial_deaths(const fs::path& data_dir, mio::Date start_date,
                                                    std::vector<int> const& vregion)
{
    // read simple casedata file to obtain the inital number of deaths
    auto path_cases      = mio::path_join(data_dir.string(), "pydata", "Germany", "cases_all_county_age_ma7.json");
    auto rki_data_status = mio::read_confirmed_cases_data(path_cases);
    auto rki_data        = rki_data_status.value();

    std::vector<std::vector<double>> initial_deaths(400, std::vector<double>(6, 0.0));
    for (auto& entry : rki_data) {
        if (entry.date == start_date) {
            auto it         = std::find_if(vregion.begin(), vregion.end(), [&entry](auto r) {
                return r == 0 || entry.county_id->get() == r;
            });
            auto region_idx = size_t(it - vregion.begin());
            initial_deaths[region_idx][static_cast<size_t>(entry.age_group)] += entry.num_deaths;
        }
    }

    double sum = 0.0;
    std::for_each(initial_deaths.begin(), initial_deaths.end(), [&](const std::vector<double>& v) {
        sum += std::accumulate(v.begin(), v.end(), 0.0);
    });
    std::cout << "Initial deaths (sum): " << sum << std::endl;

    return initial_deaths;
}

/**
 * Create the input graph for the parameter study.
 * Reads files from the data directory.
 * @param start_date start date of the simulation.
 * @param end_date end date of the simulation.
 * @param data_dir data directory.
 * @returns created graph or any io errors that happen during reading of the files.
 */
mio::IOResult<mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>>>
get_graph(mio::Date start_date, const int num_days, const fs::path& data_dir)
{
    // global parameters
    const int num_groups = 6;
    const bool long_time = true;
    auto end_date        = mio::offset_date_by_days(start_date, num_days);
    mio::osecirvvs::Parameters params(num_groups);
    params.get<mio::osecirvvs::StartDay>() = mio::get_day_in_year(start_date);
    BOOST_OUTCOME_TRY(set_covid_parameters(params, long_time));
    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, params));

    // graph of counties with populations and local parameters
    // and mobility between counties
    mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>> params_graph;
    const auto& read_function_nodes = mio::osecirvvs::read_input_data_county<mio::osecirvvs::Model<double>>;
    const auto& read_function_edges = mio::read_mobility_plain;
    const auto& node_id_function    = mio::get_node_ids;

    auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 1.0);
    auto scaling_factor_icu      = 1.0;
    auto mobile_compartments     = {mio::osecirvvs::InfectionState::SusceptibleNaive,
                                    mio::osecirvvs::InfectionState::ExposedNaive,
                                    mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive,
                                    mio::osecirvvs::InfectionState::InfectedSymptomsNaive,
                                    mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity,
                                    mio::osecirvvs::InfectionState::SusceptiblePartialImmunity,
                                    mio::osecirvvs::InfectionState::ExposedPartialImmunity,
                                    mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity,
                                    mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity,
                                    mio::osecirvvs::InfectionState::ExposedImprovedImmunity,
                                    mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity,
                                    mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity};
    auto tnt_capacity_factor     = 0.;

    const auto& set_node_function =
        mio::set_nodes<mio::osecirvvs::TestAndTraceCapacity<double>, mio::osecirvvs::ContactPatterns<double>,
                       mio::osecirvvs::Model<double>, mio::MobilityParameters<double>,
                       mio::osecirvvs::Parameters<double>, decltype(read_function_nodes), decltype(node_id_function)>;
    const auto& set_edge_function =
        mio::set_edges<ContactLocation, mio::osecirvvs::Model<double>, mio::MobilityParameters<double>,
                       mio::MobilityCoefficientGroup, mio::osecirvvs::InfectionState, decltype(read_function_edges)>;

    auto population_data_path =
        mio::path_join((data_dir / "pydata" / "Germany").string(), "county_current_population.json");

    BOOST_OUTCOME_TRY(set_node_function(params, start_date, end_date, data_dir, population_data_path, true,
                                        params_graph, read_function_nodes, node_id_function, scaling_factor_infected,
                                        scaling_factor_icu, tnt_capacity_factor, num_days, false, true));
    BOOST_OUTCOME_TRY(set_edge_function(data_dir, params_graph, mobile_compartments, contact_locations.size(),
                                        read_function_edges, std::vector<ScalarType>{0., 0., 1.0, 1.0, 0.33, 0., 0.},
                                        std::vector<std::vector<size_t>>{}));

    // add initial deaths to the DeadNaive compartment
    BOOST_OUTCOME_TRY(auto&& node_ids, mio::get_node_ids(population_data_path, true));
    auto initial_deaths_cases = get_initial_deaths(data_dir, start_date, node_ids);

    for (size_t node_indx = 0; node_indx < node_ids.size(); ++node_indx) {
        for (auto age = 0; age < num_groups; ++age) {
            params_graph.nodes()[node_indx]
                .property.populations[{(mio::AgeGroup)age, mio::osecirvvs::InfectionState::DeadNaive}] =
                initial_deaths_cases[node_indx][age];
        }
    }
    return mio::success(params_graph);
}

mio::IOResult<void> run(const int num_days_sim, mio::Date start_date, const std::string& data_dir,
                        const std::string& result_dir, bool save_single_runs, const int num_runs,
                        bool save_non_aggregated_results)
{
    //create or load graph
    mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>> params_graph;
    BOOST_OUTCOME_TRY(auto&& created, get_graph(start_date, num_days_sim, data_dir));
    params_graph = created;

    std::vector<int> county_ids(params_graph.nodes().size());
    std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });

    //run parameter study
    auto parameter_study = mio::ParameterStudy<mio::osecirvvs::Simulation<>>{
        params_graph, 0.0, static_cast<double>(num_days_sim), 0.5, static_cast<size_t>(num_runs)};

    if (mio::mpi::is_root()) {
        printf("Seeds: ");
        for (auto s : parameter_study.get_rng().get_seeds()) {
            printf("%u, ", s);
        }
        printf("\n");
    }

    auto save_single_run_result = mio::IOResult<void>(mio::success());
    auto ensemble               = parameter_study.run(
        [&](auto&& graph) {
            return draw_sample(graph, false);
        },
        [&](auto results_graph, auto&& run_idx) {
            auto interpolated_result = mio::interpolate_simulation_result(results_graph);
            auto params              = std::vector<mio::osecirvvs::Model<double>>();
            params.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(params),
                                         [](auto&& node) {
                               return node.property.get_simulation().get_model();
                           });

            if (save_single_run_result && save_single_runs) {
                save_single_run_result =
                    save_result_with_params(interpolated_result, params, county_ids, result_dir, run_idx);
            }
            std::cout << "run " << run_idx << " complete." << std::endl;
            return std::make_pair(interpolated_result, params);
        });

    if (ensemble.size() > 0) {
        auto ensemble_results = std::vector<std::vector<mio::TimeSeries<double>>>{};
        ensemble_results.reserve(ensemble.size());
        auto ensemble_params = std::vector<std::vector<mio::osecirvvs::Model<double>>>{};
        ensemble_params.reserve(ensemble.size());
        for (auto&& run : ensemble) {
            ensemble_results.emplace_back(std::move(run.first));
            ensemble_params.emplace_back(std::move(run.second));
        }

        bool save_percentiles = true;

        // BOOST_OUTCOME_TRY(save_single_run_result);
        BOOST_OUTCOME_TRY(save_results(ensemble_results, ensemble_params, county_ids, result_dir, save_single_runs,
                                       save_percentiles, size_t(num_days_sim), save_non_aggregated_results));
    }

    return mio::success();
}

int main(int argc, char** argv)
{

    mio::set_log_level(mio::LogLevel::warn);
    mio::mpi::init();

    std::string data_dir;

    mio::Date start_date             = mio::Date(2022, 6, 1);
    int num_days_sim                 = 5;
    int num_simulation_runs          = 5;
    bool save_non_aggregated_results = true;
    std::cout << argc << std::endl;

    if (argc == 1) {
        data_dir = "../../data";
    }
    else if (argc == 6) {
        data_dir     = argv[1];
        start_date   = mio::Date(std::atoi(argv[2]), std::atoi(argv[3]), std::atoi(argv[4]));
        num_days_sim = std::atoi(argv[5]);
    }
    else if (argc == 7) {
        data_dir            = argv[1];
        start_date          = mio::Date(std::atoi(argv[2]), std::atoi(argv[3]), std::atoi(argv[4]));
        num_days_sim        = std::atoi(argv[5]);
        num_simulation_runs = std::atoi(argv[6]);
    }
    else if (argc == 8) {
        data_dir                    = argv[1];
        start_date                  = mio::Date(std::atoi(argv[2]), std::atoi(argv[3]), std::atoi(argv[4]));
        num_days_sim                = std::atoi(argv[5]);
        num_simulation_runs         = std::atoi(argv[6]);
        save_non_aggregated_results = bool(argv[7]);
    }
    else {
        mio::mpi::finalize();
        return 0;
    }

    //mio::thread_local_rng().seed(
    //    {114381446, 2427727386, 806223567, 832414962, 4121923627, 1581162203}); //set seeds, e.g., for debugging
    mio::thread_local_rng().synchronize();
    if (mio::mpi::is_root()) {
        printf("Seeds: ");
        for (auto s : mio::thread_local_rng().get_seeds()) {
            printf("%u, ", s);
        }
        printf("\n");
    }

    std::string result_dir = "";

    result_dir = data_dir;
    boost::filesystem::path res_dir(result_dir);
    bool created_results = boost::filesystem::create_directories(res_dir);
    if (created_results) {
        mio::log_info("Directory '{:s}' was created.", res_dir.string());
    }
    printf("Saving results to \"%s\".\n", result_dir.c_str());

    auto result =
        run(num_days_sim, start_date, data_dir, result_dir, false, num_simulation_runs, save_non_aggregated_results);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        mio::mpi::finalize();
        return -1;
    }
    mio::mpi::finalize();
    return 0;
}
