/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Anna Wendler
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
#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"
#include "ode_secirvvs/parameter_space.h"
#include "ode_secirvvs/parameters_io.h"
#include "ode_secirvvs/parameters_io.cpp"
#include <vector>

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
    auto contact_matrices = mio::ContactMatrixGroup<double>(contact_locations.size(), size_t(params.get_num_groups()));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(
            auto&& baseline,
            mio::read_mobility_plain(
                (data_dir / "Germany/contacts" / ("baseline_" + contact_location.second + ".txt")).string()));

        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline;
        contact_matrices[size_t(contact_location.first)].get_minimum()  = Eigen::MatrixXd::Zero(6, 6);
    }
    params.get<mio::osecirvvs::ContactPatterns<double>>() = mio::UncertainContactMatrix<double>(contact_matrices);

    return mio::success();
}

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

template <typename FP, class Tag>
void set_parameters(mio::osecirvvs::Parameters<FP>& params, const std::vector<FP>& parameter_values)
{
    std::copy(parameter_values.begin(), parameter_values.end(), params.template get<Tag>());
}

template <typename FP, class Tag>
void set_parameters(mio::osecirvvs::Parameters<FP>& params, double parameter_values)
{
    params.template get<Tag>() = FP(parameter_values);
}

template <typename FP>
mio::IOResult<void> set_covid_parameters(mio::Date start_date, mio::osecirvvs::Parameters<FP>& params,
                                         const fs::path& temp_dir)
{
    BOOST_OUTCOME_TRY(auto&& parameter_list, mio::read_json((temp_dir / "parameter_list_lha.json").string()));

    std::map<std::string, std::string> id_to_name{};
    for (auto entry : parameter_list) {
        id_to_name[entry["id"].asString()] = entry["name"].asString();
    }

    // only uses group 0 for all parameters
    BOOST_OUTCOME_TRY(auto&& scenario_data_run, mio::read_json((temp_dir / "scenario_data_run_lha.json").string()));
    std::map<std::string, double> parameter_values{};
    for (auto parameter : scenario_data_run["modelParameters"]) {
        std::string parameter_name = id_to_name[parameter["parameterId"].asString()];
        parameter_values[parameter_name] =
            0.5 * (parameter["values"][0]["valueMin"].asDouble() + parameter["values"][0]["valueMax"].asDouble());
    }

    //times
    set_parameters<FP, mio::osecirvvs::TimeExposed<FP>>(params, parameter_values["TimeExposed"]);
    set_parameters<FP, mio::osecirvvs::TimeInfectedNoSymptoms<FP>>(params, parameter_values["TimeInfectedNoSymptoms"]);
    set_parameters<FP, mio::osecirvvs::TimeInfectedSymptoms<FP>>(params, parameter_values["TimeInfectedSymptoms"]);
    set_parameters<FP, mio::osecirvvs::TimeInfectedSevere<FP>>(params, parameter_values["TimeInfectedSevere"]);
    set_parameters<FP, mio::osecirvvs::TimeInfectedCritical<FP>>(params, parameter_values["TimeInfectedCritical"]);

    //probabilities
    set_parameters<FP, mio::osecirvvs::TransmissionProbabilityOnContact<FP>>(
        params, parameter_values["TransmissionProbabilityOnContact"]);
    set_parameters<FP, mio::osecirvvs::RelativeTransmissionNoSymptoms<FP>>(
        params, parameter_values["RelativeTransmissionNoSymptoms"]);
    // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
    // depends on incidence and test and trace capacity
    set_parameters<FP, mio::osecirvvs::RiskOfInfectionFromSymptomatic<FP>>(
        params, parameter_values["RiskOfInfectionFromSymptomatic"]);
    set_parameters<FP, mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<FP>>(
        params, parameter_values["MaxRiskOfInfectionFromSymptomatic"]);
    set_parameters<FP, mio::osecirvvs::RecoveredPerInfectedNoSymptoms<FP>>(
        params, parameter_values["RecoveredPerInfectedNoSymptoms"]);
    set_parameters<FP, mio::osecirvvs::SeverePerInfectedSymptoms<FP>>(params,
                                                                      parameter_values["SeverePerInfectedSymptoms"]);
    set_parameters<FP, mio::osecirvvs::CriticalPerSevere<FP>>(params, parameter_values["CriticalPerSevere"]);
    set_parameters<FP, mio::osecirvvs::DeathsPerCritical<FP>>(params, parameter_values["DeathsPerCritical"]);

    set_parameters<FP, mio::osecirvvs::ReducExposedPartialImmunity<FP>>(
        params, parameter_values["ReducedExposedPartialImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducExposedImprovedImmunity<FP>>(
        params, parameter_values["ReducedExposedImprovedImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<FP>>(
        params, parameter_values["ReducedInfectedSymptomsPartialImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<FP>>(
        params, parameter_values["ReducedInfectedSymptomsImprovedImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<FP>>(
        params, parameter_values["ReducedInfectedSevereCriticalDeadPartialImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>(
        params, parameter_values["ReducedInfectedSevereCriticalDeadImprovedImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducTimeInfectedMild<FP>>(params, parameter_values["ReducedTimeInfectedMild"]);

    params.template get<mio::osecirvvs::StartDay<FP>>() = mio::get_day_in_year(start_date);
    set_parameters<FP, mio::osecirvvs::Seasonality<FP>>(params, parameter_values["Seasonality"]);

    return mio::success();
}

mio::IOResult<mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>>>
get_graph(mio::Date start_date, std::vector<int> lha_ids, const fs::path& data_dir, const fs::path& temp_dir)
{
    // global parameters
    const int num_age_groups = 6;
    mio::osecirvvs::Parameters<double> params(num_age_groups);
    BOOST_OUTCOME_TRY(set_covid_parameters(start_date, params, temp_dir));
    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, params));

    // graph of counties with populations and local parameters
    // and mobility between counties
    mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>> params_graph;
    const auto& read_function_edges = mio::read_mobility_plain;

    auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 1.0);
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

    const auto& set_edge_function =
        mio::set_edges<double, ContactLocation, mio::osecirvvs::Model<double>, mio::MobilityParameters<double>,
                       mio::MobilityCoefficientGroup<double>, mio::osecirvvs::InfectionState,
                       decltype(read_function_edges)>;

    // auto population_data_path =
    //     mio::path_join((data_dir / "Germany" / "pydata").string(), "county_current_population.json");

    BOOST_OUTCOME_TRY(
        mio::osecirvvs::set_lha_data<double>(params, params_graph, data_dir.string(), start_date, lha_ids));
    BOOST_OUTCOME_TRY(
        set_edge_function(mio::path_join((data_dir / "Germany" / "mobility").string(), "commuter_mobility_2022.txt"),
                          params_graph, mobile_compartments, contact_locations.size(), read_function_edges,
                          std::vector<ScalarType>{0., 0., 1.0, 1.0, 0.33, 0., 0.}, std::vector<std::vector<size_t>>{}));

    return mio::success(params_graph);
}

mio::IOResult<void> run(const int num_days_sim, mio::Date start_date, std::vector<int> lha_ids,
                        const std::string& data_dir, const std::string& result_dir, const std::string& temp_dir,
                        bool save_single_runs, const int num_runs, bool save_non_aggregated_results = false)
{
    //create or load graph
    mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>> params_graph;
    BOOST_OUTCOME_TRY(auto&& created, get_graph(start_date, lha_ids, data_dir, temp_dir));
    params_graph = created;

    std::vector<int> county_ids(params_graph.nodes().size());
    std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });
    // run parameter study
    auto parameter_study =
        mio::ParameterStudy(params_graph, 0.0, static_cast<double>(num_days_sim), 0.5, static_cast<size_t>(num_runs));

    if (mio::mpi::is_root()) {
        printf("Seeds: ");
        for (auto s : parameter_study.get_rng().get_seeds()) {
            printf("%u, ", s);
        }
        printf("\n");
    }

    auto save_single_run_result = mio::IOResult<void>(mio::success());
    auto ensemble               = parameter_study.run(
        [](auto&& graph, double t0, double dt, size_t) {
            auto copy = graph;
            return mio::make_sampled_graph_simulation<double, mio::osecirvvs::Simulation<double>>(
                mio::osecirvvs::draw_sample(copy), t0, dt, dt);
        },
        [&](auto&& results_sim, auto&& run_id) {
            auto results_graph = results_sim.get_graph();

            auto interpolated_result = mio::interpolate_simulation_result(results_graph);

            auto params = std::vector<mio::osecirvvs::Model<double>>{};
            params.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(params),
                                         [](auto&& node) {
                               return node.property.get_simulation().get_model();
                           });

            auto edges = std::vector<mio::TimeSeries<double>>{};
            edges.reserve(results_graph.edges().size());
            std::transform(results_graph.edges().begin(), results_graph.edges().end(), std::back_inserter(edges),
                                         [](auto&& edge) {
                               return edge.property.get_mobility_results();
                           });
            if (save_single_run_result && save_single_runs) {
                save_single_run_result =
                    save_result_with_params(interpolated_result, params, county_ids, result_dir, run_id);
            }
            std::cout << "Run " << run_id << " done\n";
            return std::make_tuple(std::move(interpolated_result), std::move(params), std::move(edges));
        });

    if (ensemble.size() > 0) {
        auto ensemble_results = std::vector<std::vector<mio::TimeSeries<ScalarType>>>{};
        ensemble_results.reserve(ensemble.size());
        auto ensemble_params = std::vector<std::vector<mio::osecirvvs::Model<ScalarType>>>{};
        ensemble_params.reserve(ensemble.size());
        auto ensemble_edges = std::vector<std::vector<mio::TimeSeries<ScalarType>>>{};
        ensemble_edges.reserve(ensemble.size());
        for (auto&& run : ensemble) {
            ensemble_results.emplace_back(std::move(std::get<0>(run)));
            ensemble_params.emplace_back(std::move(std::get<1>(run)));
            ensemble_edges.emplace_back(std::move(std::get<2>(run)));
        }

        bool save_percentiles = true;
        BOOST_OUTCOME_TRY(save_results(ensemble_results, ensemble_params, county_ids, result_dir, save_single_runs,
                                       save_percentiles, num_days_sim, save_non_aggregated_results));
    }

    return mio::success();
}

int main(int argc, char** argv)
{
    // This example runs the simulation where one county provides LHA data and the others do not.
    // For now, we assume that there are no active interventions during the simulation period.

    mio::set_log_level(mio::LogLevel::warn);
    mio::mpi::init();

    std::string data_dir   = "../../data";
    std::string result_dir = data_dir;
    std::string temp_dir   = data_dir;
    // County ID.
    std::vector<int> lha_ids = {};
    mio::Date start_date     = mio::Date(2021, 2, 1);
    int num_days_sim         = 30;
    int num_simulation_runs  = 5;

    if (argc == 10) {
        data_dir   = argv[1];
        result_dir = argv[2];
        temp_dir   = argv[3];
        lha_ids.push_back(std::atoi(argv[4]));
        start_date          = mio::Date(std::atoi(argv[5]), std::atoi(argv[6]), std::atoi(argv[7]));
        num_days_sim        = std::atoi(argv[8]);
        num_simulation_runs = std::atoi(argv[9]);
    }

    bool save_non_aggregated_results = false;

    //mio::thread_local_rng().seed(
    //    {114381446, 2427727386, 806223567, 832414962, 4121923627, 1581162203}); //set seeds, e.g., for debugging
    mio::thread_local_rng().synchronize();
    // if (mio::mpi::is_root()) {
    //     printf("Seeds: ");
    //     for (auto s : mio::thread_local_rng().get_seeds()) {
    //         printf("%u, ", s);
    //     }
    //     printf("\n");
    // }

    boost::filesystem::path res_dir(result_dir);
    bool created_results = boost::filesystem::create_directories(res_dir);
    if (created_results) {
        mio::log_info("Directory '{:s}' was created.", res_dir.string());
    }
    printf("Saving results to \"%s\".\n", result_dir.c_str());

    auto result = run(num_days_sim, start_date, lha_ids, data_dir, result_dir, temp_dir, false, num_simulation_runs,
                      save_non_aggregated_results);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        mio::mpi::finalize();
        return -1;
    }
    mio::mpi::finalize();
    return 0;
}