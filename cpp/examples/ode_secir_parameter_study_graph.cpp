/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Henrik Zunker
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
#include "memilio/config.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/result_io.h"
#include "memilio/io/mobility_io.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/time_series.h"
#include "ode_secir/model.h"
#include "ode_secir/parameters_io.h"
#include "ode_secir/parameter_space.h"
#include <cstddef>
#include <cstdio>

/**
 * Set a value and distribution of an UncertainValue.
 * Assigns average of min and max as a value and UNIFORM(min, max) as a distribution.
 * @param p uncertain value to set.
 * @param min minimum of distribution.
 * @param max minimum of distribution.
 */
void assign_uniform_distribution(mio::UncertainValue<ScalarType>& p, ScalarType min, ScalarType max)
{
    p = mio::UncertainValue<ScalarType>(0.5 * (max + min));
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
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue<ScalarType>, mio::AgeGroup>& array,
                                       const ScalarType (&min)[N], const ScalarType (&max)[N])
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
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue<ScalarType>, mio::AgeGroup>& array,
                                       ScalarType min, ScalarType max)
{
    for (auto i = mio::AgeGroup(0); i < array.size<mio::AgeGroup>(); ++i) {
        assign_uniform_distribution(array[i], min, max);
    }
}

/**
 * Set epidemiological parameters of Sars-CoV-2 for a immune-naive
 * population and wild type variant.
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
void set_covid_parameters(mio::osecir::Parameters<ScalarType>& params)
{
    //times
    const ScalarType timeExposedMin            = 2.67;
    const ScalarType timeExposedMax            = 4.;
    const ScalarType timeInfectedNoSymptomsMin = 1.2;
    const ScalarType timeInfectedNoSymptomsMax = 2.53;

    const ScalarType timeInfectedSymptomsMin[] = {5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465};
    const ScalarType timeInfectedSymptomsMax[] = {8.427, 8.427, 8.4684, 8.3139, 8.169, 8.085};
    const ScalarType timeInfectedSevereMin[]   = {3.925, 3.925, 4.85, 6.4, 7.2, 9.};
    const ScalarType timeInfectedSevereMax[]   = {6.075, 6.075, 7., 8.7, 9.8, 13.};
    const ScalarType timeInfectedCriticalMin[] = {4.95, 4.95, 4.86, 14.14, 14.4, 10.};
    const ScalarType timeInfectedCriticalMax[] = {8.95, 8.95, 8.86, 20.58, 19.8, 13.2};

    array_assign_uniform_distribution(params.get<mio::osecir::TimeExposed<ScalarType>>(), timeExposedMin,
                                      timeExposedMax);
    array_assign_uniform_distribution(params.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>(),
                                      timeInfectedNoSymptomsMin, timeInfectedNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>(),
                                      timeInfectedSymptomsMin, timeInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecir::TimeInfectedSevere<ScalarType>>(), timeInfectedSevereMin,
                                      timeInfectedSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecir::TimeInfectedCritical<ScalarType>>(),
                                      timeInfectedCriticalMin, timeInfectedCriticalMax);

    //probabilities
    const ScalarType transmissionProbabilityOnContactMin[] = {0.02, 0.05, 0.05, 0.05, 0.08, 0.15};
    const ScalarType transmissionProbabilityOnContactMax[] = {0.04, 0.07, 0.07, 0.07, 0.10, 0.20};
    const ScalarType relativeTransmissionNoSymptomsMin     = 1;
    const ScalarType relativeTransmissionNoSymptomsMax     = 1;
    // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
    // depends on incidence and test and trace capacity
    const ScalarType riskOfInfectionFromSymptomaticMin    = 0.1;
    const ScalarType riskOfInfectionFromSymptomaticMax    = 0.3;
    const ScalarType maxRiskOfInfectionFromSymptomaticMin = 0.3;
    const ScalarType maxRiskOfInfectionFromSymptomaticMax = 0.5;
    const ScalarType recoveredPerInfectedNoSymptomsMin[]  = {0.2, 0.2, 0.15, 0.15, 0.15, 0.15};
    const ScalarType recoveredPerInfectedNoSymptomsMax[]  = {0.3, 0.3, 0.25, 0.25, 0.25, 0.25};
    const ScalarType severePerInfectedSymptomsMin[]       = {0.006, 0.006, 0.015, 0.049, 0.15, 0.20};
    const ScalarType severePerInfectedSymptomsMax[]       = {0.009, 0.009, 0.023, 0.074, 0.18, 0.25};
    const ScalarType criticalPerSevereMin[]               = {0.05, 0.05, 0.05, 0.10, 0.25, 0.35};
    const ScalarType criticalPerSevereMax[]               = {0.10, 0.10, 0.10, 0.20, 0.35, 0.45};
    const ScalarType deathsPerCriticalMin[]               = {0.00, 0.00, 0.10, 0.10, 0.30, 0.5};
    const ScalarType deathsPerCriticalMax[]               = {0.10, 0.10, 0.18, 0.18, 0.50, 0.7};

    array_assign_uniform_distribution(params.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>(),
                                      transmissionProbabilityOnContactMin, transmissionProbabilityOnContactMax);
    array_assign_uniform_distribution(params.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>(),
                                      relativeTransmissionNoSymptomsMin, relativeTransmissionNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>(),
                                      riskOfInfectionFromSymptomaticMin, riskOfInfectionFromSymptomaticMax);
    array_assign_uniform_distribution(params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<ScalarType>>(),
                                      maxRiskOfInfectionFromSymptomaticMin, maxRiskOfInfectionFromSymptomaticMax);
    array_assign_uniform_distribution(params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>(),
                                      recoveredPerInfectedNoSymptomsMin, recoveredPerInfectedNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>(),
                                      severePerInfectedSymptomsMin, severePerInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecir::CriticalPerSevere<ScalarType>>(), criticalPerSevereMin,
                                      criticalPerSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecir::DeathsPerCritical<ScalarType>>(), deathsPerCriticalMin,
                                      deathsPerCriticalMax);

    //sasonality
    const ScalarType seasonality_min = 0.1;
    const ScalarType seasonality_max = 0.3;

    assign_uniform_distribution(params.get<mio::osecir::Seasonality<ScalarType>>(), seasonality_min, seasonality_max);

    params.set<mio::osecir::StartDay<ScalarType>>(0);
}

/**
 * Set synthetic population data for testing.
 * Same total populaton but different spread of infection in each county.
 * @param counties parameters for each county.
 */
void set_synthetic_population_data(mio::osecir::Model<ScalarType>& model)
{
    ScalarType nb_total_t0 = 10000, nb_exp_t0 = 2, nb_inf_t0 = 0, nb_car_t0 = 0, nb_hosp_t0 = 0, nb_icu_t0 = 0,
               nb_rec_t0 = 0, nb_dead_t0 = 0;

    for (mio::AgeGroup i = 0; i < model.parameters.get_num_groups(); i++) {
        model.populations[{i, mio::osecir::InfectionState::Exposed}]                     = nb_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}]          = nb_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]            = nb_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]              = nb_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]            = nb_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]                   = nb_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]                        = nb_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         nb_total_t0);
    }
}

std::vector<std::vector<size_t>> get_indices_of_symptomatic_and_nonsymptomatic(mio::osecir::Model<ScalarType>& model)
{
    std::vector<std::vector<size_t>> indices_save_edges(2);
    const auto num_groups = static_cast<size_t>(model.parameters.get_num_groups());

    // Reserve Space. The multiplication by 2 is necessary because we have the
    // base and the confirmed compartments for each age group.
    for (auto& vec : indices_save_edges) {
        vec.reserve(2 * num_groups);
    }

    // get indices and write them to the vector
    for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(num_groups); ++i) {
        indices_save_edges[0].emplace_back(
            model.populations.get_flat_index({i, mio::osecir::InfectionState::InfectedNoSymptoms}));
        indices_save_edges[0].emplace_back(
            model.populations.get_flat_index({i, mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}));
        indices_save_edges[1].emplace_back(
            model.populations.get_flat_index({i, mio::osecir::InfectionState::InfectedSymptoms}));
        indices_save_edges[1].emplace_back(
            model.populations.get_flat_index({i, mio::osecir::InfectionState::InfectedSymptomsConfirmed}));
    }
    return indices_save_edges;
}

/**
 * Run the parameter study.
 * Load a previously stored graph or create a new one from data.
 * The graph is the input for the parameter study.
 * A newly created graph is saved and can be reused.
 * @param mode Mode for running the parameter study.
 * @param data_dir data directory. Not used if mode is RunMode::Load.
 * @param save_dir directory where the graph is loaded from if mode is RunMOde::Load or save to if mode is RunMode::Save.
 * @param result_dir directory where all results of the parameter study will be stored.
 * @param save_single_runs [Default: true] Defines if single run results are written to the disk.
 * @returns any io error that occurs during reading or writing of files.
 */
int main()
{
    mio::set_log_level(mio::LogLevel::warn);
    const auto num_days_sim = 30.0;
    const auto num_runs     = 10;

    mio::Graph<mio::osecir::Model<ScalarType>, mio::MobilityParameters<ScalarType>> params_graph;

    const int num_age_groups = 6;
    mio::osecir::Model<ScalarType> model(num_age_groups);

    // set parameters
    set_covid_parameters(model.parameters);

    // set contact matrix
    const auto cont_freq = 10.0;
    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        model.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(
        (size_t)num_age_groups, (size_t)num_age_groups, (1. / num_age_groups) * cont_freq));

    // set population data
    set_synthetic_population_data(model);
    auto indices_save_edges = get_indices_of_symptomatic_and_nonsymptomatic(model);

    params_graph.add_node(1001, model);
    params_graph.add_node(1002, model);
    params_graph.add_node(1003, model);

    params_graph.add_edge(
        0, 1, Eigen::VectorX<ScalarType>::Constant(num_age_groups * (size_t)mio::osecir::InfectionState::Count, 0.05),
        indices_save_edges);
    params_graph.add_edge(
        1, 0, Eigen::VectorX<ScalarType>::Constant(num_age_groups * (size_t)mio::osecir::InfectionState::Count, 0.1),
        indices_save_edges);
    params_graph.add_edge(
        1, 2, Eigen::VectorX<ScalarType>::Constant(num_age_groups * (size_t)mio::osecir::InfectionState::Count, 0.15),
        indices_save_edges);
    params_graph.add_edge(
        2, 1, Eigen::VectorX<ScalarType>::Constant(num_age_groups * (size_t)mio::osecir::InfectionState::Count, 0.2),
        indices_save_edges);

    mio::ParameterStudy parameter_study(params_graph, 0.0, num_days_sim, 0.5, size_t(num_runs));

    if (mio::mpi::is_root()) {
        printf("Seeds: ");
        for (auto s : parameter_study.get_rng().get_seeds()) {
            printf("%u, ", s);
        }
        printf("\n");
    }

    auto save_single_run_result = mio::IOResult<void>(mio::success());
    auto ensemble               = parameter_study.run(
        [](auto&& graph, ScalarType t0, ScalarType dt, size_t) {
            auto copy = graph;
            return mio::make_sampled_graph_simulation<ScalarType, mio::osecir::Simulation<ScalarType>>(
                mio::osecir::draw_sample(copy), t0, dt, dt);
        },
        [&](auto&& results_sim, auto&& run_id) {
            auto results_graph       = results_sim.get_graph();
            auto interpolated_result = mio::interpolate_simulation_result(results_graph);

            auto params = std::vector<mio::osecir::Model<ScalarType>>{};
            params.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(params),
                                         [](auto&& node) {
                               return node.property.get_simulation().get_model();
                           });

            auto edges = std::vector<mio::TimeSeries<ScalarType>>{};
            edges.reserve(results_graph.edges().size());
            std::transform(results_graph.edges().begin(), results_graph.edges().end(), std::back_inserter(edges),
                                         [](auto&& edge) {
                               return edge.property.get_mobility_results();
                           });

            std::cout << "Run " << run_id << " done\n";
            return std::make_tuple(std::move(interpolated_result), std::move(params), std::move(edges));
        });

    if (ensemble.size() > 0) {
        auto ensemble_results = std::vector<std::vector<mio::TimeSeries<ScalarType>>>{};
        ensemble_results.reserve(ensemble.size());
        auto ensemble_params = std::vector<std::vector<mio::osecir::Model<ScalarType>>>{};
        ensemble_params.reserve(ensemble.size());
        auto ensemble_edges = std::vector<std::vector<mio::TimeSeries<ScalarType>>>{};
        ensemble_edges.reserve(ensemble.size());
        for (auto&& run : ensemble) {
            ensemble_results.emplace_back(std::move(std::get<0>(run)));
            ensemble_params.emplace_back(std::move(std::get<1>(run)));
            ensemble_edges.emplace_back(std::move(std::get<2>(run)));
        }
        // create directory for results.
        boost::filesystem::path results_dir("test_results");
        bool created = boost::filesystem::create_directories(results_dir);
        if (created) {
            mio::log_info("Directory '{}' was created.", results_dir.string());
        }

        auto county_ids          = std::vector<int>{1001, 1002, 1003};
        auto save_results_status = save_results(ensemble_results, ensemble_params, county_ids, results_dir, false);
        auto pairs_edges         = std::vector<std::pair<int, int>>{};
        for (auto& edge : params_graph.edges()) {
            pairs_edges.push_back({county_ids[edge.start_node_idx], county_ids[edge.end_node_idx]});
        }
        auto save_edges_status = save_edges(ensemble_edges, pairs_edges, "test_results", false, true);
    }
}
