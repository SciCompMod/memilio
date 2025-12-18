/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Wadim Koslow, Henrik Zunker
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
#include "load_test_data.h"
#include "matchers.h"
#include "temp_file_register.h"
#include "memilio/compartments/parameter_studies.h"
#include "memilio/compartments/simulation.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/time_series.h"
#include "ode_secir/parameter_space.h"
#include "ode_secir/parameters_io.h"
#include "utils.h"
#include <gtest/gtest.h>
#include <string>

TEST(TestSaveResult, save_result)
{
    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double cont_freq = 10, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model<double> model(1);
    auto& params            = model.parameters;
    mio::AgeGroup nb_groups = params.get_num_groups();
    ;

    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        params.get<mio::osecir::TimeExposed<double>>()[i]            = 3.2;
        params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i] = 2;
        params.get<mio::osecir::TimeInfectedSymptoms<double>>()[i]   = 5.;
        params.get<mio::osecir::TimeInfectedSevere<double>>()[i]     = 10.;
        params.get<mio::osecir::TimeInfectedCritical<double>>()[i]   = 8.;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]            = nb_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}] = nb_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]   = nb_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]     = nb_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]   = nb_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]          = nb_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]               = nb_dead_t0;
        model.populations.set_difference_from_total({i, mio::osecir::InfectionState::Susceptible}, nb_total_t0);

        params.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i] = 0.06;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[i]   = 0.67;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[i]   = alpha;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i]   = beta;
        params.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i]        = rho;
        params.get<mio::osecir::CriticalPerSevere<double>>()[i]                = theta;
        params.get<mio::osecir::DeathsPerCritical<double>>()[i]                = delta;
    }

    mio::ContactMatrixGroup<double>& contact_matrix = params.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0] =
        mio::ContactMatrix<double>(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, cont_freq));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime<double>(30.));

    auto result_from_sim                                  = mio::simulate<double>(t0, tmax, dt, model);
    std::vector<mio::TimeSeries<double>> results_from_sim = {result_from_sim, result_from_sim};
    std::vector<int> ids                                  = {1, 2};

    TempFileRegister file_register;
    auto results_file_path  = file_register.get_unique_path("test_result-%%%%-%%%%.h5");
    auto save_result_status = mio::save_result(results_from_sim, ids, (int)(size_t)nb_groups, results_file_path);
    ASSERT_TRUE(save_result_status);

    auto results_from_file = mio::read_result(results_file_path);
    ASSERT_TRUE(results_from_file);
    auto result_from_file = results_from_file.value()[0];

    ASSERT_EQ(result_from_file.get_groups().get_num_time_points(), result_from_sim.get_num_time_points());
    ASSERT_EQ(result_from_file.get_totals().get_num_time_points(), result_from_sim.get_num_time_points());
    for (Eigen::Index i = 0; i < result_from_sim.get_num_time_points(); i++) {
        ASSERT_EQ(result_from_file.get_groups().get_num_elements(), result_from_sim.get_num_elements())
            << "at row " << i;
        ASSERT_EQ(result_from_file.get_totals().get_num_elements(),
                  result_from_sim.get_num_elements() / static_cast<Eigen::Index>((size_t)nb_groups))
            << "at row " << i;
        ASSERT_NEAR(result_from_sim.get_time(i), result_from_file.get_groups().get_time(i), 1e-10) << "at row " << i;
        ASSERT_NEAR(result_from_sim.get_time(i), result_from_file.get_totals().get_time(i), 1e-10) << "at row " << i;
        for (Eigen::Index l = 0; l < result_from_file.get_totals().get_num_elements(); l++) {
            double total = 0.0;
            for (Eigen::Index j = 0; j < Eigen::Index((size_t)nb_groups); j++) {
                total += result_from_sim[i][j * (size_t)mio::osecir::InfectionState::Count + l];
                EXPECT_NEAR(result_from_file.get_groups()[i][j * (size_t)mio::osecir::InfectionState::Count + l],
                            result_from_sim[i][j * (size_t)mio::osecir::InfectionState::Count + l], 1e-10)
                    << " at row " << i << " at row " << l << " at Group " << j;
            }
            EXPECT_NEAR(result_from_file.get_totals()[i][l], total, 1e-10) << " at row " << i << " at row " << l;
        }
    }
}

TEST(TestSaveResult, save_result_with_params)
{
    mio::set_log_level(mio::LogLevel::err); // suppress warnings because this is inevatible here
    // set up parameter study
    double t0           = 0;
    double tmax         = 100;
    double cont_freq    = 10; // see Polymod study
    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    size_t num_groups = 3;
    mio::osecir::Model<double> model((int)num_groups);
    double fact = 1.0 / (double)num_groups;

    auto& params = model.parameters;
    for (auto i = mio::Index<mio::AgeGroup>(0); i.get() < (size_t)num_groups; i++) {
        params.get<mio::osecir::TimeExposed<double>>()[i]            = 2.6;
        params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i] = 2.6;
        params.get<mio::osecir::TimeInfectedSymptoms<double>>()[i]   = 5.;
        params.get<mio::osecir::TimeInfectedSevere<double>>()[i]     = 10.;
        params.get<mio::osecir::TimeInfectedCritical<double>>()[i]   = 8.;

        params.get<mio::osecir::Seasonality<double>>()          = 0.0;
        params.get<mio::osecir::ICUCapacity<double>>()          = 0.09;
        params.get<mio::osecir::TestAndTraceCapacity<double>>() = 0.09;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]            = fact * num_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}] = fact * num_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]   = fact * num_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]     = fact * num_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]   = fact * num_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]          = fact * num_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]               = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * num_total_t0);

        params.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i] = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[i]   = 0.67;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[i]   = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i]   = 0.25;
        params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()   = 0.85;
        params.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i]        = 0.2;
        params.get<mio::osecir::CriticalPerSevere<double>>()[i]                = 0.25;
        params.get<mio::osecir::DeathsPerCritical<double>>()[i]                = 0.3;
    }

    mio::ContactMatrixGroup<double>& contact_matrix = params.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0] = mio::ContactMatrix<double>(Eigen::MatrixXd::Constant(num_groups, num_groups, fact * cont_freq));

    mio::osecir::set_params_distributions_normal(model, t0, tmax, 0.);

    auto graph = mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>>();
    graph.add_node(0, model);
    graph.add_node(1, model);
    graph.add_edge(0, 1,
                   mio::MobilityParameters<double>(Eigen::VectorXd::Constant(Eigen::Index(num_groups * 10), 1.0)));

    auto num_runs = 3;
    mio::ParameterStudy parameter_study(graph, 0.0, 2.0, 0.5, num_runs);
    mio::log_rng_seeds(parameter_study.get_rng(), mio::LogLevel::warn);

    TempFileRegister tmp_file_register;
    std::string tmp_results_dir = tmp_file_register.get_unique_path();
    ASSERT_THAT(mio::create_directory(tmp_results_dir), IsSuccess());

    auto ensemble_results = std::vector<std::vector<mio::TimeSeries<double>>>{};
    ensemble_results.reserve(size_t(num_runs));
    auto ensemble_params = std::vector<std::vector<mio::osecir::Model<double>>>{};
    ensemble_params.reserve(size_t(num_runs));
    auto save_result_status = mio::IOResult<void>(mio::success());
    parameter_study.run(
        [](auto&& g, auto t0_, auto dt_, auto) {
            auto copy = g;
            return mio::make_sampled_graph_simulation<double, mio::osecir::Simulation<double>>(draw_sample(copy), t0_,
                                                                                               dt_, dt_);
        },
        [&](auto&& results, auto run_idx) {
            auto results_graph = results.get_graph();
            ensemble_results.push_back(mio::interpolate_simulation_result(results_graph));

            ensemble_params.emplace_back();
            ensemble_params.back().reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(),
                           std::back_inserter(ensemble_params.back()), [](auto&& node) {
                               return node.property.get_simulation().get_model();
                           });

            save_result_status = save_result_with_params(ensemble_results.back(), ensemble_params.back(), {0, 1},
                                                         tmp_results_dir, run_idx);

            return 0; //function needs to return something
        });
    ASSERT_TRUE(save_result_status);
    auto results_from_file = mio::read_result(tmp_results_dir + "/run0/Result.h5");
    ASSERT_TRUE(results_from_file);
    auto result_from_file = results_from_file.value()[0];
    EXPECT_EQ(ensemble_results.back().back().get_num_elements(), result_from_file.get_groups().get_num_elements());
    EXPECT_EQ(ensemble_results.back().back().get_num_time_points(),
              result_from_file.get_groups().get_num_time_points());

    auto read_graph = mio::read_graph<double, mio::osecir::Model<double>>(tmp_results_dir + "/run0",
                                                                          mio::IOF_OmitDistributions, false);

    EXPECT_EQ(read_graph.value()
                  .nodes()[0]
                  .property.parameters
                  .get<mio::osecir::TransmissionProbabilityOnContact<double>>()[mio::Index<mio::AgeGroup>(0)],
              params.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[mio::Index<mio::AgeGroup>(0)]);

    EXPECT_EQ(read_graph.value()
                  .nodes()[0]
                  .property.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[mio::Index<mio::AgeGroup>(1)],
              params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[mio::Index<mio::AgeGroup>(1)]);
}

TEST(TestSaveResult, save_result_order)
{
    std::vector<mio::TimeSeries<double>> results{mio::TimeSeries<double>(0, Eigen::VectorX<double>::Constant(1, 0)),
                                                 mio::TimeSeries<double>(0, Eigen::VectorX<double>::Constant(1, 1)),
                                                 mio::TimeSeries<double>(0, Eigen::VectorX<double>::Constant(1, 2))};

    // case: check order of results, where lexical ordering would rearrange the results;
    // expect: order follows the ids
    std::vector<int> ids = {1, 2, 10};

    TempFileRegister file_register;
    auto results_file_path  = file_register.get_unique_path("test_result-%%%%-%%%%.h5");
    auto save_result_status = mio::save_result(results, ids, 1, results_file_path);
    ASSERT_TRUE(save_result_status);

    auto results_from_file = mio::read_result(results_file_path);
    ASSERT_TRUE(results_from_file);

    ASSERT_DOUBLE_EQ(results_from_file.value()[0].get_groups()[Eigen::Index(0)][Eigen::Index(0)], 0);
    ASSERT_DOUBLE_EQ(results_from_file.value()[1].get_groups()[Eigen::Index(0)][Eigen::Index(0)], 1);
    ASSERT_DOUBLE_EQ(results_from_file.value()[2].get_groups()[Eigen::Index(0)][Eigen::Index(0)], 2);

    // case: check order of results;
    // expect: order is changed due to ids not increasing
    ids = {1, 10, 2};

    results_file_path  = file_register.get_unique_path("test_result-%%%%-%%%%.h5");
    save_result_status = mio::save_result(results, ids, 1, results_file_path);
    ASSERT_TRUE(save_result_status);

    results_from_file = mio::read_result(results_file_path);
    ASSERT_TRUE(results_from_file);

    ASSERT_DOUBLE_EQ(results_from_file.value()[0].get_groups()[Eigen::Index(0)][Eigen::Index(0)], 0);
    ASSERT_DOUBLE_EQ(results_from_file.value()[1].get_groups()[Eigen::Index(0)][Eigen::Index(0)], 2);
    ASSERT_DOUBLE_EQ(results_from_file.value()[2].get_groups()[Eigen::Index(0)][Eigen::Index(0)], 1);
}

TEST(TestSaveResult, save_percentiles_and_sums)
{
    // set up parameter study
    double t0           = 0;
    double tmax         = 100;
    double cont_freq    = 10; // see Polymod study
    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    const size_t num_groups = 3;
    mio::osecir::Model<double> model((int)num_groups);
    double fact = 1.0 / (double)num_groups;

    auto& params = model.parameters;
    for (auto i = mio::Index<mio::AgeGroup>(0); i.get() < (size_t)num_groups; i++) {
        params.get<mio::osecir::TimeExposed<double>>()[i]            = 2.4;
        params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i] = 2.8;
        params.get<mio::osecir::TimeInfectedSymptoms<double>>()[i]   = 5.;
        params.get<mio::osecir::TimeInfectedSevere<double>>()[i]     = 10.;
        params.get<mio::osecir::TimeInfectedCritical<double>>()[i]   = 8.;

        params.get<mio::osecir::Seasonality<double>>()          = 0.0;
        params.get<mio::osecir::ICUCapacity<double>>()          = 100.0;
        params.get<mio::osecir::TestAndTraceCapacity<double>>() = 10.0;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]            = fact * num_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}] = fact * num_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]   = fact * num_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]     = fact * num_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]   = fact * num_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]          = fact * num_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]               = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * num_total_t0);

        params.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i] = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[i]   = 0.67;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[i]   = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i]   = 0.25;
        params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()   = 0.85;
        params.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i]        = 0.2;
        params.get<mio::osecir::CriticalPerSevere<double>>()[i]                = 0.25;
        params.get<mio::osecir::DeathsPerCritical<double>>()[i]                = 0.3;
    }

    mio::ContactMatrixGroup<double>& contact_matrix = params.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0] = mio::ContactMatrix<double>(Eigen::MatrixXd::Constant(num_groups, num_groups, fact * cont_freq));

    // get indices of INS and ISy compartments.
    std::vector<std::vector<size_t>> indices_save_edges(2);

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

    mio::osecir::set_params_distributions_normal(model, t0, tmax, 0.2);

    auto graph = mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>>();
    graph.add_node(0, model);
    graph.add_node(1, model);
    graph.add_edge(0, 1,
                   mio::MobilityParameters<double>(Eigen::VectorXd::Constant(Eigen::Index(num_groups * 10), 1.0),
                                                   indices_save_edges));

    auto num_runs = 3;
    mio::ParameterStudy parameter_study(graph, 0.0, 2.0, 0.5, num_runs);

    TempFileRegister tmp_file_register;
    std::string tmp_results_dir = tmp_file_register.get_unique_path();
    ASSERT_THAT(mio::create_directory(tmp_results_dir), IsSuccess());

    auto ensemble_results = std::vector<std::vector<mio::TimeSeries<double>>>{};
    ensemble_results.reserve(size_t(num_runs));
    auto ensemble_params = std::vector<std::vector<mio::osecir::Model<double>>>{};
    ensemble_params.reserve(size_t(num_runs));
    auto ensemble_edges = std::vector<std::vector<mio::TimeSeries<double>>>{};
    ensemble_edges.reserve(size_t(num_runs));
    parameter_study.run(
        [](auto&& g, auto t0_, auto dt_, auto) {
            mio::LogLevelOverride llo(mio::LogLevel::off);
            auto copy = g;
            return mio::make_sampled_graph_simulation<double, mio::osecir::Simulation<double>>(draw_sample(copy), t0_,
                                                                                               dt_, dt_);
        },
        [&](auto&& results, auto /*run_idx*/) {
            auto results_graph = results.get_graph();
            ensemble_results.push_back(mio::interpolate_simulation_result(results_graph));

            ensemble_params.emplace_back();
            ensemble_params.back().reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(),
                           std::back_inserter(ensemble_params.back()), [](auto&& node) {
                               return node.property.get_simulation().get_model();
                           });

            ensemble_edges.emplace_back();
            ensemble_edges.back().reserve(results_graph.edges().size());
            std::transform(results_graph.edges().begin(), results_graph.edges().end(),
                           std::back_inserter(ensemble_edges.back()), [](auto&& edge) {
                               return edge.property.get_mobility_results();
                           });

            return 0; //function needs to return something
        });

    auto save_results_status = save_results(ensemble_results, ensemble_params, {0, 1}, tmp_results_dir);
    ASSERT_TRUE(save_results_status);

    // test percentiles
    auto results_from_file_p05 = mio::read_result(tmp_results_dir + "/p05/Results.h5");
    ASSERT_TRUE(results_from_file_p05);
    auto results_from_file_p25 = mio::read_result(tmp_results_dir + "/p25/Results.h5");
    ASSERT_TRUE(results_from_file_p25);
    auto results_from_file_p50 = mio::read_result(tmp_results_dir + "/p50/Results.h5");
    ASSERT_TRUE(results_from_file_p50);
    auto results_from_file_p75 = mio::read_result(tmp_results_dir + "/p75/Results.h5");
    ASSERT_TRUE(results_from_file_p75);
    auto results_from_file_p95 = mio::read_result(tmp_results_dir + "/p95/Results.h5");
    ASSERT_TRUE(results_from_file_p95);

    auto result_from_file = results_from_file_p25.value()[0];
    EXPECT_EQ(ensemble_results.back().back().get_num_elements(), result_from_file.get_groups().get_num_elements());
    EXPECT_EQ(ensemble_results.back().back().get_num_time_points(),
              result_from_file.get_groups().get_num_time_points());

    // results_run
    auto results_run0 = mio::read_result(tmp_results_dir + "/results_run0.h5");
    ASSERT_TRUE(results_run0);
    auto results_run0_sum = mio::read_result(tmp_results_dir + "/results_run0_sum.h5");
    ASSERT_TRUE(results_run0_sum);
    auto results_run1 = mio::read_result(tmp_results_dir + "/results_run1.h5");
    ASSERT_TRUE(results_run1);
    auto results_run1_sum = mio::read_result(tmp_results_dir + "/results_run1_sum.h5");
    ASSERT_TRUE(results_run1_sum);
    auto results_run2 = mio::read_result(tmp_results_dir + "/results_run2.h5");
    ASSERT_TRUE(results_run2);
    auto results_run2_sum = mio::read_result(tmp_results_dir + "/results_run2_sum.h5");
    ASSERT_TRUE(results_run2_sum);

    // test save edges (percentiles and results from single runs)
    std::vector<std::pair<int, int>> pairs_edges = {{0, 1}};

    auto save_edges_status = save_edges(ensemble_edges, pairs_edges, tmp_results_dir, true, true);
    ASSERT_TRUE(save_edges_status);

    // percentiles
    auto results_edges_from_file_p05 = mio::read_result(tmp_results_dir + "/p05/Edges.h5");
    ASSERT_TRUE(results_edges_from_file_p05);
    auto results_edges_from_file_p25 = mio::read_result(tmp_results_dir + "/p25/Edges.h5");
    ASSERT_TRUE(results_edges_from_file_p25);
    auto results_edges_from_file_p50 = mio::read_result(tmp_results_dir + "/p50/Edges.h5");
    ASSERT_TRUE(results_edges_from_file_p50);
    auto results_edges_from_file_p75 = mio::read_result(tmp_results_dir + "/p75/Edges.h5");
    ASSERT_TRUE(results_edges_from_file_p75);
    auto results_edges_from_file_p95 = mio::read_result(tmp_results_dir + "/p95/Edges.h5");
    ASSERT_TRUE(results_edges_from_file_p95);

    auto result_edges_from_file = results_edges_from_file_p25.value()[0];
    EXPECT_EQ(ensemble_edges.back().back().get_num_elements(), result_edges_from_file.get_groups().get_num_elements());
    EXPECT_EQ(ensemble_edges.back().back().get_num_time_points(),
              result_edges_from_file.get_groups().get_num_time_points());

    // single runs
    auto results_edges_run0 = mio::read_result(tmp_results_dir + "/Edges_run0.h5");
    ASSERT_TRUE(results_edges_run0);
    auto results_edges_run1 = mio::read_result(tmp_results_dir + "/Edges_run1.h5");
    ASSERT_TRUE(results_edges_run1);
    auto results_edges_run2 = mio::read_result(tmp_results_dir + "/Edges_run2.h5");
    ASSERT_TRUE(results_edges_run2);
}

TEST(TestSaveResult, save_edges)
{
    // create some results and pairs_edges
    const auto n = Eigen::Index(3);
    std::vector<mio::TimeSeries<double>> results_edges(3, mio::TimeSeries<double>(n));
    results_edges[0].add_time_point(0.0, Eigen::VectorXd::Constant(n, 0));
    results_edges[0].add_time_point(1.0, Eigen::VectorXd::Constant(n, 1));
    results_edges[0].add_time_point(2.0, Eigen::VectorXd::Constant(n, 1));

    results_edges[1].add_time_point(0.0, Eigen::VectorXd::Constant(n, 2));
    results_edges[1].add_time_point(1.0, Eigen::VectorXd::Constant(n, 3));
    results_edges[1].add_time_point(2.0, Eigen::VectorXd::Constant(n, 5));

    results_edges[2].add_time_point(0.0, Eigen::VectorXd::Constant(n, 3));
    results_edges[2].add_time_point(1.0, Eigen::VectorXd::Constant(n, 4));
    results_edges[2].add_time_point(2.0, Eigen::VectorXd::Constant(n, 7));

    const std::vector<std::pair<int, int>> pairs_edges = {{0, 1}, {0, 2}, {1, 2}};

    // save the results to a file
    TempFileRegister file_register;
    auto results_file_path = file_register.get_unique_path("test_result-%%%%-%%%%.h5");
    auto save_edges_status = mio::save_edges(results_edges, pairs_edges, results_file_path);
    ASSERT_TRUE(save_edges_status);

    // read the results back in and check that they are correct.
    auto results_from_file = mio::read_result(results_file_path);
    ASSERT_TRUE(results_from_file);

    // group 0
    auto result_from_file_group0 = results_from_file.value()[0];
    EXPECT_EQ(result_from_file_group0.get_groups().get_num_time_points(), 3);
    EXPECT_EQ(result_from_file_group0.get_groups().get_num_elements(), 6);
    EXPECT_EQ(result_from_file_group0.get_groups().get_value(0), (Eigen::VectorXd(6) << 0, 0, 0, 2, 2, 2).finished());
    EXPECT_EQ(result_from_file_group0.get_groups().get_value(1), (Eigen::VectorXd(6) << 1, 1, 1, 3, 3, 3).finished());
    EXPECT_EQ(result_from_file_group0.get_groups().get_value(2), (Eigen::VectorXd(6) << 1, 1, 1, 5, 5, 5).finished());
    EXPECT_EQ(result_from_file_group0.get_groups().get_time(0), 0.0);
    EXPECT_EQ(result_from_file_group0.get_groups().get_time(1), 1.0);
    EXPECT_EQ(result_from_file_group0.get_groups().get_time(2), 2.0);

    EXPECT_EQ(result_from_file_group0.get_totals().get_num_time_points(), 3);
    EXPECT_EQ(result_from_file_group0.get_totals().get_num_elements(), 3);
    EXPECT_EQ(result_from_file_group0.get_totals().get_value(0), Eigen::VectorXd::Constant(3, 2));
    EXPECT_EQ(result_from_file_group0.get_totals().get_value(1), Eigen::VectorXd::Constant(3, 4));
    EXPECT_EQ(result_from_file_group0.get_totals().get_value(2), Eigen::VectorXd::Constant(3, 6));
    EXPECT_EQ(result_from_file_group0.get_totals().get_time(0), 0.0);
    EXPECT_EQ(result_from_file_group0.get_totals().get_time(1), 1.0);
    EXPECT_EQ(result_from_file_group0.get_totals().get_time(2), 2.0);

    // group 1
    auto result_from_file_group1 = results_from_file.value()[1];
    EXPECT_EQ(result_from_file_group1.get_groups().get_num_elements(), 3);
    EXPECT_EQ(result_from_file_group1.get_groups().get_value(0), Eigen::VectorXd::Constant(3, 3));
    EXPECT_EQ(result_from_file_group1.get_groups().get_value(1), Eigen::VectorXd::Constant(3, 4));
    EXPECT_EQ(result_from_file_group1.get_groups().get_value(2), Eigen::VectorXd::Constant(3, 7));
    EXPECT_EQ(result_from_file_group1.get_groups().get_time(0), 0.0);
    EXPECT_EQ(result_from_file_group1.get_groups().get_time(1), 1.0);
    EXPECT_EQ(result_from_file_group1.get_groups().get_time(2), 2.0);

    EXPECT_EQ(result_from_file_group1.get_totals().get_num_elements(), 3);
    EXPECT_EQ(result_from_file_group1.get_totals().get_value(0), Eigen::VectorXd::Constant(3, 3));
    EXPECT_EQ(result_from_file_group1.get_totals().get_value(1), Eigen::VectorXd::Constant(3, 4));
    EXPECT_EQ(result_from_file_group1.get_totals().get_value(2), Eigen::VectorXd::Constant(3, 7));
}

TEST(TestSaveEdges, save_edges_empty_ts)
{
    mio::set_log_level(mio::LogLevel::off);
    std::vector<mio::TimeSeries<double>> results;

    const auto num_elements = 2;

    // Add filled TimeSeries to the results vector
    mio::TimeSeries<double> ts_1(num_elements);
    ts_1.add_time_point(0.0, Eigen::VectorXd::Constant(num_elements, 1));
    ts_1.add_time_point(1.0, Eigen::VectorXd::Constant(num_elements, 2));
    results.push_back(ts_1);

    // Add an empty TimeSeries and add it to the results vector
    mio::TimeSeries<double> ts_2(num_elements);
    results.push_back(ts_2);

    const std::vector<std::pair<int, int>> pairs_edges = {{0, 1}, {1, 2}};

    // Create a TempFile for HDF5 output
    TempFileRegister file_register;
    auto results_file_path = file_register.get_unique_path("TestEdges-%%%%-%%%%.h5");

    // Call the save_edges function and check if it returns a failure
    ASSERT_THAT(save_edges(results, pairs_edges, results_file_path), IsFailure(mio::StatusCode::InvalidValue));
}
