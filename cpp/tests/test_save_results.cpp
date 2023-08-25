/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

    mio::osecir::Model model(1);
    auto& params            = model.parameters;
    mio::AgeGroup nb_groups = params.get_num_groups();
    ;

    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        params.get<mio::osecir::IncubationTime>()[i]       = 5.2;
        params.get<mio::osecir::TimeInfectedSymptoms>()[i] = 5.;
        params.get<mio::osecir::SerialInterval>()[i]       = 4.2;
        params.get<mio::osecir::TimeInfectedSevere>()[i]   = 10.;
        params.get<mio::osecir::TimeInfectedCritical>()[i] = 8.;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]            = nb_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}] = nb_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]   = nb_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]     = nb_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]   = nb_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]          = nb_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]               = nb_dead_t0;
        model.populations.set_difference_from_total({i, mio::osecir::InfectionState::Susceptible}, nb_total_t0);

        params.get<mio::osecir::TransmissionProbabilityOnContact>()[i] = 0.06;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms>()[i]   = 0.67;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i]   = alpha;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i]   = beta;
        params.get<mio::osecir::SeverePerInfectedSymptoms>()[i]        = rho;
        params.get<mio::osecir::CriticalPerSevere>()[i]                = theta;
        params.get<mio::osecir::DeathsPerCritical>()[i]                = delta;
    }

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns>();
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, cont_freq));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    auto result_from_sim                                  = simulate(t0, tmax, dt, model);
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
    //rng needs to be reseeded right before using parallel parameterstudies
    //to keep the rest of the tests independent, we install a temporary RNG for this test
    auto rng                = mio::thread_local_rng();
    mio::thread_local_rng() = mio::RandomNumberGenerator();
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    // set up parameter study
    double t0           = 0;
    double tmax         = 100;
    double cont_freq    = 10; // see Polymod study
    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    size_t num_groups = 3;
    mio::osecir::Model model((int)num_groups);
    double fact = 1.0 / (double)num_groups;

    auto& params = model.parameters;
    for (auto i = mio::Index<mio::AgeGroup>(0); i.get() < num_groups; i++) {
        params.get<mio::osecir::IncubationTime>()[i]       = 5.2;
        params.get<mio::osecir::TimeInfectedSymptoms>()[i] = 5.;
        params.get<mio::osecir::SerialInterval>()[i]       = 4.2;
        params.get<mio::osecir::TimeInfectedSevere>()[i]   = 10.;
        params.get<mio::osecir::TimeInfectedCritical>()[i] = 8.;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]                     = fact * num_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}]          = fact * num_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]            = fact * num_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]              = fact * num_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]            = fact * num_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]                   = fact * num_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]                        = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * num_total_t0);

        params.get<mio::osecir::TransmissionProbabilityOnContact>()[i] = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms>()[i]   = 0.67;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i]   = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i]   = 0.25;
        params.get<mio::osecir::SeverePerInfectedSymptoms>()[i]        = 0.2;
        params.get<mio::osecir::CriticalPerSevere>()[i]                = 0.25;
        params.get<mio::osecir::DeathsPerCritical>()[i]                = 0.3;
    }

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns>();
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_groups, num_groups, fact * cont_freq));

    mio::osecir::set_params_distributions_normal(model, t0, tmax, 0.);

    auto graph = mio::Graph<mio::osecir::Model, mio::MigrationParameters>();
    graph.add_node(0, model);
    graph.add_node(1, model);
    graph.add_edge(0, 1, mio::MigrationParameters(Eigen::VectorXd::Constant(Eigen::Index(num_groups * 10), 1.0)));

    auto num_runs        = 3;
    auto parameter_study = mio::ParameterStudy<mio::osecir::Simulation<>>(graph, 0.0, 2.0, 0.5, num_runs);

    TempFileRegister tmp_file_register;
    std::string tmp_results_dir = tmp_file_register.get_unique_path();
    ASSERT_THAT(mio::create_directory(tmp_results_dir), IsSuccess());

    auto ensemble_results = std::vector<std::vector<mio::TimeSeries<double>>>{};
    ensemble_results.reserve(size_t(num_runs));
    auto ensemble_params = std::vector<std::vector<mio::osecir::Model>>{};
    ensemble_params.reserve(size_t(num_runs));
    auto save_result_status = mio::IOResult<void>(mio::success());
    parameter_study.run(
        [](auto&& g) {
            return draw_sample(g);
        },
        [&](auto results_graph, auto run_idx) {
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

    auto read_graph = mio::read_graph<mio::osecir::Model>(tmp_results_dir + "/run0", mio::IOF_OmitDistributions, false);

    EXPECT_EQ(
        read_graph.value()
            .nodes()[0]
            .property.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[mio::Index<mio::AgeGroup>(0)],
        params.get<mio::osecir::TransmissionProbabilityOnContact>()[mio::Index<mio::AgeGroup>(0)]);

    EXPECT_EQ(read_graph.value()
                  .nodes()[0]
                  .property.parameters.get<mio::osecir::CriticalPerSevere>()[mio::Index<mio::AgeGroup>(0)],
              params.get<mio::osecir::CriticalPerSevere>()[mio::Index<mio::AgeGroup>(0)]);

    EXPECT_EQ(read_graph.value()
                  .nodes()[0]
                  .property.parameters.get<mio::osecir::SerialInterval>()[mio::Index<mio::AgeGroup>(1)],
              params.get<mio::osecir::SerialInterval>()[mio::Index<mio::AgeGroup>(1)]);

    mio::thread_local_rng() = rng;
}

TEST(TestSaveResult, save_percentiles_and_sums)
{
    //rng needs to be reseeded right before using parallel parameterstudies
    //to keep the rest of the tests independent, we install a temporary RNG for this test
    auto prev_rng           = mio::thread_local_rng();
    mio::thread_local_rng() = mio::RandomNumberGenerator();
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    // set up parameter study
    double t0           = 0;
    double tmax         = 100;
    double cont_freq    = 10; // see Polymod study
    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    size_t num_groups = 3;
    mio::osecir::Model model((int)num_groups);
    double fact = 1.0 / (double)num_groups;

    auto& params = model.parameters;
    for (auto i = mio::Index<mio::AgeGroup>(0); i.get() < num_groups; i++) {
        params.get<mio::osecir::IncubationTime>()[i]       = 5.2;
        params.get<mio::osecir::TimeInfectedSymptoms>()[i] = 5.;
        params.get<mio::osecir::SerialInterval>()[i]       = 4.2;
        params.get<mio::osecir::TimeInfectedSevere>()[i]   = 10.;
        params.get<mio::osecir::TimeInfectedCritical>()[i] = 8.;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]                     = fact * num_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}]          = fact * num_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]            = fact * num_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]              = fact * num_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]            = fact * num_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]                   = fact * num_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]                        = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * num_total_t0);

        params.get<mio::osecir::TransmissionProbabilityOnContact>()[i] = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms>()[i]   = 0.67;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i]   = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i]   = 0.25;
        params.get<mio::osecir::SeverePerInfectedSymptoms>()[i]        = 0.2;
        params.get<mio::osecir::CriticalPerSevere>()[i]                = 0.25;
        params.get<mio::osecir::DeathsPerCritical>()[i]                = 0.3;
    }

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns>();
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_groups, num_groups, fact * cont_freq));

    mio::osecir::set_params_distributions_normal(model, t0, tmax, 0.2);

    auto graph = mio::Graph<mio::osecir::Model, mio::MigrationParameters>();
    graph.add_node(0, model);
    graph.add_node(1, model);
    graph.add_edge(0, 1, mio::MigrationParameters(Eigen::VectorXd::Constant(Eigen::Index(num_groups * 10), 1.0)));

    auto num_runs        = 3;
    auto parameter_study = mio::ParameterStudy<mio::osecir::Simulation<>>(graph, 0.0, 2.0, 0.5, num_runs);

    TempFileRegister tmp_file_register;
    std::string tmp_results_dir = tmp_file_register.get_unique_path();
    ASSERT_THAT(mio::create_directory(tmp_results_dir), IsSuccess());

    auto ensemble_results = std::vector<std::vector<mio::TimeSeries<double>>>{};
    ensemble_results.reserve(size_t(num_runs));
    auto ensemble_params = std::vector<std::vector<mio::osecir::Model>>{};
    ensemble_params.reserve(size_t(num_runs));
    parameter_study.run(
        [](auto&& g) {
            return draw_sample(g);
        },
        [&](auto results_graph, auto /*run_idx*/) {
            ensemble_results.push_back(mio::interpolate_simulation_result(results_graph));

            ensemble_params.emplace_back();
            ensemble_params.back().reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(),
                           std::back_inserter(ensemble_params.back()), [](auto&& node) {
                               return node.property.get_simulation().get_model();
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

    mio::thread_local_rng() = prev_rng;
}
