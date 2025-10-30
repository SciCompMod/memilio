/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#include "memilio/config.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/parameter_distributions.h"
#include "ode_secir/model.h"
#include "ode_secir/parameter_space.h"
#include "memilio/compartments/parameter_studies.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/random_number_generator.h"
#include <cstddef>
#include <gtest/gtest.h>
#include <numeric>
#include <stdio.h>

TEST(ParameterStudies, sample_from_secir_params)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    double t0   = 0;
    double tmax = 100;

    double cont_freq = 10; // see Polymod study

    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    mio::osecir::Model<double> model(3);
    mio::AgeGroup num_groups = model.parameters.get_num_groups();
    double fact              = 1.0 / (double)(size_t)num_groups;

    auto& params = model.parameters;
    for (auto i = mio::Index<mio::AgeGroup>(0); i.get() < (size_t)num_groups; i++) {
        params.get<mio::osecir::TimeExposed<double>>()[i]            = 2.6;
        params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i] = 2.6;
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
    contact_matrix[0] =
        mio::ContactMatrix<double>(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, fact * cont_freq));

    mio::osecir::set_params_distributions_normal(model, t0, tmax, 0.2);

    draw_sample(model);

    for (auto i = mio::AgeGroup(0); i < params.get_num_groups(); i++) {
        ASSERT_EQ(params.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[mio::AgeGroup(0)].value(),
                  params.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[i].value());
        ASSERT_EQ(params.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[mio::AgeGroup(0)].value(),
                  params.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i].value());
        ASSERT_EQ(params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[mio::AgeGroup(0)].value(),
                  params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[i].value());

        EXPECT_GE(model.populations.get_group_total(i), 0);

        EXPECT_NEAR(model.populations.get_group_total(i), fact * num_total_t0, 1e-6);

        EXPECT_GE(params.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i], 0);
    }

    mio::ContactMatrixGroup<double>& contact_matrix_sample = params.get<mio::osecir::ContactPatterns<double>>();
    EXPECT_EQ(contact_matrix_sample[0].get_dampings().size(), 1);
}

TEST(ParameterStudies, sample_graph)
{
    double t0   = 0;
    double tmax = 100;

    double cont_freq = 10; // see Polymod study

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

    mio::osecir::set_params_distributions_normal(model, t0, tmax, 0.2);

    auto graph = mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>>();
    graph.add_node(0, model);
    graph.add_node(1, model);
    graph.add_edge(0, 1, mio::MobilityParameters<double>(Eigen::VectorXd::Constant(Eigen::Index(num_groups * 8), 1.0)));

    mio::ParameterStudy study(graph, 0.0, 0.0, 0.5, 1);
    mio::log_rng_seeds(study.get_rng(), mio::LogLevel::warn);
    auto ensemble_results = study.run_serial([](auto&& g, auto t0_, auto dt_, auto) {
        auto copy = g;
        return mio::make_sampled_graph_simulation<double, mio::osecir::Simulation<ScalarType>>(draw_sample(copy), t0_,
                                                                                               dt_, dt_);
    });

    auto& results = ensemble_results.at(0);
    EXPECT_EQ(results.get_graph().edges()[0].property.get_parameters().get_coefficients()[0].get_dampings().size(), 1);
    for (auto& node : results.get_graph().nodes()) {
        auto& result_model = node.property.get_simulation().get_model();
        EXPECT_EQ(result_model.parameters.get<mio::osecir::ContactPatterns<double>>()
                      .get_cont_freq_mat()[0]
                      .get_dampings()
                      .size(),
                  1);
    }
}

TEST(ParameterStudies, test_normal_distribution)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    mio::ParameterDistributionNormal parameter_dist_normal_1;

    // check if standard deviation is reduced if between too narrow interval [min,max] has to be sampled.
    parameter_dist_normal_1.set_upper_bound(1);
    parameter_dist_normal_1.set_lower_bound(-1);
    parameter_dist_normal_1.log_stddev_changes(false); // only avoid warning output in tests

    double std_dev_demanded = parameter_dist_normal_1.get_standard_dev();
    parameter_dist_normal_1.get_sample(mio::thread_local_rng());

    EXPECT_GE(std_dev_demanded, parameter_dist_normal_1.get_standard_dev());

    // check if full argument constructor works correctly
    mio::ParameterDistributionNormal parameter_dist_normal_2(-1.0, 1.0, 0, parameter_dist_normal_1.get_standard_dev(),
                                                             2.5758);

    EXPECT_EQ(parameter_dist_normal_1.get_lower_bound(), parameter_dist_normal_2.get_lower_bound());
    EXPECT_EQ(parameter_dist_normal_1.get_upper_bound(), parameter_dist_normal_2.get_upper_bound());
    EXPECT_EQ(parameter_dist_normal_1.get_mean(), parameter_dist_normal_2.get_mean());
    EXPECT_EQ(parameter_dist_normal_1.get_standard_dev(), parameter_dist_normal_2.get_standard_dev());

    // check if std_dev is not changed if boundaries are far enough away such that 99% of the density fits into the interval
    parameter_dist_normal_2.set_mean(5);
    parameter_dist_normal_2.set_standard_dev(1.5);
    parameter_dist_normal_2.set_lower_bound(1);
    parameter_dist_normal_2.set_upper_bound(10);
    parameter_dist_normal_2.log_stddev_changes(false); // only avoid warning output in tests
    std_dev_demanded = parameter_dist_normal_2.get_standard_dev();

    parameter_dist_normal_2.check_quantiles();
    EXPECT_EQ(std_dev_demanded, parameter_dist_normal_2.get_standard_dev());

    // check that sampling only occurs in boundaries
    for (int i = 0; i < 1000; i++) {
        double val = parameter_dist_normal_2.get_sample(mio::thread_local_rng());
        EXPECT_GE(parameter_dist_normal_2.get_upper_bound() + 1e-10, val);
        EXPECT_LE(parameter_dist_normal_2.get_lower_bound() - 1e-10, val);
    }

    //degenerate case: ub == lb //For MSVC the normal distribution cannot have a value of 0.0 for sigma
    mio::ParameterDistributionNormal dist3(0.999999999 * 3.0, 1.000000001 * 3.0, 3.0, 0.00000001);
    EXPECT_NEAR(dist3.get_sample(mio::thread_local_rng()), 3.0, 1e-07);
}

TEST(ParameterStudies, test_uniform_distribution)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    // check if full argument constructor works correctly
    mio::ParameterDistributionUniform parameter_dist_unif(1.0, 10.0);

    EXPECT_EQ(parameter_dist_unif.get_lower_bound(), 1.0);
    EXPECT_EQ(parameter_dist_unif.get_upper_bound(), 10.0);

    // check that sampling only occurs in boundaries
    for (int i = 0; i < 1000; i++) {
        double val = parameter_dist_unif.get_sample(mio::thread_local_rng());
        EXPECT_GE(parameter_dist_unif.get_upper_bound() + 1e-10, val);
        EXPECT_LE(parameter_dist_unif.get_lower_bound() - 1e-10, val);
    }
}

TEST(ParameterStudies, test_lognormal_distribution)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    // check if full argument constructor works correctly
    mio::ParameterDistributionLogNormal parameter_dist_lognorm(0.0, 0.25);

    EXPECT_EQ(parameter_dist_lognorm.get_log_mean(), 0.0);
    EXPECT_EQ(parameter_dist_lognorm.get_log_stddev(), 0.25);
}

TEST(ParameterStudies, test_exponential_distribution)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    // check if full argument constructor works correctly
    mio::ParameterDistributionExponential parameter_dist_exponential(1.);

    EXPECT_EQ(parameter_dist_exponential.get_rate(), 1.0);
}

TEST(ParameterStudies, test_constant_distribution)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    // check if full argument constructor works correctly
    mio::ParameterDistributionConstant parameter_dist_constant(2.);

    EXPECT_EQ(parameter_dist_constant.get_constant(), 2.0);
}

TEST(ParameterStudies, test_predefined_samples)
{
    mio::ParameterDistributionUniform parameter_dist_unif(1.0, 10.0);

    mio::ParameterDistributionNormal parameter_dist_normal(-1.0, 1.0, 0, 0.1);

    // set predefined sample (can be out of [min,max]) and get it
    parameter_dist_unif.add_predefined_sample(2);
    double var = parameter_dist_unif.get_sample(mio::thread_local_rng());
    EXPECT_EQ(var, 2);

    // predefined sample was deleted, get real sample which cannot be 2 due to [min,max]
    var = parameter_dist_unif.get_sample(mio::thread_local_rng());
    EXPECT_NE(var, 2);

    // set predefined sample (can be out of [min,max]) and get it
    parameter_dist_normal.add_predefined_sample(2);
    var = parameter_dist_normal.get_sample(mio::thread_local_rng());
    EXPECT_EQ(var, 2);

    // predefined sample was deleted, get real sample which cannot be 2 due to [min,max]
    var = parameter_dist_normal.get_sample(mio::thread_local_rng());
    EXPECT_NE(var, 2);

    mio::ParameterDistributionLogNormal parameter_dist_lognorm(0.0, 0.25);
    //set predefined sample
    parameter_dist_lognorm.add_predefined_sample(-5.);
    var = parameter_dist_lognorm.get_sample(mio::thread_local_rng());
    EXPECT_EQ(var, -5.);

    mio::ParameterDistributionExponential parameter_dist_exponential(1.);
    //set predefined sample
    parameter_dist_exponential.add_predefined_sample(-2.);
    var = parameter_dist_exponential.get_sample(mio::thread_local_rng());
    EXPECT_EQ(var, -2.);

    mio::ParameterDistributionConstant parameter_dist_constant(3.);
    //set predefined sample
    parameter_dist_constant.add_predefined_sample(1.);
    var = parameter_dist_constant.get_sample(mio::thread_local_rng());
    EXPECT_EQ(var, 1.);
    //get another sample with should be the constant
    var = parameter_dist_constant.get_sample(mio::thread_local_rng());
    EXPECT_EQ(var, 3.);
}

TEST(ParameterStudies, check_ensemble_run_result)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    double t0   = 0;
    double tmax = 50;

    double cont_freq = 10; // see Polymod study

    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    mio::osecir::Model<double> model(1);
    mio::AgeGroup num_groups = model.parameters.get_num_groups();
    double fact              = 1.0 / (double)(size_t)num_groups;

    auto& params = model.parameters;

    for (auto i = mio::Index<mio::AgeGroup>(0); i.get() < (size_t)num_groups; i++) {
        params.get<mio::osecir::TimeExposed<double>>()[i]            = 2.6;
        params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i] = 2.6;
        params.get<mio::osecir::TimeInfectedSymptoms<double>>()[i]   = 5.;
        params.get<mio::osecir::TimeInfectedSevere<double>>()[i]     = 10.;
        params.get<mio::osecir::TimeInfectedCritical<double>>()[i]   = 8.;
        params.get<mio::osecir::Seasonality<double>>()               = 0.0;
        params.get<mio::osecir::ICUCapacity<double>>()               = 100.0;
        params.get<mio::osecir::TestAndTraceCapacity<double>>()      = 10.0;

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
    contact_matrix[0] =
        mio::ContactMatrix<double>(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, fact * cont_freq));

    mio::osecir::set_params_distributions_normal(model, t0, tmax, 0.2);
    mio::ParameterStudy parameter_study(model, t0, tmax, 0.1, 1);
    mio::log_rng_seeds(parameter_study.get_rng(), mio::LogLevel::warn);

    // Run parameter study
    auto ensemble_results = parameter_study.run_serial([](auto&& model_, auto t0_, auto dt_, auto) {
        auto copy = model_;
        draw_sample(copy);
        return mio::osecir::Simulation<double>(copy, t0_, dt_);
    });

    const mio::TimeSeries<double>& results = ensemble_results.at(0).get_result();

    for (Eigen::Index i = 0; i < results.get_num_time_points(); i++) {
        std::vector<double> total_at_ti((size_t)mio::osecir::InfectionState::Count, 0);

        for (Eigen::Index j = 0; j < results[i].size(); j++) { // number of compartments per time step
            EXPECT_GE(results[i][j], 0.0) << " day " << results.get_time(i) << " group " << j;
            total_at_ti[static_cast<size_t>(j) / (size_t)mio::osecir::InfectionState::Count] += results[i][j];
        }

        for (auto j = mio::AgeGroup(0); j < params.get_num_groups(); j++) {
            EXPECT_NEAR(total_at_ti[(size_t)j], model.populations.get_group_total(j), 1e-3)
                << " day " << i << " group " << j;
        }
    }
}

namespace
{

struct MockStudyParams {
    const int init, run;
};

struct MockStudySim {
    MockStudySim(const MockStudyParams& p_, double t0_, double dt_)
        : p(p_)
        , t0(t0_)
        , dt(dt_)
    {
    }
    void advance(double t)
    {
        tmax = t;
    }

    MockStudyParams p;
    double t0, dt;
    double tmax = 0;
};

} // namespace

TEST(ParameterStudies, mocked_run)
{
    // run a very simple study, that works with mpi
    const double t0 = 20, tmax = 21, dt = 22;
    const MockStudyParams params{23, -1};
    const size_t num_runs = 5; // enough to notice MPI effects
    const auto make_sim   = [&](auto&& params_, auto t0_, auto dt_, auto i_) {
        MockStudyParams cp{params_.init, (int)i_};
        return MockStudySim(cp, t0_, dt_);
    };
    const auto process_sim = [&](MockStudySim&& s, size_t i) {
        return s.tmax + i;
    };
    const double process_sim_result = (num_runs * tmax) + num_runs * (num_runs - 1) / 2.;
    mio::ParameterStudy study(params, t0, tmax, dt, num_runs);
    // case: run_serial without processing; expect created simulations in order
    auto result_serial = study.run_serial(make_sim);
    EXPECT_EQ(result_serial.size(), num_runs);
    for (int i = 0; const auto& sim : result_serial) {
        EXPECT_EQ(sim.t0, t0);
        EXPECT_EQ(sim.dt, dt);
        EXPECT_EQ(sim.tmax, tmax);
        EXPECT_EQ(sim.p.init, params.init);
        EXPECT_EQ(sim.p.run, i++);
    }
    // case: run and run_serial with processing; expect the same (unordered) result for both, on all ranks
    // Note: currently the tests are not make use of MPI, so we expect the same result from each rank
    auto result_serial_processed = study.run_serial(make_sim, process_sim);
    auto result_parallel         = study.run(make_sim, process_sim);
    for (const auto& result : {result_serial_processed, result_parallel}) {
        EXPECT_EQ(result.size(), num_runs);
        EXPECT_EQ(std::accumulate(result.begin(), result.end(), 0.0), process_sim_result);
    }
}
