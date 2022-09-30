
/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#include "secir/secir.h"
#include "secir/parameter_space.h"
#include "memilio/compartments/parameter_studies.h"
#include "memilio/mobility/mobility.h"
#include "memilio/utils/random_number_generator.h"
#include <gtest/gtest.h>
#include <stdio.h>

TEST(ParameterStudies, sample_from_secir_params)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    double t0   = 0;
    double tmax = 100;

    double tinc    = 5.2, 
        tinf   = 6,
        tserint    = 4.2, // 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        tsevere = 12, 
        ticu2home  = 8, // 5-16 (=R8^(-1))
        ticu2death = 5; // 3.5-7 (=R5^(-1))

    double cont_freq = 10, // see Polymod study
        inf_prob = 0.05, carr_infec = 0.67,
           alpha = 0.09, // 0.01-0.16
        beta     = 0.25, // 0.05-0.5
        delta    = 0.3, // 0.15-0.77
        rho      = 0.2, // 0.1-0.35
        theta    = 0.25; // 0.15-0.4

    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    mio::SecirModel model(3);
    mio::AgeGroup num_groups = model.parameters.get_num_groups();
    double fact              = 1.0 / (double)(size_t)num_groups;

    auto& params = model.parameters;
    for (auto i = mio::AgeGroup(0); i < num_groups; i++) {
        params.get<mio::IncubationTime>()[i]             = tinc;
        params.get<mio::TimeInfectedSymptoms>()[i]         = tinf;
        params.get<mio::SerialInterval>()[i]             = tserint;
        params.get<mio::TimeInfectedSevere>()[i]     = tsevere;
        params.get<mio::ICUToHomeTime>()[i]              = ticu2home;
        params.get<mio::ICUToDeathTime>()[i]             = ticu2death;

        model.populations[{i, mio::InfectionState::Exposed}]      = fact * num_exp_t0;
        model.populations[{i, mio::InfectionState::Carrier}]      = fact * num_car_t0;
        model.populations[{i, mio::InfectionState::Infected}]     = fact * num_inf_t0;
        model.populations[{i, mio::InfectionState::Hospitalized}] = fact * num_hosp_t0;
        model.populations[{i, mio::InfectionState::ICU}]          = fact * num_icu_t0;
        model.populations[{i, mio::InfectionState::Recovered}]    = fact * num_rec_t0;
        model.populations[{i, mio::InfectionState::Dead}]         = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::InfectionState::Susceptible},
                                                                         fact * num_total_t0);

        params.get<mio::InfectionProbabilityFromContact>()[i] = inf_prob;
        params.get<mio::RelativeCarrierInfectability>()[i]    = carr_infec;
        params.get<mio::AsymptomaticCasesPerInfectious>()[i]    = alpha;
        params.get<mio::RiskOfInfectionFromSymptomatic>()[i]   = beta;
        params.get<mio::HospitalizedCasesPerInfectious>()[i]  = rho;
        params.get<mio::ICUCasesPerHospitalized>()[i]         = theta;
        params.get<mio::DeathsPerICU>()[i]                    = delta;
    }

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, fact * cont_freq));

    mio::set_params_distributions_normal(model, t0, tmax, 0.2);

    draw_sample(model);

    for (auto i = mio::AgeGroup(0); i < params.get_num_groups(); i++) {
        ASSERT_EQ(params.get<mio::IncubationTime>()[mio::AgeGroup(0)].value(),
                  params.get<mio::IncubationTime>()[i].value());
        ASSERT_EQ(params.get<mio::SerialInterval>()[mio::AgeGroup(0)].value(),
                  params.get<mio::SerialInterval>()[i].value());
        ASSERT_EQ(params.get<mio::TimeInfectedSymptoms>()[mio::AgeGroup(0)].value(),
                  params.get<mio::TimeInfectedSymptoms>()[i].value());
        ASSERT_EQ(params.get<mio::RelativeCarrierInfectability>()[mio::AgeGroup(0)].value(),
                  params.get<mio::RelativeCarrierInfectability>()[i].value());
        ASSERT_EQ(params.get<mio::RiskOfInfectionFromSymptomatic>()[mio::AgeGroup(0)].value(),
                  params.get<mio::RiskOfInfectionFromSymptomatic>()[i].value());
        ASSERT_EQ(params.get<mio::MaxRiskOfInfectionFromSymptomatic>()[mio::AgeGroup(0)].value(),
                  params.get<mio::MaxRiskOfInfectionFromSymptomatic>()[i].value());

        EXPECT_GE(model.populations.get_group_total(i), 0);

        EXPECT_NEAR(model.populations.get_group_total(i), fact * num_total_t0, 1e-6);

        EXPECT_GE(params.get<mio::IncubationTime>()[i], 0);

        EXPECT_GE(params.get<mio::InfectionProbabilityFromContact>()[i], 0);
    }

    mio::ContactMatrixGroup& contact_matrix_sample = params.get<mio::ContactPatterns>();
    EXPECT_EQ(contact_matrix_sample[0].get_dampings().size(), 1);
}

TEST(ParameterStudies, sample_graph)
{
    double t0   = 0;
    double tmax = 100;

    double tinc    = 5.2, 
        tinf   = 6, 
        tserint    = 4.2, // 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        tsevere = 12, 
        ticu2home  = 8, // 5-16 (=R8^(-1))
        ticu2death = 5; // 3.5-7 (=R5^(-1))

    double cont_freq = 10, // see Polymod study
        inf_prob = 0.05, carr_infec = 0.67,
           alpha = 0.09, // 0.01-0.16
        beta     = 0.25, // 0.05-0.5
        delta    = 0.3, // 0.15-0.77
        rho      = 0.2, // 0.1-0.35
        theta    = 0.25; // 0.15-0.4

    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    size_t num_groups = 3;
    mio::SecirModel model((int)num_groups);
    double fact = 1.0 / (double)num_groups;

    auto& params = model.parameters;
    for (auto i = mio::Index<mio::AgeGroup>(0); i.get() < num_groups; i++) {
        params.get<mio::IncubationTime>()[i]             = tinc;
        params.get<mio::TimeInfectedSymptoms>()[i]         = tinf;
        params.get<mio::SerialInterval>()[i]             = tserint;
        params.get<mio::TimeInfectedSevere>()[i]     = tsevere;
        params.get<mio::ICUToHomeTime>()[i]              = ticu2home;
        params.get<mio::ICUToDeathTime>()[i]             = ticu2death;

        model.populations[{i, mio::InfectionState::Exposed}]      = fact * num_exp_t0;
        model.populations[{i, mio::InfectionState::Carrier}]      = fact * num_car_t0;
        model.populations[{i, mio::InfectionState::Infected}]     = fact * num_inf_t0;
        model.populations[{i, mio::InfectionState::Hospitalized}] = fact * num_hosp_t0;
        model.populations[{i, mio::InfectionState::ICU}]          = fact * num_icu_t0;
        model.populations[{i, mio::InfectionState::Recovered}]    = fact * num_rec_t0;
        model.populations[{i, mio::InfectionState::Dead}]         = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::InfectionState::Susceptible},
                                                                         fact * num_total_t0);

        params.get<mio::InfectionProbabilityFromContact>()[i] = inf_prob;
        params.get<mio::RelativeCarrierInfectability>()[i]    = carr_infec;
        params.get<mio::AsymptomaticCasesPerInfectious>()[i]    = alpha;
        params.get<mio::RiskOfInfectionFromSymptomatic>()[i]   = beta;
        params.get<mio::HospitalizedCasesPerInfectious>()[i]  = rho;
        params.get<mio::ICUCasesPerHospitalized>()[i]         = theta;
        params.get<mio::DeathsPerICU>()[i]                    = delta;
    }

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::ContactPatterns>();
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_groups, num_groups, fact * cont_freq));

    mio::set_params_distributions_normal(model, t0, tmax, 0.2);

    auto graph = mio::Graph<mio::SecirModel, mio::MigrationParameters>();
    graph.add_node(0, model);
    graph.add_node(1, model);
    graph.add_edge(0, 1, mio::MigrationParameters(Eigen::VectorXd::Constant(Eigen::Index(num_groups * 8), 1.0)));

    auto study   = mio::ParameterStudy<mio::SecirSimulation<>>(graph, 0.0, 0.0, 0.5, 1);
    auto results = study.run([](auto&& g) {
        return draw_sample(g);
    });

    EXPECT_EQ(results[0].edges()[0].property.get_parameters().get_coefficients()[0].get_dampings().size(), 1);
    for (auto& node : results[0].nodes()) {
        auto& result_model = node.property.get_simulation().get_model();
        EXPECT_EQ(result_model.parameters.get<mio::ContactPatterns>().get_cont_freq_mat()[0].get_dampings().size(), 1);
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
    parameter_dist_normal_1.get_sample();

    EXPECT_GE(std_dev_demanded, parameter_dist_normal_1.get_standard_dev());

    // check if full argument constructor works correctly
    mio::ParameterDistributionNormal parameter_dist_normal_2(-1.0, 1.0, 0, parameter_dist_normal_1.get_standard_dev());

    EXPECT_EQ(parameter_dist_normal_1.get_lower_bound(), parameter_dist_normal_2.get_lower_bound());
    EXPECT_EQ(parameter_dist_normal_1.get_upper_bound(), parameter_dist_normal_2.get_upper_bound());
    EXPECT_EQ(parameter_dist_normal_1.get_mean(), parameter_dist_normal_2.get_mean());
    EXPECT_EQ(parameter_dist_normal_1.get_standard_dev(), parameter_dist_normal_2.get_standard_dev());

    // check if std_dev is not changed if boundaries are far enough away such that 99% of the density fits into the interval
    parameter_dist_normal_2.set_mean(5);
    parameter_dist_normal_2.set_standard_dev(1.5);
    parameter_dist_normal_2.set_lower_bound(1);
    parameter_dist_normal_2.set_upper_bound(10);
    std_dev_demanded = parameter_dist_normal_2.get_standard_dev();

    parameter_dist_normal_2.check_quantiles();
    EXPECT_EQ(std_dev_demanded, parameter_dist_normal_2.get_standard_dev());

    // check that sampling only occurs in boundaries
    for (int i = 0; i < 1000; i++) {
        double val = parameter_dist_normal_2.get_sample();
        EXPECT_GE(parameter_dist_normal_2.get_upper_bound() + 1e-10, val);
        EXPECT_LE(parameter_dist_normal_2.get_lower_bound() - 1e-10, val);
    }

    //degenerate case: ub == lb
    mio::ParameterDistributionNormal dist3(3.0, 3.0, 3.0, 3.0);
    EXPECT_EQ(dist3.get_sample(), 3.0);
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
        double val = parameter_dist_unif.get_sample();
        EXPECT_GE(parameter_dist_unif.get_upper_bound() + 1e-10, val);
        EXPECT_LE(parameter_dist_unif.get_lower_bound() - 1e-10, val);
    }
}

TEST(ParameterStudies, test_predefined_samples)
{
    mio::ParameterDistributionUniform parameter_dist_unif(1.0, 10.0);

    mio::ParameterDistributionNormal parameter_dist_normal(-1.0, 1.0, 0, 0.1);

    // set predefined sample (can be out of [min,max]) and get it
    parameter_dist_unif.add_predefined_sample(2);
    double var = parameter_dist_unif.get_sample();
    EXPECT_EQ(var, 2);

    // predefined sample was deleted, get real sample which cannot be 2 due to [min,max]
    var = parameter_dist_unif.get_sample();
    EXPECT_NE(var, 2);

    // set predefined sample (can be out of [min,max]) and get it
    parameter_dist_normal.add_predefined_sample(2);
    var = parameter_dist_normal.get_sample();
    EXPECT_EQ(var, 2);

    // predefined sample was deleted, get real sample which cannot be 2 due to [min,max]
    var = parameter_dist_normal.get_sample();
    EXPECT_NE(var, 2);
}

TEST(ParameterStudies, check_ensemble_run_result)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    double t0   = 0;
    double tmax = 50;

    double tinc    = 5.2, 
        tinf   = 6, 
        tserint    = 4.2, // 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        tsevere = 12, 
        ticu2home  = 8, // 5-16 (=R8^(-1))
        ticu2death = 5; // 3.5-7 (=R5^(-1))

    double cont_freq = 10, // see Polymod study
        inf_prob = 0.05, carr_infec = 0.67,
           alpha = 0.09, // 0.01-0.16
        beta     = 0.25, // 0.05-0.5
        delta    = 0.3, // 0.15-0.77
        rho      = 0.2, // 0.1-0.35
        theta    = 0.25; // 0.15-0.4

    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    mio::SecirModel model(1);
    mio::AgeGroup num_groups = model.parameters.get_num_groups();
    double fact              = 1.0 / (double)(size_t)num_groups;

    auto& params = model.parameters;

    for (auto i = mio::AgeGroup(0); i < num_groups; i++) {
        params.get<mio::IncubationTime>()[i]             = tinc;
        params.get<mio::TimeInfectedSymptoms>()[i]         = tinf;
        params.get<mio::SerialInterval>()[i]             = tserint;
        params.get<mio::TimeInfectedSevere>()[i]     = tsevere;
        params.get<mio::ICUToHomeTime>()[i]              = ticu2home;
        params.get<mio::ICUToDeathTime>()[i]             = ticu2death;

        model.populations.set_total(num_total_t0);
        model.populations[{i, mio::InfectionState::Exposed}]      = num_exp_t0;
        model.populations[{i, mio::InfectionState::Carrier}]      = num_car_t0;
        model.populations[{i, mio::InfectionState::Infected}]     = num_inf_t0;
        model.populations[{i, mio::InfectionState::Hospitalized}] = num_hosp_t0;
        model.populations[{i, mio::InfectionState::ICU}]          = num_icu_t0;
        model.populations[{i, mio::InfectionState::Recovered}]    = num_rec_t0;
        model.populations[{i, mio::InfectionState::Dead}]         = num_dead_t0;
        model.populations.set_difference_from_total({i, mio::InfectionState::Susceptible}, num_total_t0);

        params.get<mio::InfectionProbabilityFromContact>()[i] = inf_prob;
        params.get<mio::RelativeCarrierInfectability>()[i]    = carr_infec;
        params.get<mio::AsymptomaticCasesPerInfectious>()[i]    = alpha;
        params.get<mio::RiskOfInfectionFromSymptomatic>()[i]   = beta;
        params.get<mio::HospitalizedCasesPerInfectious>()[i]  = rho;
        params.get<mio::ICUCasesPerHospitalized>()[i]         = theta;
        params.get<mio::DeathsPerICU>()[i]                    = delta;
    }

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, fact * cont_freq));

    mio::set_params_distributions_normal(model, t0, tmax, 0.2);
    mio::ParameterStudy<mio::SecirSimulation<>> parameter_study(model, t0, tmax, 1);

    // Run parameter study
    parameter_study.set_num_runs(1);
    auto graph_results = parameter_study.run([](auto&& g) {
        return draw_sample(g);
    });

    std::vector<mio::TimeSeries<double>> results;
    for (size_t i = 0; i < graph_results.size(); i++) {
        results.push_back(std::move(graph_results[i].nodes()[0].property.get_result()));
    }

    for (Eigen::Index i = 0; i < results[0].get_num_time_points(); i++) {
        std::vector<double> total_at_ti((size_t)mio::InfectionState::Count, 0);

        for (Eigen::Index j = 0; j < results[0][i].size(); j++) { // number of compartments per time step
            EXPECT_GE(results[0][i][j], 0.0) << " day " << results[0].get_time(i) << " group " << j;
            total_at_ti[static_cast<size_t>(j) / (size_t)mio::InfectionState::Count] += results[0][i][j];
        }

        for (auto j = mio::AgeGroup(0); j < params.get_num_groups(); j++) {
            EXPECT_NEAR(total_at_ti[(size_t)j], model.populations.get_group_total(j), 1e-3)
                << " day " << i << " group " << j;
        }
    }
}
