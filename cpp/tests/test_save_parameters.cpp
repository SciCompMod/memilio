/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Wadim Koslow
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
#include "test_data_dir.h"
#include "ode_secir/model.h"
#include "ode_secir/parameter_space.h"
#include "ode_secir/parameters_io.h"
#include "ode_secirvvs/model.h"
#include "ode_secirvvs/parameters_io.h"
#include "memilio/io/mobility_io.h"
#include "memilio/io/result_io.h"
#include <distributions_helpers.h>
#include <matchers.h>
#include "temp_file_register.h"
#include "memilio/utils/date.h"
#include <gtest/gtest.h>

TEST(TestSaveParameters, json_single_sim_write_read_compare)
{
    double t0   = 0.0;
    double tmax = 50.5;

    double cont_freq = 10, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    mio::osecir::Model<double> model(2);
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

        model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i] = 0.06;
        model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[i]   = 0.67;
        model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[i]   = alpha;
        model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i]   = beta;
        model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()   = 0.85;
        model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i]        = rho;
        model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[i]                = theta;
        model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[i]                = delta;
    }

    mio::ContactMatrixGroup<double>& contact_matrix = params.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0] =
        mio::ContactMatrix<double>(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, fact * cont_freq));
    contact_matrix.add_damping(0.7, mio::SimulationTime<double>(30.));
    auto damping2  = Eigen::MatrixXd::Zero((size_t)num_groups, (size_t)num_groups).eval();
    damping2(0, 0) = 0.8;
    contact_matrix.add_damping(damping2, mio::SimulationTime<double>(35));

    mio::osecir::set_params_distributions_normal(model, t0, tmax, 0.2);

    params.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0].get_distribution()->add_predefined_sample(4711.0);

    TempFileRegister file_register;
    auto filename     = file_register.get_unique_path("TestParameters-%%%%-%%%%.json");
    auto write_status = mio::write_json(filename, model);
    ASSERT_THAT(print_wrap(write_status), IsSuccess());

    auto read_result = mio::read_json(filename, mio::Tag<mio::osecir::Model<double>>{});
    ASSERT_THAT(print_wrap(read_result), IsSuccess());
    auto& read_model = read_result.value();

    const mio::UncertainContactMatrix<double>& contact = model.parameters.get<mio::osecir::ContactPatterns<double>>();
    const mio::UncertainContactMatrix<double>& read_contact =
        read_model.parameters.get<mio::osecir::ContactPatterns<double>>();

    num_groups           = model.parameters.get_num_groups();
    auto num_groups_read = read_model.parameters.get_num_groups();
    ASSERT_EQ(num_groups, num_groups_read);

    for (auto i = mio::AgeGroup(0); i < num_groups; i++) {
        ASSERT_EQ((model.populations[{i, mio::osecir::InfectionState::Dead}]),
                  (read_model.populations[{i, mio::osecir::InfectionState::Dead}]));
        ASSERT_EQ((model.populations.get_group_total(i)), (read_model.populations.get_group_total(i)));
        ASSERT_EQ((model.populations[{i, mio::osecir::InfectionState::Exposed}]),
                  (read_model.populations[{i, mio::osecir::InfectionState::Exposed}]));
        ASSERT_EQ((model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}]),
                  (read_model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}]));
        ASSERT_EQ((model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]),
                  (read_model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]));
        ASSERT_EQ((model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]),
                  (read_model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]));
        ASSERT_EQ((model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]),
                  (read_model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]));
        ASSERT_EQ((model.populations[{i, mio::osecir::InfectionState::Recovered}]),
                  (read_model.populations[{i, mio::osecir::InfectionState::Recovered}]));

        check_distribution(*model.populations[{i, mio::osecir::InfectionState::Exposed}].get_distribution(),
                           *read_model.populations[{i, mio::osecir::InfectionState::Exposed}].get_distribution());
        check_distribution(
            *model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}].get_distribution(),
            *read_model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}].get_distribution());
        check_distribution(
            *model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}].get_distribution(),
            *read_model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}].get_distribution());
        check_distribution(
            *model.populations[{i, mio::osecir::InfectionState::InfectedSevere}].get_distribution(),
            *read_model.populations[{i, mio::osecir::InfectionState::InfectedSevere}].get_distribution());
        check_distribution(
            *model.populations[{i, mio::osecir::InfectionState::InfectedCritical}].get_distribution(),
            *read_model.populations[{i, mio::osecir::InfectionState::InfectedCritical}].get_distribution());
        check_distribution(*model.populations[{i, mio::osecir::InfectionState::Recovered}].get_distribution(),
                           *read_model.populations[{i, mio::osecir::InfectionState::Recovered}].get_distribution());

        ASSERT_EQ(model.parameters.get<mio::osecir::TimeExposed<double>>()[i],
                  read_model.parameters.get<mio::osecir::TimeExposed<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i],
                  read_model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[i],
                  read_model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[i],
                  read_model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[i],
                  read_model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[i]);

        check_distribution(*model.parameters.get<mio::osecir::TimeExposed<double>>()[i].get_distribution(),
                           *read_model.parameters.get<mio::osecir::TimeExposed<double>>()[i].get_distribution());
        check_distribution(
            *model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i].get_distribution(),
            *read_model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i].get_distribution());
        check_distribution(
            *model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[i].get_distribution(),
            *read_model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[i].get_distribution(),
                           *read_model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[i].get_distribution());
        check_distribution(
            *model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[i].get_distribution(),
            *read_model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[i].get_distribution());

        ASSERT_EQ(model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i],
                  read_model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i],
                  read_model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[i],
                  read_model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[i],
                  read_model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i],
                  read_model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[i],
                  read_model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[i]);

        check_distribution(
            *model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i].get_distribution(),
            *read_model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i].get_distribution());
        check_distribution(
            *model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i].get_distribution(),
            *read_model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i].get_distribution());
        check_distribution(
            *model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[i].get_distribution(),
            *read_model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[i].get_distribution(),
                           *read_model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[i].get_distribution());
        check_distribution(
            *model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i].get_distribution(),
            *read_model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[i].get_distribution(),
                           *read_model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[i].get_distribution());

        ASSERT_THAT(contact.get_cont_freq_mat(), testing::ContainerEq(read_contact.get_cont_freq_mat()));
        ASSERT_EQ(contact.get_dampings(), read_contact.get_dampings());
    }
}

TEST(TestSaveParameters, read_graph_without_edges)
{
    // set up parameter study
    double cont_freq    = 10; // see Polymod study
    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    size_t num_groups = 3;
    mio::osecir::Model<double> model((int)num_groups);
    double fact = 1.0 / (double)num_groups;

    auto& params = model.parameters;
    for (auto i = mio::Index<mio::AgeGroup>(0); i.get() < num_groups; i++) {
        params.get<mio::osecir::TimeExposed<double>>()[i]            = 3.2;
        params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i] = 2.;
        params.get<mio::osecir::TimeInfectedSymptoms<double>>()[i]   = 5.;
        params.get<mio::osecir::TimeInfectedSevere<double>>()[i]     = 10.;
        params.get<mio::osecir::TimeInfectedCritical<double>>()[i]   = 8.;

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
        params.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i]        = 0.2;
        params.get<mio::osecir::CriticalPerSevere<double>>()[i]                = 0.25;
        params.get<mio::osecir::DeathsPerCritical<double>>()[i]                = 0.3;
    }

    mio::ContactMatrixGroup<double>& contact_matrix = params.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0] = mio::ContactMatrix<double>(Eigen::MatrixXd::Constant(num_groups, num_groups, fact * cont_freq));

    auto graph = mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>>();
    graph.add_node(0, model);
    graph.add_node(1, model);
    graph.add_edge(0, 1, mio::MobilityParameters<double>(Eigen::VectorXd::Constant(Eigen::Index(num_groups * 8), 1.0)));

    TempFileRegister tmp_file_register;
    std::string tmp_results_dir = tmp_file_register.get_unique_path();
    ASSERT_THAT(mio::create_directory(tmp_results_dir), IsSuccess());

    std::vector<mio::osecir::Model<double>> models = {model, model};
    std::vector<int> ids                           = {0, 1};
    auto graph_no_edges =
        mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>>(ids, models);
    auto write_status = mio::write_graph(graph_no_edges, tmp_results_dir);
    ASSERT_THAT(print_wrap(write_status), IsSuccess());

    auto read_graph =
        mio::read_graph<double, mio::osecir::Model<double>>(tmp_results_dir, mio::IOF_OmitDistributions, false);

    for (auto i = mio::AgeGroup(0); i < params.get_num_groups(); i++) {
        EXPECT_EQ(read_graph.value().nodes()[0].property.parameters.get<mio::osecir::TimeExposed<double>>()[i],
                  params.get<mio::osecir::TimeExposed<double>>()[i]);

        EXPECT_EQ(read_graph.value().nodes()[0].property.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[i],
                  params.get<mio::osecir::TimeInfectedSymptoms<double>>()[i]);

        EXPECT_EQ(
            read_graph.value().nodes()[0].property.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i],
            params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i]);

        EXPECT_EQ(read_graph.value().nodes()[0].property.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[i],
                  params.get<mio::osecir::TimeInfectedSevere<double>>()[i]);

        EXPECT_EQ(read_graph.value().nodes()[0].property.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[i],
                  params.get<mio::osecir::TimeInfectedCritical<double>>()[i]);

        EXPECT_EQ(read_graph.value()
                      .nodes()[0]
                      .property.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i],
                  params.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i]);

        EXPECT_EQ(read_graph.value()
                      .nodes()[0]
                      .property.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[i],
                  params.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[i]);

        EXPECT_EQ(read_graph.value()
                      .nodes()[0]
                      .property.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i],
                  params.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i]);

        EXPECT_EQ(
            read_graph.value().nodes()[0].property.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i],
            params.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i]);

        EXPECT_EQ(read_graph.value().nodes()[0].property.parameters.get<mio::osecir::DeathsPerCritical<double>>()[i],
                  params.get<mio::osecir::DeathsPerCritical<double>>()[i]);
    }
}

TEST(TestSaveParameters, json_uncertain_matrix_write_read_compare)
{
    enum class InterventionLevelMock
    {
        level,
    };
    enum class InterventionMock
    {
        intervention,
    };
    enum class ContactLocationMock
    {
        location,
    };

    const auto start_date = mio::Date(2020, 12, 12);
    auto damping_time1    = mio::SimulationTime<double>(mio::get_offset_in_days(mio::Date(2020, 12, 1), start_date));
    auto damping_time2    = mio::SimulationTime<double>(mio::get_offset_in_days(mio::Date(2020, 12, 24), start_date));
    mio::osecir::Model<double> model(2);
    auto& contacts         = model.parameters.get<mio::osecir::ContactPatterns<double>>();
    auto& contact_dampings = contacts.get_dampings();

    auto group_weights = Eigen::VectorXd::Constant(size_t(model.parameters.get_num_groups()), 1.0);

    //add two dampings
    auto v1 = mio::UncertainValue<double>(0.5 * (0.6 + 0.4));
    v1.set_distribution(mio::ParameterDistributionUniform(0.4, 0.6));
    auto v2 = mio::UncertainValue<double>(0.5 * (0.3 + 0.2));
    v1.set_distribution(mio::ParameterDistributionUniform(0.2, 0.3));
    contact_dampings.push_back(mio::DampingSampling<double>(
        v1, mio::DampingLevel(int(InterventionLevelMock::level)), mio::DampingType(int(InterventionMock::intervention)),
        damping_time1, {size_t(ContactLocationMock::location)}, group_weights));
    contact_dampings.push_back(mio::DampingSampling<double>(
        v2, mio::DampingLevel(int(InterventionLevelMock::level)), mio::DampingType(int(InterventionMock::intervention)),
        damping_time2, {size_t(ContactLocationMock::location)}, group_weights));
    //add school_holiday_damping
    auto school_holiday_value = mio::UncertainValue<double>(1);
    school_holiday_value.set_distribution(mio::ParameterDistributionUniform(1, 1));
    contacts.get_school_holiday_damping() = mio::DampingSampling<double>(
        school_holiday_value, mio::DampingLevel(int(InterventionLevelMock::level)),
        mio::DampingType(int(InterventionMock::intervention)), mio::SimulationTime<double>(0.0),
        {size_t(ContactLocationMock::location)}, group_weights);

    //add school holidays
    auto holiday_start_time = mio::SimulationTime<double>(mio::get_offset_in_days(mio::Date(2020, 12, 23), start_date));
    auto holiday_end_time   = mio::SimulationTime<double>(mio::get_offset_in_days(mio::Date(2021, 1, 6), start_date));
    contacts.get_school_holidays() = {std::make_pair(holiday_start_time, holiday_end_time)};

    //write json
    TempFileRegister file_register;
    auto filename     = file_register.get_unique_path("TestParameters-%%%%-%%%%.json");
    auto write_status = mio::write_json(filename, model);
    ASSERT_THAT(print_wrap(write_status), IsSuccess());

    //read json
    auto read_result = mio::read_json(filename, mio::Tag<mio::osecir::Model<double>>{});
    ASSERT_THAT(print_wrap(read_result), IsSuccess());
    auto& read_model = read_result.value();

    auto& read_contacts = read_model.parameters.get<mio::osecir::ContactPatterns<double>>();

    //check parameters
    ASSERT_EQ(contacts.get_dampings(), read_contacts.get_dampings());
    ASSERT_EQ(contacts.get_school_holiday_damping(), read_contacts.get_school_holiday_damping());
    ASSERT_EQ(contacts.get_school_holidays(), read_contacts.get_school_holidays());
}

TEST(TestSaveParameters, json_graphs_write_read_compare)
{
    double t0   = 0.0;
    double tmax = 50.5;

    double cont_freq = 10, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    mio::osecir::Model<double> model(2);
    mio::AgeGroup num_groups = model.parameters.get_num_groups();
    double fact              = 1.0 / (double)(size_t)num_groups;

    model.parameters.set<mio::osecir::TestAndTraceCapacity<double>>(30);
    auto& params = model.parameters;
    for (auto i = mio::Index<mio::AgeGroup>(0); i.get() < (size_t)num_groups; i++) {
        params.get<mio::osecir::TimeExposed<double>>()[i]            = 2.6;
        params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i] = 2.6;
        params.get<mio::osecir::TimeInfectedSymptoms<double>>()[i]   = 5.;
        params.get<mio::osecir::TimeInfectedSevere<double>>()[i]     = 10.;
        params.get<mio::osecir::TimeInfectedCritical<double>>()[i]   = 8.;

        params.get<mio::osecir::Seasonality<double>>() = 0.0;
        params.get<mio::osecir::ICUCapacity<double>>() = 100.0;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]            = fact * num_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}] = fact * num_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]   = fact * num_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]     = fact * num_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]   = fact * num_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]          = fact * num_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]               = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * num_total_t0);

        model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i]  = 0.06;
        model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[i]    = 0.67;
        model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[i]    = alpha;
        model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i]    = beta;
        model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[i] = beta * 3;
        model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i]         = rho;
        model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[i]                 = theta;
        model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[i]                 = delta;
    }

    mio::ContactMatrixGroup<double>& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0] =
        mio::ContactMatrix<double>(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, fact * cont_freq));
    Eigen::MatrixXd m =
        Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, 0.7).triangularView<Eigen::Upper>();
    contact_matrix.add_damping(m, mio::SimulationTime<double>(30.));

    mio::osecir::set_params_distributions_normal(model, t0, tmax, 0.15);

    mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>> graph;
    graph.add_node(0, model);
    graph.add_node(1, model);
    graph.add_edge(0, 1, Eigen::VectorXd::Constant(model.populations.get_num_compartments(), 0.01));
    graph.add_edge(1, 0, Eigen::VectorXd::Constant(model.populations.get_num_compartments(), 0.01));

    TempFileRegister file_register;
    auto graph_dir    = file_register.get_unique_path("graph_parameters-%%%%-%%%%");
    auto write_status = mio::write_graph(graph, graph_dir);
    ASSERT_THAT(print_wrap(write_status), IsSuccess());

    auto read_result = mio::read_graph<double, mio::osecir::Model<double>>(graph_dir);
    ASSERT_THAT(print_wrap(read_result), IsSuccess());

    auto& graph_read = read_result.value();
    auto num_nodes   = graph.nodes().size();
    auto num_edges   = graph.edges().size();

    ASSERT_EQ(num_nodes, graph_read.nodes().size());
    ASSERT_EQ(num_edges, graph_read.edges().size());

    for (size_t node = 0; node < num_nodes; node++) {
        mio::osecir::Model<double> graph_model = graph.nodes()[0].property;
        mio::ContactMatrixGroup<double>& graph_cont_matrix =
            graph_model.parameters.get<mio::osecir::ContactPatterns<double>>();

        mio::osecir::Model<double> graph_read_model = graph_read.nodes()[0].property;
        mio::ContactMatrixGroup<double>& graph_read_cont_matrix =
            graph_read_model.parameters.get<mio::osecir::ContactPatterns<double>>();

        ASSERT_EQ(graph_read_cont_matrix.get_num_groups(), static_cast<Eigen::Index>((size_t)num_groups));
        ASSERT_EQ(graph_read_cont_matrix, graph_cont_matrix);
        ASSERT_EQ(graph_model.populations.get_num_compartments(), graph_read_model.populations.get_num_compartments());
        ASSERT_EQ(graph.nodes()[node].id, graph_read.nodes()[node].id);
        EXPECT_THAT(graph_read_model.parameters.get<mio::osecir::TestAndTraceCapacity<double>>().value(),
                    FloatingPointEqual(graph_model.parameters.get<mio::osecir::TestAndTraceCapacity<double>>().value(),
                                       1e-12, 1e-12));
        check_distribution(
            *graph_model.parameters.get<mio::osecir::TestAndTraceCapacity<double>>().get_distribution().get(),
            *graph_read_model.parameters.get<mio::osecir::TestAndTraceCapacity<double>>().get_distribution().get());

        for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(num_groups); group++) {
            ASSERT_EQ((graph_model.populations[{group, mio::osecir::InfectionState::Dead}]),
                      (graph_read_model.populations[{group, mio::osecir::InfectionState::Dead}]));
            ASSERT_EQ(graph_model.populations.get_total(), graph_read_model.populations.get_total());
            check_distribution(
                *graph_model.populations[{group, mio::osecir::InfectionState::Exposed}].get_distribution().get(),
                *graph_read_model.populations[{group, mio::osecir::InfectionState::Exposed}].get_distribution().get());
            check_distribution(*graph_model.populations[{group, mio::osecir::InfectionState::InfectedNoSymptoms}]
                                    .get_distribution()
                                    .get(),
                               *graph_read_model.populations[{group, mio::osecir::InfectionState::InfectedNoSymptoms}]
                                    .get_distribution()
                                    .get());
            check_distribution(*graph_model.populations[{group, mio::osecir::InfectionState::InfectedSymptoms}]
                                    .get_distribution()
                                    .get(),
                               *graph_read_model.populations[{group, mio::osecir::InfectionState::InfectedSymptoms}]
                                    .get_distribution()
                                    .get());
            check_distribution(
                *graph_model.populations[{group, mio::osecir::InfectionState::InfectedSevere}].get_distribution().get(),
                *graph_read_model.populations[{group, mio::osecir::InfectionState::InfectedSevere}]
                     .get_distribution()
                     .get());
            check_distribution(*graph_model.populations[{group, mio::osecir::InfectionState::InfectedCritical}]
                                    .get_distribution()
                                    .get(),
                               *graph_read_model.populations[{group, mio::osecir::InfectionState::InfectedCritical}]
                                    .get_distribution()
                                    .get());
            check_distribution(
                *graph_model.populations[{group, mio::osecir::InfectionState::Recovered}].get_distribution().get(),
                *graph_read_model.populations[{group, mio::osecir::InfectionState::Recovered}]
                     .get_distribution()
                     .get());
            check_distribution(
                *graph_model.populations[{group, mio::osecir::InfectionState::Exposed}].get_distribution().get(),
                *graph_read_model.populations[{group, mio::osecir::InfectionState::Exposed}].get_distribution().get());

            ASSERT_EQ(graph_model.parameters.get<mio::osecir::TimeExposed<double>>()[group],
                      graph_read_model.parameters.get<mio::osecir::TimeExposed<double>>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[group],
                      graph_read_model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[group],
                      graph_read_model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[group],
                      graph_read_model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[group],
                      graph_read_model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[group]);

            ASSERT_EQ(graph_model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[group],
                      graph_read_model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[group],
                      graph_read_model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[group],
                      graph_read_model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[group],
                      graph_read_model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[group],
                      graph_read_model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[group],
                      graph_read_model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[group],
                      graph_read_model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[group]);

            check_distribution(
                *graph_model.parameters.get<mio::osecir::TimeExposed<double>>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::osecir::TimeExposed<double>>()[group].get_distribution().get());
            check_distribution(*graph_model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[group]
                                    .get_distribution()
                                    .get(),
                               *graph_read_model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[group]
                                    .get_distribution()
                                    .get());
            check_distribution(*graph_model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[group]
                                    .get_distribution()
                                    .get(),
                               *graph_read_model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[group]
                                    .get_distribution()
                                    .get());
            check_distribution(
                *graph_model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[group]
                     .get_distribution()
                     .get());
            check_distribution(*graph_model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[group]
                                    .get_distribution()
                                    .get(),
                               *graph_read_model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[group]
                                    .get_distribution()
                                    .get());

            check_distribution(*graph_model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[group]
                                    .get_distribution()
                                    .get(),
                               *graph_read_model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[group]
                                    .get_distribution()
                                    .get());
            check_distribution(
                *graph_model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[group]
                     .get_distribution()
                     .get());
            check_distribution(
                *graph_model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[group]
                     .get_distribution()
                     .get(),
                *graph_read_model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[group]
                     .get_distribution()
                     .get());
            check_distribution(
                *graph_model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[group]
                     .get_distribution()
                     .get());
            check_distribution(*graph_model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[group]
                                    .get_distribution()
                                    .get(),
                               *graph_read_model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[group]
                                    .get_distribution()
                                    .get());

            check_distribution(*graph_model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[group]
                                    .get_distribution()
                                    .get(),
                               *graph_read_model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[group]
                                    .get_distribution()
                                    .get());
            check_distribution(
                *graph_model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[group]
                     .get_distribution()
                     .get());
            check_distribution(
                *graph_model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[group]
                     .get_distribution()
                     .get(),
                *graph_read_model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[group]
                     .get_distribution()
                     .get());
            check_distribution(
                *graph_model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[group]
                     .get_distribution()
                     .get());
            check_distribution(
                *graph_model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[group]
                     .get_distribution()
                     .get());

            ASSERT_EQ(graph_model.parameters.get<mio::osecir::ContactPatterns<double>>().get_dampings(),
                      graph_read_model.parameters.get<mio::osecir::ContactPatterns<double>>().get_dampings());
        }

        ASSERT_THAT(graph_read.edges(), testing::ElementsAreArray(graph.edges()));
    }
}

TEST(TestSaveParameters, json_write_read_parameters_secirvvs)
{
    mio::osecirvvs::Model<double> model(2);
    mio::AgeGroup num_groups = model.parameters.get_num_groups();
    double fact              = 1.0 / (double)(size_t)num_groups;

    auto& params = model.parameters;

    for (auto i = mio::AgeGroup(0); i < num_groups; i++) {
        params.get<mio::osecirvvs::TimeExposed<double>>()[i]            = 3.2;
        params.get<mio::osecirvvs::TimeInfectedNoSymptoms<double>>()[i] = 2.;
        params.get<mio::osecirvvs::TimeInfectedSymptoms<double>>()[i]   = 5.;
        params.get<mio::osecirvvs::TimeInfectedSevere<double>>()[i]     = 10.;
        params.get<mio::osecirvvs::TimeInfectedCritical<double>>()[i]   = 8.;

        model.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>()[i] = 0.06;
        model.parameters.get<mio::osecirvvs::RelativeTransmissionNoSymptoms<double>>()[i]   = 0.67;
        model.parameters.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<double>>()[i]   = 0.09;
        model.parameters.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic<double>>()[i]   = 0.25;
        model.parameters.get<mio::osecirvvs::SeverePerInfectedSymptoms<double>>()[i]        = 0.2;
        model.parameters.get<mio::osecirvvs::CriticalPerSevere<double>>()[i]                = 0.25;
        model.parameters.get<mio::osecirvvs::DeathsPerCritical<double>>()[i]                = 0.3;

        model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[i]  = 0.8;
        model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[i] = 0.7;
        model.parameters.get<mio::osecirvvs::ReducExposedPartialImmunity<double>>()[i]                     = 0.8;
        model.parameters.get<mio::osecirvvs::ReducExposedImprovedImmunity<double>>()[i]                    = 0.7;
        model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<double>>()[i]            = 0.8;
        model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<double>>()[i]           = 0.7;
        model.parameters.get<mio::osecirvvs::ReducTimeInfectedMild<double>>()[i]                           = 0.5;
    }

    mio::ContactMatrixGroup<double>& contact_matrix = params.get<mio::osecirvvs::ContactPatterns<double>>();
    contact_matrix[0] =
        mio::ContactMatrix<double>(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, fact * 10));
    contact_matrix.add_damping(0.7, mio::SimulationTime<double>(30.));
    auto damping2  = Eigen::MatrixXd::Zero((size_t)num_groups, (size_t)num_groups).eval();
    damping2(0, 0) = 0.8;
    contact_matrix.add_damping(damping2, mio::SimulationTime<double>(35));

    TempFileRegister file_register;
    auto filename     = file_register.get_unique_path("TestParameters-%%%%-%%%%.json");
    auto write_status = mio::write_json(filename, model);
    ASSERT_THAT(print_wrap(write_status), IsSuccess());

    auto read_result = mio::read_json(filename, mio::Tag<mio::osecirvvs::Model<double>>{});
    ASSERT_THAT(print_wrap(read_result), IsSuccess());
    auto& read_model = read_result.value();

    const mio::UncertainContactMatrix<double>& contact =
        model.parameters.get<mio::osecirvvs::ContactPatterns<double>>();
    const mio::UncertainContactMatrix<double>& read_contact =
        read_model.parameters.get<mio::osecirvvs::ContactPatterns<double>>();

    num_groups           = model.parameters.get_num_groups();
    auto num_groups_read = read_model.parameters.get_num_groups();
    ASSERT_EQ(num_groups, num_groups_read);

    for (auto i = mio::AgeGroup(0); i < num_groups; i++) {
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::TimeExposed<double>>()[i],
                  read_model.parameters.get<mio::osecirvvs::TimeExposed<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedNoSymptoms<double>>()[i],
                  read_model.parameters.get<mio::osecirvvs::TimeInfectedNoSymptoms<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedSymptoms<double>>()[i],
                  read_model.parameters.get<mio::osecirvvs::TimeInfectedSymptoms<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedSevere<double>>()[i],
                  read_model.parameters.get<mio::osecirvvs::TimeInfectedSevere<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedCritical<double>>()[i],
                  read_model.parameters.get<mio::osecirvvs::TimeInfectedCritical<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>()[i],
                  read_model.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic<double>>()[i],
                  read_model.parameters.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<double>>()[i],
                  read_model.parameters.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::DeathsPerCritical<double>>()[i],
                  read_model.parameters.get<mio::osecirvvs::DeathsPerCritical<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::SeverePerInfectedSymptoms<double>>()[i],
                  read_model.parameters.get<mio::osecirvvs::SeverePerInfectedSymptoms<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::CriticalPerSevere<double>>()[i],
                  read_model.parameters.get<mio::osecirvvs::CriticalPerSevere<double>>()[i]);
        ASSERT_EQ(
            model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[i],
            read_model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[i]);
        ASSERT_EQ(
            model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[i],
            read_model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::ReducExposedPartialImmunity<double>>()[i],
                  read_model.parameters.get<mio::osecirvvs::ReducExposedPartialImmunity<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::ReducExposedImprovedImmunity<double>>()[i],
                  read_model.parameters.get<mio::osecirvvs::ReducExposedImprovedImmunity<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<double>>()[i],
                  read_model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<double>>()[i],
                  read_model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<double>>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::ReducTimeInfectedMild<double>>()[i],
                  read_model.parameters.get<mio::osecirvvs::ReducTimeInfectedMild<double>>()[i]);

        ASSERT_THAT(contact.get_cont_freq_mat(), testing::ContainerEq(read_contact.get_cont_freq_mat()));
        ASSERT_EQ(contact.get_dampings(), read_contact.get_dampings());
    }
}
