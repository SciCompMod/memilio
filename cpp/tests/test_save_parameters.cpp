/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

    mio::osecir::Model model(2);
    mio::AgeGroup num_groups = model.parameters.get_num_groups();
    double fact              = 1.0 / (double)(size_t)num_groups;

    auto& params = model.parameters;

    for (auto i = mio::AgeGroup(0); i < num_groups; i++) {
        params.get<mio::osecir::IncubationTime>()[i]       = 5.2;
        params.get<mio::osecir::TimeInfectedSymptoms>()[i] = 5.;
        params.get<mio::osecir::SerialInterval>()[i]       = 4.2;
        params.get<mio::osecir::TimeInfectedSevere>()[i]   = 10.;
        params.get<mio::osecir::TimeInfectedCritical>()[i] = 8.;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]            = fact * num_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}] = fact * num_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]   = fact * num_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]     = fact * num_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]   = fact * num_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]          = fact * num_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]               = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * num_total_t0);

        model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[i] = 0.06;
        model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[i]   = 0.67;
        model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i]   = alpha;
        model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i]   = beta;
        model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[i]        = rho;
        model.parameters.get<mio::osecir::CriticalPerSevere>()[i]                = theta;
        model.parameters.get<mio::osecir::DeathsPerCritical>()[i]                = delta;
    }

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, fact * cont_freq));
    contact_matrix.add_damping(0.7, mio::SimulationTime(30.));
    auto damping2  = Eigen::MatrixXd::Zero((size_t)num_groups, (size_t)num_groups).eval();
    damping2(0, 0) = 0.8;
    contact_matrix.add_damping(damping2, mio::SimulationTime(35));

    mio::osecir::set_params_distributions_normal(model, t0, tmax, 0.2);

    params.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0].get_distribution()->add_predefined_sample(4711.0);

    TempFileRegister file_register;
    auto filename     = file_register.get_unique_path("TestParameters-%%%%-%%%%.json");
    auto write_status = mio::write_json(filename, model);
    ASSERT_THAT(print_wrap(write_status), IsSuccess());

    auto read_result = mio::read_json(filename, mio::Tag<mio::osecir::Model>{});
    ASSERT_THAT(print_wrap(read_result), IsSuccess());
    auto& read_model = read_result.value();

    const mio::UncertainContactMatrix& contact      = model.parameters.get<mio::osecir::ContactPatterns>();
    const mio::UncertainContactMatrix& read_contact = read_model.parameters.get<mio::osecir::ContactPatterns>();

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

        ASSERT_EQ(model.parameters.get<mio::osecir::IncubationTime>()[i],
                  read_model.parameters.get<mio::osecir::IncubationTime>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[i],
                  read_model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::SerialInterval>()[i],
                  read_model.parameters.get<mio::osecir::SerialInterval>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::TimeInfectedSevere>()[i],
                  read_model.parameters.get<mio::osecir::TimeInfectedSevere>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::TimeInfectedCritical>()[i],
                  read_model.parameters.get<mio::osecir::TimeInfectedCritical>()[i]);

        check_distribution(*model.parameters.get<mio::osecir::IncubationTime>()[i].get_distribution(),
                           *read_model.parameters.get<mio::osecir::IncubationTime>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[i].get_distribution(),
                           *read_model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::osecir::SerialInterval>()[i].get_distribution(),
                           *read_model.parameters.get<mio::osecir::SerialInterval>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::osecir::TimeInfectedSevere>()[i].get_distribution(),
                           *read_model.parameters.get<mio::osecir::TimeInfectedSevere>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::osecir::TimeInfectedCritical>()[i].get_distribution(),
                           *read_model.parameters.get<mio::osecir::TimeInfectedCritical>()[i].get_distribution());

        ASSERT_EQ(model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[i],
                  read_model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i],
                  read_model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i],
                  read_model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::DeathsPerCritical>()[i],
                  read_model.parameters.get<mio::osecir::DeathsPerCritical>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[i],
                  read_model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecir::CriticalPerSevere>()[i],
                  read_model.parameters.get<mio::osecir::CriticalPerSevere>()[i]);

        check_distribution(
            *model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[i].get_distribution(),
            *read_model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[i].get_distribution());
        check_distribution(
            *model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i].get_distribution(),
            *read_model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i].get_distribution());
        check_distribution(
            *model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i].get_distribution(),
            *read_model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::osecir::DeathsPerCritical>()[i].get_distribution(),
                           *read_model.parameters.get<mio::osecir::DeathsPerCritical>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[i].get_distribution(),
                           *read_model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::osecir::CriticalPerSevere>()[i].get_distribution(),
                           *read_model.parameters.get<mio::osecir::CriticalPerSevere>()[i].get_distribution());

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
    mio::osecir::Model model((int)num_groups);
    double fact = 1.0 / (double)num_groups;

    auto& params = model.parameters;
    for (auto i = mio::Index<mio::AgeGroup>(0); i.get() < num_groups; i++) {
        params.get<mio::osecir::IncubationTime>()[i]       = 5.2;
        params.get<mio::osecir::TimeInfectedSymptoms>()[i] = 5.;
        params.get<mio::osecir::SerialInterval>()[i]       = 4.2;
        params.get<mio::osecir::TimeInfectedSevere>()[i]   = 10.;
        params.get<mio::osecir::TimeInfectedCritical>()[i] = 8.;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]            = fact * num_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}] = fact * num_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]   = fact * num_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]     = fact * num_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]   = fact * num_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]          = fact * num_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]               = fact * num_dead_t0;
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

    auto graph = mio::Graph<mio::osecir::Model, mio::MigrationParameters>();
    graph.add_node(0, model);
    graph.add_node(1, model);
    graph.add_edge(0, 1, mio::MigrationParameters(Eigen::VectorXd::Constant(Eigen::Index(num_groups * 8), 1.0)));

    TempFileRegister tmp_file_register;
    std::string tmp_results_dir = tmp_file_register.get_unique_path();
    ASSERT_THAT(mio::create_directory(tmp_results_dir), IsSuccess());

    std::vector<mio::osecir::Model> models = {model, model};
    std::vector<int> ids                   = {0, 1};
    auto graph_no_edges = mio::create_graph_without_edges<mio::osecir::Model, mio::MigrationParameters>(models, ids);
    auto write_status   = mio::write_graph(graph_no_edges, tmp_results_dir);
    ASSERT_THAT(print_wrap(write_status), IsSuccess());

    auto read_graph = mio::read_graph<mio::osecir::Model>(tmp_results_dir, mio::IOF_OmitDistributions, false);

    for (auto i = mio::AgeGroup(0); i < params.get_num_groups(); i++) {
        EXPECT_EQ(read_graph.value().nodes()[0].property.parameters.get<mio::osecir::IncubationTime>()[i],
                  params.get<mio::osecir::IncubationTime>()[i]);

        EXPECT_EQ(read_graph.value().nodes()[0].property.parameters.get<mio::osecir::TimeInfectedSymptoms>()[i],
                  params.get<mio::osecir::TimeInfectedSymptoms>()[i]);

        EXPECT_EQ(read_graph.value().nodes()[0].property.parameters.get<mio::osecir::SerialInterval>()[i],
                  params.get<mio::osecir::SerialInterval>()[i]);

        EXPECT_EQ(read_graph.value().nodes()[0].property.parameters.get<mio::osecir::TimeInfectedSevere>()[i],
                  params.get<mio::osecir::TimeInfectedSevere>()[i]);

        EXPECT_EQ(read_graph.value().nodes()[0].property.parameters.get<mio::osecir::TimeInfectedCritical>()[i],
                  params.get<mio::osecir::TimeInfectedCritical>()[i]);

        EXPECT_EQ(
            read_graph.value().nodes()[0].property.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[i],
            params.get<mio::osecir::TransmissionProbabilityOnContact>()[i]);

        EXPECT_EQ(
            read_graph.value().nodes()[0].property.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[i],
            params.get<mio::osecir::RelativeTransmissionNoSymptoms>()[i]);

        EXPECT_EQ(
            read_graph.value().nodes()[0].property.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i],
            params.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i]);

        EXPECT_EQ(read_graph.value().nodes()[0].property.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[i],
                  params.get<mio::osecir::SeverePerInfectedSymptoms>()[i]);

        EXPECT_EQ(read_graph.value().nodes()[0].property.parameters.get<mio::osecir::DeathsPerCritical>()[i],
                  params.get<mio::osecir::DeathsPerCritical>()[i]);
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
    auto damping_time1    = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 12, 1), start_date));
    auto damping_time2    = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 12, 24), start_date));
    mio::osecir::Model model(2);
    auto& contacts         = model.parameters.get<mio::osecir::ContactPatterns>();
    auto& contact_dampings = contacts.get_dampings();

    auto group_weights = Eigen::VectorXd::Constant(size_t(model.parameters.get_num_groups()), 1.0);

    //add two dampings
    auto v1 = mio::UncertainValue(0.5 * (0.6 + 0.4));
    v1.set_distribution(mio::ParameterDistributionUniform(0.4, 0.6));
    auto v2 = mio::UncertainValue(0.5 * (0.3 + 0.2));
    v1.set_distribution(mio::ParameterDistributionUniform(0.2, 0.3));
    contact_dampings.push_back(mio::DampingSampling(
        v1, mio::DampingLevel(int(InterventionLevelMock::level)), mio::DampingType(int(InterventionMock::intervention)),
        damping_time1, {size_t(ContactLocationMock::location)}, group_weights));
    contact_dampings.push_back(mio::DampingSampling(
        v2, mio::DampingLevel(int(InterventionLevelMock::level)), mio::DampingType(int(InterventionMock::intervention)),
        damping_time2, {size_t(ContactLocationMock::location)}, group_weights));
    //add school_holiday_damping
    auto school_holiday_value = mio::UncertainValue(1);
    school_holiday_value.set_distribution(mio::ParameterDistributionUniform(1, 1));
    contacts.get_school_holiday_damping() =
        mio::DampingSampling(school_holiday_value, mio::DampingLevel(int(InterventionLevelMock::level)),
                             mio::DampingType(int(InterventionMock::intervention)), mio::SimulationTime(0.0),
                             {size_t(ContactLocationMock::location)}, group_weights);

    //add school holidays
    auto holiday_start_time        = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 12, 23), start_date));
    auto holiday_end_time          = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2021, 1, 6), start_date));
    contacts.get_school_holidays() = {std::make_pair(holiday_start_time, holiday_end_time)};

    //write json
    TempFileRegister file_register;
    auto filename     = file_register.get_unique_path("TestParameters-%%%%-%%%%.json");
    auto write_status = mio::write_json(filename, model);
    ASSERT_THAT(print_wrap(write_status), IsSuccess());

    //read json
    auto read_result = mio::read_json(filename, mio::Tag<mio::osecir::Model>{});
    ASSERT_THAT(print_wrap(read_result), IsSuccess());
    auto& read_model = read_result.value();

    auto& read_contacts = read_model.parameters.get<mio::osecir::ContactPatterns>();

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

    mio::osecir::Model model(2);
    mio::AgeGroup num_groups = model.parameters.get_num_groups();
    double fact              = 1.0 / (double)(size_t)num_groups;

    model.parameters.set<mio::osecir::TestAndTraceCapacity>(30);

    for (auto i = mio::AgeGroup(0); i < num_groups; i++) {
        model.parameters.get<mio::osecir::IncubationTime>()[i]       = 5.2;
        model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[i] = 5.;
        model.parameters.get<mio::osecir::SerialInterval>()[i]       = 4.2;
        model.parameters.get<mio::osecir::TimeInfectedSevere>()[i]   = 10.;
        model.parameters.get<mio::osecir::TimeInfectedCritical>()[i] = 8.;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]            = fact * num_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}] = fact * num_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]   = fact * num_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]     = fact * num_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]   = fact * num_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]          = fact * num_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]               = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * num_total_t0);

        model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[i]  = 0.06;
        model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[i]    = 0.67;
        model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i]    = alpha;
        model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i]    = beta;
        model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>()[i] = beta * 3;
        model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[i]         = rho;
        model.parameters.get<mio::osecir::CriticalPerSevere>()[i]                 = theta;
        model.parameters.get<mio::osecir::DeathsPerCritical>()[i]                 = delta;
    }

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, fact * cont_freq));
    Eigen::MatrixXd m =
        Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, 0.7).triangularView<Eigen::Upper>();
    contact_matrix.add_damping(m, mio::SimulationTime(30.));

    mio::osecir::set_params_distributions_normal(model, t0, tmax, 0.15);

    mio::Graph<mio::osecir::Model, mio::MigrationParameters> graph;
    graph.add_node(0, model);
    graph.add_node(1, model);
    graph.add_edge(0, 1, Eigen::VectorXd::Constant(model.populations.get_num_compartments(), 0.01));
    graph.add_edge(1, 0, Eigen::VectorXd::Constant(model.populations.get_num_compartments(), 0.01));

    TempFileRegister file_register;
    auto graph_dir    = file_register.get_unique_path("graph_parameters-%%%%-%%%%");
    auto write_status = mio::write_graph(graph, graph_dir);
    ASSERT_THAT(print_wrap(write_status), IsSuccess());

    auto read_result = mio::read_graph<mio::osecir::Model>(graph_dir);
    ASSERT_THAT(print_wrap(read_result), IsSuccess());

    auto& graph_read = read_result.value();
    auto num_nodes   = graph.nodes().size();
    auto num_edges   = graph.edges().size();

    ASSERT_EQ(num_nodes, graph_read.nodes().size());
    ASSERT_EQ(num_edges, graph_read.edges().size());

    for (size_t node = 0; node < num_nodes; node++) {
        mio::osecir::Model graph_model             = graph.nodes()[0].property;
        mio::ContactMatrixGroup& graph_cont_matrix = graph_model.parameters.get<mio::osecir::ContactPatterns>();

        mio::osecir::Model graph_read_model = graph_read.nodes()[0].property;
        mio::ContactMatrixGroup& graph_read_cont_matrix =
            graph_read_model.parameters.get<mio::osecir::ContactPatterns>();

        ASSERT_EQ(graph_read_cont_matrix.get_num_groups(), static_cast<Eigen::Index>((size_t)num_groups));
        ASSERT_EQ(graph_read_cont_matrix, graph_cont_matrix);
        ASSERT_EQ(graph_model.populations.get_num_compartments(), graph_read_model.populations.get_num_compartments());
        ASSERT_EQ(graph.nodes()[node].id, graph_read.nodes()[node].id);
        EXPECT_THAT(
            graph_read_model.parameters.get<mio::osecir::TestAndTraceCapacity>().value(),
            FloatingPointEqual(graph_model.parameters.get<mio::osecir::TestAndTraceCapacity>().value(), 1e-12, 1e-12));
        check_distribution(
            *graph_model.parameters.get<mio::osecir::TestAndTraceCapacity>().get_distribution().get(),
            *graph_read_model.parameters.get<mio::osecir::TestAndTraceCapacity>().get_distribution().get());

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

            ASSERT_EQ(graph_model.parameters.get<mio::osecir::IncubationTime>()[group],
                      graph_read_model.parameters.get<mio::osecir::IncubationTime>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[group],
                      graph_read_model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::SerialInterval>()[group],
                      graph_read_model.parameters.get<mio::osecir::SerialInterval>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::TimeInfectedSevere>()[group],
                      graph_read_model.parameters.get<mio::osecir::TimeInfectedSevere>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::TimeInfectedCritical>()[group],
                      graph_read_model.parameters.get<mio::osecir::TimeInfectedCritical>()[group]);

            ASSERT_EQ(graph_model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[group],
                      graph_read_model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[group],
                      graph_read_model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>()[group],
                      graph_read_model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[group],
                      graph_read_model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::DeathsPerCritical>()[group],
                      graph_read_model.parameters.get<mio::osecir::DeathsPerCritical>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[group],
                      graph_read_model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::osecir::CriticalPerSevere>()[group],
                      graph_read_model.parameters.get<mio::osecir::CriticalPerSevere>()[group]);

            check_distribution(
                *graph_model.parameters.get<mio::osecir::IncubationTime>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::osecir::IncubationTime>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<mio::osecir::SerialInterval>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::osecir::SerialInterval>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<mio::osecir::TimeInfectedSevere>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::osecir::TimeInfectedSevere>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<mio::osecir::TimeInfectedCritical>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::osecir::TimeInfectedCritical>()[group].get_distribution().get());

            check_distribution(
                *graph_model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<mio::osecir::TimeInfectedSevere>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::osecir::TimeInfectedSevere>()[group].get_distribution().get());
            check_distribution(*graph_model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>()[group]
                                    .get_distribution()
                                    .get(),
                               *graph_read_model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>()[group]
                                    .get_distribution()
                                    .get());
            check_distribution(
                *graph_model.parameters.get<mio::osecir::DeathsPerCritical>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::osecir::DeathsPerCritical>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<mio::osecir::CriticalPerSevere>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::osecir::CriticalPerSevere>()[group].get_distribution().get());

            ASSERT_EQ(graph_model.parameters.get<mio::osecir::ContactPatterns>().get_dampings(),
                      graph_read_model.parameters.get<mio::osecir::ContactPatterns>().get_dampings());
        }

        ASSERT_THAT(graph_read.edges(), testing::ElementsAreArray(graph.edges()));
    }
}

TEST(TestSaveParameters, ReadPopulationDataRKIAges)
{
    std::vector<mio::osecir::Model> model(1, {6});
    model[0].apply_constraints();
    std::vector<double> scaling_factor_inf(6, 1.0);
    double scaling_factor_icu = 1.0;
    mio::Date date(2020, 12, 10);

    std::string path = TEST_DATA_DIR;

    for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(6); group++) {
        model[0].parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[group] = 0.1 * ((size_t)group + 1);
        model[0].parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[group]      = 0.11 * ((size_t)group + 1);
        model[0].parameters.get<mio::osecir::CriticalPerSevere>()[group]              = 0.12 * ((size_t)group + 1);
    }
    auto read_result = mio::osecir::read_input_data_germany(model, date, scaling_factor_inf, scaling_factor_icu, path);
    ASSERT_THAT(print_wrap(read_result), IsSuccess());

    std::vector<double> sus   = {3443857.42, 7665093.95, 18792870.93, 29503629.76, 16307262.45, 6049150.54};
    std::vector<double> exp   = {433.015, 1771.61, 8856.33, 14757.62, 7222.86, 6626.07};
    std::vector<double> car   = {434.444, 1772.14, 8724.49, 14386.90, 6995.14, 6307.14};
    std::vector<double> inf   = {375.429, 1393.43, 6007.14, 8438.71, 3377.57, 2421.57};
    std::vector<double> hosp  = {39.9614, 303.191, 1934.84, 3621.2, 1793.39, 1557.03};
    std::vector<double> icu   = {47.6813, 190.725, 429.132, 762.901, 1192.03, 1716.53};
    std::vector<double> rec   = {23557.7, 78946.3, 398585.142, 487273.71, 178660.14, 96021.9};
    std::vector<double> death = {2, 4, 48, 1137.86, 8174.14, 18528.9};

    for (size_t i = 0; i < 6; i++) {
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Susceptible}]), sus[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Exposed}]), exp[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedNoSymptoms}]), car[i],
                    1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedSymptoms}]), inf[i],
                    1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedSevere}]), hosp[i],
                    1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedCritical}]), icu[i],
                    1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Recovered}]), rec[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Dead}]), death[i], 1e-1);
    }

    EXPECT_NEAR(model[0].populations.get_total(), 83166695, 1e-6);
}

TEST(TestSaveParameters, ReadPopulationDataStateAllAges)
{
    std::vector<mio::osecir::Model> model(1, {6});
    model[0].apply_constraints();
    std::vector<double> scaling_factor_inf(6, 1.0);
    double scaling_factor_icu = 1.0;
    mio::Date date(2020, 12, 10);

    std::vector<int> state = {1};

    std::string path = TEST_DATA_DIR;

    for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(6); group++) {
        model[0].parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[group] = 0.1 * ((size_t)group + 1);
        model[0].parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[group]      = 0.11 * ((size_t)group + 1);
        model[0].parameters.get<mio::osecir::CriticalPerSevere>()[group]              = 0.12 * ((size_t)group + 1);
    }
    auto read_result =
        mio::osecir::read_input_data_state(model, date, state, scaling_factor_inf, scaling_factor_icu, path);
    ASSERT_THAT(print_wrap(read_result), IsSuccess());

    std::vector<double> sus   = {116692.2, 283912.8, 622795.86, 1042178.3, 606450.7, 212836.9};
    std::vector<double> exp   = {8.57143, 30.5357, 149.388, 228.809, 87.1429, 99.2857};
    std::vector<double> car   = {7.77778, 26.0714, 143.061, 217.143, 84.8571, 92.1429};
    std::vector<double> inf   = {7.00000, 18.7143, 97.7143, 122.000, 40.8571, 36.1429};
    std::vector<double> hosp  = {0.707143, 3.92857, 30.6429, 50.5371, 20.35, 19.9886};
    std::vector<double> icu   = {0.274725, 1.0989, 2.47253, 4.3956, 6.86813, 9.89011};
    std::vector<double> rec   = {393.143, 1216.14, 5467.86, 6543.57, 2281.29, 1045.71};
    std::vector<double> death = {0, 0, 0, 16.2857, 99.5714, 198.286};

    for (size_t i = 0; i < 6; i++) {
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Susceptible}]), sus[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Exposed}]), exp[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedNoSymptoms}]), car[i],
                    1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedSymptoms}]), inf[i],
                    1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedSevere}]), hosp[i],
                    1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedCritical}]), icu[i],
                    1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Recovered}]), rec[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Dead}]), death[i], 1e-1);
    }

    EXPECT_NEAR(model[0].populations.get_total(), 2903777, 1e-6);
}

TEST(TestSaveParameters, ReadPopulationDataCountyAllAges)
{

    std::vector<mio::osecir::Model> model1(1, {6});
    model1[0].apply_constraints();
    std::vector<mio::osecir::Model> model2(1, {6});
    model2[0].apply_constraints();
    std::vector<mio::osecir::Model> model3(1, {6});
    model3[0].apply_constraints();
    std::vector<double> scaling_factor_inf(6, 1.0);
    double scaling_factor_icu = 1.0;
    mio::Date date(2020, 12, 10);

    std::vector<int> county = {1002};

    std::string path = TEST_DATA_DIR;

    for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(6); group++) {
        model1[0].parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[group] = 0.1 * ((size_t)group + 1);
        model1[0].parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[group]      = 0.11 * ((size_t)group + 1);
        model1[0].parameters.get<mio::osecir::CriticalPerSevere>()[group]              = 0.12 * ((size_t)group + 1);
        model2[0].parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[group] = 0.1 * ((size_t)group + 1);
        model2[0].parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[group]      = 0.11 * ((size_t)group + 1);
        model2[0].parameters.get<mio::osecir::CriticalPerSevere>()[group]              = 0.12 * ((size_t)group + 1);
        model3[0].parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[group] = 0.1 * ((size_t)group + 1);
        model3[0].parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[group]      = 0.11 * ((size_t)group + 1);
        model3[0].parameters.get<mio::osecir::CriticalPerSevere>()[group]              = 0.12 * ((size_t)group + 1);
    }
    auto read_result1 =
        mio::osecir::read_input_data_county(model1, date, county, scaling_factor_inf, scaling_factor_icu, path);
    auto read_result2 =
        mio::osecir::read_input_data(model2, date, county, scaling_factor_inf, scaling_factor_icu, path);
    auto read_result_district = mio::osecir::read_input_data(
        model3, date, county, scaling_factor_inf, scaling_factor_icu, mio::path_join(path, "pydata/District"));
    ASSERT_THAT(print_wrap(read_result1), IsSuccess());
    ASSERT_THAT(print_wrap(read_result2), IsSuccess());
    ASSERT_THAT(print_wrap(read_result_district), IsSuccess());

    std::vector<double> sus   = {10284.13, 19082.86, 73783.12, 82494.81, 43725.08, 15612.70};
    std::vector<double> exp   = {0.571429, 4.82143, 20.8163, 22.1429, 4.57143, 4.64286};
    std::vector<double> car   = {0.557143, 4.46429, 22.0408, 20.7143, 4.28571, 4.64286};
    std::vector<double> inf   = {0.42857, 3.285714, 15.2857, 13.0000, 2.42857, 2.00000};
    std::vector<double> hosp  = {0.0942857, 0.691429, 4.90286, 5.34286, 1.41429, 2.45143};
    std::vector<double> icu   = {0.0769231, 0.307692, 0.692308, 1.23077, 1.92308, 2.76923};
    std::vector<double> rec   = {35, 108.571, 640.143, 573.429, 180.429, 75.5714};
    std::vector<double> death = {0, 0, 0, 0, 10, 14.4286};

    for (size_t i = 0; i < 6; i++) {
        EXPECT_NEAR((model1[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Susceptible}]), sus[i],
                    1e-1);
        EXPECT_NEAR((model2[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Susceptible}]), sus[i],
                    1e-1);
        EXPECT_NEAR((model3[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Susceptible}]), sus[i],
                    1e-1);
        EXPECT_NEAR((model1[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Exposed}]), exp[i], 1e-1);
        EXPECT_NEAR((model2[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Exposed}]), exp[i], 1e-1);
        EXPECT_NEAR((model3[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Exposed}]), exp[i], 1e-1);
        EXPECT_NEAR((model1[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedNoSymptoms}]),
                    car[i], 1e-1);
        EXPECT_NEAR((model2[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedNoSymptoms}]),
                    car[i], 1e-1);
        EXPECT_NEAR((model3[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedNoSymptoms}]),
                    car[i], 1e-1);
        EXPECT_NEAR((model1[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedSymptoms}]), inf[i],
                    1e-1);
        EXPECT_NEAR((model2[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedSymptoms}]), inf[i],
                    1e-1);
        EXPECT_NEAR((model3[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedSymptoms}]), inf[i],
                    1e-1);
        EXPECT_NEAR((model1[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedSevere}]), hosp[i],
                    1e-1);
        EXPECT_NEAR((model2[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedSevere}]), hosp[i],
                    1e-1);
        EXPECT_NEAR((model3[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedSevere}]), hosp[i],
                    1e-1);
        EXPECT_NEAR((model1[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedCritical}]), icu[i],
                    1e-1);
        EXPECT_NEAR((model2[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedCritical}]), icu[i],
                    1e-1);
        EXPECT_NEAR((model3[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::InfectedCritical}]), icu[i],
                    1e-1);
        EXPECT_NEAR((model1[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Recovered}]), rec[i], 1e-1);
        EXPECT_NEAR((model2[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Recovered}]), rec[i], 1e-1);
        EXPECT_NEAR((model3[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Recovered}]), rec[i], 1e-1);
        EXPECT_NEAR((model1[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Dead}]), death[i], 1e-1);
        EXPECT_NEAR((model2[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Dead}]), death[i], 1e-1);
        EXPECT_NEAR((model3[0].populations[{mio::AgeGroup(i), mio::osecir::InfectionState::Dead}]), death[i], 1e-1);
    }

    EXPECT_NEAR(model1[0].populations.get_total(), 246793, 1e-6);
    EXPECT_NEAR(model2[0].populations.get_total(), 246793, 1e-6);
    EXPECT_NEAR(model3[0].populations.get_total(), 246793, 1e-6);
}

TEST(TestSaveParameters, ExtrapolateRKI)
{
    std::vector<mio::osecir::Model> model{mio::osecir::Model(6)};

    model[0].apply_constraints();
    std::vector<double> scaling_factor_inf(6, 1.0);
    double scaling_factor_icu = 1.0;
    mio::Date date(2020, 12, 10);

    std::vector<int> county = {1002};

    for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(6); group++) {
        model[0].parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[group] = 0.1 * ((size_t)group + 1);
        model[0].parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[group]      = 0.11 * ((size_t)group + 1);
        model[0].parameters.get<mio::osecir::CriticalPerSevere>()[group]              = 0.12 * ((size_t)group + 1);
    }

    TempFileRegister file_register;
    auto results_dir = file_register.get_unique_path("ExtrapolateRKI-%%%%-%%%%");
    boost::filesystem::create_directory(results_dir);
    auto extrapolate_result = mio::osecir::export_input_data_county_timeseries(
        model, results_dir, county, date, scaling_factor_inf, scaling_factor_icu, 1,
        mio::path_join(TEST_DATA_DIR, "county_divi_ma7.json"),
        mio::path_join(TEST_DATA_DIR, "cases_all_county_age_ma7.json"),
        mio::path_join(TEST_DATA_DIR, "county_current_population.json"));
    ASSERT_THAT(print_wrap(extrapolate_result), IsSuccess());

    auto read_result = mio::read_result(mio::path_join(results_dir, "Results_rki.h5"));
    ASSERT_THAT(print_wrap(read_result), IsSuccess());
    auto& file_results = read_result.value();
    auto results       = file_results[0].get_groups();

    std::vector<double> sus   = {10284.1, 19082.9, 73783.1, 82494.8, 43725.1, 15612.7};
    std::vector<double> exp   = {0.571429, 4.82143, 20.8163, 22.1429, 4.57143, 4.64286};
    std::vector<double> car   = {0.557143, 4.46429, 22.0408, 20.7143, 4.28571, 4.64286};
    std::vector<double> inf   = {0.428571, 3.28571, 15.2857, 13.0000, 2.42857, 2.00000};
    std::vector<double> hosp  = {0.0942857, 0.691429, 4.90286, 5.34286, 1.41429, 2.45143};
    std::vector<double> icu   = {0.0769231, 0.307692, 0.692308, 1.23077, 1.92308, 2.76923};
    std::vector<double> rec   = {35, 108.571, 640.143, 573.429, 180.429, 75.5714};
    std::vector<double> death = {0, 0, 0, 0, 10, 14.4286};

    for (size_t i = 0; i < 6; i++) {
        EXPECT_NEAR(results[0]((size_t)mio::osecir::InfectionState::Susceptible +
                               (size_t)mio::osecir::InfectionState::Count * i),
                    sus[i], 1e-1);
        EXPECT_NEAR(
            results[0]((size_t)mio::osecir::InfectionState::Exposed + (size_t)mio::osecir::InfectionState::Count * i),
            exp[i], 1e-1);
        EXPECT_NEAR(results[0]((size_t)mio::osecir::InfectionState::InfectedNoSymptoms +
                               (size_t)mio::osecir::InfectionState::Count * i),
                    car[i], 1e-1);
        EXPECT_NEAR(results[0]((size_t)mio::osecir::InfectionState::InfectedSymptoms +
                               (size_t)mio::osecir::InfectionState::Count * i),
                    inf[i], 1e-1);
        EXPECT_NEAR(results[0]((size_t)mio::osecir::InfectionState::InfectedSevere +
                               (size_t)mio::osecir::InfectionState::Count * i),
                    hosp[i], 1e-1);
        EXPECT_NEAR(results[0]((size_t)mio::osecir::InfectionState::InfectedCritical +
                               (size_t)mio::osecir::InfectionState::Count * i),
                    icu[i], 1e-1);
        EXPECT_NEAR(
            results[0]((size_t)mio::osecir::InfectionState::Recovered + (size_t)mio::osecir::InfectionState::Count * i),
            rec[i], 1e-1);
        EXPECT_NEAR(
            results[0]((size_t)mio::osecir::InfectionState::Dead + (size_t)mio::osecir::InfectionState::Count * i),
            death[i], 1e-1);
    }
}

TEST(TestSaveParameters, json_write_read_parameters_secirvvs)
{
    mio::osecirvvs::Model model(2);
    mio::AgeGroup num_groups = model.parameters.get_num_groups();
    double fact              = 1.0 / (double)(size_t)num_groups;

    auto& params = model.parameters;

    for (auto i = mio::AgeGroup(0); i < num_groups; i++) {
        params.get<mio::osecirvvs::IncubationTime>()[i]       = 5.2;
        params.get<mio::osecirvvs::TimeInfectedSymptoms>()[i] = 5.;
        params.get<mio::osecirvvs::SerialInterval>()[i]       = 4.2;
        params.get<mio::osecirvvs::TimeInfectedSevere>()[i]   = 10.;
        params.get<mio::osecirvvs::TimeInfectedCritical>()[i] = 8.;

        model.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>()[i] = 0.06;
        model.parameters.get<mio::osecirvvs::RelativeTransmissionNoSymptoms>()[i]   = 0.67;
        model.parameters.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms>()[i]   = 0.09;
        model.parameters.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic>()[i]   = 0.25;
        model.parameters.get<mio::osecirvvs::SeverePerInfectedSymptoms>()[i]        = 0.2;
        model.parameters.get<mio::osecirvvs::CriticalPerSevere>()[i]                = 0.25;
        model.parameters.get<mio::osecirvvs::DeathsPerCritical>()[i]                = 0.3;

        model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity>()[i]  = 0.8;
        model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity>()[i] = 0.7;
        model.parameters.get<mio::osecirvvs::ReducExposedPartialImmunity>()[i]                     = 0.8;
        model.parameters.get<mio::osecirvvs::ReducExposedImprovedImmunity>()[i]                    = 0.7;
        model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity>()[i]            = 0.8;
        model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity>()[i]           = 0.7;
        model.parameters.get<mio::osecirvvs::ReducTimeInfectedMild>()[i]                           = 0.5;
    }

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecirvvs::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, fact * 10));
    contact_matrix.add_damping(0.7, mio::SimulationTime(30.));
    auto damping2  = Eigen::MatrixXd::Zero((size_t)num_groups, (size_t)num_groups).eval();
    damping2(0, 0) = 0.8;
    contact_matrix.add_damping(damping2, mio::SimulationTime(35));

    TempFileRegister file_register;
    auto filename     = file_register.get_unique_path("TestParameters-%%%%-%%%%.json");
    auto write_status = mio::write_json(filename, model);
    ASSERT_THAT(print_wrap(write_status), IsSuccess());

    auto read_result = mio::read_json(filename, mio::Tag<mio::osecirvvs::Model>{});
    ASSERT_THAT(print_wrap(read_result), IsSuccess());
    auto& read_model = read_result.value();

    const mio::UncertainContactMatrix& contact      = model.parameters.get<mio::osecirvvs::ContactPatterns>();
    const mio::UncertainContactMatrix& read_contact = read_model.parameters.get<mio::osecirvvs::ContactPatterns>();

    num_groups           = model.parameters.get_num_groups();
    auto num_groups_read = read_model.parameters.get_num_groups();
    ASSERT_EQ(num_groups, num_groups_read);

    for (auto i = mio::AgeGroup(0); i < num_groups; i++) {
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::IncubationTime>()[i],
                  read_model.parameters.get<mio::osecirvvs::IncubationTime>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedSymptoms>()[i],
                  read_model.parameters.get<mio::osecirvvs::TimeInfectedSymptoms>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::SerialInterval>()[i],
                  read_model.parameters.get<mio::osecirvvs::SerialInterval>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedSevere>()[i],
                  read_model.parameters.get<mio::osecirvvs::TimeInfectedSevere>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedCritical>()[i],
                  read_model.parameters.get<mio::osecirvvs::TimeInfectedCritical>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>()[i],
                  read_model.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic>()[i],
                  read_model.parameters.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms>()[i],
                  read_model.parameters.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::DeathsPerCritical>()[i],
                  read_model.parameters.get<mio::osecirvvs::DeathsPerCritical>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::SeverePerInfectedSymptoms>()[i],
                  read_model.parameters.get<mio::osecirvvs::SeverePerInfectedSymptoms>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::CriticalPerSevere>()[i],
                  read_model.parameters.get<mio::osecirvvs::CriticalPerSevere>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity>()[i],
                  read_model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity>()[i],
                  read_model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::ReducExposedPartialImmunity>()[i],
                  read_model.parameters.get<mio::osecirvvs::ReducExposedPartialImmunity>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::ReducExposedImprovedImmunity>()[i],
                  read_model.parameters.get<mio::osecirvvs::ReducExposedImprovedImmunity>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity>()[i],
                  read_model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity>()[i],
                  read_model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity>()[i]);
        ASSERT_EQ(model.parameters.get<mio::osecirvvs::ReducTimeInfectedMild>()[i],
                  read_model.parameters.get<mio::osecirvvs::ReducTimeInfectedMild>()[i]);

        ASSERT_THAT(contact.get_cont_freq_mat(), testing::ContainerEq(read_contact.get_cont_freq_mat()));
        ASSERT_EQ(contact.get_dampings(), read_contact.get_dampings());
    }
}
