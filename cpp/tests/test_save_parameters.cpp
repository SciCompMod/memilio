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
#include "secir/secir.h"
#include "secir/parameter_space.h"
#include "secir/secir_parameters_io.h"
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

    mio::SecirModel model(2);
    mio::AgeGroup num_groups = model.parameters.get_num_groups();
    double fact              = 1.0 / (double)(size_t)num_groups;

    auto& params = model.parameters;

    for (auto i = mio::AgeGroup(0); i < num_groups; i++) {
        params.get<mio::IncubationTime>()[i]       = 5.2;
        params.get<mio::TimeInfectedSymptoms>()[i] = 5.;
        params.get<mio::SerialInterval>()[i]       = 4.2;
        params.get<mio::TimeInfectedSevere>()[i]   = 10.;
        params.get<mio::TimeInfectedCritical>()[i] = 8.;

        model.populations[{i, mio::InfectionState::Exposed}]            = fact * num_exp_t0;
        model.populations[{i, mio::InfectionState::InfectedNoSymptoms}] = fact * num_car_t0;
        model.populations[{i, mio::InfectionState::InfectedSymptoms}]   = fact * num_inf_t0;
        model.populations[{i, mio::InfectionState::InfectedSevere}]     = fact * num_hosp_t0;
        model.populations[{i, mio::InfectionState::InfectedCritical}]   = fact * num_icu_t0;
        model.populations[{i, mio::InfectionState::Recovered}]          = fact * num_rec_t0;
        model.populations[{i, mio::InfectionState::Dead}]               = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::InfectionState::Susceptible},
                                                                         fact * num_total_t0);

        model.parameters.get<mio::TransmissionProbabilityOnContact>()[i] = 0.06;
        model.parameters.get<mio::RelativeTransmissionNoSymptoms>()[i]   = 0.67;
        model.parameters.get<mio::RecoveredPerInfectedNoSymptoms>()[i]   = alpha;
        model.parameters.get<mio::RiskOfInfectionFromSymptomatic>()[i]   = beta;
        model.parameters.get<mio::SeverePerInfectedSymptoms>()[i]        = rho;
        model.parameters.get<mio::CriticalPerSevere>()[i]                = theta;
        model.parameters.get<mio::DeathsPerCritical>()[i]                = delta;
    }

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, fact * cont_freq));
    contact_matrix.add_damping(0.7, mio::SimulationTime(30.));
    auto damping2  = Eigen::MatrixXd::Zero((size_t)num_groups, (size_t)num_groups).eval();
    damping2(0, 0) = 0.8;
    contact_matrix.add_damping(damping2, mio::SimulationTime(35));

    mio::set_params_distributions_normal(model, t0, tmax, 0.2);

    params.get<mio::IncubationTime>()[(mio::AgeGroup)0].get_distribution()->add_predefined_sample(4711.0);

    TempFileRegister file_register;
    auto filename     = file_register.get_unique_path("TestParameters-%%%%-%%%%.json");
    auto write_status = mio::write_json(filename, model);
    ASSERT_THAT(print_wrap(write_status), IsSuccess());

    auto read_result = mio::read_json(filename, mio::Tag<mio::SecirModel>{});
    ASSERT_THAT(print_wrap(read_result), IsSuccess());
    auto& read_model = read_result.value();

    const mio::UncertainContactMatrix& contact      = model.parameters.get<mio::ContactPatterns>();
    const mio::UncertainContactMatrix& read_contact = read_model.parameters.get<mio::ContactPatterns>();

    num_groups           = model.parameters.get_num_groups();
    auto num_groups_read = read_model.parameters.get_num_groups();
    ASSERT_EQ(num_groups, num_groups_read);

    for (auto i = mio::AgeGroup(0); i < num_groups; i++) {
        ASSERT_EQ((model.populations[{i, mio::InfectionState::Dead}]),
                  (read_model.populations[{i, mio::InfectionState::Dead}]));
        ASSERT_EQ((model.populations.get_group_total(i)), (read_model.populations.get_group_total(i)));
        ASSERT_EQ((model.populations[{i, mio::InfectionState::Exposed}]),
                  (read_model.populations[{i, mio::InfectionState::Exposed}]));
        ASSERT_EQ((model.populations[{i, mio::InfectionState::InfectedNoSymptoms}]),
                  (read_model.populations[{i, mio::InfectionState::InfectedNoSymptoms}]));
        ASSERT_EQ((model.populations[{i, mio::InfectionState::InfectedSymptoms}]),
                  (read_model.populations[{i, mio::InfectionState::InfectedSymptoms}]));
        ASSERT_EQ((model.populations[{i, mio::InfectionState::InfectedSevere}]),
                  (read_model.populations[{i, mio::InfectionState::InfectedSevere}]));
        ASSERT_EQ((model.populations[{i, mio::InfectionState::InfectedCritical}]),
                  (read_model.populations[{i, mio::InfectionState::InfectedCritical}]));
        ASSERT_EQ((model.populations[{i, mio::InfectionState::Recovered}]),
                  (read_model.populations[{i, mio::InfectionState::Recovered}]));

        check_distribution(*model.populations[{i, mio::InfectionState::Exposed}].get_distribution(),
                           *read_model.populations[{i, mio::InfectionState::Exposed}].get_distribution());
        check_distribution(*model.populations[{i, mio::InfectionState::InfectedNoSymptoms}].get_distribution(),
                           *read_model.populations[{i, mio::InfectionState::InfectedNoSymptoms}].get_distribution());
        check_distribution(*model.populations[{i, mio::InfectionState::InfectedSymptoms}].get_distribution(),
                           *read_model.populations[{i, mio::InfectionState::InfectedSymptoms}].get_distribution());
        check_distribution(*model.populations[{i, mio::InfectionState::InfectedSevere}].get_distribution(),
                           *read_model.populations[{i, mio::InfectionState::InfectedSevere}].get_distribution());
        check_distribution(*model.populations[{i, mio::InfectionState::InfectedCritical}].get_distribution(),
                           *read_model.populations[{i, mio::InfectionState::InfectedCritical}].get_distribution());
        check_distribution(*model.populations[{i, mio::InfectionState::Recovered}].get_distribution(),
                           *read_model.populations[{i, mio::InfectionState::Recovered}].get_distribution());

        ASSERT_EQ(model.parameters.get<mio::IncubationTime>()[i], read_model.parameters.get<mio::IncubationTime>()[i]);
        ASSERT_EQ(model.parameters.get<mio::TimeInfectedSymptoms>()[i],
                  read_model.parameters.get<mio::TimeInfectedSymptoms>()[i]);
        ASSERT_EQ(model.parameters.get<mio::SerialInterval>()[i], read_model.parameters.get<mio::SerialInterval>()[i]);
        ASSERT_EQ(model.parameters.get<mio::TimeInfectedSevere>()[i],
                  read_model.parameters.get<mio::TimeInfectedSevere>()[i]);
        ASSERT_EQ(model.parameters.get<mio::TimeInfectedCritical>()[i],
                  read_model.parameters.get<mio::TimeInfectedCritical>()[i]);

        check_distribution(*model.parameters.get<mio::IncubationTime>()[i].get_distribution(),
                           *read_model.parameters.get<mio::IncubationTime>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::TimeInfectedSymptoms>()[i].get_distribution(),
                           *read_model.parameters.get<mio::TimeInfectedSymptoms>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::SerialInterval>()[i].get_distribution(),
                           *read_model.parameters.get<mio::SerialInterval>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::TimeInfectedSevere>()[i].get_distribution(),
                           *read_model.parameters.get<mio::TimeInfectedSevere>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::TimeInfectedCritical>()[i].get_distribution(),
                           *read_model.parameters.get<mio::TimeInfectedCritical>()[i].get_distribution());

        ASSERT_EQ(model.parameters.get<mio::TransmissionProbabilityOnContact>()[i],
                  read_model.parameters.get<mio::TransmissionProbabilityOnContact>()[i]);
        ASSERT_EQ(model.parameters.get<mio::RiskOfInfectionFromSymptomatic>()[i],
                  read_model.parameters.get<mio::RiskOfInfectionFromSymptomatic>()[i]);
        ASSERT_EQ(model.parameters.get<mio::RecoveredPerInfectedNoSymptoms>()[i],
                  read_model.parameters.get<mio::RecoveredPerInfectedNoSymptoms>()[i]);
        ASSERT_EQ(model.parameters.get<mio::DeathsPerCritical>()[i],
                  read_model.parameters.get<mio::DeathsPerCritical>()[i]);
        ASSERT_EQ(model.parameters.get<mio::SeverePerInfectedSymptoms>()[i],
                  read_model.parameters.get<mio::SeverePerInfectedSymptoms>()[i]);
        ASSERT_EQ(model.parameters.get<mio::CriticalPerSevere>()[i],
                  read_model.parameters.get<mio::CriticalPerSevere>()[i]);

        check_distribution(*model.parameters.get<mio::TransmissionProbabilityOnContact>()[i].get_distribution(),
                           *read_model.parameters.get<mio::TransmissionProbabilityOnContact>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::RiskOfInfectionFromSymptomatic>()[i].get_distribution(),
                           *read_model.parameters.get<mio::RiskOfInfectionFromSymptomatic>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::RecoveredPerInfectedNoSymptoms>()[i].get_distribution(),
                           *read_model.parameters.get<mio::RecoveredPerInfectedNoSymptoms>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::DeathsPerCritical>()[i].get_distribution(),
                           *read_model.parameters.get<mio::DeathsPerCritical>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::SeverePerInfectedSymptoms>()[i].get_distribution(),
                           *read_model.parameters.get<mio::SeverePerInfectedSymptoms>()[i].get_distribution());
        check_distribution(*model.parameters.get<mio::CriticalPerSevere>()[i].get_distribution(),
                           *read_model.parameters.get<mio::CriticalPerSevere>()[i].get_distribution());

        ASSERT_THAT(contact.get_cont_freq_mat(), testing::ContainerEq(read_contact.get_cont_freq_mat()));
        ASSERT_EQ(contact.get_dampings(), read_contact.get_dampings());
    }
}

TEST(TestSaveParameters, json_uncertain_matrix_write_read_compare)
{
    const auto start_date   = mio::Date(2020, 12, 12);
    const auto end_date     = mio::offset_date_by_days(start_date, int(std::ceil(20.0)));
    auto damping_time1      = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 1, 5), start_date));
    auto damping_time2      = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 12, 24), start_date));
    mio::SecirModel model(2);
    auto& contacts = model.parameters.get<mio::ContactPatterns>();
    auto& contact_dampings = contacts.get_dampings();

    auto group_weights     = Eigen::VectorXd::Constant(size_t(model.parameters.get_num_groups()), 1.0);

    //add two dampings
    auto v1 = mio::UncertainValue(0.5 * (0.6 + 0.4));
    v1.set_distribution(mio::ParameterDistributionUniform(0.4, 0.6));
    auto v2 = mio::UncertainValue(0.5 * (0.3 + 0.2));
    v1.set_distribution(mio::ParameterDistributionUniform(0.2, 0.3));
    contact_dampings.push_back(mio::DampingSampling(v1, mio::DampingLevel(0),
                                mio::DampingType(0), damping_time1, {size_t(0)}, group_weights));
    contact_dampings.push_back(mio::DampingSampling(v2, mio::DampingLevel(0),
                                mio::DampingType(0), damping_time2, {size_t(0)}, group_weights));
    //add school_holiday_damping
    auto school_holiday_value = mio::UncertainValue(1);
    school_holiday_value.set_distribution(mio::ParameterDistributionUniform(1, 1));
    contacts.get_school_holiday_damping() =
        mio::DampingSampling(school_holiday_value, mio::DampingLevel(0),
                             mio::DampingType(0), mio::SimulationTime(0.0),
                             {size_t(0)}, group_weights);

    //add school holidays
    auto holiday_start_time = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 12, 23), start_date));
    auto holiday_end_time = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2021, 1, 6), end_date));
    contacts.get_school_holidays() = {std::make_pair(holiday_start_time,holiday_end_time)};

    //write json
    TempFileRegister file_register;
    auto filename     = file_register.get_unique_path("TestParameters-%%%%-%%%%.json");
    auto write_status = mio::write_json(filename, model);
    ASSERT_THAT(print_wrap(write_status), IsSuccess());

    //read json
    auto read_result = mio::read_json(filename, mio::Tag<mio::SecirModel>{});
    ASSERT_THAT(print_wrap(read_result), IsSuccess());
    auto& read_model = read_result.value();

    auto& read_contacts = read_model.parameters.get<mio::ContactPatterns>();

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

    mio::SecirModel model(2);
    mio::AgeGroup num_groups = model.parameters.get_num_groups();
    double fact              = 1.0 / (double)(size_t)num_groups;

    model.parameters.set<mio::TestAndTraceCapacity>(30);

    for (auto i = mio::AgeGroup(0); i < num_groups; i++) {
        model.parameters.get<mio::IncubationTime>()[i]       = 5.2;
        model.parameters.get<mio::TimeInfectedSymptoms>()[i] = 5.;
        model.parameters.get<mio::SerialInterval>()[i]       = 4.2;
        model.parameters.get<mio::TimeInfectedSevere>()[i]   = 10.;
        model.parameters.get<mio::TimeInfectedCritical>()[i] = 8.;

        model.populations[{i, mio::InfectionState::Exposed}]            = fact * num_exp_t0;
        model.populations[{i, mio::InfectionState::InfectedNoSymptoms}] = fact * num_car_t0;
        model.populations[{i, mio::InfectionState::InfectedSymptoms}]   = fact * num_inf_t0;
        model.populations[{i, mio::InfectionState::InfectedSevere}]     = fact * num_hosp_t0;
        model.populations[{i, mio::InfectionState::InfectedCritical}]   = fact * num_icu_t0;
        model.populations[{i, mio::InfectionState::Recovered}]          = fact * num_rec_t0;
        model.populations[{i, mio::InfectionState::Dead}]               = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::InfectionState::Susceptible},
                                                                         fact * num_total_t0);

        model.parameters.get<mio::TransmissionProbabilityOnContact>()[i]  = 0.06;
        model.parameters.get<mio::RelativeTransmissionNoSymptoms>()[i]    = 0.67;
        model.parameters.get<mio::RecoveredPerInfectedNoSymptoms>()[i]    = alpha;
        model.parameters.get<mio::RiskOfInfectionFromSymptomatic>()[i]    = beta;
        model.parameters.get<mio::MaxRiskOfInfectionFromSymptomatic>()[i] = beta * 3;
        model.parameters.get<mio::SeverePerInfectedSymptoms>()[i]         = rho;
        model.parameters.get<mio::CriticalPerSevere>()[i]                 = theta;
        model.parameters.get<mio::DeathsPerCritical>()[i]                 = delta;
    }

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, fact * cont_freq));
    Eigen::MatrixXd m =
        Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, 0.7).triangularView<Eigen::Upper>();
    contact_matrix.add_damping(m, mio::SimulationTime(30.));

    mio::set_params_distributions_normal(model, t0, tmax, 0.15);

    mio::Graph<mio::SecirModel, mio::MigrationParameters> graph;
    graph.add_node(0, model);
    graph.add_node(1, model);
    graph.add_edge(0, 1, Eigen::VectorXd::Constant(model.populations.get_num_compartments(), 0.01));
    graph.add_edge(1, 0, Eigen::VectorXd::Constant(model.populations.get_num_compartments(), 0.01));

    TempFileRegister file_register;
    auto graph_dir    = file_register.get_unique_path("graph_parameters-%%%%-%%%%");
    auto write_status = mio::write_graph(graph, graph_dir);
    ASSERT_THAT(print_wrap(write_status), IsSuccess());

    auto read_result = mio::read_graph<mio::SecirModel>(graph_dir);
    ASSERT_THAT(print_wrap(read_result), IsSuccess());

    auto& graph_read = read_result.value();
    auto num_nodes   = graph.nodes().size();
    auto num_edges   = graph.edges().size();

    ASSERT_EQ(num_nodes, graph_read.nodes().size());
    ASSERT_EQ(num_edges, graph_read.edges().size());

    for (size_t node = 0; node < num_nodes; node++) {
        mio::SecirModel graph_model                = graph.nodes()[0].property;
        mio::ContactMatrixGroup& graph_cont_matrix = graph_model.parameters.get<mio::ContactPatterns>();

        mio::SecirModel graph_read_model                = graph_read.nodes()[0].property;
        mio::ContactMatrixGroup& graph_read_cont_matrix = graph_read_model.parameters.get<mio::ContactPatterns>();

        ASSERT_EQ(graph_read_cont_matrix.get_num_groups(), static_cast<Eigen::Index>((size_t)num_groups));
        ASSERT_EQ(graph_read_cont_matrix, graph_cont_matrix);
        ASSERT_EQ(graph_model.populations.get_num_compartments(), graph_read_model.populations.get_num_compartments());
        ASSERT_EQ(graph.nodes()[node].id, graph_read.nodes()[node].id);
        EXPECT_THAT(graph_read_model.parameters.get<mio::TestAndTraceCapacity>().value(),
                    FloatingPointEqual(graph_model.parameters.get<mio::TestAndTraceCapacity>().value(), 1e-12, 1e-12));
        check_distribution(*graph_model.parameters.get<mio::TestAndTraceCapacity>().get_distribution().get(),
                           *graph_read_model.parameters.get<mio::TestAndTraceCapacity>().get_distribution().get());

        for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(num_groups); group++) {
            ASSERT_EQ((graph_model.populations[{group, mio::InfectionState::Dead}]),
                      (graph_read_model.populations[{group, mio::InfectionState::Dead}]));
            ASSERT_EQ(graph_model.populations.get_total(), graph_read_model.populations.get_total());
            check_distribution(
                *graph_model.populations[{group, mio::InfectionState::Exposed}].get_distribution().get(),
                *graph_read_model.populations[{group, mio::InfectionState::Exposed}].get_distribution().get());
            check_distribution(
                *graph_model.populations[{group, mio::InfectionState::InfectedNoSymptoms}].get_distribution().get(),
                *graph_read_model.populations[{group, mio::InfectionState::InfectedNoSymptoms}]
                     .get_distribution()
                     .get());
            check_distribution(
                *graph_model.populations[{group, mio::InfectionState::InfectedSymptoms}].get_distribution().get(),
                *graph_read_model.populations[{group, mio::InfectionState::InfectedSymptoms}].get_distribution().get());
            check_distribution(
                *graph_model.populations[{group, mio::InfectionState::InfectedSevere}].get_distribution().get(),
                *graph_read_model.populations[{group, mio::InfectionState::InfectedSevere}].get_distribution().get());
            check_distribution(
                *graph_model.populations[{group, mio::InfectionState::InfectedCritical}].get_distribution().get(),
                *graph_read_model.populations[{group, mio::InfectionState::InfectedCritical}].get_distribution().get());
            check_distribution(
                *graph_model.populations[{group, mio::InfectionState::Recovered}].get_distribution().get(),
                *graph_read_model.populations[{group, mio::InfectionState::Recovered}].get_distribution().get());
            check_distribution(
                *graph_model.populations[{group, mio::InfectionState::Exposed}].get_distribution().get(),
                *graph_read_model.populations[{group, mio::InfectionState::Exposed}].get_distribution().get());

            ASSERT_EQ(graph_model.parameters.get<mio::IncubationTime>()[group],
                      graph_read_model.parameters.get<mio::IncubationTime>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::TimeInfectedSymptoms>()[group],
                      graph_read_model.parameters.get<mio::TimeInfectedSymptoms>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::SerialInterval>()[group],
                      graph_read_model.parameters.get<mio::SerialInterval>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::TimeInfectedSevere>()[group],
                      graph_read_model.parameters.get<mio::TimeInfectedSevere>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::TimeInfectedCritical>()[group],
                      graph_read_model.parameters.get<mio::TimeInfectedCritical>()[group]);

            ASSERT_EQ(graph_model.parameters.get<mio::TransmissionProbabilityOnContact>()[group],
                      graph_read_model.parameters.get<mio::TransmissionProbabilityOnContact>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::RiskOfInfectionFromSymptomatic>()[group],
                      graph_read_model.parameters.get<mio::RiskOfInfectionFromSymptomatic>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::MaxRiskOfInfectionFromSymptomatic>()[group],
                      graph_read_model.parameters.get<mio::MaxRiskOfInfectionFromSymptomatic>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::RecoveredPerInfectedNoSymptoms>()[group],
                      graph_read_model.parameters.get<mio::RecoveredPerInfectedNoSymptoms>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::DeathsPerCritical>()[group],
                      graph_read_model.parameters.get<mio::DeathsPerCritical>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::SeverePerInfectedSymptoms>()[group],
                      graph_read_model.parameters.get<mio::SeverePerInfectedSymptoms>()[group]);
            ASSERT_EQ(graph_model.parameters.get<mio::CriticalPerSevere>()[group],
                      graph_read_model.parameters.get<mio::CriticalPerSevere>()[group]);

            check_distribution(*graph_model.parameters.get<mio::IncubationTime>()[group].get_distribution().get(),
                               *graph_read_model.parameters.get<mio::IncubationTime>()[group].get_distribution().get());
            check_distribution(*graph_model.parameters.get<mio::SerialInterval>()[group].get_distribution().get(),
                               *graph_read_model.parameters.get<mio::SerialInterval>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<mio::TimeInfectedSymptoms>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::TimeInfectedSymptoms>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<mio::TimeInfectedSevere>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::TimeInfectedSevere>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<mio::TimeInfectedCritical>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::TimeInfectedCritical>()[group].get_distribution().get());

            check_distribution(
                *graph_model.parameters.get<mio::TimeInfectedSymptoms>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::TimeInfectedSymptoms>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<mio::TimeInfectedSevere>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::TimeInfectedSevere>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<mio::MaxRiskOfInfectionFromSymptomatic>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::MaxRiskOfInfectionFromSymptomatic>()[group]
                     .get_distribution()
                     .get());
            check_distribution(
                *graph_model.parameters.get<mio::DeathsPerCritical>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::DeathsPerCritical>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<mio::CriticalPerSevere>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<mio::CriticalPerSevere>()[group].get_distribution().get());

            ASSERT_EQ(graph_model.parameters.get<mio::ContactPatterns>().get_dampings(),
                      graph_read_model.parameters.get<mio::ContactPatterns>().get_dampings());
        }

        ASSERT_THAT(graph_read.edges(), testing::ElementsAreArray(graph.edges()));
    }
}

TEST(TestSaveParameters, ReadPopulationDataRKIAges)
{
    std::vector<mio::SecirModel> model(1, {6});
    model[0].apply_constraints();
    std::vector<double> scaling_factor_inf(6, 1.0);
    double scaling_factor_icu = 1.0;
    mio::Date date(2020, 12, 10);

    std::string path = TEST_DATA_DIR;

    for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(6); group++) {
        model[0].parameters.get<mio::RecoveredPerInfectedNoSymptoms>()[group] = 0.1 * ((size_t)group + 1);
        model[0].parameters.get<mio::SeverePerInfectedSymptoms>()[group]      = 0.11 * ((size_t)group + 1);
        model[0].parameters.get<mio::CriticalPerSevere>()[group]              = 0.12 * ((size_t)group + 1);
    }
    auto read_result = mio::read_population_data_germany(model, date, scaling_factor_inf, scaling_factor_icu, path);
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
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::Susceptible}]), sus[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::Exposed}]), exp[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::InfectedNoSymptoms}]), car[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::InfectedSymptoms}]), inf[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::InfectedSevere}]), hosp[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::InfectedCritical}]), icu[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::Recovered}]), rec[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::Dead}]), death[i], 1e-1);
    }

    EXPECT_NEAR(model[0].populations.get_total(), 83166695, 1e-6);
}

TEST(TestSaveParameters, ReadPopulationDataStateAllAges)
{
    std::vector<mio::SecirModel> model(1, {6});
    model[0].apply_constraints();
    std::vector<double> scaling_factor_inf(6, 1.0);
    double scaling_factor_icu = 1.0;
    mio::Date date(2020, 12, 10);

    std::vector<int> state = {1};

    std::string path = TEST_DATA_DIR;

    for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(6); group++) {
        model[0].parameters.get<mio::RecoveredPerInfectedNoSymptoms>()[group] = 0.1 * ((size_t)group + 1);
        model[0].parameters.get<mio::SeverePerInfectedSymptoms>()[group]      = 0.11 * ((size_t)group + 1);
        model[0].parameters.get<mio::CriticalPerSevere>()[group]              = 0.12 * ((size_t)group + 1);
    }
    auto read_result =
        mio::read_population_data_state(model, date, state, scaling_factor_inf, scaling_factor_icu, path);
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
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::Susceptible}]), sus[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::Exposed}]), exp[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::InfectedNoSymptoms}]), car[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::InfectedSymptoms}]), inf[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::InfectedSevere}]), hosp[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::InfectedCritical}]), icu[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::Recovered}]), rec[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::Dead}]), death[i], 1e-1);
    }

    EXPECT_NEAR(model[0].populations.get_total(), 2903777, 1e-6);
}

TEST(TestSaveParameters, ReadPopulationDataCountyAllAges)
{

    std::vector<mio::SecirModel> model(1, {6});
    model[0].apply_constraints();
    std::vector<double> scaling_factor_inf(6, 1.0);
    double scaling_factor_icu = 1.0;
    mio::Date date(2020, 12, 10);

    std::vector<int> county = {1002};

    std::string path = TEST_DATA_DIR;

    for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(6); group++) {
        model[0].parameters.get<mio::RecoveredPerInfectedNoSymptoms>()[group] = 0.1 * ((size_t)group + 1);
        model[0].parameters.get<mio::SeverePerInfectedSymptoms>()[group]      = 0.11 * ((size_t)group + 1);
        model[0].parameters.get<mio::CriticalPerSevere>()[group]              = 0.12 * ((size_t)group + 1);
    }
    auto read_result =
        mio::read_population_data_county(model, date, county, scaling_factor_inf, scaling_factor_icu, path);
    ASSERT_THAT(print_wrap(read_result), IsSuccess());

    std::vector<double> sus   = {10284.13, 19082.86, 73783.12, 82494.81, 43725.08, 15612.70};
    std::vector<double> exp   = {0.571429, 4.82143, 20.8163, 22.1429, 4.57143, 4.64286};
    std::vector<double> car   = {0.557143, 4.46429, 22.0408, 20.7143, 4.28571, 4.64286};
    std::vector<double> inf   = {0.42857, 3.285714, 15.2857, 13.0000, 2.42857, 2.00000};
    std::vector<double> hosp  = {0.0942857, 0.691429, 4.90286, 5.34286, 1.41429, 2.45143};
    std::vector<double> icu   = {0.0769231, 0.307692, 0.692308, 1.23077, 1.92308, 2.76923};
    std::vector<double> rec   = {35, 108.571, 640.143, 573.429, 180.429, 75.5714};
    std::vector<double> death = {0, 0, 0, 0, 10, 14.4286};

    for (size_t i = 0; i < 6; i++) {
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::Susceptible}]), sus[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::Exposed}]), exp[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::InfectedNoSymptoms}]), car[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::InfectedSymptoms}]), inf[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::InfectedSevere}]), hosp[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::InfectedCritical}]), icu[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::Recovered}]), rec[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{mio::AgeGroup(i), mio::InfectionState::Dead}]), death[i], 1e-1);
    }

    EXPECT_NEAR(model[0].populations.get_total(), 246793, 1e-6);
}

TEST(TestSaveParameters, ExtrapolateRKI)
{
    std::vector<mio::SecirModel> model{mio::SecirModel(6)};

    model[0].apply_constraints();
    std::vector<double> scaling_factor_inf(6, 1.0);
    double scaling_factor_icu = 1.0;
    mio::Date date(2020, 12, 10);

    std::vector<int> county = {1002};

    for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(6); group++) {
        model[0].parameters.get<mio::RecoveredPerInfectedNoSymptoms>()[group] = 0.1 * ((size_t)group + 1);
        model[0].parameters.get<mio::SeverePerInfectedSymptoms>()[group]      = 0.11 * ((size_t)group + 1);
        model[0].parameters.get<mio::CriticalPerSevere>()[group]              = 0.12 * ((size_t)group + 1);
    }

    TempFileRegister file_register;
    auto results_dir = file_register.get_unique_path("ExtrapolateRKI-%%%%-%%%%");
    boost::filesystem::create_directory(results_dir);
    auto extrapolate_result = mio::export_input_data_county_timeseries(model, TEST_DATA_DIR, results_dir, county, date,
                                                                       scaling_factor_inf, scaling_factor_icu, 1);
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
        EXPECT_NEAR(results[0]((size_t)mio::InfectionState::Susceptible + (size_t)mio::InfectionState::Count * i),
                    sus[i], 1e-1);
        EXPECT_NEAR(results[0]((size_t)mio::InfectionState::Exposed + (size_t)mio::InfectionState::Count * i), exp[i],
                    1e-1);
        EXPECT_NEAR(
            results[0]((size_t)mio::InfectionState::InfectedNoSymptoms + (size_t)mio::InfectionState::Count * i),
            car[i], 1e-1);
        EXPECT_NEAR(results[0]((size_t)mio::InfectionState::InfectedSymptoms + (size_t)mio::InfectionState::Count * i),
                    inf[i], 1e-1);
        EXPECT_NEAR(results[0]((size_t)mio::InfectionState::InfectedSevere + (size_t)mio::InfectionState::Count * i),
                    hosp[i], 1e-1);
        EXPECT_NEAR(results[0]((size_t)mio::InfectionState::InfectedCritical + (size_t)mio::InfectionState::Count * i),
                    icu[i], 1e-1);
        EXPECT_NEAR(results[0]((size_t)mio::InfectionState::Recovered + (size_t)mio::InfectionState::Count * i), rec[i],
                    1e-1);
        EXPECT_NEAR(results[0]((size_t)mio::InfectionState::Dead + (size_t)mio::InfectionState::Count * i), death[i],
                    1e-1);
    }
}
