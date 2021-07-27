#include "load_test_data.h"
#include "test_data_dir.h"
#include "epidemiology/secir/secir.h"
#include <epidemiology/secir/parameter_space.h>
#include <epidemiology/secir/parameter_studies.h>
#include <epidemiology_io/secir_result_io.h>
#include <epidemiology_io/secir_parameters_io.h>
#include <epidemiology/migration/migration.h>
#include <distributions_helpers.h>
#include <matchers.h>
#include "temp_file_register.h"
#include <epidemiology/utils/date.h>
#include <gtest/gtest.h>

TEST(TestSaveParameters, json_single_sim_write_read_compare)
{
    double t0   = 0.0;
    double tmax = 50.5;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           tinfasy = 6.2, ticu2death = 5;

    double cont_freq = 10, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    epi::SecirModel model(2);
    epi::AgeGroup num_groups = model.parameters.get_num_groups();
    double fact       = 1.0 / (double)(size_t)num_groups;

    auto& params = model.parameters;

    for (auto i = epi::AgeGroup(0); i < num_groups; i++) {
        model.parameters.get<epi::IncubationTime>()[i] = tinc;
        model.parameters.get<epi::InfectiousTimeMild>()[i] = tinfmild;
        model.parameters.get<epi::SerialInterval>()[i] = tserint;
        model.parameters.get<epi::HospitalizedToHomeTime>()[i] = thosp2home;
        model.parameters.get<epi::HomeToHospitalizedTime>()[i] = thome2hosp;
        model.parameters.get<epi::HospitalizedToICUTime>()[i] = thosp2icu;
        model.parameters.get<epi::ICUToHomeTime>()[i] = ticu2home;
        model.parameters.get<epi::InfectiousTimeAsymptomatic>()[i] = tinfasy;
        model.parameters.get<epi::ICUToDeathTime>()[i] = ticu2death;

        model.populations[{i, epi::InfectionState::Exposed}] = fact * num_exp_t0;
        model.populations[{i, epi::InfectionState::Carrier}] = fact * num_car_t0;
        model.populations[{i, epi::InfectionState::Infected}] = fact * num_inf_t0;
        model.populations[{i, epi::InfectionState::Hospitalized}] = fact * num_hosp_t0;
        model.populations[{i, epi::InfectionState::ICU}] = fact * num_icu_t0;
        model.populations[{i, epi::InfectionState::Recovered}] = fact * num_rec_t0;
        model.populations[{i, epi::InfectionState::Dead}] = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<epi::AgeGroup>({i, epi::InfectionState::Susceptible}, fact * num_total_t0);

        model.parameters.get<epi::InfectionProbabilityFromContact>()[i] = 0.06;
        model.parameters.get<epi::RelativeCarrierInfectability>()[i] = 0.67;
        model.parameters.get<epi::AsymptoticCasesPerInfectious>()[i] = alpha;
        model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[i] = beta;
        model.parameters.get<epi::HospitalizedCasesPerInfectious>()[i] = rho;
        model.parameters.get<epi::ICUCasesPerHospitalized>()[i] = theta;
        model.parameters.get<epi::DeathsPerHospitalized>()[i] = delta;
    }

    epi::ContactMatrixGroup& contact_matrix = params.get<epi::ContactPatterns>();
    contact_matrix[0] = epi::ContactMatrix(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, fact * cont_freq));
    contact_matrix.add_damping(0.7, epi::SimulationTime(30.));
    auto damping2  = Eigen::MatrixXd::Zero((size_t)num_groups, (size_t)num_groups).eval();
    damping2(0, 0) = 0.8;
    contact_matrix.add_damping(damping2, epi::SimulationTime(35));

    epi::set_params_distributions_normal(model, t0, tmax, 0.2);

    params.get<epi::IncubationTime>()[(epi::AgeGroup)0].get_distribution()->add_predefined_sample(4711.0);
    
    TempFileRegister file_register;
    auto filename = file_register.get_unique_path("TestParameters-%%%%-%%%%.json");
    auto write_status = epi::write_json(filename, model);
    ASSERT_THAT(print_wrap(write_status), IsSuccess());

    auto read_result = epi::read_json(filename, epi::Tag<epi::SecirModel>{});
    ASSERT_THAT(print_wrap(read_result), IsSuccess());
    auto& read_model = read_result.value();

    const epi::UncertainContactMatrix& contact      = model.parameters.get<epi::ContactPatterns>();
    const epi::UncertainContactMatrix& read_contact = read_model.parameters.get<epi::ContactPatterns>();

    num_groups             = model.parameters.get_num_groups();
    auto num_groups_read = read_model.parameters.get_num_groups();
    ASSERT_EQ(num_groups, num_groups_read);

    for (auto i = epi::AgeGroup(0); i < num_groups; i++) {
        ASSERT_EQ((model.populations[{i, epi::InfectionState::Dead}]),
                  (read_model.populations[{i, epi::InfectionState::Dead}]));
        ASSERT_EQ((model.populations.get_group_total(i)), (read_model.populations.get_group_total(i)));
        ASSERT_EQ((model.populations[{i, epi::InfectionState::Exposed}]),
                  (read_model.populations[{i, epi::InfectionState::Exposed}]));
        ASSERT_EQ((model.populations[{i, epi::InfectionState::Carrier}]),
                  (read_model.populations[{i, epi::InfectionState::Carrier}]));
        ASSERT_EQ((model.populations[{i, epi::InfectionState::Infected}]),
                  (read_model.populations[{i, epi::InfectionState::Infected}]));
        ASSERT_EQ((model.populations[{i, epi::InfectionState::Hospitalized}]),
                  (read_model.populations[{i, epi::InfectionState::Hospitalized}]));
        ASSERT_EQ((model.populations[{i, epi::InfectionState::ICU}]),
                  (read_model.populations[{i, epi::InfectionState::ICU}]));
        ASSERT_EQ((model.populations[{i, epi::InfectionState::Recovered}]),
                  (read_model.populations[{i, epi::InfectionState::Recovered}]));

        check_distribution(*model.populations[{i, epi::InfectionState::Exposed}].get_distribution(),
                           *read_model.populations[{i, epi::InfectionState::Exposed}].get_distribution());
        check_distribution(*model.populations[{i, epi::InfectionState::Carrier}].get_distribution(),
                           *read_model.populations[{i, epi::InfectionState::Carrier}].get_distribution());
        check_distribution(*model.populations[{i, epi::InfectionState::Infected}].get_distribution(),
                           *read_model.populations[{i, epi::InfectionState::Infected}].get_distribution());
        check_distribution(*model.populations[{i, epi::InfectionState::Hospitalized}].get_distribution(),
                           *read_model.populations[{i, epi::InfectionState::Hospitalized}].get_distribution());
        check_distribution(*model.populations[{i, epi::InfectionState::ICU}].get_distribution(),
                           *read_model.populations[{i, epi::InfectionState::ICU}].get_distribution());
        check_distribution(*model.populations[{i, epi::InfectionState::Recovered}].get_distribution(),
                           *read_model.populations[{i, epi::InfectionState::Recovered}].get_distribution());

        ASSERT_EQ(model.parameters.get<epi::IncubationTime>()[i], read_model.parameters.get<epi::IncubationTime>()[i]);
        ASSERT_EQ(model.parameters.get<epi::InfectiousTimeMild>()[i],
                  read_model.parameters.get<epi::InfectiousTimeMild>()[i]);
        ASSERT_EQ(model.parameters.get<epi::SerialInterval>()[i], read_model.parameters.get<epi::SerialInterval>()[i]);
        ASSERT_EQ(model.parameters.get<epi::HospitalizedToHomeTime>()[i],
                  read_model.parameters.get<epi::HospitalizedToHomeTime>()[i]);
        ASSERT_EQ(model.parameters.get<epi::HomeToHospitalizedTime>()[i],
                  read_model.parameters.get<epi::HomeToHospitalizedTime>()[i]);
        ASSERT_EQ(model.parameters.get<epi::InfectiousTimeAsymptomatic>()[i],
                  read_model.parameters.get<epi::InfectiousTimeAsymptomatic>()[i]);
        ASSERT_EQ(model.parameters.get<epi::HospitalizedToICUTime>()[i],
                  read_model.parameters.get<epi::HospitalizedToICUTime>()[i]);
        ASSERT_EQ(model.parameters.get<epi::ICUToHomeTime>()[i], read_model.parameters.get<epi::ICUToHomeTime>()[i]);
        ASSERT_EQ(model.parameters.get<epi::ICUToDeathTime>()[i], read_model.parameters.get<epi::ICUToDeathTime>()[i]);

        check_distribution(*model.parameters.get<epi::IncubationTime>()[i].get_distribution(),
                           *read_model.parameters.get<epi::IncubationTime>()[i].get_distribution());
        check_distribution(*model.parameters.get<epi::InfectiousTimeMild>()[i].get_distribution(),
                           *read_model.parameters.get<epi::InfectiousTimeMild>()[i].get_distribution());
        check_distribution(*model.parameters.get<epi::SerialInterval>()[i].get_distribution(),
                           *read_model.parameters.get<epi::SerialInterval>()[i].get_distribution());
        check_distribution(*model.parameters.get<epi::HospitalizedToHomeTime>()[i].get_distribution(),
                           *read_model.parameters.get<epi::HospitalizedToHomeTime>()[i].get_distribution());
        check_distribution(*model.parameters.get<epi::HomeToHospitalizedTime>()[i].get_distribution(),
                           *read_model.parameters.get<epi::HomeToHospitalizedTime>()[i].get_distribution());
        check_distribution(*model.parameters.get<epi::InfectiousTimeAsymptomatic>()[i].get_distribution(),
                           *read_model.parameters.get<epi::InfectiousTimeAsymptomatic>()[i].get_distribution());
        check_distribution(*model.parameters.get<epi::HospitalizedToICUTime>()[i].get_distribution(),
                           *read_model.parameters.get<epi::HospitalizedToICUTime>()[i].get_distribution());
        check_distribution(*model.parameters.get<epi::ICUToHomeTime>()[i].get_distribution(),
                           *read_model.parameters.get<epi::ICUToHomeTime>()[i].get_distribution());
        check_distribution(*model.parameters.get<epi::ICUToDeathTime>()[i].get_distribution(),
                           *read_model.parameters.get<epi::ICUToDeathTime>()[i].get_distribution());

        ASSERT_EQ(model.parameters.get<epi::InfectionProbabilityFromContact>()[i],
                  read_model.parameters.get<epi::InfectionProbabilityFromContact>()[i]);
        ASSERT_EQ(model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[i],
                  read_model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[i]);
        ASSERT_EQ(model.parameters.get<epi::AsymptoticCasesPerInfectious>()[i],
                  read_model.parameters.get<epi::AsymptoticCasesPerInfectious>()[i]);
        ASSERT_EQ(model.parameters.get<epi::DeathsPerHospitalized>()[i],
                  read_model.parameters.get<epi::DeathsPerHospitalized>()[i]);
        ASSERT_EQ(model.parameters.get<epi::HospitalizedCasesPerInfectious>()[i],
                  read_model.parameters.get<epi::HospitalizedCasesPerInfectious>()[i]);
        ASSERT_EQ(model.parameters.get<epi::ICUCasesPerHospitalized>()[i],
                  read_model.parameters.get<epi::ICUCasesPerHospitalized>()[i]);

        check_distribution(*model.parameters.get<epi::InfectionProbabilityFromContact>()[i].get_distribution(),
                           *read_model.parameters.get<epi::InfectionProbabilityFromContact>()[i].get_distribution());
        check_distribution(*model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[i].get_distribution(),
                           *read_model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[i].get_distribution());
        check_distribution(*model.parameters.get<epi::AsymptoticCasesPerInfectious>()[i].get_distribution(),
                           *read_model.parameters.get<epi::AsymptoticCasesPerInfectious>()[i].get_distribution());
        check_distribution(*model.parameters.get<epi::DeathsPerHospitalized>()[i].get_distribution(),
                           *read_model.parameters.get<epi::DeathsPerHospitalized>()[i].get_distribution());
        check_distribution(
            *model.parameters.get<epi::HospitalizedCasesPerInfectious>()[i].get_distribution(),
            *read_model.parameters.get<epi::HospitalizedCasesPerInfectious>()[i].get_distribution());
        check_distribution(*model.parameters.get<epi::ICUCasesPerHospitalized>()[i].get_distribution(),
                           *read_model.parameters.get<epi::ICUCasesPerHospitalized>()[i].get_distribution());

        ASSERT_THAT(contact.get_cont_freq_mat(), testing::ContainerEq(read_contact.get_cont_freq_mat()));
        ASSERT_EQ(contact.get_dampings(), read_contact.get_dampings());
    }
}

TEST(TestSaveParameters, json_graphs_write_read_compare)
{
    double t0   = 0.0;
    double tmax = 50.5;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           tinfasy = 6.2, ticu2death = 5;

    double cont_freq = 10, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    epi::SecirModel model(2);
    epi::AgeGroup num_groups = model.parameters.get_num_groups();
    double fact       = 1.0 / (double)(size_t)num_groups;

    model.parameters.set<epi::TestAndTraceCapacity>(30);

    for (auto i = epi::AgeGroup(0); i < num_groups; i++) {
        model.parameters.get<epi::IncubationTime>()[i] = tinc;
        model.parameters.get<epi::InfectiousTimeMild>()[i] = tinfmild;
        model.parameters.get<epi::SerialInterval>()[i] = tserint;
        model.parameters.get<epi::HospitalizedToHomeTime>()[i] = thosp2home;
        model.parameters.get<epi::HomeToHospitalizedTime>()[i] = thome2hosp;
        model.parameters.get<epi::HospitalizedToICUTime>()[i] = thosp2icu;
        model.parameters.get<epi::ICUToHomeTime>()[i] = ticu2home;
        model.parameters.get<epi::InfectiousTimeAsymptomatic>()[i] = tinfasy;
        model.parameters.get<epi::ICUToDeathTime>()[i] = ticu2death;

        model.populations[{i, epi::InfectionState::Exposed}] = fact * num_exp_t0;
        model.populations[{i, epi::InfectionState::Carrier}] = fact * num_car_t0;
        model.populations[{i, epi::InfectionState::Infected}] = fact * num_inf_t0;
        model.populations[{i, epi::InfectionState::Hospitalized}] = fact * num_hosp_t0;
        model.populations[{i, epi::InfectionState::ICU}] = fact * num_icu_t0;
        model.populations[{i, epi::InfectionState::Recovered}] = fact * num_rec_t0;
        model.populations[{i, epi::InfectionState::Dead}] = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<epi::AgeGroup>({i, epi::InfectionState::Susceptible},
                                                                         fact * num_total_t0);

        model.parameters.get<epi::InfectionProbabilityFromContact>()[i] = 0.06;
        model.parameters.get<epi::RelativeCarrierInfectability>()[i] = 0.67;
        model.parameters.get<epi::AsymptoticCasesPerInfectious>()[i] = alpha;
        model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[i] = beta;
        model.parameters.get<epi::MaxRiskOfInfectionFromSympomatic>()[i] = beta * 3;
        model.parameters.get<epi::HospitalizedCasesPerInfectious>()[i] = rho;
        model.parameters.get<epi::ICUCasesPerHospitalized>()[i] = theta;
        model.parameters.get<epi::DeathsPerHospitalized>()[i] = delta;
    }

    epi::ContactMatrixGroup& contact_matrix = model.parameters.get<epi::ContactPatterns>();
    contact_matrix[0] = epi::ContactMatrix(Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, fact * cont_freq));
    Eigen::MatrixXd m = Eigen::MatrixXd::Constant((size_t)num_groups, (size_t)num_groups, 0.7).triangularView<Eigen::Upper>();
    contact_matrix.add_damping(m, epi::SimulationTime(30.));

    epi::set_params_distributions_normal(model, t0, tmax, 0.15);

    epi::Graph<epi::SecirModel, epi::MigrationParameters> graph;
    graph.add_node(0, model);
    graph.add_node(1, model);
    graph.add_edge(0, 1, Eigen::VectorXd::Constant(model.populations.get_num_compartments(), 0.01));
    graph.add_edge(1, 0, Eigen::VectorXd::Constant(model.populations.get_num_compartments(), 0.01));

    TempFileRegister file_register;
    auto graph_dir = file_register.get_unique_path("graph_parameters-%%%%-%%%%");
    auto write_status = epi::write_graph(graph, graph_dir);
    ASSERT_THAT(print_wrap(write_status), IsSuccess());

    auto read_result = epi::read_graph<epi::SecirModel>(graph_dir);
    ASSERT_THAT(print_wrap(read_result), IsSuccess());

    auto& graph_read = read_result.value();
    auto num_nodes   = graph.nodes().size();
    auto num_edges   = graph.edges().size();

    ASSERT_EQ(num_nodes, graph_read.nodes().size());
    ASSERT_EQ(num_edges, graph_read.edges().size());

    for (size_t node = 0; node < num_nodes; node++) {
        epi::SecirModel graph_model = graph.nodes()[0].property;
        epi::ContactMatrixGroup& graph_cont_matrix  = graph_model.parameters.get<epi::ContactPatterns>();

        epi::SecirModel graph_read_model = graph_read.nodes()[0].property;
        epi::ContactMatrixGroup& graph_read_cont_matrix  = graph_read_model.parameters.get<epi::ContactPatterns>();

        ASSERT_EQ(graph_read_cont_matrix.get_num_groups(), static_cast<Eigen::Index>((size_t)num_groups));
        ASSERT_EQ(graph_read_cont_matrix, graph_cont_matrix);
        ASSERT_EQ(graph_model.populations.get_num_compartments(), graph_read_model.populations.get_num_compartments());
        ASSERT_EQ(graph.nodes()[node].id, graph_read.nodes()[node].id);
        EXPECT_THAT(graph_read_model.parameters.get<epi::TestAndTraceCapacity>().value(),
                    FloatingPointEqual(graph_model.parameters.get<epi::TestAndTraceCapacity>().value(), 1e-12, 1e-12));
        check_distribution(*graph_model.parameters.get<epi::TestAndTraceCapacity>().get_distribution().get(),
                           *graph_read_model.parameters.get<epi::TestAndTraceCapacity>().get_distribution().get());

        for (auto group = epi::AgeGroup(0); group < epi::AgeGroup(num_groups); group++) {
            ASSERT_EQ((graph_model.populations[{group, epi::InfectionState::Dead}]),
                      (graph_read_model.populations[{group, epi::InfectionState::Dead}]));
            ASSERT_EQ(graph_model.populations.get_total(), graph_read_model.populations.get_total());
            check_distribution(
                *graph_model.populations[{group, epi::InfectionState::Exposed}].get_distribution().get(),
                *graph_read_model.populations[{group, epi::InfectionState::Exposed}].get_distribution().get());
            check_distribution(
                *graph_model.populations[{group, epi::InfectionState::Carrier}].get_distribution().get(),
                *graph_read_model.populations[{group, epi::InfectionState::Carrier}].get_distribution().get());
            check_distribution(
                *graph_model.populations[{group, epi::InfectionState::Infected}].get_distribution().get(),
                *graph_read_model.populations[{group, epi::InfectionState::Infected}].get_distribution().get());
            check_distribution(
                *graph_model.populations[{group, epi::InfectionState::Hospitalized}].get_distribution().get(),
                *graph_read_model.populations[{group, epi::InfectionState::Hospitalized}].get_distribution().get());
            check_distribution(
                *graph_model.populations[{group, epi::InfectionState::ICU}].get_distribution().get(),
                *graph_read_model.populations[{group, epi::InfectionState::ICU}].get_distribution().get());
            check_distribution(
                *graph_model.populations[{group, epi::InfectionState::Recovered}].get_distribution().get(),
                *graph_read_model.populations[{group, epi::InfectionState::Recovered}].get_distribution().get());
            check_distribution(
                *graph_model.populations[{group, epi::InfectionState::Exposed}].get_distribution().get(),
                *graph_read_model.populations[{group, epi::InfectionState::Exposed}].get_distribution().get());

            ASSERT_EQ(graph_model.parameters.get<epi::IncubationTime>()[group],
                      graph_read_model.parameters.get<epi::IncubationTime>()[group]);
            ASSERT_EQ(graph_model.parameters.get<epi::InfectiousTimeMild>()[group],
                      graph_read_model.parameters.get<epi::InfectiousTimeMild>()[group]);
            ASSERT_EQ(graph_model.parameters.get<epi::SerialInterval>()[group],
                      graph_read_model.parameters.get<epi::SerialInterval>()[group]);
            ASSERT_EQ(graph_model.parameters.get<epi::HospitalizedToHomeTime>()[group],
                      graph_read_model.parameters.get<epi::HospitalizedToHomeTime>()[group]);
            ASSERT_EQ(graph_model.parameters.get<epi::HomeToHospitalizedTime>()[group],
                      graph_read_model.parameters.get<epi::HomeToHospitalizedTime>()[group]);
            ASSERT_EQ(graph_model.parameters.get<epi::InfectiousTimeAsymptomatic>()[group],
                      graph_read_model.parameters.get<epi::InfectiousTimeAsymptomatic>()[group]);
            ASSERT_EQ(graph_model.parameters.get<epi::HospitalizedToICUTime>()[group],
                      graph_read_model.parameters.get<epi::HospitalizedToICUTime>()[group]);
            ASSERT_EQ(graph_model.parameters.get<epi::ICUToHomeTime>()[group],
                      graph_read_model.parameters.get<epi::ICUToHomeTime>()[group]);
            ASSERT_EQ(graph_model.parameters.get<epi::ICUToDeathTime>()[group],
                      graph_read_model.parameters.get<epi::ICUToDeathTime>()[group]);

            ASSERT_EQ(graph_model.parameters.get<epi::InfectionProbabilityFromContact>()[group],
                      graph_read_model.parameters.get<epi::InfectionProbabilityFromContact>()[group]);
            ASSERT_EQ(graph_model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[group],
                      graph_read_model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[group]);
            ASSERT_EQ(graph_model.parameters.get<epi::MaxRiskOfInfectionFromSympomatic>()[group],
                      graph_read_model.parameters.get<epi::MaxRiskOfInfectionFromSympomatic>()[group]);
            ASSERT_EQ(graph_model.parameters.get<epi::AsymptoticCasesPerInfectious>()[group],
                      graph_read_model.parameters.get<epi::AsymptoticCasesPerInfectious>()[group]);
            ASSERT_EQ(graph_model.parameters.get<epi::DeathsPerHospitalized>()[group],
                      graph_read_model.parameters.get<epi::DeathsPerHospitalized>()[group]);
            ASSERT_EQ(graph_model.parameters.get<epi::HospitalizedCasesPerInfectious>()[group],
                      graph_read_model.parameters.get<epi::HospitalizedCasesPerInfectious>()[group]);
            ASSERT_EQ(graph_model.parameters.get<epi::ICUCasesPerHospitalized>()[group],
                      graph_read_model.parameters.get<epi::ICUCasesPerHospitalized>()[group]);

            check_distribution(*graph_model.parameters.get<epi::IncubationTime>()[group].get_distribution().get(),
                               *graph_read_model.parameters.get<epi::IncubationTime>()[group].get_distribution().get());
            check_distribution(*graph_model.parameters.get<epi::SerialInterval>()[group].get_distribution().get(),
                               *graph_read_model.parameters.get<epi::SerialInterval>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<epi::InfectiousTimeMild>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<epi::InfectiousTimeMild>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<epi::HospitalizedToHomeTime>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<epi::HospitalizedToHomeTime>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<epi::HomeToHospitalizedTime>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<epi::HomeToHospitalizedTime>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<epi::InfectiousTimeAsymptomatic>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<epi::InfectiousTimeAsymptomatic>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<epi::HospitalizedToICUTime>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<epi::HospitalizedToICUTime>()[group].get_distribution().get());
            check_distribution(*graph_model.parameters.get<epi::ICUToHomeTime>()[group].get_distribution().get(),
                               *graph_read_model.parameters.get<epi::ICUToHomeTime>()[group].get_distribution().get());
            check_distribution(*graph_model.parameters.get<epi::ICUToDeathTime>()[group].get_distribution().get(),
                               *graph_read_model.parameters.get<epi::ICUToDeathTime>()[group].get_distribution().get());

            check_distribution(
                *graph_model.parameters.get<epi::InfectiousTimeMild>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<epi::InfectiousTimeMild>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<epi::HospitalizedToHomeTime>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<epi::HospitalizedToHomeTime>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<epi::HomeToHospitalizedTime>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<epi::HomeToHospitalizedTime>()[group].get_distribution().get());
            check_distribution(*graph_model.parameters.get<epi::MaxRiskOfInfectionFromSympomatic>()[group]
                                    .get_distribution()
                                    .get(),
                               *graph_read_model.parameters.get<epi::MaxRiskOfInfectionFromSympomatic>()[group]
                                    .get_distribution()
                                    .get());
            check_distribution(
                *graph_model.parameters.get<epi::DeathsPerHospitalized>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<epi::DeathsPerHospitalized>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<epi::InfectiousTimeAsymptomatic>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<epi::InfectiousTimeAsymptomatic>()[group].get_distribution().get());
            check_distribution(
                *graph_model.parameters.get<epi::ICUCasesPerHospitalized>()[group].get_distribution().get(),
                *graph_read_model.parameters.get<epi::ICUCasesPerHospitalized>()[group].get_distribution().get());

            ASSERT_EQ(graph_model.parameters.get<epi::ContactPatterns>().get_dampings(),
                      graph_read_model.parameters.get<epi::ContactPatterns>().get_dampings());
        }

        ASSERT_THAT(graph_read.edges(), testing::ElementsAreArray(graph.edges()));
    }
}

TEST(TestSaveParameters, ReadPopulationDataRKIAges)
{
    std::vector<epi::SecirModel> model(1, {6});
    model[0].apply_constraints();
    std::vector<double> scaling_factor_inf(6, 1.0);
    double scaling_factor_icu = 1.0;
    epi::Date date(2020, 12, 10);

    std::string path = TEST_DATA_DIR;

    for (auto group = epi::AgeGroup(0); group < epi::AgeGroup(6); group++) {
        model[0].parameters.get<epi::AsymptoticCasesPerInfectious>()[group] = 0.1 * ((size_t)group + 1);
        model[0].parameters.get<epi::HospitalizedCasesPerInfectious>()[group] = 0.11 * ((size_t)group + 1);
        model[0].parameters.get<epi::ICUCasesPerHospitalized>()[group] = 0.12 * ((size_t)group + 1);
    }
    auto read_result = epi::read_population_data_germany(model, date, scaling_factor_inf, scaling_factor_icu, path);
    ASSERT_THAT(print_wrap(read_result), IsSuccess());

    std::vector<double> sus   = {3444023.09, 7666389.350, 18801939.83, 29522450.59, 16317865.95, 6059469.35};
    std::vector<double> exp   = {389.843, 1417.37, 6171.74, 8765.6, 3554.5, 2573.89};
    std::vector<double> car   = {389.443, 1412.86, 6077.14, 8554.77, 3437.57, 2462.09};
    std::vector<double> inf   = {297.924, 811.551, 2270.16, 1442.03, 0, 0};
    std::vector<double> hosp  = {39.9614, 303.191, 1934.84, 3621.2, 1793.39, 1557.03};
    std::vector<double> icu   = {47.6813, 190.725, 429.132, 762.901, 1192.03, 1716.53};
    std::vector<double> rec   = {23557.7, 78946.3, 398585.142, 487273.71, 178660.14, 96021.9};
    std::vector<double> death = {2, 4, 48, 1137.86, 8174.14, 18528.9};

    for (size_t i = 0; i < 6; i++) {
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Susceptible}]), sus[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Exposed}]), exp[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Carrier}]), car[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Infected}]), inf[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Hospitalized}]), hosp[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::ICU}]), icu[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Recovered}]), rec[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Dead}]), death[i], 1e-1);
    }

    EXPECT_NEAR(model[0].populations.get_total(), 83166695, 1e-6);
}

TEST(TestSaveParameters, ReadPopulationDataStateAllAges)
{
    std::vector<epi::SecirModel> model(1, {6});
    model[0].apply_constraints();
    std::vector<double> scaling_factor_inf(6, 1.0);
    double scaling_factor_icu = 1.0;
    epi::Date date(2020, 12, 10);

    std::vector<int> state = {1};

    std::string path = TEST_DATA_DIR;

    for (auto group = epi::AgeGroup(0); group < epi::AgeGroup(6); group++) {
        model[0].parameters.get<epi::AsymptoticCasesPerInfectious>()[group] = 0.1 * ((size_t)group + 1);
        model[0].parameters.get<epi::HospitalizedCasesPerInfectious>()[group] = 0.11 * ((size_t)group + 1);
        model[0].parameters.get<epi::ICUCasesPerHospitalized>()[group] = 0.12 * ((size_t)group + 1);
    }
    auto read_result = epi::read_population_data_state(model, date, state, scaling_factor_inf, scaling_factor_icu, path);
    ASSERT_THAT(print_wrap(read_result), IsSuccess());

    std::vector<double> sus   = {116695.3, 283933, 622945.61, 1042462.09, 606578.8, 212990};
    std::vector<double> exp   = {7.64286, 23.7143, 103.243, 134.486, 43, 38};
    std::vector<double> car   = {7, 20.4286, 99.4143, 126.971, 41.6429, 36.4286};
    std::vector<double> inf   = {5.59286, 11.0429, 37.7571, 22.6629, 0.0785714, 0};
    std::vector<double> hosp  = {0.707143, 3.92857, 30.6429, 50.5371, 20.35, 19.9886};
    std::vector<double> icu   = {0.274725, 1.0989, 2.47253, 4.3956, 6.86813, 9.89011};
    std::vector<double> rec   = {393.143, 1216.14, 5467.86, 6543.57, 2281.29, 1045.71};
    std::vector<double> death = {0, 0, 0, 16.2857, 99.5714, 198.286};

    for (size_t i = 0; i < 6; i++) {
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Susceptible}]), sus[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Exposed}]), exp[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Carrier}]), car[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Infected}]), inf[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Hospitalized}]), hosp[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::ICU}]), icu[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Recovered}]), rec[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Dead}]), death[i], 1e-1);
    }

    EXPECT_NEAR(model[0].populations.get_total(), 2903777, 1e-6);
}

TEST(TestSaveParameters, ReadPopulationDataCountyAllAges)
{

    std::vector<epi::SecirModel> model(1, {6});
    model[0].apply_constraints();
    std::vector<double> scaling_factor_inf(6, 1.0);
    double scaling_factor_icu = 1.0;
    epi::Date date(2020, 12, 10);

    std::vector<int> county = {1002};

    std::string path = TEST_DATA_DIR;

    for (auto group = epi::AgeGroup(0); group < epi::AgeGroup(6); group++) {
        model[0].parameters.get<epi::AsymptoticCasesPerInfectious>()[group] = 0.1 * ((size_t)group + 1);
        model[0].parameters.get<epi::HospitalizedCasesPerInfectious>()[group] = 0.11 * ((size_t)group + 1);
        model[0].parameters.get<epi::ICUCasesPerHospitalized>()[group] = 0.12 * ((size_t)group + 1);
    }
    auto read_result =
        epi::read_population_data_county(model, date, county, scaling_factor_inf, scaling_factor_icu, path);
    ASSERT_THAT(print_wrap(read_result), IsSuccess());

    std::vector<double> sus   = {10284.4, 19086.2, 73805.3, 82522.6, 43731.9, 15620.2};
    std::vector<double> exp   = {0.571429, 3.8, 14.8286, 12.9429, 2.21429, 1.85714};
    std::vector<double> car   = {0.557143, 3.51429, 15.3857, 12.6571, 2.28571, 1.94286};
    std::vector<double> inf   = {0.291429, 1.93714, 5.79714, 2.45714, 0, 0};
    std::vector<double> hosp  = {0.0942857, 0.691429, 4.90286, 5.34286, 1.41429, 2.45143};
    std::vector<double> icu   = {0.0769231, 0.307692, 0.692308, 1.23077, 1.92308, 2.76923};
    std::vector<double> rec   = {35, 108.571, 640.143, 573.429, 180.429, 75.5714};
    std::vector<double> death = {0, 0, 0, 0, 10, 14.4286};

    for (size_t i = 0; i < 6; i++) {
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Susceptible}]), sus[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Exposed}]), exp[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Carrier}]), car[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Infected}]), inf[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Hospitalized}]), hosp[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::ICU}]), icu[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Recovered}]), rec[i], 1e-1);
        EXPECT_NEAR((model[0].populations[{epi::AgeGroup(i), epi::InfectionState::Dead}]), death[i], 1e-1);
    }

    EXPECT_NEAR(model[0].populations.get_total(), 246793, 1e-6);
}

TEST(TestSaveParameters, GetCountyIDs)
{
    std::vector<int> true_ids = {
        1001,  1002,  1003,  1004,  1051,  1053,  1054,  1055,  1056,  1057,  1058,  1059,  1060,  1061,  1062,  2000,
        3101,  3102,  3103,  3151,  3153,  3154,  3155,  3157,  3158,  3159,  3241,  3251,  3252,  3254,  3255,  3256,
        3257,  3351,  3352,  3353,  3354,  3355,  3356,  3357,  3358,  3359,  3360,  3361,  3401,  3402,  3403,  3404,
        3405,  3451,  3452,  3453,  3454,  3455,  3456,  3457,  3458,  3459,  3460,  3461,  3462,  4011,  4012,  5111,
        5112,  5113,  5114,  5116,  5117,  5119,  5120,  5122,  5124,  5154,  5158,  5162,  5166,  5170,  5314,  5315,
        5316,  5334,  5358,  5362,  5366,  5370,  5374,  5378,  5382,  5512,  5513,  5515,  5554,  5558,  5562,  5566,
        5570,  5711,  5754,  5758,  5762,  5766,  5770,  5774,  5911,  5913,  5914,  5915,  5916,  5954,  5958,  5962,
        5966,  5970,  5974,  5978,  6411,  6412,  6413,  6414,  6431,  6432,  6433,  6434,  6435,  6436,  6437,  6438,
        6439,  6440,  6531,  6532,  6533,  6534,  6535,  6611,  6631,  6632,  6633,  6634,  6635,  6636,  7111,  7131,
        7132,  7133,  7134,  7135,  7137,  7138,  7140,  7141,  7143,  7211,  7231,  7232,  7233,  7235,  7311,  7312,
        7313,  7314,  7315,  7316,  7317,  7318,  7319,  7320,  7331,  7332,  7333,  7334,  7335,  7336,  7337,  7338,
        7339,  7340,  8111,  8115,  8116,  8117,  8118,  8119,  8121,  8125,  8126,  8127,  8128,  8135,  8136,  8211,
        8212,  8215,  8216,  8221,  8222,  8225,  8226,  8231,  8235,  8236,  8237,  8311,  8315,  8316,  8317,  8325,
        8326,  8327,  8335,  8336,  8337,  8415,  8416,  8417,  8421,  8425,  8426,  8435,  8436,  8437,  9161,  9162,
        9163,  9171,  9172,  9173,  9174,  9175,  9176,  9177,  9178,  9179,  9180,  9181,  9182,  9183,  9184,  9185,
        9186,  9187,  9188,  9189,  9190,  9261,  9262,  9263,  9271,  9272,  9273,  9274,  9275,  9276,  9277,  9278,
        9279,  9361,  9362,  9363,  9371,  9372,  9373,  9374,  9375,  9376,  9377,  9461,  9462,  9463,  9464,  9471,
        9472,  9473,  9474,  9475,  9476,  9477,  9478,  9479,  9561,  9562,  9563,  9564,  9565,  9571,  9572,  9573,
        9574,  9575,  9576,  9577,  9661,  9662,  9663,  9671,  9672,  9673,  9674,  9675,  9676,  9677,  9678,  9679,
        9761,  9762,  9763,  9764,  9771,  9772,  9773,  9774,  9775,  9776,  9777,  9778,  9779,  9780,  10041, 10042,
        10043, 10044, 10045, 10046, 11000, 12051, 12052, 12053, 12054, 12060, 12061, 12062, 12063, 12064, 12065, 12066,
        12067, 12068, 12069, 12070, 12071, 12072, 12073, 13003, 13004, 13071, 13072, 13073, 13074, 13075, 13076, 14511,
        14521, 14522, 14523, 14524, 14612, 14625, 14626, 14627, 14628, 14713, 14729, 14730, 15001, 15002, 15003, 15081,
        15082, 15083, 15084, 15085, 15086, 15087, 15088, 15089, 15090, 15091, 16051, 16052, 16053, 16054, 16055, 16056,
        16061, 16062, 16063, 16064, 16065, 16066, 16067, 16068, 16069, 16070, 16071, 16072, 16073, 16074, 16075, 16076,
        16077};

    std::string path = TEST_DATA_DIR;
    auto read_ids    = epi::get_county_ids(path);
    ASSERT_THAT(print_wrap(read_ids), IsSuccess());

    EXPECT_THAT(read_ids.value(), testing::ElementsAreArray(true_ids));
}

TEST(TestSaveParameters, ExtrapolateRKI)
{
    std::vector<epi::SecirModel> model{epi::SecirModel(6)};

    model[0].apply_constraints();
    std::vector<double> scaling_factor_inf(6, 1.0);
    double scaling_factor_icu = 1.0;
    epi::Date date(2020, 12, 10);

    std::vector<int> county = {1002};

    for (auto group = epi::AgeGroup(0); group < epi::AgeGroup(6); group++) {
        model[0].parameters.get<epi::AsymptoticCasesPerInfectious>()[group] = 0.1 * ((size_t)group + 1);
        model[0].parameters.get<epi::HospitalizedCasesPerInfectious>()[group] = 0.11 * ((size_t)group + 1);
        model[0].parameters.get<epi::ICUCasesPerHospitalized>()[group] = 0.12 * ((size_t)group + 1);
    }

    TempFileRegister file_register;
    auto results_dir = file_register.get_unique_path("ExtrapolateRKI-%%%%-%%%%");
    boost::filesystem::create_directory(results_dir);
    auto extrapolate_result =
        epi::extrapolate_rki_results(model, TEST_DATA_DIR, results_dir, county, date, scaling_factor_inf, scaling_factor_icu, 1);
    ASSERT_THAT(print_wrap(extrapolate_result), IsSuccess());

    auto read_result = epi::read_result(epi::path_join(results_dir, "Results_rki.h5"), 6);
    ASSERT_THAT(print_wrap(read_result), IsSuccess());
    auto& file_results = read_result.value();
    auto results = file_results[0].get_groups();

    std::vector<double> sus   = {10284.4, 19086.2, 73805.3, 82522.6, 43731.9, 15620.2};
    std::vector<double> exp   = {0.571429, 3.8, 14.8286, 12.9429, 2.21429, 1.85714};
    std::vector<double> car   = {0.557143, 3.51429, 15.3857, 12.6571, 2.28571, 1.94286};
    std::vector<double> inf   = {0.291429, 1.93714, 5.79714, 2.45714, 0, 0};
    std::vector<double> hosp  = {0.0942857, 0.691429, 4.90286, 5.34286, 1.41429, 2.45143};
    std::vector<double> icu   = {0.0769231, 0.307692, 0.692308, 1.23077, 1.92308, 2.76923};
    std::vector<double> rec   = {35, 108.571, 640.143, 573.429, 180.429, 75.5714};
    std::vector<double> death = {0, 0, 0, 0, 10, 14.4286};

    for (size_t i = 0; i < 6; i++) {
        EXPECT_NEAR(results[0]((size_t)epi::InfectionState::Susceptible + (size_t)epi::InfectionState::Count * i),
                    sus[i], 1e-1);
        EXPECT_NEAR(results[0]((size_t)epi::InfectionState::Exposed + (size_t)epi::InfectionState::Count * i), exp[i],
                    1e-1);
        EXPECT_NEAR(results[0]((size_t)epi::InfectionState::Carrier + (size_t)epi::InfectionState::Count * i), car[i],
                    1e-1);
        EXPECT_NEAR(results[0]((size_t)epi::InfectionState::Infected + (size_t)epi::InfectionState::Count * i), inf[i],
                    1e-1);
        EXPECT_NEAR(results[0]((size_t)epi::InfectionState::Hospitalized + (size_t)epi::InfectionState::Count * i),
                    hosp[i], 1e-1);
        EXPECT_NEAR(results[0]((size_t)epi::InfectionState::ICU + (size_t)epi::InfectionState::Count * i), icu[i],
                    1e-1);
        EXPECT_NEAR(results[0]((size_t)epi::InfectionState::Recovered + (size_t)epi::InfectionState::Count * i), rec[i],
                    1e-1);
        EXPECT_NEAR(results[0]((size_t)epi::InfectionState::Dead + (size_t)epi::InfectionState::Count * i), death[i],
                    1e-1);
    }
}
