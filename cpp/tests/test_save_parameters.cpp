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
#include <epidemiology/utils/date.h>
#include <gtest/gtest.h>

TEST(TestSaveParameters, compareParameterStudy)
{
    double t0   = 0.0;
    double tmax = 50.5;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           tinfasy = 6.2, ticu2death = 5;

    double cont_freq = 10, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    epi::SecirModel<epi::AgeGroup2> model;
    size_t num_groups = model.parameters.get_num_groups();
    double fact       = 1.0 / (double)num_groups;

    auto& params = model.parameters;

    for (size_t i = 0; i < num_groups; i++) {
        model.parameters.times[i].set_incubation(tinc);
        model.parameters.times[i].set_infectious_mild(tinfmild);
        model.parameters.times[i].set_serialinterval(tserint);
        model.parameters.times[i].set_hospitalized_to_home(thosp2home);
        model.parameters.times[i].set_home_to_hospitalized(thome2hosp);
        model.parameters.times[i].set_hospitalized_to_icu(thosp2icu);
        model.parameters.times[i].set_icu_to_home(ticu2home);
        model.parameters.times[i].set_infectious_asymp(tinfasy);
        model.parameters.times[i].set_icu_to_death(ticu2death);

        model.populations.set(fact * num_exp_t0, (epi::AgeGroup2)i, epi::InfectionType::E);
        model.populations.set(fact * num_car_t0, (epi::AgeGroup2)i, epi::InfectionType::C);
        model.populations.set(fact * num_inf_t0, (epi::AgeGroup2)i, epi::InfectionType::I);
        model.populations.set(fact * num_hosp_t0, (epi::AgeGroup2)i, epi::InfectionType::H);
        model.populations.set(fact * num_icu_t0, (epi::AgeGroup2)i, epi::InfectionType::U);
        model.populations.set(fact * num_rec_t0, (epi::AgeGroup2)i, epi::InfectionType::R);
        model.populations.set(fact * num_dead_t0, (epi::AgeGroup2)i, epi::InfectionType::D);
        model.populations.set_difference_from_group_total(fact * num_total_t0, (epi::AgeGroup2)i, (epi::AgeGroup2)i,
                                                          epi::InfectionType::S);

        model.parameters.probabilities[i].set_infection_from_contact(0.06);
        model.parameters.probabilities[i].set_carrier_infectability(0.67);
        model.parameters.probabilities[i].set_asymp_per_infectious(alpha);
        model.parameters.probabilities[i].set_risk_from_symptomatic(beta);
        model.parameters.probabilities[i].set_hospitalized_per_infectious(rho);
        model.parameters.probabilities[i].set_icu_per_hospitalized(theta);
        model.parameters.probabilities[i].set_dead_per_icu(delta);
    }

    epi::ContactMatrixGroup& contact_matrix = params.get_contact_patterns();
    contact_matrix[0] = epi::ContactMatrix(Eigen::MatrixXd::Constant(num_groups, num_groups, fact * cont_freq));
    contact_matrix.add_damping(0.7, epi::SimulationTime(30.));
    auto damping2  = Eigen::MatrixXd::Zero(num_groups, num_groups).eval();
    damping2(0, 0) = 0.8;
    contact_matrix.add_damping(damping2, epi::SimulationTime(35));

    int num_runs     = 5;
    std::string path = "/Parameters";
    TixiDocumentHandle handle;

    epi::set_params_distributions_normal(model, t0, tmax, 0.2);

    params.times[0].get_incubation().get_distribution()->add_predefined_sample(4711.0);

    params.get_contact_patterns().get_distribution_damp_days()->add_predefined_sample(4711.0);
    tixiCreateDocument("Parameters", &handle);
    epi::ParameterStudy<epi::SecirModel<epi::AgeGroup2>> study(model, t0, tmax, num_runs);

    epi::write_parameter_study(handle, path, study);
    tixiSaveDocument(handle, "TestParameters.xml");
    tixiCloseDocument(handle);

    tixiOpenDocument("TestParameters.xml", &handle);
    epi::ParameterStudy<epi::SecirModel<epi::AgeGroup2>> read_study =
        epi::read_parameter_study<epi::AgeGroup2>(handle, path);
    tixiCloseDocument(handle);

    ASSERT_EQ(study.get_num_runs(), read_study.get_num_runs());
    ASSERT_EQ(study.get_t0(), read_study.get_t0());
    ASSERT_EQ(study.get_tmax(), read_study.get_tmax());

    const auto& read_model = read_study.get_model();

    const epi::UncertainContactMatrix& contact      = study.get_model().parameters.get_contact_patterns();
    const epi::UncertainContactMatrix& read_contact = read_model.parameters.get_contact_patterns();

    num_groups             = study.get_model().parameters.get_num_groups();
    size_t num_groups_read = read_model.parameters.get_num_groups();
    ASSERT_EQ(num_groups, num_groups_read);

    for (size_t i = 0; i < num_groups; i++) {
        ASSERT_EQ(model.populations.get((epi::AgeGroup2)i, epi::InfectionType::D),
                  read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::D));
        ASSERT_EQ(model.populations.get_group_total((epi::AgeGroup2)i),
                  read_model.populations.get_group_total((epi::AgeGroup2)i));
        ASSERT_EQ(model.populations.get((epi::AgeGroup2)i, epi::InfectionType::E),
                  read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::E));
        ASSERT_EQ(model.populations.get((epi::AgeGroup2)i, epi::InfectionType::C),
                  read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::C));
        ASSERT_EQ(model.populations.get((epi::AgeGroup2)i, epi::InfectionType::I),
                  read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::I));
        ASSERT_EQ(model.populations.get((epi::AgeGroup2)i, epi::InfectionType::H),
                  read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::H));
        ASSERT_EQ(model.populations.get((epi::AgeGroup2)i, epi::InfectionType::U),
                  read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::U));
        ASSERT_EQ(model.populations.get((epi::AgeGroup2)i, epi::InfectionType::R),
                  read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::R));

        check_distribution(*model.populations.get((epi::AgeGroup2)i, epi::InfectionType::E).get_distribution(),
                           *read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::E).get_distribution());
        check_distribution(*model.populations.get((epi::AgeGroup2)i, epi::InfectionType::C).get_distribution(),
                           *read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::C).get_distribution());
        check_distribution(*model.populations.get((epi::AgeGroup2)i, epi::InfectionType::I).get_distribution(),
                           *read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::I).get_distribution());
        check_distribution(*model.populations.get((epi::AgeGroup2)i, epi::InfectionType::H).get_distribution(),
                           *read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::H).get_distribution());
        check_distribution(*model.populations.get((epi::AgeGroup2)i, epi::InfectionType::U).get_distribution(),
                           *read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::U).get_distribution());
        check_distribution(*model.populations.get((epi::AgeGroup2)i, epi::InfectionType::R).get_distribution(),
                           *read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::R).get_distribution());

        ASSERT_EQ(model.parameters.times[i].get_incubation(), read_model.parameters.times[i].get_incubation());
        ASSERT_EQ(model.parameters.times[i].get_infectious_mild(),
                  read_model.parameters.times[i].get_infectious_mild());
        ASSERT_EQ(model.parameters.times[i].get_serialinterval(), read_model.parameters.times[i].get_serialinterval());
        ASSERT_EQ(model.parameters.times[i].get_hospitalized_to_home(),
                  read_model.parameters.times[i].get_hospitalized_to_home());
        ASSERT_EQ(model.parameters.times[i].get_home_to_hospitalized(),
                  read_model.parameters.times[i].get_home_to_hospitalized());
        ASSERT_EQ(model.parameters.times[i].get_infectious_asymp(),
                  read_model.parameters.times[i].get_infectious_asymp());
        ASSERT_EQ(model.parameters.times[i].get_hospitalized_to_icu(),
                  read_model.parameters.times[i].get_hospitalized_to_icu());
        ASSERT_EQ(model.parameters.times[i].get_icu_to_home(), read_model.parameters.times[i].get_icu_to_home());
        ASSERT_EQ(model.parameters.times[i].get_icu_to_dead(), read_model.parameters.times[i].get_icu_to_dead());

        check_distribution(*model.parameters.times[i].get_incubation().get_distribution(),
                           *read_model.parameters.times[i].get_incubation().get_distribution());
        check_distribution(*model.parameters.times[i].get_infectious_mild().get_distribution(),
                           *read_model.parameters.times[i].get_infectious_mild().get_distribution());
        check_distribution(*model.parameters.times[i].get_serialinterval().get_distribution(),
                           *read_model.parameters.times[i].get_serialinterval().get_distribution());
        check_distribution(*model.parameters.times[i].get_hospitalized_to_home().get_distribution(),
                           *read_model.parameters.times[i].get_hospitalized_to_home().get_distribution());
        check_distribution(*model.parameters.times[i].get_home_to_hospitalized().get_distribution(),
                           *read_model.parameters.times[i].get_home_to_hospitalized().get_distribution());
        check_distribution(*model.parameters.times[i].get_infectious_asymp().get_distribution(),
                           *read_model.parameters.times[i].get_infectious_asymp().get_distribution());
        check_distribution(*model.parameters.times[i].get_hospitalized_to_icu().get_distribution(),
                           *read_model.parameters.times[i].get_hospitalized_to_icu().get_distribution());
        check_distribution(*model.parameters.times[i].get_icu_to_home().get_distribution(),
                           *read_model.parameters.times[i].get_icu_to_home().get_distribution());
        check_distribution(*model.parameters.times[i].get_icu_to_dead().get_distribution(),
                           *read_model.parameters.times[i].get_icu_to_dead().get_distribution());

        ASSERT_EQ(model.parameters.probabilities[i].get_infection_from_contact(),
                  read_model.parameters.probabilities[i].get_infection_from_contact());
        ASSERT_EQ(model.parameters.probabilities[i].get_risk_from_symptomatic(),
                  read_model.parameters.probabilities[i].get_risk_from_symptomatic());
        ASSERT_EQ(model.parameters.probabilities[i].get_asymp_per_infectious(),
                  read_model.parameters.probabilities[i].get_asymp_per_infectious());
        ASSERT_EQ(model.parameters.probabilities[i].get_dead_per_icu(),
                  read_model.parameters.probabilities[i].get_dead_per_icu());
        ASSERT_EQ(model.parameters.probabilities[i].get_hospitalized_per_infectious(),
                  read_model.parameters.probabilities[i].get_hospitalized_per_infectious());
        ASSERT_EQ(model.parameters.probabilities[i].get_icu_per_hospitalized(),
                  read_model.parameters.probabilities[i].get_icu_per_hospitalized());

        check_distribution(*model.parameters.probabilities[i].get_infection_from_contact().get_distribution(),
                           *read_model.parameters.probabilities[i].get_infection_from_contact().get_distribution());
        check_distribution(*model.parameters.probabilities[i].get_risk_from_symptomatic().get_distribution(),
                           *read_model.parameters.probabilities[i].get_risk_from_symptomatic().get_distribution());
        check_distribution(*model.parameters.probabilities[i].get_asymp_per_infectious().get_distribution(),
                           *read_model.parameters.probabilities[i].get_asymp_per_infectious().get_distribution());
        check_distribution(*model.parameters.probabilities[i].get_dead_per_icu().get_distribution(),
                           *read_model.parameters.probabilities[i].get_dead_per_icu().get_distribution());
        check_distribution(
            *model.parameters.probabilities[i].get_hospitalized_per_infectious().get_distribution(),
            *read_model.parameters.probabilities[i].get_hospitalized_per_infectious().get_distribution());
        check_distribution(*model.parameters.probabilities[i].get_icu_per_hospitalized().get_distribution(),
                           *read_model.parameters.probabilities[i].get_icu_per_hospitalized().get_distribution());

        ASSERT_THAT(contact.get_cont_freq_mat(), testing::ContainerEq(read_contact.get_cont_freq_mat()));
    }

    check_distribution(*contact.get_distribution_damp_nb().get(), *read_contact.get_distribution_damp_nb().get());
    check_distribution(*contact.get_distribution_damp_days().get(), *read_contact.get_distribution_damp_days().get());
    check_distribution(*contact.get_distribution_damp_diag_base().get(),
                       *read_contact.get_distribution_damp_diag_base().get());
    check_distribution(*contact.get_distribution_damp_diag_rel().get(),
                       *read_contact.get_distribution_damp_diag_rel().get());
    check_distribution(*contact.get_distribution_damp_offdiag_rel().get(),
                       *read_contact.get_distribution_damp_offdiag_rel().get());
}

TEST(TestSaveParameters, compareSingleRun)
{
    double t0   = 0.0;
    double tmax = 50.5;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           tinfasy = 6.2, ticu2death = 5;

    double cont_freq = 10, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    epi::SecirModel<epi::AgeGroup2> model;
    auto num_groups = size_t(epi::AgeGroup2::Count);
    double fact     = 1.0 / (double)num_groups;

    for (size_t i = 0; i < num_groups; i++) {
        model.parameters.times[i].set_incubation(tinc);
        model.parameters.times[i].set_infectious_mild(tinfmild);
        model.parameters.times[i].set_serialinterval(tserint);
        model.parameters.times[i].set_hospitalized_to_home(thosp2home);
        model.parameters.times[i].set_home_to_hospitalized(thome2hosp);
        model.parameters.times[i].set_hospitalized_to_icu(thosp2icu);
        model.parameters.times[i].set_icu_to_home(ticu2home);
        model.parameters.times[i].set_infectious_asymp(tinfasy);
        model.parameters.times[i].set_icu_to_death(ticu2death);

        model.populations.set(fact * num_exp_t0, (epi::AgeGroup2)i, epi::InfectionType::E);
        model.populations.set(fact * num_car_t0, (epi::AgeGroup2)i, epi::InfectionType::C);
        model.populations.set(fact * num_inf_t0, (epi::AgeGroup2)i, epi::InfectionType::I);
        model.populations.set(fact * num_hosp_t0, (epi::AgeGroup2)i, epi::InfectionType::H);
        model.populations.set(fact * num_icu_t0, (epi::AgeGroup2)i, epi::InfectionType::U);
        model.populations.set(fact * num_rec_t0, (epi::AgeGroup2)i, epi::InfectionType::R);
        model.populations.set(fact * num_dead_t0, (epi::AgeGroup2)i, epi::InfectionType::D);
        model.populations.set_difference_from_group_total(fact * num_total_t0, (epi::AgeGroup2)i, (epi::AgeGroup2)i,
                                                          epi::InfectionType::S);

        model.parameters.probabilities[i].set_infection_from_contact(0.06);
        model.parameters.probabilities[i].set_carrier_infectability(0.67);
        model.parameters.probabilities[i].set_asymp_per_infectious(alpha);
        model.parameters.probabilities[i].set_risk_from_symptomatic(beta);
        model.parameters.probabilities[i].set_hospitalized_per_infectious(rho);
        model.parameters.probabilities[i].set_icu_per_hospitalized(theta);
        model.parameters.probabilities[i].set_dead_per_icu(delta);
    }

    epi::ContactMatrixGroup& contact_matrix = model.parameters.get_contact_patterns();
    contact_matrix[0] = epi::ContactMatrix(Eigen::MatrixXd::Constant(num_groups, num_groups, fact * cont_freq));
    contact_matrix.add_damping(0.7, epi::SimulationTime(30.));
    auto damping2  = Eigen::MatrixXd::Zero(num_groups, num_groups).eval();
    damping2(0, 0) = 0.8;
    contact_matrix.add_damping(damping2, epi::SimulationTime(35));

    int num_runs     = 5;
    std::string path = "/Parameters";
    TixiDocumentHandle handle;

    tixiCreateDocument("Parameters", &handle);
    epi::set_params_distributions_normal(model, t0, tmax, 0.0);
    epi::ParameterStudy<epi::SecirModel<epi::AgeGroup2>> study(model, t0, tmax, num_runs);

    epi::write_parameter_study(handle, path, study, 0);
    tixiSaveDocument(handle, "TestParameterValues.xml");
    tixiCloseDocument(handle);

    tixiOpenDocument("TestParameterValues.xml", &handle);
    epi::ParameterStudy<epi::SecirModel<epi::AgeGroup2>> read_study =
        epi::read_parameter_study<epi::AgeGroup2>(handle, path);
    tixiCloseDocument(handle);

    ASSERT_EQ(study.get_num_runs(), read_study.get_num_runs());
    ASSERT_EQ(study.get_t0(), read_study.get_t0());
    ASSERT_EQ(study.get_tmax(), read_study.get_tmax());

    const epi::SecirModel<epi::AgeGroup2>& read_model = read_study.get_model();

    const epi::UncertainContactMatrix& contact      = study.get_model().parameters.get_contact_patterns();
    const epi::UncertainContactMatrix& read_contact = read_model.parameters.get_contact_patterns();

    num_groups             = study.get_model().parameters.get_num_groups();
    size_t num_groups_read = read_model.parameters.get_num_groups();

    ASSERT_EQ(num_groups, num_groups_read);

    for (size_t i = 0; i < num_groups; i++) {
        ASSERT_EQ(model.populations.get((epi::AgeGroup2)i, epi::InfectionType::D),
                  read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::D));
        ASSERT_EQ(model.populations.get_group_total((epi::AgeGroup2)i),
                  read_model.populations.get_group_total((epi::AgeGroup2)i));
        ASSERT_EQ(model.populations.get((epi::AgeGroup2)i, epi::InfectionType::E),
                  read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::E));
        ASSERT_EQ(model.populations.get((epi::AgeGroup2)i, epi::InfectionType::C),
                  read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::C));
        ASSERT_EQ(model.populations.get((epi::AgeGroup2)i, epi::InfectionType::I),
                  read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::I));
        ASSERT_EQ(model.populations.get((epi::AgeGroup2)i, epi::InfectionType::H),
                  read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::H));
        ASSERT_EQ(model.populations.get((epi::AgeGroup2)i, epi::InfectionType::U),
                  read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::U));
        ASSERT_EQ(model.populations.get((epi::AgeGroup2)i, epi::InfectionType::R),
                  read_model.populations.get((epi::AgeGroup2)i, epi::InfectionType::R));

        ASSERT_EQ(model.parameters.times[i].get_incubation(), read_model.parameters.times[i].get_incubation());
        ASSERT_EQ(model.parameters.times[i].get_infectious_mild(),
                  read_model.parameters.times[i].get_infectious_mild());
        ASSERT_EQ(model.parameters.times[i].get_serialinterval(), read_model.parameters.times[i].get_serialinterval());
        ASSERT_EQ(model.parameters.times[i].get_hospitalized_to_home(),
                  read_model.parameters.times[i].get_hospitalized_to_home());
        ASSERT_EQ(model.parameters.times[i].get_home_to_hospitalized(),
                  read_model.parameters.times[i].get_home_to_hospitalized());
        ASSERT_EQ(model.parameters.times[i].get_infectious_asymp(),
                  read_model.parameters.times[i].get_infectious_asymp());
        ASSERT_EQ(model.parameters.times[i].get_hospitalized_to_icu(),
                  read_model.parameters.times[i].get_hospitalized_to_icu());
        ASSERT_EQ(model.parameters.times[i].get_icu_to_home(), read_model.parameters.times[i].get_icu_to_home());
        ASSERT_EQ(model.parameters.times[i].get_icu_to_dead(), read_model.parameters.times[i].get_icu_to_dead());

        ASSERT_EQ(model.parameters.probabilities[i].get_infection_from_contact(),
                  read_model.parameters.probabilities[i].get_infection_from_contact());
        ASSERT_EQ(model.parameters.probabilities[i].get_risk_from_symptomatic(),
                  read_model.parameters.probabilities[i].get_risk_from_symptomatic());
        ASSERT_EQ(model.parameters.probabilities[i].get_asymp_per_infectious(),
                  read_model.parameters.probabilities[i].get_asymp_per_infectious());
        ASSERT_EQ(model.parameters.probabilities[i].get_dead_per_icu(),
                  read_model.parameters.probabilities[i].get_dead_per_icu());
        ASSERT_EQ(model.parameters.probabilities[i].get_hospitalized_per_infectious(),
                  read_model.parameters.probabilities[i].get_hospitalized_per_infectious());
        ASSERT_EQ(model.parameters.probabilities[i].get_icu_per_hospitalized(),
                  read_model.parameters.probabilities[i].get_icu_per_hospitalized());

        ASSERT_EQ(contact.get_cont_freq_mat(), read_contact.get_cont_freq_mat());
    }
}

TEST(TestSaveParameters, compareGraphs)
{
    double t0   = 0.0;
    double tmax = 50.5;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           tinfasy = 6.2, ticu2death = 5;

    double cont_freq = 10, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    epi::SecirModel<epi::AgeGroup2> model;
    size_t num_groups = size_t(epi::AgeGroup2::Count);
    double fact       = 1.0 / (double)num_groups;

    model.parameters.set_test_and_trace_capacity(30);

    for (size_t i = 0; i < num_groups; i++) {
        model.parameters.times[i].set_incubation(tinc);
        model.parameters.times[i].set_infectious_mild(tinfmild);
        model.parameters.times[i].set_serialinterval(tserint);
        model.parameters.times[i].set_hospitalized_to_home(thosp2home);
        model.parameters.times[i].set_home_to_hospitalized(thome2hosp);
        model.parameters.times[i].set_hospitalized_to_icu(thosp2icu);
        model.parameters.times[i].set_icu_to_home(ticu2home);
        model.parameters.times[i].set_infectious_asymp(tinfasy);
        model.parameters.times[i].set_icu_to_death(ticu2death);

        model.populations.set(fact * num_exp_t0, (epi::AgeGroup2)i, epi::InfectionType::E);
        model.populations.set(fact * num_car_t0, (epi::AgeGroup2)i, epi::InfectionType::C);
        model.populations.set(fact * num_inf_t0, (epi::AgeGroup2)i, epi::InfectionType::I);
        model.populations.set(fact * num_hosp_t0, (epi::AgeGroup2)i, epi::InfectionType::H);
        model.populations.set(fact * num_icu_t0, (epi::AgeGroup2)i, epi::InfectionType::U);
        model.populations.set(fact * num_rec_t0, (epi::AgeGroup2)i, epi::InfectionType::R);
        model.populations.set(fact * num_dead_t0, (epi::AgeGroup2)i, epi::InfectionType::D);
        model.populations.set_difference_from_group_total(fact * num_total_t0, (epi::AgeGroup2)i, (epi::AgeGroup2)i,
                                                          epi::InfectionType::S);

        model.parameters.probabilities[i].set_infection_from_contact(0.06);
        model.parameters.probabilities[i].set_carrier_infectability(0.67);
        model.parameters.probabilities[i].set_asymp_per_infectious(alpha);
        model.parameters.probabilities[i].set_risk_from_symptomatic(beta);
        model.parameters.probabilities[i].set_test_and_trace_max_risk_from_symptomatic(beta * 3);
        model.parameters.probabilities[i].set_hospitalized_per_infectious(rho);
        model.parameters.probabilities[i].set_icu_per_hospitalized(theta);
        model.parameters.probabilities[i].set_dead_per_icu(delta);
    }

    epi::ContactMatrixGroup& contact_matrix = model.parameters.get_contact_patterns();
    contact_matrix[0] = epi::ContactMatrix(Eigen::MatrixXd::Constant(num_groups, num_groups, fact * cont_freq));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant(num_groups, num_groups, 0.7).triangularView<Eigen::Upper>(),
                               epi::SimulationTime(30.));

    epi::set_params_distributions_normal(model, t0, tmax, 0.15);

    epi::Graph<epi::SecirModel<epi::AgeGroup2>, epi::MigrationEdge> graph;
    graph.add_node(0, model);
    graph.add_node(1, model);
    graph.add_edge(0, 1, Eigen::VectorXd::Constant(model.populations.get_num_compartments(), 0.01));
    graph.add_edge(1, 0, Eigen::VectorXd::Constant(model.populations.get_num_compartments(), 0.01));

    epi::write_graph(graph, "graph_parameters");

    epi::Graph<epi::SecirModel<epi::AgeGroup2>, epi::MigrationEdge> graph_read =
        epi::read_graph<epi::SecirModel<epi::AgeGroup2>>("graph_parameters");

    auto num_nodes = graph.nodes().size();
    auto num_edges = graph.edges().size();

    ASSERT_EQ(num_nodes, graph_read.nodes().size());
    ASSERT_EQ(num_edges, graph_read.edges().size());

    for (size_t node = 0; node < num_nodes; node++) {
        epi::SecirModel<epi::AgeGroup2> graph_model = graph.nodes()[0].property;
        epi::ContactMatrixGroup& graph_cont_matrix  = graph_model.parameters.get_contact_patterns();

        epi::SecirModel<epi::AgeGroup2> graph_read_model = graph_read.nodes()[0].property;
        epi::ContactMatrixGroup& graph_read_cont_matrix  = graph_read_model.parameters.get_contact_patterns();

        ASSERT_EQ(graph_read_cont_matrix.get_num_groups(), static_cast<Eigen::Index>(num_groups));
        ASSERT_EQ(graph_read_cont_matrix, graph_cont_matrix);
        ASSERT_EQ(graph_model.populations.get_num_compartments(), graph_read_model.populations.get_num_compartments());
        ASSERT_EQ(graph.nodes()[node].id, graph_read.nodes()[node].id);
        EXPECT_THAT(graph_read_model.parameters.get_test_and_trace_capacity().value(),
                    FloatingPointEqual(graph_model.parameters.get_test_and_trace_capacity().value(), 1e-12, 1e-12));
        check_distribution(*graph_model.parameters.get_test_and_trace_capacity().get_distribution().get(),
                           *graph_read_model.parameters.get_test_and_trace_capacity().get_distribution().get());

        for (size_t group = 0; group < num_groups; group++) {
            ASSERT_EQ(graph_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::D),
                      graph_read_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::D));
            ASSERT_EQ(graph_model.populations.get_total(), graph_read_model.populations.get_total());
            check_distribution(
                *graph_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::E).get_distribution().get(),
                *graph_read_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::E)
                     .get_distribution()
                     .get());
            check_distribution(
                *graph_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::C).get_distribution().get(),
                *graph_read_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::C)
                     .get_distribution()
                     .get());
            check_distribution(
                *graph_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::I).get_distribution().get(),
                *graph_read_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::I)
                     .get_distribution()
                     .get());
            check_distribution(
                *graph_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::H).get_distribution().get(),
                *graph_read_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::H)
                     .get_distribution()
                     .get());
            check_distribution(
                *graph_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::U).get_distribution().get(),
                *graph_read_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::U)
                     .get_distribution()
                     .get());
            check_distribution(
                *graph_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::R).get_distribution().get(),
                *graph_read_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::R)
                     .get_distribution()
                     .get());
            check_distribution(
                *graph_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::E).get_distribution().get(),
                *graph_read_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::E)
                     .get_distribution()
                     .get());

            ASSERT_EQ(graph_model.parameters.times[group].get_incubation(),
                      graph_read_model.parameters.times[group].get_incubation());
            ASSERT_EQ(graph_model.parameters.times[group].get_infectious_mild(),
                      graph_read_model.parameters.times[group].get_infectious_mild());
            ASSERT_EQ(graph_model.parameters.times[group].get_serialinterval(),
                      graph_read_model.parameters.times[group].get_serialinterval());
            ASSERT_EQ(graph_model.parameters.times[group].get_hospitalized_to_home(),
                      graph_read_model.parameters.times[group].get_hospitalized_to_home());
            ASSERT_EQ(graph_model.parameters.times[group].get_home_to_hospitalized(),
                      graph_read_model.parameters.times[group].get_home_to_hospitalized());
            ASSERT_EQ(graph_model.parameters.times[group].get_infectious_asymp(),
                      graph_read_model.parameters.times[group].get_infectious_asymp());
            ASSERT_EQ(graph_model.parameters.times[group].get_hospitalized_to_icu(),
                      graph_read_model.parameters.times[group].get_hospitalized_to_icu());
            ASSERT_EQ(graph_model.parameters.times[group].get_icu_to_home(),
                      graph_read_model.parameters.times[group].get_icu_to_home());
            ASSERT_EQ(graph_model.parameters.times[group].get_icu_to_dead(),
                      graph_read_model.parameters.times[group].get_icu_to_dead());

            ASSERT_EQ(graph_model.parameters.probabilities[group].get_infection_from_contact(),
                      graph_read_model.parameters.probabilities[group].get_infection_from_contact());
            ASSERT_EQ(graph_model.parameters.probabilities[group].get_risk_from_symptomatic(),
                      graph_read_model.parameters.probabilities[group].get_risk_from_symptomatic());
            ASSERT_EQ(graph_model.parameters.probabilities[group].get_test_and_trace_max_risk_from_symptomatic(),
                      graph_read_model.parameters.probabilities[group].get_test_and_trace_max_risk_from_symptomatic());
            ASSERT_EQ(graph_model.parameters.probabilities[group].get_asymp_per_infectious(),
                      graph_read_model.parameters.probabilities[group].get_asymp_per_infectious());
            ASSERT_EQ(graph_model.parameters.probabilities[group].get_dead_per_icu(),
                      graph_read_model.parameters.probabilities[group].get_dead_per_icu());
            ASSERT_EQ(graph_model.parameters.probabilities[group].get_hospitalized_per_infectious(),
                      graph_read_model.parameters.probabilities[group].get_hospitalized_per_infectious());
            ASSERT_EQ(graph_model.parameters.probabilities[group].get_icu_per_hospitalized(),
                      graph_read_model.parameters.probabilities[group].get_icu_per_hospitalized());

            check_distribution(*graph_model.parameters.times[group].get_incubation().get_distribution().get(),
                               *graph_read_model.parameters.times[group].get_incubation().get_distribution().get());
            check_distribution(*graph_model.parameters.times[group].get_serialinterval().get_distribution().get(),
                               *graph_read_model.parameters.times[group].get_serialinterval().get_distribution().get());
            check_distribution(
                *graph_model.parameters.times[group].get_infectious_mild().get_distribution().get(),
                *graph_read_model.parameters.times[group].get_infectious_mild().get_distribution().get());
            check_distribution(
                *graph_model.parameters.times[group].get_hospitalized_to_home().get_distribution().get(),
                *graph_read_model.parameters.times[group].get_hospitalized_to_home().get_distribution().get());
            check_distribution(
                *graph_model.parameters.times[group].get_home_to_hospitalized().get_distribution().get(),
                *graph_read_model.parameters.times[group].get_home_to_hospitalized().get_distribution().get());
            check_distribution(
                *graph_model.parameters.times[group].get_infectious_asymp().get_distribution().get(),
                *graph_read_model.parameters.times[group].get_infectious_asymp().get_distribution().get());
            check_distribution(
                *graph_model.parameters.times[group].get_hospitalized_to_icu().get_distribution().get(),
                *graph_read_model.parameters.times[group].get_hospitalized_to_icu().get_distribution().get());
            check_distribution(*graph_model.parameters.times[group].get_icu_to_home().get_distribution().get(),
                               *graph_read_model.parameters.times[group].get_icu_to_home().get_distribution().get());
            check_distribution(*graph_model.parameters.times[group].get_icu_to_dead().get_distribution().get(),
                               *graph_read_model.parameters.times[group].get_icu_to_dead().get_distribution().get());

            check_distribution(
                *graph_model.parameters.times[group].get_infectious_mild().get_distribution().get(),
                *graph_read_model.parameters.times[group].get_infectious_mild().get_distribution().get());
            check_distribution(
                *graph_model.parameters.times[group].get_hospitalized_to_home().get_distribution().get(),
                *graph_read_model.parameters.times[group].get_hospitalized_to_home().get_distribution().get());
            check_distribution(
                *graph_model.parameters.times[group].get_home_to_hospitalized().get_distribution().get(),
                *graph_read_model.parameters.times[group].get_home_to_hospitalized().get_distribution().get());
            check_distribution(*graph_model.parameters.probabilities[group]
                                    .get_test_and_trace_max_risk_from_symptomatic()
                                    .get_distribution()
                                    .get(),
                               *graph_read_model.parameters.probabilities[group]
                                    .get_test_and_trace_max_risk_from_symptomatic()
                                    .get_distribution()
                                    .get());
            check_distribution(
                *graph_model.parameters.probabilities[group].get_dead_per_icu().get_distribution().get(),
                *graph_read_model.parameters.probabilities[group].get_dead_per_icu().get_distribution().get());
            check_distribution(
                *graph_model.parameters.times[group].get_infectious_asymp().get_distribution().get(),
                *graph_read_model.parameters.times[group].get_infectious_asymp().get_distribution().get());
            check_distribution(
                *graph_model.parameters.probabilities[group].get_icu_per_hospitalized().get_distribution().get(),
                *graph_read_model.parameters.probabilities[group].get_icu_per_hospitalized().get_distribution().get());

            check_distribution(*graph_model.parameters.get_contact_patterns().get_distribution_damp_nb().get(),
                               *graph_read_model.parameters.get_contact_patterns().get_distribution_damp_nb().get());
            check_distribution(*graph_model.parameters.get_contact_patterns().get_distribution_damp_days().get(),
                               *graph_read_model.parameters.get_contact_patterns().get_distribution_damp_days().get());
            check_distribution(
                *graph_model.parameters.get_contact_patterns().get_distribution_damp_diag_base().get(),
                *graph_read_model.parameters.get_contact_patterns().get_distribution_damp_diag_base().get());
            check_distribution(
                *graph_model.parameters.get_contact_patterns().get_distribution_damp_diag_rel().get(),
                *graph_read_model.parameters.get_contact_patterns().get_distribution_damp_diag_rel().get());
            check_distribution(
                *graph_model.parameters.get_contact_patterns().get_distribution_damp_offdiag_rel().get(),
                *graph_read_model.parameters.get_contact_patterns().get_distribution_damp_offdiag_rel().get());
        }

        ASSERT_THAT(graph_read.edges(), testing::ElementsAreArray(graph.edges()));
    }
}

TEST(TestSaveParameters, compareGraphWithFile)
{
    double t0   = 0.0;
    double tmax = 50.5;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           tinfasy = 6.2, ticu2death = 5;

    double cont_freq = 10, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    size_t num_groups = 2;
    double fact       = 1.0 / (double)num_groups;

    using Model = epi::SecirModel<epi::AgeGroup2>;
    Model model;
    auto& params = model.parameters;

    for (size_t i = 0; i < num_groups; i++) {
        params.times[i].set_incubation(tinc);
        params.times[i].set_infectious_mild(tinfmild);
        params.times[i].set_serialinterval(tserint);
        params.times[i].set_hospitalized_to_home(thosp2home);
        params.times[i].set_home_to_hospitalized(thome2hosp);
        params.times[i].set_hospitalized_to_icu(thosp2icu);
        params.times[i].set_icu_to_home(ticu2home);
        params.times[i].set_infectious_asymp(tinfasy);
        params.times[i].set_icu_to_death(ticu2death);

        model.populations.set(fact * num_exp_t0, (epi::AgeGroup2)i, epi::InfectionType::E);
        model.populations.set(fact * num_car_t0, (epi::AgeGroup2)i, epi::InfectionType::C);
        model.populations.set(fact * num_inf_t0, (epi::AgeGroup2)i, epi::InfectionType::I);
        model.populations.set(fact * num_hosp_t0, (epi::AgeGroup2)i, epi::InfectionType::H);
        model.populations.set(fact * num_icu_t0, (epi::AgeGroup2)i, epi::InfectionType::U);
        model.populations.set(fact * num_rec_t0, (epi::AgeGroup2)i, epi::InfectionType::R);
        model.populations.set(fact * num_dead_t0, (epi::AgeGroup2)i, epi::InfectionType::D);
        model.populations.set_difference_from_group_total(fact * num_total_t0, (epi::AgeGroup2)i, (epi::AgeGroup2)i,
                                                          epi::InfectionType::S);

        params.probabilities[i].set_infection_from_contact(0.06);
        params.probabilities[i].set_carrier_infectability(0.67);
        params.probabilities[i].set_asymp_per_infectious(alpha);
        params.probabilities[i].set_risk_from_symptomatic(beta);
        params.probabilities[i].set_hospitalized_per_infectious(rho);
        params.probabilities[i].set_icu_per_hospitalized(theta);
        params.probabilities[i].set_dead_per_icu(delta);
    }

    epi::ContactMatrixGroup& contact_matrix = params.get_contact_patterns();
    contact_matrix[0] = epi::ContactMatrix(Eigen::MatrixXd::Constant(num_groups, num_groups, fact * cont_freq));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant(num_groups, num_groups, 0.7).triangularView<Eigen::Upper>(),
                               epi::SimulationTime(30.));

    epi::set_params_distributions_normal(model, t0, tmax, 0.15);

    epi::Graph<Model, epi::MigrationEdge> graph;
    graph.add_node(0, model);
    graph.add_node(1, model);
    graph.add_edge(0, 1, Eigen::VectorXd::Constant(model.populations.get_num_compartments(), 0.01));
    graph.add_edge(1, 0, Eigen::VectorXd::Constant(model.populations.get_num_compartments(), 0.01));

    epi::Graph<Model, epi::MigrationEdge> graph_read = epi::read_graph<Model>(TEST_DATA_DIR);

    auto num_nodes = graph.nodes().size();
    auto num_edges = graph.edges().size();

    ASSERT_EQ(num_nodes, graph_read.nodes().size());
    ASSERT_EQ(num_edges, graph_read.edges().size());

    for (size_t node = 0; node < num_nodes; node++) {
        Model graph_model                       = graph.nodes()[0].property;
        epi::ContactMatrixGroup graph_cont_freq = graph_model.parameters.get_contact_patterns();

        Model graph_read_model                       = graph_read.nodes()[0].property;
        epi::ContactMatrixGroup graph_read_cont_freq = graph_read_model.parameters.get_contact_patterns();

        ASSERT_EQ(num_groups, static_cast<size_t>(graph_read_cont_freq.get_num_groups()));
        ASSERT_EQ(graph_model.populations.get_num_compartments(), graph_read_model.populations.get_num_compartments());

        for (size_t group = 0; group < num_groups; group++) {
            ASSERT_EQ(graph_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::D),
                      graph_read_model.populations.get((epi::AgeGroup2)group, epi::InfectionType::D));
            ASSERT_EQ(graph_model.populations.get_total(), graph_read_model.populations.get_total());
            check_distribution(
                *graph_model.parameters.probabilities[group].get_infection_from_contact().get_distribution().get(),
                *graph_read_model.parameters.probabilities[group]
                     .get_infection_from_contact()
                     .get_distribution()
                     .get());
            check_distribution(
                *graph_model.parameters.probabilities[group].get_asymp_per_infectious().get_distribution().get(),
                *graph_read_model.parameters.probabilities[group].get_asymp_per_infectious().get_distribution().get());
            check_distribution(
                *graph_model.parameters.probabilities[group].get_risk_from_symptomatic().get_distribution().get(),
                *graph_read_model.parameters.probabilities[group].get_risk_from_symptomatic().get_distribution().get());
            check_distribution(
                *graph_model.parameters.probabilities[group].get_dead_per_icu().get_distribution().get(),
                *graph_read_model.parameters.probabilities[group].get_dead_per_icu().get_distribution().get());
            check_distribution(
                *graph_model.parameters.probabilities[group].get_hospitalized_per_infectious().get_distribution().get(),
                *graph_read_model.parameters.probabilities[group]
                     .get_hospitalized_per_infectious()
                     .get_distribution()
                     .get());
            check_distribution(
                *graph_model.parameters.probabilities[group].get_icu_per_hospitalized().get_distribution().get(),
                *graph_read_model.parameters.probabilities[group].get_icu_per_hospitalized().get_distribution().get());

            ASSERT_EQ(graph_cont_freq, graph_read_cont_freq);

            check_distribution(*graph_model.parameters.get_contact_patterns().get_distribution_damp_nb().get(),
                               *graph_read_model.parameters.get_contact_patterns().get_distribution_damp_nb().get());
            check_distribution(*graph_model.parameters.get_contact_patterns().get_distribution_damp_days().get(),
                               *graph_read_model.parameters.get_contact_patterns().get_distribution_damp_days().get());
            check_distribution(
                *graph_model.parameters.get_contact_patterns().get_distribution_damp_diag_base().get(),
                *graph_read_model.parameters.get_contact_patterns().get_distribution_damp_diag_base().get());
            check_distribution(
                *graph_model.parameters.get_contact_patterns().get_distribution_damp_diag_rel().get(),
                *graph_read_model.parameters.get_contact_patterns().get_distribution_damp_diag_rel().get());
            check_distribution(
                *graph_model.parameters.get_contact_patterns().get_distribution_damp_offdiag_rel().get(),
                *graph_read_model.parameters.get_contact_patterns().get_distribution_damp_offdiag_rel().get());
        }

        ASSERT_THAT(graph_read.edges(), testing::ElementsAreArray(graph.edges()));
    }
}

TEST(TestSaveParameters, ReadPopulationDataRKIAges)
{
    std::vector<epi::SecirModel<epi::AgeGroup6>> model(1);
    model[0].apply_constraints();
    std::vector<double> scaling_factor_inf(6, 1.0);
    double scaling_factor_icu = 1.0;
    epi::Date date(2020, 12, 10);

    std::string path = TEST_DATA_DIR;

    for (size_t group = 0; group < 6; group++) {
        model[0].parameters.probabilities[group].set_asymp_per_infectious(0.1 * (group + 1));
        model[0].parameters.probabilities[group].set_hospitalized_per_infectious(0.11 * (group + 1));
        model[0].parameters.probabilities[group].set_icu_per_hospitalized(0.12 * (group + 1));
    }
    epi::read_population_data_germany(model, date, scaling_factor_inf, scaling_factor_icu, path);

    std::vector<double> sus   = {3444023.09, 7666389.350, 18801939.83, 29522450.59, 16317865.95, 6059469.35};
    std::vector<double> exp   = {389.843, 1417.37, 6171.74, 8765.6, 3554.5, 2573.89};
    std::vector<double> car   = {389.443, 1412.86, 6077.14, 8554.77, 3437.57, 2462.09};
    std::vector<double> inf   = {297.924, 811.551, 2270.16, 1442.03, 0, 0};
    std::vector<double> hosp  = {39.9614, 303.191, 1934.84, 3621.2, 1793.39, 1557.03};
    std::vector<double> icu   = {47.6813, 190.725, 429.132, 762.901, 1192.03, 1716.53};
    std::vector<double> rec   = {23557.7, 78946.3, 398585.142, 487273.71, 178660.14, 96021.9};
    std::vector<double> death = {2, 4, 48, 1137.86, 8174.14, 18528.9};

    for (size_t i = 0; i < 6; i++) {
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::S), sus[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::E), exp[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::C), car[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::I), inf[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::H), hosp[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::U), icu[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::R), rec[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::D), death[i], 1e-1);
    }

    EXPECT_NEAR(model[0].populations.get_total(), 83166695, 1e-6);
}

TEST(TestSaveParameters, ReadPopulationDataStateAllAges)
{
    std::vector<epi::SecirModel<epi::AgeGroup6>> model(1);
    model[0].apply_constraints();
    std::vector<double> scaling_factor_inf(6, 1.0);
    double scaling_factor_icu = 1.0;
    epi::Date date(2020, 12, 10);

    std::vector<int> state = {1};

    std::string path = TEST_DATA_DIR;

    for (size_t group = 0; group < 6; group++) {
        model[0].parameters.probabilities[group].set_asymp_per_infectious(0.1 * (group + 1));
        model[0].parameters.probabilities[group].set_hospitalized_per_infectious(0.11 * (group + 1));
        model[0].parameters.probabilities[group].set_icu_per_hospitalized(0.12 * (group + 1));
    }
    epi::read_population_data_state(model, date, state, scaling_factor_inf, scaling_factor_icu, path);

    std::vector<double> sus   = {116695.3, 283933, 622945.61, 1042462.09, 606578.8, 212990};
    std::vector<double> exp   = {7.64286, 23.7143, 103.243, 134.486, 43, 38};
    std::vector<double> car   = {7, 20.4286, 99.4143, 126.971, 41.6429, 36.4286};
    std::vector<double> inf   = {5.59286, 11.0429, 37.7571, 22.6629, 0.0785714, 0};
    std::vector<double> hosp  = {0.707143, 3.92857, 30.6429, 50.5371, 20.35, 19.9886};
    std::vector<double> icu   = {0.274725, 1.0989, 2.47253, 4.3956, 6.86813, 9.89011};
    std::vector<double> rec   = {393.143, 1216.14, 5467.86, 6543.57, 2281.29, 1045.71};
    std::vector<double> death = {0, 0, 0, 16.2857, 99.5714, 198.286};

    for (size_t i = 0; i < 6; i++) {
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::S), sus[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::E), exp[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::C), car[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::I), inf[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::H), hosp[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::U), icu[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::R), rec[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::D), death[i], 1e-1);
    }

    EXPECT_NEAR(model[0].populations.get_total(), 2903777, 1e-6);
}

TEST(TestSaveParameters, ReadPopulationDataCountyAllAges)
{

    std::vector<epi::SecirModel<epi::AgeGroup6>> model(1);
    model[0].apply_constraints();
    std::vector<double> scaling_factor_inf(6, 1.0);
    double scaling_factor_icu = 1.0;
    epi::Date date(2020, 12, 10);

    std::vector<int> county = {1002};

    std::string path = TEST_DATA_DIR;

    for (size_t group = 0; group < 6; group++) {
        model[0].parameters.probabilities[group].set_asymp_per_infectious(0.1 * (group + 1));
        model[0].parameters.probabilities[group].set_hospitalized_per_infectious(0.11 * (group + 1));
        model[0].parameters.probabilities[group].set_icu_per_hospitalized(0.12 * (group + 1));
    }
    epi::read_population_data_county(model, date, county, scaling_factor_inf, scaling_factor_icu, path);

    std::vector<double> sus   = {10284.4, 19086.2, 73805.3, 82522.6, 43731.9, 15620.2};
    std::vector<double> exp   = {0.571429, 3.8, 14.8286, 12.9429, 2.21429, 1.85714};
    std::vector<double> car   = {0.557143, 3.51429, 15.3857, 12.6571, 2.28571, 1.94286};
    std::vector<double> inf   = {0.291429, 1.93714, 5.79714, 2.45714, 0, 0};
    std::vector<double> hosp  = {0.0942857, 0.691429, 4.90286, 5.34286, 1.41429, 2.45143};
    std::vector<double> icu   = {0.0769231, 0.307692, 0.692308, 1.23077, 1.92308, 2.76923};
    std::vector<double> rec   = {35, 108.571, 640.143, 573.429, 180.429, 75.5714};
    std::vector<double> death = {0, 0, 0, 0, 10, 14.4286};

    for (size_t i = 0; i < 6; i++) {
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::S), sus[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::E), exp[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::C), car[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::I), inf[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::H), hosp[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::U), icu[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::R), rec[i], 1e-1);
        EXPECT_NEAR(model[0].populations.get((epi::AgeGroup6)i, epi::InfectionType::D), death[i], 1e-1);
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

    EXPECT_THAT(read_ids, testing::ElementsAreArray(true_ids));
}
