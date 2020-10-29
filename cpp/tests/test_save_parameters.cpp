#include "load_test_data.h"
#include "test_data_dir.h"
#include "epidemiology/secir/secir.h"
#include <epidemiology/secir/parameter_space.h>
#include <epidemiology/secir/parameter_studies.h>
#include <epidemiology_io/secir_result_io.h>
#include <epidemiology_io/secir_parameters_io.h>
#include <epidemiology/migration/migration.h>
#include <distributions_helpers.h>
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

    epi::ContactFrequencyMatrix& cont_freq_matrix = params.get_contact_patterns();
    epi::Damping dummy(30., 0.3);
    for (int i = 0; i < static_cast<int>(num_groups); i++) {
        for (int j = i; j < static_cast<int>(num_groups); j++) {
            cont_freq_matrix.set_cont_freq(fact * cont_freq, i, j);
            cont_freq_matrix.add_damping(dummy, i, j);
        }
    }
    cont_freq_matrix.add_damping(epi::Damping{35, 0.2}, 0, 0);
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

        for (size_t j = 0; j < num_groups; j++) {
            ASSERT_EQ(contact.get_cont_freq_mat().get_cont_freq(static_cast<int>(i), static_cast<int>(j)),
                      read_contact.get_cont_freq_mat().get_cont_freq(static_cast<int>(i), static_cast<int>(j)));

            auto& dampings_vector = contact.get_cont_freq_mat()
                                        .get_dampings(static_cast<int>(i), static_cast<int>(j))
                                        .get_dampings_vector();
            auto& cmp_dampings_vector = read_contact.get_cont_freq_mat()
                                            .get_dampings(static_cast<int>(i), static_cast<int>(j))
                                            .get_dampings_vector();
            ASSERT_THAT(dampings_vector, testing::ContainerEq(cmp_dampings_vector));
        }

        for (size_t j = 0; j < num_groups; j++) {
            ASSERT_EQ(contact.get_cont_freq_mat().get_cont_freq(static_cast<int>(i), static_cast<int>(j)),
                      read_contact.get_cont_freq_mat().get_cont_freq(static_cast<int>(i), static_cast<int>(j)));

            auto& dampings_vector = contact.get_cont_freq_mat()
                                        .get_dampings(static_cast<int>(i), static_cast<int>(j))
                                        .get_dampings_vector();
            auto& cmp_dampings_vector = read_contact.get_cont_freq_mat()
                                            .get_dampings(static_cast<int>(i), static_cast<int>(j))
                                            .get_dampings_vector();
            ASSERT_THAT(dampings_vector, testing::ContainerEq(cmp_dampings_vector));
        }
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

    epi::ContactFrequencyMatrix& cont_freq_matrix = model.parameters.get_contact_patterns();
    epi::Damping dummy(30., 0.3);
    for (int i = 0; i < static_cast<int>(num_groups); i++) {
        for (int j = i; j < static_cast<int>(num_groups); j++) {
            cont_freq_matrix.set_cont_freq(fact * cont_freq, i, j);
            cont_freq_matrix.add_damping(dummy, i, j);
        }
    }
    cont_freq_matrix.add_damping(epi::Damping{35, 0.2}, 0, 0);
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

        for (size_t j = 0; j < num_groups; j++) {
            ASSERT_EQ(contact.get_cont_freq_mat().get_cont_freq(static_cast<int>(i), static_cast<int>(j)),
                      read_contact.get_cont_freq_mat().get_cont_freq(static_cast<int>(i), static_cast<int>(j)));

            auto& dampings_vector = contact.get_cont_freq_mat()
                                        .get_dampings(static_cast<int>(i), static_cast<int>(j))
                                        .get_dampings_vector();
            auto& cmp_dampings_vector = read_contact.get_cont_freq_mat()
                                            .get_dampings(static_cast<int>(i), static_cast<int>(j))
                                            .get_dampings_vector();
            ASSERT_THAT(dampings_vector, testing::ContainerEq(cmp_dampings_vector));
        }
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

    epi::ContactFrequencyMatrix& cont_freq_matrix = model.parameters.get_contact_patterns();
    epi::Damping dummy(30., 0.3);
    for (int i = 0; i < static_cast<int>(num_groups); i++) {
        for (int j = i; j < static_cast<int>(num_groups); j++) {
            cont_freq_matrix.set_cont_freq(fact * cont_freq, i, j);
            cont_freq_matrix.add_damping(dummy, i, j);
            if (j > i) {
                cont_freq_matrix.set_cont_freq(fact * cont_freq, j, i);
            }
        }
    }

    epi::set_params_distributions_normal(model, t0, tmax, 0.15);

    epi::Graph<epi::SecirModel<epi::AgeGroup2>, epi::MigrationEdge> graph;
    graph.add_node(model);
    graph.add_node(model);
    graph.add_edge(0, 1, Eigen::VectorXd::Constant(model.populations.get_num_compartments(), 0.01));
    graph.add_edge(1, 0, Eigen::VectorXd::Constant(model.populations.get_num_compartments(), 0.01));

    epi::write_graph(graph);

    epi::Graph<epi::SecirModel<epi::AgeGroup2>, epi::MigrationEdge> graph_read =
        epi::read_graph<epi::SecirModel<epi::AgeGroup2>>();

    auto num_nodes = graph.nodes().size();
    auto num_edges = graph.edges().size();

    ASSERT_EQ(num_nodes, graph_read.nodes().size());
    ASSERT_EQ(num_edges, graph_read.edges().size());

    for (size_t node = 0; node < num_nodes; node++) {
        epi::SecirModel<epi::AgeGroup2> graph_model = graph.nodes()[0];
        epi::ContactFrequencyMatrix graph_cont_freq = graph_model.parameters.get_contact_patterns();

        epi::SecirModel<epi::AgeGroup2> graph_read_model = graph_read.nodes()[0];
        epi::ContactFrequencyMatrix graph_read_cont_freq = graph_read_model.parameters.get_contact_patterns();

        ASSERT_EQ(num_groups, graph_read_cont_freq.get_size());
        ASSERT_EQ(graph_model.populations.get_num_compartments(), graph_read_model.populations.get_num_compartments());

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

            for (size_t contact_group = 0; contact_group < num_groups; contact_group++) {
                ASSERT_EQ(graph_cont_freq.get_cont_freq(static_cast<int>(group), static_cast<int>(contact_group)),
                          graph_read_cont_freq.get_cont_freq(static_cast<int>(group), static_cast<int>(contact_group)));

                auto& dampings_v =
                    graph_cont_freq.get_dampings(static_cast<int>(group), static_cast<int>(contact_group))
                        .get_dampings_vector();
                auto& cmp_dampings_v =
                    graph_read_cont_freq.get_dampings(static_cast<int>(group), static_cast<int>(contact_group))
                        .get_dampings_vector();
                ASSERT_THAT(dampings_v, testing::ContainerEq(cmp_dampings_v));
            }

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

TEST(TestSaveParameters, ReadPopulationDataAllAges)
{
    epi::SecirModel<epi::AgeGroup1> model;
    std::vector<double> ranges = {100};

    std::string path = TEST_DATA_DIR;
    epi::read_population_data_germany(model, ranges, 5, 5, path);

    ASSERT_EQ(model.populations.get((epi::AgeGroup1)0, epi::InfectionType::I), 0);
    ASSERT_EQ(model.populations.get((epi::AgeGroup1)0, epi::InfectionType::D), 8626);
    ASSERT_EQ(model.populations.get((epi::AgeGroup1)0, epi::InfectionType::R), 160148);
    ASSERT_EQ(model.populations.get((epi::AgeGroup1)0, epi::InfectionType::U), 1937);
}

TEST(TestSaveParameters, ReadPopulationDataRKIAges)
{
    epi::SecirModel<epi::AgeGroup6> model;
    std::vector<double> ranges = {5., 10., 20., 25., 20., 20.};

    std::string path = TEST_DATA_DIR;
    epi::read_population_data_germany(model, ranges, 5, 5, path);

    std::vector<std::string> age_names = {"A00-A04", "A05-A14", "A15-A34", "A35-A59", "A60-A79", "A80+", "unknown"};

    std::vector<double> infected  = {0, 0, 0, 0, 0, 0};
    std::vector<double> deaths    = {1, 0, 18, 391, 2791, 5425};
    std::vector<double> recovered = {1516, 3656, 41947, 70301, 29224, 13504};

    for (size_t i = 0; i < ranges.size(); i++) {
        ASSERT_EQ(model.populations.get((epi::AgeGroup6)i, epi::InfectionType::I), infected[i]);
        ASSERT_EQ(model.populations.get((epi::AgeGroup6)i, epi::InfectionType::D), deaths[i]);
        ASSERT_EQ(model.populations.get((epi::AgeGroup6)i, epi::InfectionType::R), recovered[i]);
        ASSERT_EQ(model.populations.get((epi::AgeGroup6)i, epi::InfectionType::U), 1937 / (double)ranges.size());
    }
}

TEST(TestSaveParameters, ReadPopulationDataMultipleAges)
{
    epi::SecirModel<epi::AgeGroup8> model;
    std::vector<double> ranges = {5., 15., 15., 10., 10., 10., 10., 25.};

    std::string path = TEST_DATA_DIR;
    epi::read_population_data_germany(model, ranges, 5, 5, path);

    double infected_param  = 0.;
    double deaths_param    = 0.;
    double recovered_param = 0.;
    double icu_param       = 0.;

    for (size_t i = 0; i < ranges.size(); i++) {
        infected_param += model.populations.get((epi::AgeGroup8)i, epi::InfectionType::I);
        deaths_param += model.populations.get((epi::AgeGroup8)i, epi::InfectionType::D);
        recovered_param += model.populations.get((epi::AgeGroup8)i, epi::InfectionType::R);
        icu_param += model.populations.get((epi::AgeGroup8)i, epi::InfectionType::U);
    }

    std::vector<double> infected  = {0, 0, 0, 0, 0, 0};
    std::vector<double> deaths    = {1, 0, 18, 391, 2791, 5425};
    std::vector<double> recovered = {1516, 3656, 41947, 70301, 29224, 13504};

    ASSERT_EQ(model.populations.get((epi::AgeGroup8)0, epi::InfectionType::I), infected[0]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)0, epi::InfectionType::D), deaths[0]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)0, epi::InfectionType::R), recovered[0]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)0, epi::InfectionType::U), 1937 / (double)ranges.size());

    ASSERT_EQ(model.populations.get((epi::AgeGroup8)1, epi::InfectionType::I), infected[1] + (1.0 / 4.0) * infected[2]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)1, epi::InfectionType::D), deaths[1] + (1.0 / 4.0) * deaths[2]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)1, epi::InfectionType::R), recovered[1] + (1.0 / 4.0) * recovered[2]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)1, epi::InfectionType::U), 1937 / (double)ranges.size());

    ASSERT_EQ(model.populations.get((epi::AgeGroup8)2, epi::InfectionType::I), (3.0 / 4.0) * infected[2]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)2, epi::InfectionType::D), (3.0 / 4.0) * deaths[2]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)2, epi::InfectionType::R), (3.0 / 4.0) * recovered[2]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)2, epi::InfectionType::U), 1937 / (double)ranges.size());

    ASSERT_EQ(model.populations.get((epi::AgeGroup8)3, epi::InfectionType::I), (2.0 / 5.0) * infected[3]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)3, epi::InfectionType::D), (2.0 / 5.0) * deaths[3]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)3, epi::InfectionType::R), (2.0 / 5.0) * recovered[3]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)3, epi::InfectionType::U), 1937 / (double)ranges.size());

    ASSERT_EQ(model.populations.get((epi::AgeGroup8)4, epi::InfectionType::I), (2.0 / 5.0) * infected[3]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)4, epi::InfectionType::D), (2.0 / 5.0) * deaths[3]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)4, epi::InfectionType::R), (2.0 / 5.0) * recovered[3]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)4, epi::InfectionType::U), 1937 / (double)ranges.size());

    ASSERT_EQ(model.populations.get((epi::AgeGroup8)5, epi::InfectionType::I),
              (1.0 / 5.0) * infected[3] + (1.0 / 4.0) * infected[4]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)5, epi::InfectionType::D),
              (1.0 / 5.0) * deaths[3] + (1.0 / 4.0) * deaths[4]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)5, epi::InfectionType::R),
              (1.0 / 5.0) * recovered[3] + (1.0 / 4.0) * recovered[4]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)5, epi::InfectionType::U), 1937 / (double)ranges.size());

    ASSERT_EQ(model.populations.get((epi::AgeGroup8)6, epi::InfectionType::I), (2.0 / 4.0) * infected[4]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)6, epi::InfectionType::D), (2.0 / 4.0) * deaths[4]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)6, epi::InfectionType::R), (2.0 / 4.0) * recovered[4]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)6, epi::InfectionType::U), 1937 / (double)ranges.size());

    ASSERT_EQ(model.populations.get((epi::AgeGroup8)7, epi::InfectionType::I), (1.0 / 4.0) * infected[4] + infected[5]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)7, epi::InfectionType::D), (1.0 / 4.0) * deaths[4] + deaths[5]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)7, epi::InfectionType::R), (1.0 / 4.0) * recovered[4] + recovered[5]);
    ASSERT_EQ(model.populations.get((epi::AgeGroup8)7, epi::InfectionType::U), 1937 / (double)ranges.size());

    EXPECT_NEAR(0, infected_param, 1e-6);
    EXPECT_NEAR(8626, deaths_param, 1e-6);
    EXPECT_NEAR(160148, recovered_param, 1e-6);
    EXPECT_NEAR(1937, icu_param, 1e-6);
}
