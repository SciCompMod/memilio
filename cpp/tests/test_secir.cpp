#include "load_test_data.h"
#include "epidemiology/secir/secir.h"
#include "epidemiology/secir/analyze_result.h"
#include <distributions_helpers.h>
#include <gtest/gtest.h>

TEST(TestSecir, compareWithPreviousRun)
{
    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           ticu2death = 5;

    double cont_freq = 10, inf_prob = 0.05, carr_infec = 1, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2,
           theta = 0.25;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    epi::SecirModel<epi::AgeGroup1> model;

    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs

    model.parameters.times[0].set_incubation(tinc);
    model.parameters.times[0].set_infectious_mild(tinfmild);
    model.parameters.times[0].set_serialinterval(tserint);
    model.parameters.times[0].set_hospitalized_to_home(thosp2home);
    model.parameters.times[0].set_home_to_hospitalized(thome2hosp);
    model.parameters.times[0].set_hospitalized_to_icu(thosp2icu);
    model.parameters.times[0].set_icu_to_home(ticu2home);
    model.parameters.times[0].set_icu_to_death(ticu2death);

    epi::ContactFrequencyMatrix& cont_freq_matrix = model.parameters.get_contact_patterns();
    cont_freq_matrix.set_cont_freq(cont_freq, 0, 0);
    epi::Damping dummy(30., 0.3);
    cont_freq_matrix.add_damping(dummy, 0, 0);

    model.populations.set_total(nb_total_t0);
    model.populations.set(nb_exp_t0, (epi::AgeGroup1)0, epi::InfectionType::E);
    model.populations.set(nb_car_t0, (epi::AgeGroup1)0, epi::InfectionType::C);
    model.populations.set(nb_inf_t0, (epi::AgeGroup1)0, epi::InfectionType::I);
    model.populations.set(nb_hosp_t0, (epi::AgeGroup1)0, epi::InfectionType::H);
    model.populations.set(nb_icu_t0, (epi::AgeGroup1)0, epi::InfectionType::U);
    model.populations.set(nb_rec_t0, (epi::AgeGroup1)0, epi::InfectionType::R);
    model.populations.set(nb_dead_t0, (epi::AgeGroup1)0, epi::InfectionType::D);
    model.populations.set_difference_from_total(nb_total_t0, (epi::AgeGroup1)0, epi::InfectionType::S);

    model.parameters.probabilities[0].set_infection_from_contact(inf_prob);
    model.parameters.probabilities[0].set_carrier_infectability(carr_infec);
    model.parameters.probabilities[0].set_asymp_per_infectious(alpha);
    model.parameters.probabilities[0].set_risk_from_symptomatic(beta);
    model.parameters.probabilities[0].set_hospitalized_per_infectious(rho);
    model.parameters.probabilities[0].set_icu_per_hospitalized(theta);
    model.parameters.probabilities[0].set_dead_per_icu(delta);

    model.apply_constraints();

    epi::TimeSeries<double> secihurd = simulate(t0, tmax, dt, model);

    auto compare = load_test_data_csv<double>("secihurd-compare.csv");

    ASSERT_EQ(compare.size(), static_cast<size_t>(secihurd.get_num_time_points()));
    for (size_t i = 0; i < compare.size(); i++) {
        ASSERT_EQ(compare[i].size(), static_cast<size_t>(secihurd.get_num_elements()) + 1) << "at row " << i;
        EXPECT_NEAR(secihurd.get_time(i), compare[i][0], 1e-10) << "at row " << i;
        for (size_t j = 1; j < compare[i].size(); j++) {
            EXPECT_NEAR(secihurd.get_value(i)[j - 1], compare[i][j], 1e-10) << " at row " << i;
        }
    }
}

TEST(TestSecir, testParamConstructors)
{
    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           ticu2death = 5;

    double cont_freq = 10, inf_prob = 0.05, carr_infec = 0.67, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2,
           theta = 0.24;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 54, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 11, nb_dead_t0 = 0;

    double icu_cap   = 4444;
    double start_day = 30, seasonality = 0.3;

    epi::SecirModel<epi::AgeGroup1> model;

    model.parameters.set_icu_capacity(icu_cap);

    model.parameters.set_start_day(start_day);
    model.parameters.set_seasonality(seasonality);

    model.parameters.times[0].set_incubation(tinc);
    model.parameters.times[0].set_infectious_mild(tinfmild);
    model.parameters.times[0].set_serialinterval(tserint);
    model.parameters.times[0].set_hospitalized_to_home(thosp2home);
    model.parameters.times[0].set_home_to_hospitalized(thome2hosp);
    model.parameters.times[0].set_hospitalized_to_icu(thosp2icu);
    model.parameters.times[0].set_icu_to_home(ticu2home);
    model.parameters.times[0].set_icu_to_death(ticu2death);

    model.populations.set_total(nb_total_t0);
    model.populations.set(nb_exp_t0, (epi::AgeGroup1)0, epi::InfectionType::E);
    model.populations.set(nb_car_t0, (epi::AgeGroup1)0, epi::InfectionType::C);
    model.populations.set(nb_inf_t0, (epi::AgeGroup1)0, epi::InfectionType::I);
    model.populations.set(nb_hosp_t0, (epi::AgeGroup1)0, epi::InfectionType::H);
    model.populations.set(nb_icu_t0, (epi::AgeGroup1)0, epi::InfectionType::U);
    model.populations.set(nb_rec_t0, (epi::AgeGroup1)0, epi::InfectionType::R);
    model.populations.set(nb_dead_t0, (epi::AgeGroup1)0, epi::InfectionType::D);
    model.populations.set_difference_from_total(nb_total_t0, (epi::AgeGroup1)0, epi::InfectionType::S);

    model.parameters.probabilities[0].set_infection_from_contact(inf_prob);
    model.parameters.probabilities[0].set_carrier_infectability(carr_infec);
    model.parameters.probabilities[0].set_asymp_per_infectious(alpha);
    model.parameters.probabilities[0].set_risk_from_symptomatic(beta);
    model.parameters.probabilities[0].set_hospitalized_per_infectious(rho);
    model.parameters.probabilities[0].set_icu_per_hospitalized(theta);
    model.parameters.probabilities[0].set_dead_per_icu(delta);

    epi::ContactFrequencyMatrix& cont_freq_matrix = model.parameters.get_contact_patterns();
    cont_freq_matrix.set_cont_freq(cont_freq, 0, 0);
    epi::Damping dummy(30., 0.3);
    cont_freq_matrix.add_damping(dummy, 0, 0);

    epi::SecirModel<epi::AgeGroup1> model2{model}; // copy constructor

    EXPECT_EQ(model.parameters.get_icu_capacity(), model2.parameters.get_icu_capacity());
    EXPECT_EQ(model.parameters.get_start_day(), model2.parameters.get_start_day());
    EXPECT_EQ(model.parameters.get_seasonality(), model2.parameters.get_seasonality());

    EXPECT_EQ(model.populations.get_total(), model2.populations.get_total());
    EXPECT_EQ(model.populations.get((epi::AgeGroup1)0, epi::InfectionType::S),
              model.populations.get((epi::AgeGroup1)0, epi::InfectionType::S));
    EXPECT_EQ(model.populations.get((epi::AgeGroup1)0, epi::InfectionType::E),
              model.populations.get((epi::AgeGroup1)0, epi::InfectionType::E));
    EXPECT_EQ(model.populations.get((epi::AgeGroup1)0, epi::InfectionType::C),
              model.populations.get((epi::AgeGroup1)0, epi::InfectionType::C));
    EXPECT_EQ(model.populations.get((epi::AgeGroup1)0, epi::InfectionType::I),
              model.populations.get((epi::AgeGroup1)0, epi::InfectionType::I));
    EXPECT_EQ(model.populations.get((epi::AgeGroup1)0, epi::InfectionType::H),
              model.populations.get((epi::AgeGroup1)0, epi::InfectionType::H));
    EXPECT_EQ(model.populations.get((epi::AgeGroup1)0, epi::InfectionType::U),
              model.populations.get((epi::AgeGroup1)0, epi::InfectionType::U));
    EXPECT_EQ(model.populations.get((epi::AgeGroup1)0, epi::InfectionType::R),
              model.populations.get((epi::AgeGroup1)0, epi::InfectionType::R));
    EXPECT_EQ(model.populations.get((epi::AgeGroup1)0, epi::InfectionType::D),
              model.populations.get((epi::AgeGroup1)0, epi::InfectionType::D));

    EXPECT_EQ(model.parameters.times[0].get_incubation(), model2.parameters.times[0].get_incubation());
    EXPECT_EQ(model.parameters.times[0].get_serialinterval(), model2.parameters.times[0].get_serialinterval());
    EXPECT_EQ(model.parameters.times[0].get_infectious_mild(), model2.parameters.times[0].get_infectious_mild());
    EXPECT_EQ(model.parameters.times[0].get_infectious_asymp(), model2.parameters.times[0].get_infectious_asymp());
    EXPECT_EQ(model.parameters.times[0].get_home_to_hospitalized(),
              model2.parameters.times[0].get_home_to_hospitalized());
    EXPECT_EQ(model.parameters.times[0].get_hospitalized_to_home(),
              model2.parameters.times[0].get_hospitalized_to_home());
    EXPECT_EQ(model.parameters.times[0].get_hospitalized_to_icu(),
              model2.parameters.times[0].get_hospitalized_to_icu());
    EXPECT_EQ(model.parameters.times[0].get_icu_to_dead(), model2.parameters.times[0].get_icu_to_dead());
    EXPECT_EQ(model.parameters.times[0].get_icu_to_home(), model2.parameters.times[0].get_icu_to_home());
    EXPECT_EQ(model.parameters.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35),
              model2.parameters.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35));

    EXPECT_EQ(model.parameters.probabilities[0].get_infection_from_contact(),
              model2.parameters.probabilities[0].get_infection_from_contact());
    EXPECT_EQ(model.parameters.probabilities[0].get_carrier_infectability(),
              model2.parameters.probabilities[0].get_carrier_infectability());
    EXPECT_EQ(model.parameters.probabilities[0].get_risk_from_symptomatic(),
              model2.parameters.probabilities[0].get_risk_from_symptomatic());
    EXPECT_EQ(model.parameters.probabilities[0].get_asymp_per_infectious(),
              model2.parameters.probabilities[0].get_asymp_per_infectious());
    EXPECT_EQ(model.parameters.probabilities[0].get_hospitalized_per_infectious(),
              model2.parameters.probabilities[0].get_hospitalized_per_infectious());
    EXPECT_EQ(model.parameters.probabilities[0].get_icu_per_hospitalized(),
              model2.parameters.probabilities[0].get_icu_per_hospitalized());
    EXPECT_EQ(model.parameters.probabilities[0].get_dead_per_icu(),
              model2.parameters.probabilities[0].get_dead_per_icu());
    EXPECT_EQ(model.parameters.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35),
              model2.parameters.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35));

    epi::SecirModel<epi::AgeGroup1> model3 = std::move(model2); // move constructor

    EXPECT_EQ(model.parameters.get_icu_capacity(), model3.parameters.get_icu_capacity());
    EXPECT_EQ(model.parameters.get_start_day(), model3.parameters.get_start_day());
    EXPECT_EQ(model.parameters.get_seasonality(), model3.parameters.get_seasonality());

    EXPECT_EQ(model3.populations.get_total(), model.populations.get_total());
    EXPECT_EQ(model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::S),
              model.populations.get((epi::AgeGroup1)0, epi::InfectionType::S));
    EXPECT_EQ(model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::E),
              model.populations.get((epi::AgeGroup1)0, epi::InfectionType::E));
    EXPECT_EQ(model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::C),
              model.populations.get((epi::AgeGroup1)0, epi::InfectionType::C));
    EXPECT_EQ(model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::I),
              model.populations.get((epi::AgeGroup1)0, epi::InfectionType::I));
    EXPECT_EQ(model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::H),
              model.populations.get((epi::AgeGroup1)0, epi::InfectionType::H));
    EXPECT_EQ(model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::U),
              model.populations.get((epi::AgeGroup1)0, epi::InfectionType::U));
    EXPECT_EQ(model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::R),
              model.populations.get((epi::AgeGroup1)0, epi::InfectionType::R));
    EXPECT_EQ(model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::D),
              model.populations.get((epi::AgeGroup1)0, epi::InfectionType::D));

    EXPECT_EQ(model3.parameters.times[0].get_incubation(), model.parameters.times[0].get_incubation());
    EXPECT_EQ(model3.parameters.times[0].get_serialinterval(), model.parameters.times[0].get_serialinterval());
    EXPECT_EQ(model3.parameters.times[0].get_infectious_mild(), model.parameters.times[0].get_infectious_mild());
    EXPECT_EQ(model3.parameters.times[0].get_infectious_asymp(), model.parameters.times[0].get_infectious_asymp());
    EXPECT_EQ(model3.parameters.times[0].get_home_to_hospitalized(),
              model.parameters.times[0].get_home_to_hospitalized());
    EXPECT_EQ(model3.parameters.times[0].get_hospitalized_to_home(),
              model.parameters.times[0].get_hospitalized_to_home());
    EXPECT_EQ(model3.parameters.times[0].get_hospitalized_to_icu(),
              model.parameters.times[0].get_hospitalized_to_icu());
    EXPECT_EQ(model3.parameters.times[0].get_icu_to_dead(), model.parameters.times[0].get_icu_to_dead());
    EXPECT_EQ(model3.parameters.times[0].get_icu_to_home(), model.parameters.times[0].get_icu_to_home());

    EXPECT_EQ(model3.parameters.probabilities[0].get_infection_from_contact(),
              model.parameters.probabilities[0].get_infection_from_contact());
    EXPECT_EQ(model3.parameters.probabilities[0].get_carrier_infectability(),
              model.parameters.probabilities[0].get_carrier_infectability());
    EXPECT_EQ(model3.parameters.probabilities[0].get_risk_from_symptomatic(),
              model.parameters.probabilities[0].get_risk_from_symptomatic());
    EXPECT_EQ(model3.parameters.probabilities[0].get_asymp_per_infectious(),
              model.parameters.probabilities[0].get_asymp_per_infectious());
    EXPECT_EQ(model3.parameters.probabilities[0].get_hospitalized_per_infectious(),
              model.parameters.probabilities[0].get_hospitalized_per_infectious());
    EXPECT_EQ(model3.parameters.probabilities[0].get_icu_per_hospitalized(),
              model.parameters.probabilities[0].get_icu_per_hospitalized());
    EXPECT_EQ(model3.parameters.probabilities[0].get_dead_per_icu(),
              model.parameters.probabilities[0].get_dead_per_icu());

    EXPECT_EQ(model.parameters.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35),
              model3.parameters.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35));

    epi::SecirModel<epi::AgeGroup1> model4 = model3; // copy assignment constructor

    EXPECT_EQ(model4.parameters.get_icu_capacity(), model3.parameters.get_icu_capacity());
    EXPECT_EQ(model4.parameters.get_start_day(), model3.parameters.get_start_day());
    EXPECT_EQ(model4.parameters.get_seasonality(), model3.parameters.get_seasonality());

    EXPECT_EQ(model3.populations.get_total(), model4.populations.get_total());
    EXPECT_EQ(model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::S),
              model4.populations.get((epi::AgeGroup1)0, epi::InfectionType::S));
    EXPECT_EQ(model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::E),
              model4.populations.get((epi::AgeGroup1)0, epi::InfectionType::E));
    EXPECT_EQ(model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::C),
              model4.populations.get((epi::AgeGroup1)0, epi::InfectionType::C));
    EXPECT_EQ(model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::I),
              model4.populations.get((epi::AgeGroup1)0, epi::InfectionType::I));
    EXPECT_EQ(model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::H),
              model4.populations.get((epi::AgeGroup1)0, epi::InfectionType::H));
    EXPECT_EQ(model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::U),
              model4.populations.get((epi::AgeGroup1)0, epi::InfectionType::U));
    EXPECT_EQ(model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::R),
              model4.populations.get((epi::AgeGroup1)0, epi::InfectionType::R));
    EXPECT_EQ(model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::D),
              model4.populations.get((epi::AgeGroup1)0, epi::InfectionType::D));

    EXPECT_EQ(model3.parameters.times[0].get_incubation(), model4.parameters.times[0].get_incubation());
    EXPECT_EQ(model3.parameters.times[0].get_serialinterval(), model4.parameters.times[0].get_serialinterval());
    EXPECT_EQ(model3.parameters.times[0].get_infectious_mild(), model4.parameters.times[0].get_infectious_mild());
    EXPECT_EQ(model3.parameters.times[0].get_infectious_asymp(), model4.parameters.times[0].get_infectious_asymp());
    EXPECT_EQ(model3.parameters.times[0].get_home_to_hospitalized(),
              model4.parameters.times[0].get_home_to_hospitalized());
    EXPECT_EQ(model3.parameters.times[0].get_hospitalized_to_home(),
              model4.parameters.times[0].get_hospitalized_to_home());
    EXPECT_EQ(model3.parameters.times[0].get_hospitalized_to_icu(),
              model4.parameters.times[0].get_hospitalized_to_icu());
    EXPECT_EQ(model3.parameters.times[0].get_icu_to_dead(), model4.parameters.times[0].get_icu_to_dead());
    EXPECT_EQ(model3.parameters.times[0].get_icu_to_home(), model4.parameters.times[0].get_icu_to_home());

    EXPECT_EQ(model3.parameters.probabilities[0].get_infection_from_contact(),
              model4.parameters.probabilities[0].get_infection_from_contact());
    EXPECT_EQ(model3.parameters.probabilities[0].get_carrier_infectability(),
              model4.parameters.probabilities[0].get_carrier_infectability());
    EXPECT_EQ(model3.parameters.probabilities[0].get_risk_from_symptomatic(),
              model4.parameters.probabilities[0].get_risk_from_symptomatic());
    EXPECT_EQ(model3.parameters.probabilities[0].get_asymp_per_infectious(),
              model4.parameters.probabilities[0].get_asymp_per_infectious());
    EXPECT_EQ(model3.parameters.probabilities[0].get_hospitalized_per_infectious(),
              model4.parameters.probabilities[0].get_hospitalized_per_infectious());
    EXPECT_EQ(model3.parameters.probabilities[0].get_icu_per_hospitalized(),
              model4.parameters.probabilities[0].get_icu_per_hospitalized());
    EXPECT_EQ(model3.parameters.probabilities[0].get_dead_per_icu(),
              model4.parameters.probabilities[0].get_dead_per_icu());

    EXPECT_EQ(model4.parameters.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35),
              model3.parameters.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35));

    epi::SecirModel<epi::AgeGroup1> model5 = std::move(model4); // move assignment constructor

    EXPECT_EQ(model5.parameters.get_icu_capacity(), model3.parameters.get_icu_capacity());
    EXPECT_EQ(model5.parameters.get_start_day(), model3.parameters.get_start_day());
    EXPECT_EQ(model5.parameters.get_seasonality(), model3.parameters.get_seasonality());

    EXPECT_EQ(model5.populations.get_total(), model3.populations.get_total());
    EXPECT_EQ(model5.populations.get((epi::AgeGroup1)0, epi::InfectionType::S),
              model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::S));
    EXPECT_EQ(model5.populations.get((epi::AgeGroup1)0, epi::InfectionType::E),
              model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::E));
    EXPECT_EQ(model5.populations.get((epi::AgeGroup1)0, epi::InfectionType::C),
              model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::C));
    EXPECT_EQ(model5.populations.get((epi::AgeGroup1)0, epi::InfectionType::I),
              model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::I));
    EXPECT_EQ(model5.populations.get((epi::AgeGroup1)0, epi::InfectionType::H),
              model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::H));
    EXPECT_EQ(model5.populations.get((epi::AgeGroup1)0, epi::InfectionType::U),
              model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::U));
    EXPECT_EQ(model5.populations.get((epi::AgeGroup1)0, epi::InfectionType::R),
              model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::R));
    EXPECT_EQ(model5.populations.get((epi::AgeGroup1)0, epi::InfectionType::D),
              model3.populations.get((epi::AgeGroup1)0, epi::InfectionType::D));

    EXPECT_EQ(model5.parameters.times[0].get_incubation(), model3.parameters.times[0].get_incubation());
    EXPECT_EQ(model5.parameters.times[0].get_serialinterval(), model3.parameters.times[0].get_serialinterval());
    EXPECT_EQ(model5.parameters.times[0].get_infectious_mild(), model3.parameters.times[0].get_infectious_mild());
    EXPECT_EQ(model5.parameters.times[0].get_infectious_asymp(), model3.parameters.times[0].get_infectious_asymp());
    EXPECT_EQ(model5.parameters.times[0].get_home_to_hospitalized(),
              model3.parameters.times[0].get_home_to_hospitalized());
    EXPECT_EQ(model5.parameters.times[0].get_hospitalized_to_home(),
              model3.parameters.times[0].get_hospitalized_to_home());
    EXPECT_EQ(model5.parameters.times[0].get_hospitalized_to_icu(),
              model3.parameters.times[0].get_hospitalized_to_icu());
    EXPECT_EQ(model5.parameters.times[0].get_icu_to_dead(), model3.parameters.times[0].get_icu_to_dead());
    EXPECT_EQ(model5.parameters.times[0].get_icu_to_home(), model3.parameters.times[0].get_icu_to_home());

    EXPECT_EQ(model5.parameters.probabilities[0].get_infection_from_contact(),
              model3.parameters.probabilities[0].get_infection_from_contact());
    EXPECT_EQ(model5.parameters.probabilities[0].get_carrier_infectability(),
              model3.parameters.probabilities[0].get_carrier_infectability());
    EXPECT_EQ(model5.parameters.probabilities[0].get_risk_from_symptomatic(),
              model3.parameters.probabilities[0].get_risk_from_symptomatic());
    EXPECT_EQ(model5.parameters.probabilities[0].get_asymp_per_infectious(),
              model3.parameters.probabilities[0].get_asymp_per_infectious());
    EXPECT_EQ(model5.parameters.probabilities[0].get_hospitalized_per_infectious(),
              model3.parameters.probabilities[0].get_hospitalized_per_infectious());
    EXPECT_EQ(model5.parameters.probabilities[0].get_icu_per_hospitalized(),
              model3.parameters.probabilities[0].get_icu_per_hospitalized());
    EXPECT_EQ(model5.parameters.probabilities[0].get_dead_per_icu(),
              model3.parameters.probabilities[0].get_dead_per_icu());

    EXPECT_EQ(model5.parameters.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35),
              model3.parameters.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35));

    epi::ContactFrequencyMatrix contact_freq_matrix{2};

    for (int i = 0; i < 2; i++) {
        for (int j = i; j < 2; j++) {
            contact_freq_matrix.set_cont_freq(0.5 * cont_freq, i, j);
        }
    }

    epi::ContactFrequencyMatrix contact_freq_matrix2{contact_freq_matrix};

    for (int i = 0; i < 2; i++) {
        for (int j = i; j < 2; j++) {
            EXPECT_EQ(contact_freq_matrix.get_cont_freq(i, j), contact_freq_matrix2.get_cont_freq(i, j));
        }
    }

    epi::ContactFrequencyMatrix contact_freq_matrix3{std::move(contact_freq_matrix2)};

    for (int i = 0; i < 2; i++) {
        for (int j = i; j < 2; j++) {
            EXPECT_EQ(contact_freq_matrix3.get_cont_freq(i, j), contact_freq_matrix.get_cont_freq(i, j));
        }
    }

    epi::ContactFrequencyMatrix contact_freq_matrix4 = std::move(contact_freq_matrix3);

    for (int i = 0; i < 2; i++) {
        for (int j = i; j < 2; j++) {
            EXPECT_EQ(contact_freq_matrix4.get_cont_freq(i, j), contact_freq_matrix.get_cont_freq(i, j));
        }
    }
}

TEST(TestSecir, testSettersAndGetters)
{
    std::vector<epi::UncertainValue> vec;

    for (int i = 0; i < 26; i++) {
        epi::UncertainValue val = epi::UncertainValue(i);
        val.set_distribution(epi::ParameterDistributionNormal(i, 10 * i, 5 * i, i / 10.0));
        vec.push_back(val);
    }

    epi::SecirModel<epi::AgeGroup1> model;

    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs

    EXPECT_EQ(model.parameters.times[0].get_incubation().get_distribution().get(), nullptr);

    model.parameters.set_icu_capacity(vec[0]);

    model.parameters.times[0].set_incubation(vec[1]);
    model.parameters.times[0].set_infectious_mild(vec[2]);
    model.parameters.times[0].set_serialinterval(vec[3]);
    model.parameters.times[0].set_hospitalized_to_home(vec[4]);
    model.parameters.times[0].set_home_to_hospitalized(vec[5]);
    model.parameters.times[0].set_hospitalized_to_icu(vec[6]);
    model.parameters.times[0].set_icu_to_home(vec[7]);
    model.parameters.times[0].set_infectious_asymp(vec[8]);
    model.parameters.times[0].set_icu_to_death(vec[9]);

    model.populations.set(vec[10], (epi::AgeGroup1)0, epi::InfectionType::E);
    model.populations.set(vec[11], (epi::AgeGroup1)0, epi::InfectionType::C);
    model.populations.set(vec[12], (epi::AgeGroup1)0, epi::InfectionType::I);
    model.populations.set(vec[13], (epi::AgeGroup1)0, epi::InfectionType::H);
    model.populations.set(vec[14], (epi::AgeGroup1)0, epi::InfectionType::U);
    model.populations.set(vec[15], (epi::AgeGroup1)0, epi::InfectionType::R);
    model.populations.set(vec[16], (epi::AgeGroup1)0, epi::InfectionType::D);

    model.parameters.probabilities[0].set_infection_from_contact(vec[17]);
    model.parameters.probabilities[0].set_carrier_infectability(vec[18]);
    model.parameters.probabilities[0].set_asymp_per_infectious(vec[19]);
    model.parameters.probabilities[0].set_risk_from_symptomatic(vec[20]);
    model.parameters.probabilities[0].set_hospitalized_per_infectious(vec[21]);
    model.parameters.probabilities[0].set_icu_per_hospitalized(vec[22]);
    model.parameters.probabilities[0].set_dead_per_icu(vec[23]);

    EXPECT_NE(model.parameters.times[0].get_incubation().get_distribution().get(), nullptr);

    check_distribution(*vec[0].get_distribution(), *model.parameters.get_icu_capacity().get_distribution());

    model.parameters.set_start_day(vec[24]);
    model.parameters.set_seasonality(vec[25]);

    EXPECT_NE(model.parameters.times[0].get_incubation().get_distribution().get(), nullptr);

    check_distribution(*vec[1].get_distribution(), *model.parameters.times[0].get_incubation().get_distribution());
    check_distribution(*vec[2].get_distribution(), *model.parameters.times[0].get_infectious_mild().get_distribution());
    check_distribution(*vec[3].get_distribution(), *model.parameters.times[0].get_serialinterval().get_distribution());
    check_distribution(*vec[4].get_distribution(),
                       *model.parameters.times[0].get_hospitalized_to_home().get_distribution());
    check_distribution(*vec[5].get_distribution(),
                       *model.parameters.times[0].get_home_to_hospitalized().get_distribution());
    check_distribution(*vec[6].get_distribution(),
                       *model.parameters.times[0].get_hospitalized_to_icu().get_distribution());
    check_distribution(*vec[7].get_distribution(), *model.parameters.times[0].get_icu_to_home().get_distribution());
    check_distribution(*vec[8].get_distribution(),
                       *model.parameters.times[0].get_infectious_asymp().get_distribution());
    check_distribution(*vec[9].get_distribution(), *model.parameters.times[0].get_icu_to_dead().get_distribution());
    check_distribution(*vec[10].get_distribution(),
                       *model.populations.get((epi::AgeGroup1)0, epi::InfectionType::E).get_distribution());
    check_distribution(*vec[11].get_distribution(),
                       *model.populations.get((epi::AgeGroup1)0, epi::InfectionType::C).get_distribution());
    check_distribution(*vec[12].get_distribution(),
                       *model.populations.get((epi::AgeGroup1)0, epi::InfectionType::I).get_distribution());
    check_distribution(*vec[13].get_distribution(),
                       *model.populations.get((epi::AgeGroup1)0, epi::InfectionType::H).get_distribution());
    check_distribution(*vec[14].get_distribution(),
                       *model.populations.get((epi::AgeGroup1)0, epi::InfectionType::U).get_distribution());
    check_distribution(*vec[15].get_distribution(),
                       *model.populations.get((epi::AgeGroup1)0, epi::InfectionType::R).get_distribution());
    check_distribution(*vec[16].get_distribution(),
                       *model.populations.get((epi::AgeGroup1)0, epi::InfectionType::D).get_distribution());
    check_distribution(*vec[17].get_distribution(),
                       *model.parameters.probabilities[0].get_infection_from_contact().get_distribution());
    check_distribution(*vec[18].get_distribution(),
                       *model.parameters.probabilities[0].get_carrier_infectability().get_distribution());
    check_distribution(*vec[19].get_distribution(),
                       *model.parameters.probabilities[0].get_asymp_per_infectious().get_distribution());
    check_distribution(*vec[20].get_distribution(),
                       *model.parameters.probabilities[0].get_risk_from_symptomatic().get_distribution());
    check_distribution(*vec[21].get_distribution(),
                       *model.parameters.probabilities[0].get_hospitalized_per_infectious().get_distribution());
    check_distribution(*vec[22].get_distribution(),
                       *model.parameters.probabilities[0].get_icu_per_hospitalized().get_distribution());
    check_distribution(*vec[23].get_distribution(),
                       *model.parameters.probabilities[0].get_dead_per_icu().get_distribution());
    // no dist for start day
    check_distribution(*vec[25].get_distribution(), *model.parameters.get_seasonality().get_distribution());

    EXPECT_EQ(vec[0], model.parameters.get_icu_capacity());
    EXPECT_EQ(vec[1], model.parameters.times[0].get_incubation());
    EXPECT_EQ(vec[2], model.parameters.times[0].get_infectious_mild());
    EXPECT_EQ(vec[3], model.parameters.times[0].get_serialinterval());
    EXPECT_EQ(vec[4], model.parameters.times[0].get_hospitalized_to_home());
    EXPECT_EQ(vec[5], model.parameters.times[0].get_home_to_hospitalized());
    EXPECT_EQ(vec[6], model.parameters.times[0].get_hospitalized_to_icu());
    EXPECT_EQ(vec[7], model.parameters.times[0].get_icu_to_home());
    EXPECT_EQ(vec[8], model.parameters.times[0].get_infectious_asymp());
    EXPECT_EQ(vec[9], model.parameters.times[0].get_icu_to_dead());
    EXPECT_EQ(vec[10], model.populations.get((epi::AgeGroup1)0, epi::InfectionType::E));
    EXPECT_EQ(vec[11], model.populations.get((epi::AgeGroup1)0, epi::InfectionType::C));
    EXPECT_EQ(vec[12], model.populations.get((epi::AgeGroup1)0, epi::InfectionType::I));
    EXPECT_EQ(vec[13], model.populations.get((epi::AgeGroup1)0, epi::InfectionType::H));
    EXPECT_EQ(vec[14], model.populations.get((epi::AgeGroup1)0, epi::InfectionType::U));
    EXPECT_EQ(vec[15], model.populations.get((epi::AgeGroup1)0, epi::InfectionType::R));
    EXPECT_EQ(vec[16], model.populations.get((epi::AgeGroup1)0, epi::InfectionType::D));
    EXPECT_EQ(vec[17], model.parameters.probabilities[0].get_infection_from_contact());
    EXPECT_EQ(vec[18], model.parameters.probabilities[0].get_carrier_infectability());
    EXPECT_EQ(vec[19], model.parameters.probabilities[0].get_asymp_per_infectious());
    EXPECT_EQ(vec[20], model.parameters.probabilities[0].get_risk_from_symptomatic());
    EXPECT_EQ(vec[21], model.parameters.probabilities[0].get_hospitalized_per_infectious());
    EXPECT_EQ(vec[22], model.parameters.probabilities[0].get_icu_per_hospitalized());
    EXPECT_EQ(vec[23], model.parameters.probabilities[0].get_dead_per_icu());
    EXPECT_EQ(vec[24], model.parameters.get_start_day());
    EXPECT_EQ(vec[25], model.parameters.get_seasonality());
}

TEST(TestSecir, testValueConstraints)
{
    double tinc    = 5.1, // R_2^(-1)+R_3^(-1)
        tinfmild   = 5.86642, // 4-14  (=R4^(-1))
        tserint    = 5.08993, // 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        thosp2home = 11.6138, // 7-16 (=R5^(-1))
        thome2hosp = 4.45361, // 2.5-7 (=R6^(-1))
        thosp2icu  = 2.15791, // 1-3.5 (=R7^(-1))
        ticu2home  = 9.16291, // 5-16 (=R8^(-1))
        ticu2death = 5.90264; // 3.5-7 (=R5^(-1))

    double cont_freq = 10, // 0.2-0.75
        inf_prob = 0.064519, carr_infec = 0.56758,
           alpha = 2.124921, // 0.01-0.16
        beta     = 0.190609, // 0.05-0.5
        delta    = 0.245801, // 0.15-0.77
        rho      = 0.183693, // 0.1-0.35
        theta    = 0.185556; // 0.15-0.4

    double nb_total_t0 = 10000, nb_exp_t0 = -91, nb_inf_t0 = 39, nb_car_t0 = 36, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 8, nb_dead_t0 = 0;

    epi::SecirModel<epi::AgeGroup1> model;

    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs

    model.parameters.times[0].set_incubation(tinc);
    model.parameters.times[0].set_infectious_mild(tinfmild);
    model.parameters.times[0].set_serialinterval(tserint);
    model.parameters.times[0].set_hospitalized_to_home(thosp2home);
    model.parameters.times[0].set_home_to_hospitalized(thome2hosp);
    model.parameters.times[0].set_hospitalized_to_icu(thosp2icu);
    model.parameters.times[0].set_icu_to_home(ticu2home);
    model.parameters.times[0].set_icu_to_death(ticu2death);

    epi::ContactFrequencyMatrix& cont_freq_matrix = model.parameters.get_contact_patterns();
    cont_freq_matrix.set_cont_freq(cont_freq, 0, 0);
    epi::Damping dummy(30., 0.3);
    cont_freq_matrix.add_damping(dummy, 0, 0);

    model.populations.set_total(nb_total_t0);
    model.populations.set(nb_exp_t0, (epi::AgeGroup1)0, epi::InfectionType::E);
    model.populations.set(nb_car_t0, (epi::AgeGroup1)0, epi::InfectionType::C);
    model.populations.set(nb_inf_t0, (epi::AgeGroup1)0, epi::InfectionType::I);
    model.populations.set(nb_hosp_t0, (epi::AgeGroup1)0, epi::InfectionType::H);
    model.populations.set(nb_icu_t0, (epi::AgeGroup1)0, epi::InfectionType::U);
    model.populations.set(nb_rec_t0, (epi::AgeGroup1)0, epi::InfectionType::R);
    model.populations.set(nb_dead_t0, (epi::AgeGroup1)0, epi::InfectionType::D);
    model.populations.set_difference_from_total(nb_total_t0, (epi::AgeGroup1)0, epi::InfectionType::S);

    model.parameters.probabilities[0].set_infection_from_contact(inf_prob);
    model.parameters.probabilities[0].set_carrier_infectability(carr_infec);
    model.parameters.probabilities[0].set_asymp_per_infectious(alpha);
    model.parameters.probabilities[0].set_risk_from_symptomatic(beta);
    model.parameters.probabilities[0].set_hospitalized_per_infectious(rho);
    model.parameters.probabilities[0].set_icu_per_hospitalized(theta);
    model.parameters.probabilities[0].set_dead_per_icu(delta);

    epi::set_log_level(epi::LogLevel::off);
    model.parameters.check_constraints();

    EXPECT_EQ(-91, model.populations.get((epi::AgeGroup1)0, epi::InfectionType::E));
    EXPECT_EQ(2.124921, model.parameters.probabilities[0].get_asymp_per_infectious().value());
    EXPECT_NEAR(5.08993, model.parameters.times[0].get_serialinterval(), 1e-14);

    model.apply_constraints();

    EXPECT_EQ(0.0, model.populations.get((epi::AgeGroup1)0, epi::InfectionType::E));
    EXPECT_EQ(0.0, model.parameters.probabilities[0].get_asymp_per_infectious().value());
    EXPECT_NEAR(4.6, model.parameters.times[0].get_serialinterval(), 1e-14);
}

TEST(TestSecir, testModelConstraints)
{
    double t0   = 0;
    double tmax = 60; // after 60 days with cont_freq 10 and winter, the virus would already decline
    double dt   = 0.1;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           ticu2death = 5;

    double cont_freq = 10, inf_prob = 0.05, carr_infec = 1, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2,
           theta = 0.25;

    double nb_total_t0 = 1000000, nb_exp_t0 = 10000, nb_inf_t0 = 5000, nb_car_t0 = 500, nb_hosp_t0 = 20, nb_icu_t0 = 0,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    epi::SecirModel1 model;

    model.parameters.times[0].set_incubation(tinc);
    model.parameters.times[0].set_infectious_mild(tinfmild);
    model.parameters.times[0].set_serialinterval(tserint);
    model.parameters.times[0].set_hospitalized_to_home(thosp2home);
    model.parameters.times[0].set_home_to_hospitalized(thome2hosp);
    model.parameters.times[0].set_hospitalized_to_icu(thosp2icu);
    model.parameters.times[0].set_icu_to_home(ticu2home);
    model.parameters.times[0].set_icu_to_death(ticu2death);

    model.populations.set(nb_exp_t0, (epi::AgeGroup1)0, epi::InfectionType::E);
    model.populations.set(nb_car_t0, (epi::AgeGroup1)0, epi::InfectionType::C);
    model.populations.set(nb_inf_t0, (epi::AgeGroup1)0, epi::InfectionType::I);
    model.populations.set(nb_hosp_t0, (epi::AgeGroup1)0, epi::InfectionType::H);
    model.populations.set(nb_icu_t0, (epi::AgeGroup1)0, epi::InfectionType::U);
    model.populations.set(nb_rec_t0, (epi::AgeGroup1)0, epi::InfectionType::R);
    model.populations.set(nb_dead_t0, (epi::AgeGroup1)0, epi::InfectionType::D);
    model.populations.set_difference_from_total(nb_total_t0, (epi::AgeGroup1)0, epi::InfectionType::S);

    model.parameters.probabilities[0].set_infection_from_contact(inf_prob);
    model.parameters.probabilities[0].set_carrier_infectability(carr_infec);
    model.parameters.probabilities[0].set_asymp_per_infectious(alpha);
    model.parameters.probabilities[0].set_risk_from_symptomatic(beta);
    model.parameters.probabilities[0].set_hospitalized_per_infectious(rho);
    model.parameters.probabilities[0].set_icu_per_hospitalized(theta);
    model.parameters.probabilities[0].set_dead_per_icu(delta);

    epi::ContactFrequencyMatrix& cont_freq_matrix = model.parameters.get_contact_patterns();
    cont_freq_matrix.set_cont_freq(cont_freq, 0, 0);

    model.apply_constraints();

    epi::TimeSeries<double> secihurd = simulate(t0, tmax, dt, model);
    double max_icu_cap               = 0;
    for (Eigen::Index i = 0; i < secihurd.get_num_time_points(); i++) {
        if (secihurd.get_value(i)[5] > max_icu_cap) {
            max_icu_cap = secihurd.get_value(i)[5];
        }
    }

    epi::TimeSeries<double> secihurd_interp = epi::interpolate_simulation_result(secihurd);

    model.parameters.set_start_day(100);
    model.parameters.set_seasonality(0.5);

    epi::TimeSeries<double> secihurd_season        = simulate(t0, tmax, dt, model);
    epi::TimeSeries<double> secihurd_season_interp = epi::interpolate_simulation_result(secihurd_season);

    for (Eigen::Index i = 0; i < secihurd_interp.get_num_time_points(); i++) {
        EXPECT_LE(secihurd_season_interp.get_value(i)[3], secihurd_interp.get_value(i)[3]) << " at row " << i;
    }

    model.parameters.set_start_day(280);

    epi::TimeSeries<double> secihurd_season2        = simulate(t0, tmax, dt, model);
    epi::TimeSeries<double> secihurd_season2_interp = epi::interpolate_simulation_result(secihurd_season2);

    for (Eigen::Index i = 0; i < secihurd_interp.get_num_time_points(); i++) {
        EXPECT_GE(secihurd_season2_interp.get_value(i)[3], secihurd_interp.get_value(i)[3]) << " at row " << i;
    }

    // params.set_icu_capacity(max_icu_cap - 3);

    // secihurd = simulate(t0, tmax, dt, params);
    // for (Eigen::Index i = 0; i < secihurd.get_num_time_points(); i++) {
    //     EXPECT_LE(secihurd.get_value(i)[5], max_icu_cap - 2.5) << " at row " << i;
    // }

    // temporary test for random variables
    set_params_distributions_normal(model, t0, tmax, 0.2);

    for (size_t j = 0; j < 10; j++) {
        draw_sample(model);
        model.parameters.set_icu_capacity(8000);
        secihurd = simulate(t0, tmax, dt, model);
        // max_icu_cap = 0;
        // for (Eigen::Index i = 0; i < secihurd.get_num_time_points(); i++) {
        //     if (secihurd.get_value(i)[5] > max_icu_cap) {
        //         max_icu_cap = secihurd.get_value(i)[5];
        //     }
        // }
        // printf("\n max cap: %.4f ", max_icu_cap);
        for (Eigen::Index i = 0; i < secihurd.get_num_time_points(); i++) {
            EXPECT_LE(secihurd.get_value(i)[5], 9000) << " at row " << i;
        }
    }
}
