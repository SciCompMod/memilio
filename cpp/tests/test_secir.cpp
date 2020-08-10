#include "load_test_data.h"
#include "epidemiology/secir.h"
#include <gtest/gtest.h>

TEST(TestSecir, compareWithPreviousRun)
{
    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           tinfasy = 6.2, ticu2death = 5;

    double cont_freq = 0.5, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    epi::SecirParams params;

    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs

    params.times[0].set_incubation(tinc);
    params.times[0].set_infectious_mild(tinfmild);
    params.times[0].set_serialinterval(tserint);
    params.times[0].set_hospitalized_to_home(thosp2home);
    params.times[0].set_home_to_hospitalized(thome2hosp);
    params.times[0].set_hospitalized_to_icu(thosp2icu);
    params.times[0].set_icu_to_home(ticu2home);
    params.times[0].set_infectious_asymp(tinfasy);
    params.times[0].set_icu_to_death(ticu2death);

    params.get_cont_freq_matrix().set_cont_freq(cont_freq, 0, 0);
    epi::Damping dummy(30., 0.3);
    params.get_cont_freq_matrix().add_damping(dummy, 0, 0);

    params.populations.set({0, epi::SecirCompartments::E}, nb_exp_t0);
    params.populations.set({0, epi::SecirCompartments::C}, nb_car_t0);
    params.populations.set({0, epi::SecirCompartments::I}, nb_inf_t0);
    params.populations.set({0, epi::SecirCompartments::H}, nb_hosp_t0);
    params.populations.set({0, epi::SecirCompartments::U}, nb_icu_t0);
    params.populations.set({0, epi::SecirCompartments::R}, nb_rec_t0);
    params.populations.set({0, epi::SecirCompartments::D}, nb_dead_t0);
    params.populations.set_difference_from_total({0, epi::SecirCompartments::S}, nb_total_t0);

    params.probabilities[0].set_asymp_per_infectious(alpha);
    params.probabilities[0].set_risk_from_symptomatic(beta);
    params.probabilities[0].set_hospitalized_per_infectious(rho);
    params.probabilities[0].set_icu_per_hospitalized(theta);
    params.probabilities[0].set_dead_per_icu(delta);

    std::vector<Eigen::VectorXd> secihurd(0);
    auto t = simulate(t0, tmax, dt, params, secihurd);

    auto compare = load_test_data_csv<double>("secihurd-compare.csv");

    ASSERT_EQ(compare.size(), t.size());
    ASSERT_EQ(compare.size(), secihurd.size());
    for (size_t i = 0; i < compare.size(); i++) {
        ASSERT_EQ(compare[i].size(), secihurd[i].size() + 1) << "at row " << i;
        EXPECT_NEAR(t[i], compare[i][0], 1e-10) << "at row " << i;
        for (size_t j = 1; j < compare[i].size(); j++) {
            EXPECT_NEAR(secihurd[i][j - 1], compare[i][j], 1e-10) << " at row " << i;
        }
    }
}

TEST(TestSecir, testParamConstructors)
{
    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           tinfasy = 6.2, ticu2death = 5;

    double cont_freq = 0.5, inf_cont = 0.9, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.24;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 54, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 11, nb_dead_t0 = 0;

    epi::SecirParams params;

    params.times[0].set_incubation(tinc);
    params.times[0].set_infectious_mild(tinfmild);
    params.times[0].set_serialinterval(tserint);
    params.times[0].set_hospitalized_to_home(thosp2home);
    params.times[0].set_home_to_hospitalized(thome2hosp);
    params.times[0].set_hospitalized_to_icu(thosp2icu);
    params.times[0].set_icu_to_home(ticu2home);
    params.times[0].set_infectious_asymp(tinfasy);
    params.times[0].set_icu_to_death(ticu2death);

    params.populations.set({0, epi::SecirCompartments::E}, nb_exp_t0);
    params.populations.set({0, epi::SecirCompartments::C}, nb_car_t0);
    params.populations.set({0, epi::SecirCompartments::I}, nb_inf_t0);
    params.populations.set({0, epi::SecirCompartments::H}, nb_hosp_t0);
    params.populations.set({0, epi::SecirCompartments::U}, nb_icu_t0);
    params.populations.set({0, epi::SecirCompartments::R}, nb_rec_t0);
    params.populations.set({0, epi::SecirCompartments::D}, nb_dead_t0);
    params.populations.set_difference_from_total({0, epi::SecirCompartments::S}, nb_total_t0);

    params.probabilities[0].set_infection_from_contact(inf_cont);
    params.probabilities[0].set_asymp_per_infectious(alpha);
    params.probabilities[0].set_risk_from_symptomatic(beta);
    params.probabilities[0].set_hospitalized_per_infectious(rho);
    params.probabilities[0].set_icu_per_hospitalized(theta);
    params.probabilities[0].set_dead_per_icu(delta);

    epi::SecirParams params2{params}; // copy constructor

    EXPECT_EQ(params.populations.get_total(), params2.populations.get_total());
    EXPECT_EQ(params.populations.get({0, epi::SecirCompartments::S}),
              params2.populations.get({0, epi::SecirCompartments::S}));
    EXPECT_EQ(params.populations.get({0, epi::SecirCompartments::E}),
              params2.populations.get({0, epi::SecirCompartments::E}));
    EXPECT_EQ(params.populations.get({0, epi::SecirCompartments::C}),
              params2.populations.get({0, epi::SecirCompartments::C}));
    EXPECT_EQ(params.populations.get({0, epi::SecirCompartments::I}),
              params2.populations.get({0, epi::SecirCompartments::I}));
    EXPECT_EQ(params.populations.get({0, epi::SecirCompartments::H}),
              params2.populations.get({0, epi::SecirCompartments::H}));
    EXPECT_EQ(params.populations.get({0, epi::SecirCompartments::U}),
              params2.populations.get({0, epi::SecirCompartments::U}));
    EXPECT_EQ(params.populations.get({0, epi::SecirCompartments::R}),
              params2.populations.get({0, epi::SecirCompartments::R}));
    EXPECT_EQ(params.populations.get({0, epi::SecirCompartments::D}),
              params2.populations.get({0, epi::SecirCompartments::D}));

    EXPECT_EQ(params.times[0].get_incubation(), params2.times[0].get_incubation());
    EXPECT_EQ(params.times[0].get_serialinterval(), params2.times[0].get_serialinterval());
    EXPECT_EQ(params.times[0].get_infectious_mild(), params2.times[0].get_infectious_mild());
    EXPECT_EQ(params.times[0].get_infectious_asymp(), params2.times[0].get_infectious_asymp());
    EXPECT_EQ(params.times[0].get_home_to_hospitalized(), params2.times[0].get_home_to_hospitalized());
    EXPECT_EQ(params.times[0].get_hospitalized_to_home(), params2.times[0].get_hospitalized_to_home());
    EXPECT_EQ(params.times[0].get_hospitalized_to_icu(), params2.times[0].get_hospitalized_to_icu());
    EXPECT_EQ(params.times[0].get_icu_to_dead(), params2.times[0].get_icu_to_dead());
    EXPECT_EQ(params.times[0].get_icu_to_home(), params2.times[0].get_icu_to_home());

    EXPECT_EQ(params.probabilities[0].get_infection_from_contact(),
              params2.probabilities[0].get_infection_from_contact());
    EXPECT_EQ(params.probabilities[0].get_risk_from_symptomatic(),
              params2.probabilities[0].get_risk_from_symptomatic());
    EXPECT_EQ(params.probabilities[0].get_asymp_per_infectious(), params2.probabilities[0].get_asymp_per_infectious());
    EXPECT_EQ(params.probabilities[0].get_hospitalized_per_infectious(),
              params2.probabilities[0].get_hospitalized_per_infectious());
    EXPECT_EQ(params.probabilities[0].get_icu_per_hospitalized(), params2.probabilities[0].get_icu_per_hospitalized());
    EXPECT_EQ(params.probabilities[0].get_dead_per_icu(), params2.probabilities[0].get_dead_per_icu());

    epi::SecirParams params3 = std::move(params2); // move constructor

    EXPECT_EQ(params3.populations.get_total(), params.populations.get_total());
    EXPECT_EQ(params3.populations.get({0, epi::SecirCompartments::S}),
              params.populations.get({0, epi::SecirCompartments::S}));
    EXPECT_EQ(params3.populations.get({0, epi::SecirCompartments::E}),
              params.populations.get({0, epi::SecirCompartments::E}));
    EXPECT_EQ(params3.populations.get({0, epi::SecirCompartments::C}),
              params.populations.get({0, epi::SecirCompartments::C}));
    EXPECT_EQ(params3.populations.get({0, epi::SecirCompartments::I}),
              params.populations.get({0, epi::SecirCompartments::I}));
    EXPECT_EQ(params3.populations.get({0, epi::SecirCompartments::H}),
              params.populations.get({0, epi::SecirCompartments::H}));
    EXPECT_EQ(params3.populations.get({0, epi::SecirCompartments::U}),
              params.populations.get({0, epi::SecirCompartments::U}));
    EXPECT_EQ(params3.populations.get({0, epi::SecirCompartments::R}),
              params.populations.get({0, epi::SecirCompartments::R}));
    EXPECT_EQ(params3.populations.get({0, epi::SecirCompartments::D}),
              params.populations.get({0, epi::SecirCompartments::D}));

    EXPECT_EQ(params3.times[0].get_incubation(), params.times[0].get_incubation());
    EXPECT_EQ(params3.times[0].get_serialinterval(), params.times[0].get_serialinterval());
    EXPECT_EQ(params3.times[0].get_infectious_mild(), params.times[0].get_infectious_mild());
    EXPECT_EQ(params3.times[0].get_infectious_asymp(), params.times[0].get_infectious_asymp());
    EXPECT_EQ(params3.times[0].get_home_to_hospitalized(), params.times[0].get_home_to_hospitalized());
    EXPECT_EQ(params3.times[0].get_hospitalized_to_home(), params.times[0].get_hospitalized_to_home());
    EXPECT_EQ(params3.times[0].get_hospitalized_to_icu(), params.times[0].get_hospitalized_to_icu());
    EXPECT_EQ(params3.times[0].get_icu_to_dead(), params.times[0].get_icu_to_dead());
    EXPECT_EQ(params3.times[0].get_icu_to_home(), params.times[0].get_icu_to_home());

    EXPECT_EQ(params3.probabilities[0].get_infection_from_contact(),
              params.probabilities[0].get_infection_from_contact());
    EXPECT_EQ(params3.probabilities[0].get_risk_from_symptomatic(),
              params.probabilities[0].get_risk_from_symptomatic());
    EXPECT_EQ(params3.probabilities[0].get_asymp_per_infectious(), params.probabilities[0].get_asymp_per_infectious());
    EXPECT_EQ(params3.probabilities[0].get_hospitalized_per_infectious(),
              params.probabilities[0].get_hospitalized_per_infectious());
    EXPECT_EQ(params3.probabilities[0].get_icu_per_hospitalized(), params.probabilities[0].get_icu_per_hospitalized());
    EXPECT_EQ(params3.probabilities[0].get_dead_per_icu(), params.probabilities[0].get_dead_per_icu());

    epi::SecirParams params4 = params3; // copy assignment constructor

    EXPECT_EQ(params3.populations.get_total(), params4.populations.get_total());
    EXPECT_EQ(params3.populations.get({0, epi::SecirCompartments::S}),
              params4.populations.get({0, epi::SecirCompartments::S}));
    EXPECT_EQ(params3.populations.get({0, epi::SecirCompartments::E}),
              params4.populations.get({0, epi::SecirCompartments::E}));
    EXPECT_EQ(params3.populations.get({0, epi::SecirCompartments::C}),
              params4.populations.get({0, epi::SecirCompartments::C}));
    EXPECT_EQ(params3.populations.get({0, epi::SecirCompartments::I}),
              params4.populations.get({0, epi::SecirCompartments::I}));
    EXPECT_EQ(params3.populations.get({0, epi::SecirCompartments::H}),
              params4.populations.get({0, epi::SecirCompartments::H}));
    EXPECT_EQ(params3.populations.get({0, epi::SecirCompartments::U}),
              params4.populations.get({0, epi::SecirCompartments::U}));
    EXPECT_EQ(params3.populations.get({0, epi::SecirCompartments::R}),
              params4.populations.get({0, epi::SecirCompartments::R}));
    EXPECT_EQ(params3.populations.get({0, epi::SecirCompartments::D}),
              params4.populations.get({0, epi::SecirCompartments::D}));

    EXPECT_EQ(params3.times[0].get_incubation(), params4.times[0].get_incubation());
    EXPECT_EQ(params3.times[0].get_serialinterval(), params4.times[0].get_serialinterval());
    EXPECT_EQ(params3.times[0].get_infectious_mild(), params4.times[0].get_infectious_mild());
    EXPECT_EQ(params3.times[0].get_infectious_asymp(), params4.times[0].get_infectious_asymp());
    EXPECT_EQ(params3.times[0].get_home_to_hospitalized(), params4.times[0].get_home_to_hospitalized());
    EXPECT_EQ(params3.times[0].get_hospitalized_to_home(), params4.times[0].get_hospitalized_to_home());
    EXPECT_EQ(params3.times[0].get_hospitalized_to_icu(), params4.times[0].get_hospitalized_to_icu());
    EXPECT_EQ(params3.times[0].get_icu_to_dead(), params4.times[0].get_icu_to_dead());
    EXPECT_EQ(params3.times[0].get_icu_to_home(), params4.times[0].get_icu_to_home());

    EXPECT_EQ(params3.probabilities[0].get_infection_from_contact(),
              params4.probabilities[0].get_infection_from_contact());
    EXPECT_EQ(params3.probabilities[0].get_risk_from_symptomatic(),
              params4.probabilities[0].get_risk_from_symptomatic());
    EXPECT_EQ(params3.probabilities[0].get_asymp_per_infectious(), params4.probabilities[0].get_asymp_per_infectious());
    EXPECT_EQ(params3.probabilities[0].get_hospitalized_per_infectious(),
              params4.probabilities[0].get_hospitalized_per_infectious());
    EXPECT_EQ(params3.probabilities[0].get_icu_per_hospitalized(), params4.probabilities[0].get_icu_per_hospitalized());
    EXPECT_EQ(params3.probabilities[0].get_dead_per_icu(), params4.probabilities[0].get_dead_per_icu());

    epi::SecirParams params5 = std::move(params4); // move assignment constructor

    EXPECT_EQ(params5.populations.get_total(), params3.populations.get_total());
    EXPECT_EQ(params5.populations.get({0, epi::SecirCompartments::S}),
              params3.populations.get({0, epi::SecirCompartments::S}));
    EXPECT_EQ(params5.populations.get({0, epi::SecirCompartments::E}),
              params3.populations.get({0, epi::SecirCompartments::E}));
    EXPECT_EQ(params5.populations.get({0, epi::SecirCompartments::C}),
              params3.populations.get({0, epi::SecirCompartments::C}));
    EXPECT_EQ(params5.populations.get({0, epi::SecirCompartments::I}),
              params3.populations.get({0, epi::SecirCompartments::I}));
    EXPECT_EQ(params5.populations.get({0, epi::SecirCompartments::H}),
              params3.populations.get({0, epi::SecirCompartments::H}));
    EXPECT_EQ(params5.populations.get({0, epi::SecirCompartments::U}),
              params3.populations.get({0, epi::SecirCompartments::U}));
    EXPECT_EQ(params5.populations.get({0, epi::SecirCompartments::R}),
              params3.populations.get({0, epi::SecirCompartments::R}));
    EXPECT_EQ(params5.populations.get({0, epi::SecirCompartments::D}),
              params3.populations.get({0, epi::SecirCompartments::D}));

    EXPECT_EQ(params5.times[0].get_incubation(), params3.times[0].get_incubation());
    EXPECT_EQ(params5.times[0].get_serialinterval(), params3.times[0].get_serialinterval());
    EXPECT_EQ(params5.times[0].get_infectious_mild(), params3.times[0].get_infectious_mild());
    EXPECT_EQ(params5.times[0].get_infectious_asymp(), params3.times[0].get_infectious_asymp());
    EXPECT_EQ(params5.times[0].get_home_to_hospitalized(), params3.times[0].get_home_to_hospitalized());
    EXPECT_EQ(params5.times[0].get_hospitalized_to_home(), params3.times[0].get_hospitalized_to_home());
    EXPECT_EQ(params5.times[0].get_hospitalized_to_icu(), params3.times[0].get_hospitalized_to_icu());
    EXPECT_EQ(params5.times[0].get_icu_to_dead(), params3.times[0].get_icu_to_dead());
    EXPECT_EQ(params5.times[0].get_icu_to_home(), params3.times[0].get_icu_to_home());

    EXPECT_EQ(params5.probabilities[0].get_infection_from_contact(),
              params3.probabilities[0].get_infection_from_contact());
    EXPECT_EQ(params5.probabilities[0].get_risk_from_symptomatic(),
              params3.probabilities[0].get_risk_from_symptomatic());
    EXPECT_EQ(params5.probabilities[0].get_asymp_per_infectious(), params3.probabilities[0].get_asymp_per_infectious());
    EXPECT_EQ(params5.probabilities[0].get_hospitalized_per_infectious(),
              params3.probabilities[0].get_hospitalized_per_infectious());
    EXPECT_EQ(params5.probabilities[0].get_icu_per_hospitalized(), params3.probabilities[0].get_icu_per_hospitalized());
    EXPECT_EQ(params5.probabilities[0].get_dead_per_icu(), params3.probabilities[0].get_dead_per_icu());

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
