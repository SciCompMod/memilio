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

    std::vector<epi::SecirParams> params{epi::SecirParams{}};
    epi::ContactFrequencyMatrix contact_freq_matrix{};

    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs

    params[0].times.set_incubation(tinc);
    params[0].times.set_infectious_mild(tinfmild);
    params[0].times.set_serialinterval(tserint);
    params[0].times.set_hospitalized_to_home(thosp2home);
    params[0].times.set_home_to_hospitalized(thome2hosp);
    params[0].times.set_hospitalized_to_icu(thosp2icu);
    params[0].times.set_icu_to_home(ticu2home);
    params[0].times.set_infectious_asymp(tinfasy);
    params[0].times.set_icu_to_death(ticu2death);

    contact_freq_matrix.set_cont_freq(cont_freq, 0, 0);
    epi::Damping dummy(30., 0.3);
    contact_freq_matrix.add_damping(dummy, 0, 0);

    params[0].populations.set_total_t0(nb_total_t0);
    params[0].populations.set_exposed_t0(nb_exp_t0);
    params[0].populations.set_carrier_t0(nb_car_t0);
    params[0].populations.set_infectious_t0(nb_inf_t0);
    params[0].populations.set_hospital_t0(nb_hosp_t0);
    params[0].populations.set_icu_t0(nb_icu_t0);
    params[0].populations.set_recovered_t0(nb_rec_t0);
    params[0].populations.set_dead_t0(nb_dead_t0);

    params[0].probabilities.set_asymp_per_infectious(alpha);
    params[0].probabilities.set_risk_from_symptomatic(beta);
    params[0].probabilities.set_hospitalized_per_infectious(rho);
    params[0].probabilities.set_icu_per_hospitalized(theta);
    params[0].probabilities.set_dead_per_icu(delta);

    std::vector<Eigen::VectorXd> secihurd(0);
    auto t = simulate(t0, tmax, dt, contact_freq_matrix, params, secihurd);

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

    epi::SecirParams params{epi::SecirParams{}};

    params.times.set_incubation(tinc);
    params.times.set_infectious_mild(tinfmild);
    params.times.set_serialinterval(tserint);
    params.times.set_hospitalized_to_home(thosp2home);
    params.times.set_home_to_hospitalized(thome2hosp);
    params.times.set_hospitalized_to_icu(thosp2icu);
    params.times.set_icu_to_home(ticu2home);
    params.times.set_infectious_asymp(tinfasy);
    params.times.set_icu_to_death(ticu2death);

    params.populations.set_total_t0(nb_total_t0);
    params.populations.set_exposed_t0(nb_exp_t0);
    params.populations.set_carrier_t0(nb_car_t0);
    params.populations.set_infectious_t0(nb_inf_t0);
    params.populations.set_hospital_t0(nb_hosp_t0);
    params.populations.set_icu_t0(nb_icu_t0);
    params.populations.set_recovered_t0(nb_rec_t0);
    params.populations.set_dead_t0(nb_dead_t0);

    params.probabilities.set_infection_from_contact(inf_cont);
    params.probabilities.set_asymp_per_infectious(alpha);
    params.probabilities.set_risk_from_symptomatic(beta);
    params.probabilities.set_hospitalized_per_infectious(rho);
    params.probabilities.set_icu_per_hospitalized(theta);
    params.probabilities.set_dead_per_icu(delta);

    epi::SecirParams params2{params}; // copy constructor

    EXPECT_EQ(params.populations.get_total_t0(), params2.populations.get_total_t0());
    EXPECT_EQ(params.populations.get_suscetible_t0(), params2.populations.get_suscetible_t0());
    EXPECT_EQ(params.populations.get_exposed_t0(), params2.populations.get_exposed_t0());
    EXPECT_EQ(params.populations.get_carrier_t0(), params2.populations.get_carrier_t0());
    EXPECT_EQ(params.populations.get_infectious_t0(), params2.populations.get_infectious_t0());
    EXPECT_EQ(params.populations.get_hospitalized_t0(), params2.populations.get_hospitalized_t0());
    EXPECT_EQ(params.populations.get_icu_t0(), params2.populations.get_icu_t0());
    EXPECT_EQ(params.populations.get_recovered_t0(), params2.populations.get_recovered_t0());
    EXPECT_EQ(params.populations.get_dead_t0(), params2.populations.get_dead_t0());

    EXPECT_EQ(params.times.get_incubation_inv(), params2.times.get_incubation_inv());
    EXPECT_EQ(params.times.get_serialinterval_inv(), params2.times.get_serialinterval_inv());
    EXPECT_EQ(params.times.get_infectious_mild_inv(), params2.times.get_infectious_mild_inv());
    EXPECT_EQ(params.times.get_infectious_asymp_inv(), params2.times.get_infectious_asymp_inv());
    EXPECT_EQ(params.times.get_home_to_hospitalized_inv(), params2.times.get_home_to_hospitalized_inv());
    EXPECT_EQ(params.times.get_hospitalized_to_home_inv(), params2.times.get_hospitalized_to_home_inv());
    EXPECT_EQ(params.times.get_hospitalized_to_icu_inv(), params2.times.get_hospitalized_to_icu_inv());
    EXPECT_EQ(params.times.get_icu_to_dead_inv(), params2.times.get_icu_to_dead_inv());
    EXPECT_EQ(params.times.get_icu_to_home_inv(), params2.times.get_icu_to_home_inv());

    EXPECT_EQ(params.probabilities.get_infection_from_contact(), params2.probabilities.get_infection_from_contact());
    EXPECT_EQ(params.probabilities.get_risk_from_symptomatic(), params2.probabilities.get_risk_from_symptomatic());
    EXPECT_EQ(params.probabilities.get_asymp_per_infectious(), params2.probabilities.get_asymp_per_infectious());
    EXPECT_EQ(params.probabilities.get_hospitalized_per_infectious(),
              params2.probabilities.get_hospitalized_per_infectious());
    EXPECT_EQ(params.probabilities.get_icu_per_hospitalized(), params2.probabilities.get_icu_per_hospitalized());
    EXPECT_EQ(params.probabilities.get_dead_per_icu(), params2.probabilities.get_dead_per_icu());

    epi::SecirParams params3{std::move(params2)}; // move constructor

    EXPECT_EQ(params3.populations.get_total_t0(), params2.populations.get_total_t0());
    EXPECT_EQ(params3.populations.get_suscetible_t0(), params2.populations.get_suscetible_t0());
    EXPECT_EQ(params3.populations.get_exposed_t0(), params2.populations.get_exposed_t0());
    EXPECT_EQ(params3.populations.get_carrier_t0(), params2.populations.get_carrier_t0());
    EXPECT_EQ(params3.populations.get_infectious_t0(), params2.populations.get_infectious_t0());
    EXPECT_EQ(params3.populations.get_hospitalized_t0(), params2.populations.get_hospitalized_t0());
    EXPECT_EQ(params3.populations.get_icu_t0(), params2.populations.get_icu_t0());
    EXPECT_EQ(params3.populations.get_recovered_t0(), params2.populations.get_recovered_t0());
    EXPECT_EQ(params3.populations.get_dead_t0(), params2.populations.get_dead_t0());

    EXPECT_EQ(params3.times.get_incubation_inv(), params2.times.get_incubation_inv());
    EXPECT_EQ(params3.times.get_serialinterval_inv(), params2.times.get_serialinterval_inv());
    EXPECT_EQ(params3.times.get_infectious_mild_inv(), params2.times.get_infectious_mild_inv());
    EXPECT_EQ(params3.times.get_infectious_asymp_inv(), params2.times.get_infectious_asymp_inv());
    EXPECT_EQ(params3.times.get_home_to_hospitalized_inv(), params2.times.get_home_to_hospitalized_inv());
    EXPECT_EQ(params3.times.get_hospitalized_to_home_inv(), params2.times.get_hospitalized_to_home_inv());
    EXPECT_EQ(params3.times.get_hospitalized_to_icu_inv(), params2.times.get_hospitalized_to_icu_inv());
    EXPECT_EQ(params3.times.get_icu_to_dead_inv(), params2.times.get_icu_to_dead_inv());
    EXPECT_EQ(params3.times.get_icu_to_home_inv(), params2.times.get_icu_to_home_inv());

    EXPECT_EQ(params3.probabilities.get_infection_from_contact(), params2.probabilities.get_infection_from_contact());
    EXPECT_EQ(params3.probabilities.get_risk_from_symptomatic(), params2.probabilities.get_risk_from_symptomatic());
    EXPECT_EQ(params3.probabilities.get_asymp_per_infectious(), params2.probabilities.get_asymp_per_infectious());
    EXPECT_EQ(params3.probabilities.get_hospitalized_per_infectious(),
              params2.probabilities.get_hospitalized_per_infectious());
    EXPECT_EQ(params3.probabilities.get_icu_per_hospitalized(), params2.probabilities.get_icu_per_hospitalized());
    EXPECT_EQ(params3.probabilities.get_dead_per_icu(), params2.probabilities.get_dead_per_icu());

    epi::SecirParams params4 = params3; // copy assignment constructor

    EXPECT_EQ(params3.populations.get_total_t0(), params4.populations.get_total_t0());
    EXPECT_EQ(params3.populations.get_suscetible_t0(), params4.populations.get_suscetible_t0());
    EXPECT_EQ(params3.populations.get_exposed_t0(), params4.populations.get_exposed_t0());
    EXPECT_EQ(params3.populations.get_carrier_t0(), params4.populations.get_carrier_t0());
    EXPECT_EQ(params3.populations.get_infectious_t0(), params4.populations.get_infectious_t0());
    EXPECT_EQ(params3.populations.get_hospitalized_t0(), params4.populations.get_hospitalized_t0());
    EXPECT_EQ(params3.populations.get_icu_t0(), params4.populations.get_icu_t0());
    EXPECT_EQ(params3.populations.get_recovered_t0(), params4.populations.get_recovered_t0());
    EXPECT_EQ(params3.populations.get_dead_t0(), params4.populations.get_dead_t0());

    EXPECT_EQ(params3.times.get_incubation_inv(), params4.times.get_incubation_inv());
    EXPECT_EQ(params3.times.get_serialinterval_inv(), params4.times.get_serialinterval_inv());
    EXPECT_EQ(params3.times.get_infectious_mild_inv(), params4.times.get_infectious_mild_inv());
    EXPECT_EQ(params3.times.get_infectious_asymp_inv(), params4.times.get_infectious_asymp_inv());
    EXPECT_EQ(params3.times.get_home_to_hospitalized_inv(), params4.times.get_home_to_hospitalized_inv());
    EXPECT_EQ(params3.times.get_hospitalized_to_home_inv(), params4.times.get_hospitalized_to_home_inv());
    EXPECT_EQ(params3.times.get_hospitalized_to_icu_inv(), params4.times.get_hospitalized_to_icu_inv());
    EXPECT_EQ(params3.times.get_icu_to_dead_inv(), params4.times.get_icu_to_dead_inv());
    EXPECT_EQ(params3.times.get_icu_to_home_inv(), params4.times.get_icu_to_home_inv());

    EXPECT_EQ(params3.probabilities.get_infection_from_contact(), params4.probabilities.get_infection_from_contact());
    EXPECT_EQ(params3.probabilities.get_risk_from_symptomatic(), params4.probabilities.get_risk_from_symptomatic());
    EXPECT_EQ(params3.probabilities.get_asymp_per_infectious(), params4.probabilities.get_asymp_per_infectious());
    EXPECT_EQ(params3.probabilities.get_hospitalized_per_infectious(),
              params4.probabilities.get_hospitalized_per_infectious());
    EXPECT_EQ(params3.probabilities.get_icu_per_hospitalized(), params4.probabilities.get_icu_per_hospitalized());
    EXPECT_EQ(params3.probabilities.get_dead_per_icu(), params4.probabilities.get_dead_per_icu());

    epi::SecirParams params5 = std::move(params4); // move assignment constructor

    EXPECT_EQ(params5.populations.get_total_t0(), params4.populations.get_total_t0());
    EXPECT_EQ(params5.populations.get_suscetible_t0(), params4.populations.get_suscetible_t0());
    EXPECT_EQ(params5.populations.get_exposed_t0(), params4.populations.get_exposed_t0());
    EXPECT_EQ(params5.populations.get_carrier_t0(), params4.populations.get_carrier_t0());
    EXPECT_EQ(params5.populations.get_infectious_t0(), params4.populations.get_infectious_t0());
    EXPECT_EQ(params5.populations.get_hospitalized_t0(), params4.populations.get_hospitalized_t0());
    EXPECT_EQ(params5.populations.get_icu_t0(), params4.populations.get_icu_t0());
    EXPECT_EQ(params5.populations.get_recovered_t0(), params4.populations.get_recovered_t0());
    EXPECT_EQ(params5.populations.get_dead_t0(), params4.populations.get_dead_t0());

    EXPECT_EQ(params5.times.get_incubation_inv(), params4.times.get_incubation_inv());
    EXPECT_EQ(params5.times.get_serialinterval_inv(), params4.times.get_serialinterval_inv());
    EXPECT_EQ(params5.times.get_infectious_mild_inv(), params4.times.get_infectious_mild_inv());
    EXPECT_EQ(params5.times.get_infectious_asymp_inv(), params4.times.get_infectious_asymp_inv());
    EXPECT_EQ(params5.times.get_home_to_hospitalized_inv(), params4.times.get_home_to_hospitalized_inv());
    EXPECT_EQ(params5.times.get_hospitalized_to_home_inv(), params4.times.get_hospitalized_to_home_inv());
    EXPECT_EQ(params5.times.get_hospitalized_to_icu_inv(), params4.times.get_hospitalized_to_icu_inv());
    EXPECT_EQ(params5.times.get_icu_to_dead_inv(), params4.times.get_icu_to_dead_inv());
    EXPECT_EQ(params5.times.get_icu_to_home_inv(), params4.times.get_icu_to_home_inv());

    EXPECT_EQ(params5.probabilities.get_infection_from_contact(), params4.probabilities.get_infection_from_contact());
    EXPECT_EQ(params5.probabilities.get_risk_from_symptomatic(), params4.probabilities.get_risk_from_symptomatic());
    EXPECT_EQ(params5.probabilities.get_asymp_per_infectious(), params4.probabilities.get_asymp_per_infectious());
    EXPECT_EQ(params5.probabilities.get_hospitalized_per_infectious(),
              params4.probabilities.get_hospitalized_per_infectious());
    EXPECT_EQ(params5.probabilities.get_icu_per_hospitalized(), params4.probabilities.get_icu_per_hospitalized());
    EXPECT_EQ(params5.probabilities.get_dead_per_icu(), params4.probabilities.get_dead_per_icu());

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