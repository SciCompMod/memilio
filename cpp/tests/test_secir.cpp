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
    params.times[0].set_icu_to_death(ticu2death);

    epi::ContactFrequencyMatrix& cont_freq_matrix = params.get_contact_patterns();
    cont_freq_matrix.set_cont_freq(cont_freq, 0, 0);
    epi::Damping dummy(30., 0.3);
    cont_freq_matrix.add_damping(dummy, 0, 0);

    params.populations.set({0, epi::SecirCompartments::E}, nb_exp_t0);
    params.populations.set({0, epi::SecirCompartments::C}, nb_car_t0);
    params.populations.set({0, epi::SecirCompartments::I}, nb_inf_t0);
    params.populations.set({0, epi::SecirCompartments::H}, nb_hosp_t0);
    params.populations.set({0, epi::SecirCompartments::U}, nb_icu_t0);
    params.populations.set({0, epi::SecirCompartments::R}, nb_rec_t0);
    params.populations.set({0, epi::SecirCompartments::D}, nb_dead_t0);
    params.populations.set_difference_from_total({0, epi::SecirCompartments::S}, nb_total_t0);

    params.probabilities[0].set_infection_from_contact(inf_prob);
    params.probabilities[0].set_carrier_infectability(carr_infec);
    params.probabilities[0].set_asymp_per_infectious(alpha);
    params.probabilities[0].set_risk_from_symptomatic(beta);
    params.probabilities[0].set_hospitalized_per_infectious(rho);
    params.probabilities[0].set_icu_per_hospitalized(theta);
    params.probabilities[0].set_dead_per_icu(delta);

    params.apply_constraints();

    epi::TimeSeries<double> secihurd = simulate(t0, tmax, dt, params);

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

    double start_day = 30, seasonality = 0.3;

    epi::SecirParams params;

    params.set_start_day(start_day);
    params.set_seasonality(seasonality);

    params.times[0].set_incubation(tinc);
    params.times[0].set_infectious_mild(tinfmild);
    params.times[0].set_serialinterval(tserint);
    params.times[0].set_hospitalized_to_home(thosp2home);
    params.times[0].set_home_to_hospitalized(thome2hosp);
    params.times[0].set_hospitalized_to_icu(thosp2icu);
    params.times[0].set_icu_to_home(ticu2home);
    params.times[0].set_icu_to_death(ticu2death);

    params.populations.set({0, epi::SecirCompartments::E}, nb_exp_t0);
    params.populations.set({0, epi::SecirCompartments::C}, nb_car_t0);
    params.populations.set({0, epi::SecirCompartments::I}, nb_inf_t0);
    params.populations.set({0, epi::SecirCompartments::H}, nb_hosp_t0);
    params.populations.set({0, epi::SecirCompartments::U}, nb_icu_t0);
    params.populations.set({0, epi::SecirCompartments::R}, nb_rec_t0);
    params.populations.set({0, epi::SecirCompartments::D}, nb_dead_t0);
    params.populations.set_difference_from_total({0, epi::SecirCompartments::S}, nb_total_t0);

    params.probabilities[0].set_infection_from_contact(inf_prob);
    params.probabilities[0].set_carrier_infectability(carr_infec);
    params.probabilities[0].set_asymp_per_infectious(alpha);
    params.probabilities[0].set_risk_from_symptomatic(beta);
    params.probabilities[0].set_hospitalized_per_infectious(rho);
    params.probabilities[0].set_icu_per_hospitalized(theta);
    params.probabilities[0].set_dead_per_icu(delta);

    epi::ContactFrequencyMatrix& cont_freq_matrix = params.get_contact_patterns();
    cont_freq_matrix.set_cont_freq(cont_freq, 0, 0);
    epi::Damping dummy(30., 0.3);
    cont_freq_matrix.add_damping(dummy, 0, 0);

    epi::SecirParams params2{params}; // copy constructor

    EXPECT_EQ(params.get_start_day(), params2.get_start_day());
    EXPECT_EQ(params.get_seasonality(), params2.get_seasonality());

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
    EXPECT_EQ(params.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35),
              params2.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35));

    EXPECT_EQ(params.probabilities[0].get_infection_from_contact(),
              params2.probabilities[0].get_infection_from_contact());
    EXPECT_EQ(params.probabilities[0].get_carrier_infectability(),
              params2.probabilities[0].get_carrier_infectability());
    EXPECT_EQ(params.probabilities[0].get_risk_from_symptomatic(),
              params2.probabilities[0].get_risk_from_symptomatic());
    EXPECT_EQ(params.probabilities[0].get_asymp_per_infectious(), params2.probabilities[0].get_asymp_per_infectious());
    EXPECT_EQ(params.probabilities[0].get_hospitalized_per_infectious(),
              params2.probabilities[0].get_hospitalized_per_infectious());
    EXPECT_EQ(params.probabilities[0].get_icu_per_hospitalized(), params2.probabilities[0].get_icu_per_hospitalized());
    EXPECT_EQ(params.probabilities[0].get_dead_per_icu(), params2.probabilities[0].get_dead_per_icu());
    EXPECT_EQ(params.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35),
              params2.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35));

    epi::SecirParams params3 = std::move(params2); // move constructor

    EXPECT_EQ(params.get_start_day(), params3.get_start_day());
    EXPECT_EQ(params.get_seasonality(), params3.get_seasonality());

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
    EXPECT_EQ(params3.probabilities[0].get_carrier_infectability(),
              params.probabilities[0].get_carrier_infectability());
    EXPECT_EQ(params3.probabilities[0].get_risk_from_symptomatic(),
              params.probabilities[0].get_risk_from_symptomatic());
    EXPECT_EQ(params3.probabilities[0].get_asymp_per_infectious(), params.probabilities[0].get_asymp_per_infectious());
    EXPECT_EQ(params3.probabilities[0].get_hospitalized_per_infectious(),
              params.probabilities[0].get_hospitalized_per_infectious());
    EXPECT_EQ(params3.probabilities[0].get_icu_per_hospitalized(), params.probabilities[0].get_icu_per_hospitalized());
    EXPECT_EQ(params3.probabilities[0].get_dead_per_icu(), params.probabilities[0].get_dead_per_icu());

    EXPECT_EQ(params.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35),
              params3.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35));

    epi::SecirParams params4 = params3; // copy assignment constructor

    EXPECT_EQ(params4.get_start_day(), params3.get_start_day());
    EXPECT_EQ(params4.get_seasonality(), params3.get_seasonality());

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
    EXPECT_EQ(params3.probabilities[0].get_carrier_infectability(),
              params4.probabilities[0].get_carrier_infectability());
    EXPECT_EQ(params3.probabilities[0].get_risk_from_symptomatic(),
              params4.probabilities[0].get_risk_from_symptomatic());
    EXPECT_EQ(params3.probabilities[0].get_asymp_per_infectious(), params4.probabilities[0].get_asymp_per_infectious());
    EXPECT_EQ(params3.probabilities[0].get_hospitalized_per_infectious(),
              params4.probabilities[0].get_hospitalized_per_infectious());
    EXPECT_EQ(params3.probabilities[0].get_icu_per_hospitalized(), params4.probabilities[0].get_icu_per_hospitalized());
    EXPECT_EQ(params3.probabilities[0].get_dead_per_icu(), params4.probabilities[0].get_dead_per_icu());

    EXPECT_EQ(params4.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35),
              params3.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35));

    epi::SecirParams params5 = std::move(params4); // move assignment constructor

    EXPECT_EQ(params5.get_start_day(), params3.get_start_day());
    EXPECT_EQ(params5.get_seasonality(), params3.get_seasonality());

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
    EXPECT_EQ(params5.probabilities[0].get_carrier_infectability(),
              params3.probabilities[0].get_carrier_infectability());
    EXPECT_EQ(params5.probabilities[0].get_risk_from_symptomatic(),
              params3.probabilities[0].get_risk_from_symptomatic());
    EXPECT_EQ(params5.probabilities[0].get_asymp_per_infectious(), params3.probabilities[0].get_asymp_per_infectious());
    EXPECT_EQ(params5.probabilities[0].get_hospitalized_per_infectious(),
              params3.probabilities[0].get_hospitalized_per_infectious());
    EXPECT_EQ(params5.probabilities[0].get_icu_per_hospitalized(), params3.probabilities[0].get_icu_per_hospitalized());
    EXPECT_EQ(params5.probabilities[0].get_dead_per_icu(), params3.probabilities[0].get_dead_per_icu());

    EXPECT_EQ(params5.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35),
              params3.get_contact_patterns().get_cont_freq_mat().get_dampings(0, 0).get_factor(35));

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

    epi::SecirParams params;

    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs

    EXPECT_EQ(params.times[0].get_incubation().get_distribution().get(), nullptr);

    params.times[0].set_incubation(vec[1]);
    params.times[0].set_infectious_mild(vec[2]);
    params.times[0].set_serialinterval(vec[3]);
    params.times[0].set_hospitalized_to_home(vec[4]);
    params.times[0].set_home_to_hospitalized(vec[5]);
    params.times[0].set_hospitalized_to_icu(vec[6]);
    params.times[0].set_icu_to_home(vec[7]);
    params.times[0].set_infectious_asymp(vec[8]);
    params.times[0].set_icu_to_death(vec[9]);

    params.populations.set({0, epi::SecirCompartments::E}, vec[10]);
    params.populations.set({0, epi::SecirCompartments::C}, vec[11]);
    params.populations.set({0, epi::SecirCompartments::I}, vec[12]);
    params.populations.set({0, epi::SecirCompartments::H}, vec[13]);
    params.populations.set({0, epi::SecirCompartments::U}, vec[14]);
    params.populations.set({0, epi::SecirCompartments::R}, vec[15]);
    params.populations.set({0, epi::SecirCompartments::D}, vec[16]);

    params.probabilities[0].set_infection_from_contact(vec[17]);
    params.probabilities[0].set_carrier_infectability(vec[18]);
    params.probabilities[0].set_asymp_per_infectious(vec[19]);
    params.probabilities[0].set_risk_from_symptomatic(vec[20]);
    params.probabilities[0].set_hospitalized_per_infectious(vec[21]);
    params.probabilities[0].set_icu_per_hospitalized(vec[22]);
    params.probabilities[0].set_dead_per_icu(vec[23]);

    params.set_start_day(vec[24]);
    params.set_seasonality(vec[25]);

    EXPECT_NE(params.times[0].get_incubation().get_distribution().get(), nullptr);

    check_distribution(*vec[1].get_distribution(), *params.times[0].get_incubation().get_distribution());
    check_distribution(*vec[2].get_distribution(), *params.times[0].get_infectious_mild().get_distribution());
    check_distribution(*vec[3].get_distribution(), *params.times[0].get_serialinterval().get_distribution());
    check_distribution(*vec[4].get_distribution(), *params.times[0].get_hospitalized_to_home().get_distribution());
    check_distribution(*vec[5].get_distribution(), *params.times[0].get_home_to_hospitalized().get_distribution());
    check_distribution(*vec[6].get_distribution(), *params.times[0].get_hospitalized_to_icu().get_distribution());
    check_distribution(*vec[7].get_distribution(), *params.times[0].get_icu_to_home().get_distribution());
    check_distribution(*vec[8].get_distribution(), *params.times[0].get_infectious_asymp().get_distribution());
    check_distribution(*vec[9].get_distribution(), *params.times[0].get_icu_to_dead().get_distribution());
    check_distribution(*vec[10].get_distribution(),
                       *params.populations.get({0, epi::SecirCompartments::E}).get_distribution());
    check_distribution(*vec[11].get_distribution(),
                       *params.populations.get({0, epi::SecirCompartments::C}).get_distribution());
    check_distribution(*vec[12].get_distribution(),
                       *params.populations.get({0, epi::SecirCompartments::I}).get_distribution());
    check_distribution(*vec[13].get_distribution(),
                       *params.populations.get({0, epi::SecirCompartments::H}).get_distribution());
    check_distribution(*vec[14].get_distribution(),
                       *params.populations.get({0, epi::SecirCompartments::U}).get_distribution());
    check_distribution(*vec[15].get_distribution(),
                       *params.populations.get({0, epi::SecirCompartments::R}).get_distribution());
    check_distribution(*vec[16].get_distribution(),
                       *params.populations.get({0, epi::SecirCompartments::D}).get_distribution());
    check_distribution(*vec[17].get_distribution(),
                       *params.probabilities[0].get_infection_from_contact().get_distribution());
    check_distribution(*vec[18].get_distribution(),
                       *params.probabilities[0].get_carrier_infectability().get_distribution());
    check_distribution(*vec[19].get_distribution(),
                       *params.probabilities[0].get_asymp_per_infectious().get_distribution());
    check_distribution(*vec[20].get_distribution(),
                       *params.probabilities[0].get_risk_from_symptomatic().get_distribution());
    check_distribution(*vec[21].get_distribution(),
                       *params.probabilities[0].get_hospitalized_per_infectious().get_distribution());
    check_distribution(*vec[22].get_distribution(),
                       *params.probabilities[0].get_icu_per_hospitalized().get_distribution());
    check_distribution(*vec[23].get_distribution(), *params.probabilities[0].get_dead_per_icu().get_distribution());

    check_distribution(*vec[25].get_distribution(), *params.get_seasonality().get_distribution());

    EXPECT_EQ(vec[1], params.times[0].get_incubation());
    EXPECT_EQ(vec[2], params.times[0].get_infectious_mild());
    EXPECT_EQ(vec[3], params.times[0].get_serialinterval());
    EXPECT_EQ(vec[4], params.times[0].get_hospitalized_to_home());
    EXPECT_EQ(vec[5], params.times[0].get_home_to_hospitalized());
    EXPECT_EQ(vec[6], params.times[0].get_hospitalized_to_icu());
    EXPECT_EQ(vec[7], params.times[0].get_icu_to_home());
    EXPECT_EQ(vec[8], params.times[0].get_infectious_asymp());
    EXPECT_EQ(vec[9], params.times[0].get_icu_to_dead());
    EXPECT_EQ(vec[10], params.populations.get({0, epi::SecirCompartments::E}));
    EXPECT_EQ(vec[11], params.populations.get({0, epi::SecirCompartments::C}));
    EXPECT_EQ(vec[12], params.populations.get({0, epi::SecirCompartments::I}));
    EXPECT_EQ(vec[13], params.populations.get({0, epi::SecirCompartments::H}));
    EXPECT_EQ(vec[14], params.populations.get({0, epi::SecirCompartments::U}));
    EXPECT_EQ(vec[15], params.populations.get({0, epi::SecirCompartments::R}));
    EXPECT_EQ(vec[16], params.populations.get({0, epi::SecirCompartments::D}));
    EXPECT_EQ(vec[17], params.probabilities[0].get_infection_from_contact());
    EXPECT_EQ(vec[18], params.probabilities[0].get_carrier_infectability());
    EXPECT_EQ(vec[19], params.probabilities[0].get_asymp_per_infectious());
    EXPECT_EQ(vec[20], params.probabilities[0].get_risk_from_symptomatic());
    EXPECT_EQ(vec[21], params.probabilities[0].get_hospitalized_per_infectious());
    EXPECT_EQ(vec[22], params.probabilities[0].get_icu_per_hospitalized());
    EXPECT_EQ(vec[23], params.probabilities[0].get_dead_per_icu());

    EXPECT_EQ(vec[24], params.get_start_day());
    EXPECT_EQ(vec[25], params.get_seasonality());
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
    params.times[0].set_icu_to_death(ticu2death);

    epi::ContactFrequencyMatrix& cont_freq_matrix = params.get_contact_patterns();
    cont_freq_matrix.set_cont_freq(cont_freq, 0, 0);
    epi::Damping dummy(30., 0.3);
    cont_freq_matrix.add_damping(dummy, 0, 0);

    params.populations.set({0, epi::SecirCompartments::E}, nb_exp_t0);
    params.populations.set({0, epi::SecirCompartments::C}, nb_car_t0);
    params.populations.set({0, epi::SecirCompartments::I}, nb_inf_t0);
    params.populations.set({0, epi::SecirCompartments::H}, nb_hosp_t0);
    params.populations.set({0, epi::SecirCompartments::U}, nb_icu_t0);
    params.populations.set({0, epi::SecirCompartments::R}, nb_rec_t0);
    params.populations.set({0, epi::SecirCompartments::D}, nb_dead_t0);
    params.populations.set_difference_from_total({0, epi::SecirCompartments::S}, nb_total_t0);

    params.probabilities[0].set_infection_from_contact(inf_prob);
    params.probabilities[0].set_carrier_infectability(carr_infec);
    params.probabilities[0].set_asymp_per_infectious(alpha);
    params.probabilities[0].set_risk_from_symptomatic(beta);
    params.probabilities[0].set_hospitalized_per_infectious(rho);
    params.probabilities[0].set_icu_per_hospitalized(theta);
    params.probabilities[0].set_dead_per_icu(delta);

    epi::set_log_level(epi::LogLevel::off);
    params.check_constraints();

    EXPECT_EQ(-91, params.populations.get({0, epi::SecirCompartments::E}));
    EXPECT_EQ(2.124921, params.probabilities[0].get_asymp_per_infectious().value());
    EXPECT_NEAR(5.08993, params.times[0].get_serialinterval(), 1e-14);

    params.apply_constraints();

    EXPECT_EQ(0.0, params.populations.get({0, epi::SecirCompartments::E}));
    EXPECT_EQ(0.0, params.probabilities[0].get_asymp_per_infectious().value());
    EXPECT_NEAR(4.6, params.times[0].get_serialinterval(), 1e-14);
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

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 0,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    epi::SecirParams params;

    params.times[0].set_incubation(tinc);
    params.times[0].set_infectious_mild(tinfmild);
    params.times[0].set_serialinterval(tserint);
    params.times[0].set_hospitalized_to_home(thosp2home);
    params.times[0].set_home_to_hospitalized(thome2hosp);
    params.times[0].set_hospitalized_to_icu(thosp2icu);
    params.times[0].set_icu_to_home(ticu2home);
    params.times[0].set_icu_to_death(ticu2death);

    params.populations.set({0, epi::SecirCompartments::E}, nb_exp_t0);
    params.populations.set({0, epi::SecirCompartments::C}, nb_car_t0);
    params.populations.set({0, epi::SecirCompartments::I}, nb_inf_t0);
    params.populations.set({0, epi::SecirCompartments::H}, nb_hosp_t0);
    params.populations.set({0, epi::SecirCompartments::U}, nb_icu_t0);
    params.populations.set({0, epi::SecirCompartments::R}, nb_rec_t0);
    params.populations.set({0, epi::SecirCompartments::D}, nb_dead_t0);
    params.populations.set_difference_from_total({0, epi::SecirCompartments::S}, nb_total_t0);

    params.probabilities[0].set_infection_from_contact(inf_prob);
    params.probabilities[0].set_carrier_infectability(carr_infec);
    params.probabilities[0].set_asymp_per_infectious(alpha);
    params.probabilities[0].set_risk_from_symptomatic(beta);
    params.probabilities[0].set_hospitalized_per_infectious(rho);
    params.probabilities[0].set_icu_per_hospitalized(theta);
    params.probabilities[0].set_dead_per_icu(delta);

    epi::ContactFrequencyMatrix& cont_freq_matrix = params.get_contact_patterns();
    cont_freq_matrix.set_cont_freq(cont_freq, 0, 0);

    params.apply_constraints();

    epi::TimeSeries<double> secihurd = simulate(t0, tmax, dt, params);

    epi::TimeSeries<double> secihurd_interp = epi::interpolate_simulation_result(secihurd);

    params.set_start_day(100);
    params.set_seasonality(0.5);

    epi::TimeSeries<double> secihurd_season        = simulate(t0, tmax, dt, params);
    epi::TimeSeries<double> secihurd_season_interp = epi::interpolate_simulation_result(secihurd_season);

    for (Eigen::Index i = 0; i < secihurd_interp.get_num_time_points(); i++) {
        EXPECT_LE(secihurd_season_interp.get_value(i)[3], secihurd_interp.get_value(i)[3]) << " at row " << i;
    }

    params.set_start_day(280);

    epi::TimeSeries<double> secihurd_season2        = simulate(t0, tmax, dt, params);
    epi::TimeSeries<double> secihurd_season2_interp = epi::interpolate_simulation_result(secihurd_season2);

    for (Eigen::Index i = 0; i < secihurd_interp.get_num_time_points(); i++) {
        EXPECT_GE(secihurd_season2_interp.get_value(i)[3], secihurd_interp.get_value(i)[3]) << " at row " << i;
    }
}