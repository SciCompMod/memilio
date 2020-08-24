#include "load_test_data.h"
#include "epidemiology/secir.h"
#include <epidemiology_io/secir_result_io.h>
#include <gtest/gtest.h>

TEST(TestSaveResult, compareResultWithH5)
{
    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           tinfasy = 6.2, ticu2death = 5;

    double cont_freq = 0.5, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    int nb_groups = 1;
    double fact   = 1.0 / (double)nb_groups;

    epi::SecirParams params(nb_groups);

    for (size_t i = 0; i < nb_groups; i++) {
        params.times[i].set_incubation(tinc);
        params.times[i].set_infectious_mild(tinfmild);
        params.times[i].set_serialinterval(tserint);
        params.times[i].set_hospitalized_to_home(thosp2home);
        params.times[i].set_home_to_hospitalized(thome2hosp);
        params.times[i].set_hospitalized_to_icu(thosp2icu);
        params.times[i].set_icu_to_home(ticu2home);
        params.times[i].set_infectious_asymp(tinfasy);
        params.times[i].set_icu_to_death(ticu2death);

        params.populations.set({i, epi::SecirCompartments::E}, fact * nb_exp_t0);
        params.populations.set({i, epi::SecirCompartments::C}, fact * nb_car_t0);
        params.populations.set({i, epi::SecirCompartments::I}, fact * nb_inf_t0);
        params.populations.set({i, epi::SecirCompartments::H}, fact * nb_hosp_t0);
        params.populations.set({i, epi::SecirCompartments::U}, fact * nb_icu_t0);
        params.populations.set({i, epi::SecirCompartments::R}, fact * nb_rec_t0);
        params.populations.set({i, epi::SecirCompartments::D}, fact * nb_dead_t0);
        params.populations.set_difference_from_group_total({i, epi::SecirCompartments::S}, epi::SecirCategory::AgeGroup,
                                                           i, fact * nb_total_t0);

        params.probabilities[i].set_asymp_per_infectious(alpha);
        params.probabilities[i].set_risk_from_symptomatic(beta);
        params.probabilities[i].set_hospitalized_per_infectious(rho);
        params.probabilities[i].set_icu_per_hospitalized(theta);
        params.probabilities[i].set_dead_per_icu(delta);
    }

    epi::ContactFrequencyMatrix& cont_freq_matrix = params.get_contact_patterns();
    epi::Damping dummy(30., 0.3);
    for (int i = 0; i < nb_groups; i++) {
        for (int j = i; j < nb_groups; j++) {
            cont_freq_matrix.set_cont_freq(fact * cont_freq, i, j);
            cont_freq_matrix.add_damping(dummy, i, j);
        }
    }

    std::vector<Eigen::VectorXd> secihurd(0);
    auto t = simulate(t0, tmax, dt, params, secihurd);

    epi::save_result(t, secihurd, "test_result.h5");

    epi::SecirSimulationResult test_result{epi::read_result("test_result.h5", nb_groups)};

    ASSERT_EQ(test_result.get_time_vector().size(), t.size());
    ASSERT_EQ(test_result.get_groups_vectors().size(), secihurd.size());
    for (size_t i = 0; i < test_result.get_time_vector().size(); i++) {
        ASSERT_EQ(test_result.get_groups_vectors()[i].size(), secihurd[i].size()) << "at row " << i;
        ASSERT_NEAR(t[i], test_result.get_time_vector()[i], 1e-10) << "at row " << i;
        for (size_t l = 0; l < test_result.get_groups_vectors()[i].size() / nb_groups; l++) {
            double dummy = 0.0;
            for (size_t j = 0; j < nb_groups; j++) {
                dummy += secihurd[i][j * epi::SecirCompartments::SecirCount + l];
                EXPECT_NEAR(test_result.get_groups_vectors()[i][j * epi::SecirCompartments::SecirCount + l],
                            secihurd[i][j * epi::SecirCompartments::SecirCount + l], 1e-10)
                    << " at row " << i << " at row " << l << " at Group " << j;
            }
            EXPECT_NEAR(test_result.get_totals_vector()[i][l], dummy, 1e-10) << " at row " << i << " at row " << l;
        }
    }
}
