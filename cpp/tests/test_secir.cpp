#include "load_test_data.h"
#include "epidemiology/secir.h"
#include <gtest/gtest.h>

TEST(TestSecir, compareWithPreviousRun)
{
    double t0   = 0;
    double tmax = 5;
    double dt   = 0.1;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           tinfasy = 6.2, ticu2death = 5;

    double cont_freq = 0.5, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    epi::SecirParams params(tinc, tinfmild, tserint, thosp2home, thome2hosp, thosp2icu, ticu2home, tinfasy, ticu2death,
                            cont_freq, alpha, beta, delta, rho, theta, nb_total_t0, nb_exp_t0, nb_car_t0, nb_inf_t0,
                            nb_hosp_t0, nb_icu_t0, nb_rec_t0, nb_dead_t0);

    params.dampings.add(epi::Damping(30., 0.3));

    std::vector<std::vector<double>> secihurd(0);
    auto t = simulate(t0, tmax, dt, params, secihurd);

    auto compare = load_test_data_csv<double>("data/secihurd-compare.csv");

    ASSERT_EQ(compare.size(), t.size());
    ASSERT_EQ(compare.size(), secihurd.size());
    for (size_t i = 0; i < compare.size(); i++) {
        ASSERT_EQ(compare[i].size(), secihurd[i].size() + 1) << "at row " << i;
        ASSERT_FLOAT_EQ(t[i], compare[i][0]) << "at row " << i;
        for (size_t j = 1; j < compare[i].size(); j++) {
            EXPECT_FLOAT_EQ(secihurd[i][j - 1], compare[i][j]) << " at row " << i;
        }
    }
}