#include "load_test_data.h"
#include "epidemiology/secir/secir.h"
#include "epidemiology/math/euler.h"
#include <gtest/gtest.h>

TEST(TestImplicitEuler, compareOneTimeStep)
{
    double t0 = 0;
    double dt = 0.1;

    // working_params
    double tinc    = 5.2, // R_2^(-1)+R_3^(-1)
        tinfmild   = 6, // 4-14  (=R4^(-1))
        tserint    = 4.2, // 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        thosp2home = 12, // 7-16 (=R5^(-1))
        thome2hosp = 5, // 2.5-7 (=R6^(-1))
        thosp2icu  = 2, // 1-3.5 (=R7^(-1))
        ticu2home  = 8, // 5-16 (=R8^(-1))
        tinfasy    = 6.2, // (=R9^(-1)=R_3^(-1)+0.5*R_4^(-1))
        ticu2death = 5; // 3.5-7 (=R5^(-1))

    double cont_freq = 0.5, // 0.2-0.75
        alpha        = 0.09, // 0.01-0.16
        beta         = 0.25, // 0.05-0.5
        delta        = 0.3, // 0.15-0.77
        rho          = 0.2, // 0.1-0.35
        theta        = 0.25; // 0.15-0.4

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

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

    epi::ContactMatrixGroup& contact_matrix = params.get_contact_patterns();
    contact_matrix[0]                       = epi::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix.add_damping(0.7, epi::SimulationTime(30.));

    params.populations.set_total(nb_total_t0);
    params.populations.set({0, epi::SecirCompartments::E}, nb_exp_t0);
    params.populations.set({0, epi::SecirCompartments::C}, nb_car_t0);
    params.populations.set({0, epi::SecirCompartments::I}, nb_inf_t0);
    params.populations.set({0, epi::SecirCompartments::H}, nb_hosp_t0);
    params.populations.set({0, epi::SecirCompartments::U}, nb_icu_t0);
    params.populations.set({0, epi::SecirCompartments::R}, nb_rec_t0);
    params.populations.set({0, epi::SecirCompartments::D}, nb_dead_t0);
    params.populations.set_difference_from_total({0, epi::SecirCompartments::S}, nb_total_t0);

    params.probabilities[0].set_infection_from_contact(1.0);
    params.probabilities[0].set_asymp_per_infectious(alpha);
    params.probabilities[0].set_risk_from_symptomatic(beta);
    params.probabilities[0].set_hospitalized_per_infectious(rho);
    params.probabilities[0].set_icu_per_hospitalized(theta);
    params.probabilities[0].set_dead_per_icu(delta);

    epi::DerivFunction dummy_f; // only required for explicit time integrators

    auto y0 = params.populations.get_compartments();
    Eigen::VectorXd y1(y0.size()); // solution at time t=\Delta t=0.1

    epi::ImplicitEulerIntegratorCore(params)
        .step(dummy_f, y0, t0, dt, y1); // just one iteration of implicit Euler scheme

    EXPECT_NEAR(y1[0], 9756.897564185243, 1e-10);
    EXPECT_NEAR(y1[1], 99.97811957794602, 1e-10);
    EXPECT_NEAR(y1[2], 50.74190209182218, 1e-10);
    EXPECT_NEAR(y1[3], 51.41751953982101, 1e-10);
    EXPECT_NEAR(y1[4], 19.833786579788253, 1e-10);
    EXPECT_NEAR(y1[5], 10.098962633404634, 1e-10);
    EXPECT_NEAR(y1[6], 10.97155161617429, 1e-10);
    EXPECT_NEAR(y1[7], 0.0605937758004278, 1e-10);
}
