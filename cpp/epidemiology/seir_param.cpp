#include <epidemiology/seir_param.h>

Damping::Damping(double day_in, double factor_in)
    : day(day_in)
    , factor(factor_in)
{
}

void printSeirParams(const SeirParams& params)
{
    if (params.model == 0) {
        printf("\n SEIR model set.\n Parameters:\n\t Time incubation:\t %.4f \n\t Time infectious:\t %.4f \n\t b:\t "
               "%.4f \n\t N:\t %d \n\t E0:\t %d \n\t I0:\t %d \n\t R0:\t %d\n",
               1.0 / params.tinc_inv, 1.0 / params.tinfmild_inv, params.b, (int)params.nb_total_t0,
               (int)params.nb_exp_t0, (int)params.nb_inf_t0, (int)params.nb_rec_t0);
    }
    else {
        printf("\n SECIR (SECIHURD) model set.\n Parameters:\n\t Time incubation:\t %.4f \n\t Time infectious "
               "(mild):\t %.4f \n\t Serial interval:\t %.4f \n\t Time hosp.->home:\t %.4f \n\t Time home->hosp.:\t "
               "%.4f \n\t Time hosp.->icu:\t %.4f \n\t Time infectious (asymp.):\t %.4f \n\t Time icu->death:\t\t "
               "%.4f\n\t contact freq.:\t %.4f \n\t alpha:\t %.4f \n\t beta:\t %.4f \n\t delta:\t %.4f \n\t rho:\t "
               "%.4f \n\t theta:\t %.4f \n\t N0:\t %d \n\t E0:\t %d \n\t C0:\t %d\n\t I0:\t %d \n\t H0:\t %d \n\t "
               "U0:\t %d \n\t R0:\t %d \n\t D0:\t %d\n\t Calculated R_0: %.4f\n",
               1.0 / params.tinc_inv, 1.0 / params.tinfmild_inv, 1.0 / params.tserint_inv, 1.0 / params.thosp2home_inv,
               1.0 / params.thome2hosp_inv, 1.0 / params.thosp2icu_inv, 1.0 / params.tinfasy_inv,
               1.0 / params.ticu2death_inv, params.cont_freq, params.alpha, params.beta, params.delta, params.rho,
               params.theta, (int)params.nb_total_t0, (int)params.nb_exp_t0, (int)params.nb_car_t0,
               (int)params.nb_inf_t0, (int)params.nb_hosp_t0, (int)params.nb_icu_t0, (int)params.nb_rec_t0,
               (int)params.nb_dead_t0, params.base_reprod);
    }
}

SeirParams::SeirParams()
{
    // assume an incubation period of 5.2 days;
    // an infectious period of (nonhospitalized) people (after disease) of 6 days
    // and a basis reproduction number (R0) of 2.7
    tinc_inv     = 1.0 / 5.2; // 1.0/5.2 (in JS version)
    tinfmild_inv = 1.0 / 6.0; // 1.0/2.0 (in JS version)
    base_reprod  = 2.7; // 3.5 (in JS version)

    // contact rate beta
    b = base_reprod * tinfmild_inv;

    model = 0; // standard is SEIR

    // Initial Population size
    nb_total_t0 = 10000;
    // Initial Number of exposed
    nb_exp_t0 = 50.0;
    // Intial Number of infected
    nb_inf_t0 = 10.0;
    // Initial Number of recovered
    nb_rec_t0 = 0.0;

    nb_sus_t0 = nb_total_t0 - nb_exp_t0 - nb_inf_t0 - nb_rec_t0;

    // additional numbers for SECIR (SECIHURD) model
    nb_car_t0  = 0;
    nb_hosp_t0 = 0;
    nb_icu_t0  = 0;
    nb_dead_t0 = 0;
    // List of damping initially empty
    dampings.push_back(Damping(0.0, 1.0));
}

SeirParams::SeirParams(double tinc, double tinfmild, double base_reprod_in, double nb_total_t0_in, double nb_exp_t0_in,
                       double nb_inf_t0_in, double nb_rec_t0_in)
{
    tinc_inv     = 1.0 / tinc;
    tinfmild_inv = 1.0 / tinfmild;
    base_reprod  = base_reprod_in;

    // contact rate beta
    b = base_reprod * tinfmild_inv;

    model = 0; // standard is SEIR

    // Initial Population size
    nb_total_t0 = nb_total_t0_in;
    // Initial Number of exposed
    nb_exp_t0 = nb_exp_t0_in;
    // Intial Number of infected
    nb_inf_t0 = nb_inf_t0_in;
    // Initial Number of recovered
    nb_rec_t0 = nb_rec_t0_in;
    // List of damping initially empty
    dampings.push_back(Damping(0.0, 1.0));
}

SeirParams::SeirParams(double tinc, double tinfmild, double tserint, double thosp2home, double thome2hosp,
                       double thosp2icu, double ticu2home, double tinfasy, double ticu2death, double cont_freq_in,
                       double alpha_in, double beta_in, double delta_in, double rho_in, double theta_in,
                       double nb_total_t0_in, double nb_exp_t0_in, double nb_car_t0_in, double nb_inf_t0_in,
                       double nb_hosp_t0_in, double nb_icu_t0_in, double nb_rec_t0_in, double nb_dead_t0_in)
{

    tinc_inv       = 1.0 / tinc; // inverse incubation period
    tinfmild_inv   = 1.0 / tinfmild; // inverse infectious period (of nonhospitalized cases)
    tserint_inv    = 1.0 / tserint; // inverse serial interval
    thosp2home_inv = 1.0 / thosp2home; // inverse hospital to home time
    thome2hosp_inv = 1.0 / thome2hosp; // inverse home to hospital time
    thosp2icu_inv  = 1.0 / thosp2icu; // inverse hospital to ICU time
    ticu2home_inv  = 1.0 / ticu2home; // inverse icu to home time
    tinfasy_inv    = 1.0 / tinfasy; // inverse infectious period for asymptomatic cases
    ticu2death_inv = 1.0 / ticu2death; // inverse ICU to death time

    cont_freq = cont_freq_in; // contact frequency
    alpha     = alpha_in; // percentage of asymptomatic cases
    beta      = beta_in; // risk of infection from the infected symptomatic patients
    rho       = rho_in; // hospitalized per infected
    theta     = theta_in; // icu per hospitalized
    delta     = delta_in; // deaths per ICUs

    model = 1;

    // Initial Population size
    nb_total_t0 = nb_total_t0_in;
    // Initial Number of exposed
    nb_exp_t0 = nb_exp_t0_in;
    // Intial Number of infected
    nb_inf_t0 = nb_inf_t0_in;
    // Initial Number of recovered
    nb_rec_t0 = nb_rec_t0_in;

    // Initial Number of carriers
    nb_car_t0 = nb_car_t0_in;
    // Initial Number of hospitalized
    nb_hosp_t0 = nb_hosp_t0_in;
    // Initial Number of ICU
    nb_icu_t0 = nb_icu_t0_in;
    // Initial Number of deaths
    nb_dead_t0 = nb_dead_t0_in;
    // List of damping initially empty
    dampings.push_back(Damping(0.0, 1.0));

    auto dummy_R3 = 0.5 / (tinfmild - tserint);

    nb_sus_t0 = nb_total_t0 - nb_exp_t0 - nb_car_t0 - nb_inf_t0 - nb_hosp_t0 - nb_icu_t0 - nb_rec_t0 - nb_dead_t0;

    // only for output... (not explicitly used as single parameter in SECIR)
    base_reprod = cont_freq * ((1 - rho) * tinfmild_inv + dummy_R3 * beta * (1 - alpha) + rho * thome2hosp) /
                  ((dummy_R3 * (1 - alpha) + alpha * tinfasy_inv) * (tinfmild_inv * (1 - rho) + rho * thome2hosp_inv)) *
                  nb_sus_t0 / nb_total_t0;
}

void SeirParams::add_damping(const Damping& d)
{
    dampings.push_back(d);
}
