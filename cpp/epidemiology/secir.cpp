#include <epidemiology/secir.h>

//#include <epidemiology/euler.h>
#include <epidemiology/adapt_rk.h>

#include <cmath>
#include <cassert>
#include <cstdio>
#include <algorithm>

namespace epi
{

void print_secir_params(const SecirParams& params)
{
    printf("\n SECIR (SECIHURD) model set.\n Parameters:\n\t Time incubation:\t %.4f \n\t Time infectious "
           "(mild):\t %.4f \n\t Serial interval:\t %.4f \n\t Time hosp.->home:\t %.4f \n\t Time home->hosp.:\t "
           "%.4f \n\t Time hosp.->icu:\t %.4f \n\t Time infectious (asymp.):\t %.4f \n\t Time icu->death:\t\t "
           "%.4f\n\t contact freq.:\t %.4f \n\t alpha:\t %.4f \n\t beta:\t %.4f \n\t delta:\t %.4f \n\t rho:\t "
           "%.4f \n\t theta:\t %.4f \n\t N0:\t %d \n\t E0:\t %d \n\t C0:\t %d\n\t I0:\t %d \n\t H0:\t %d \n\t "
           "U0:\t %d \n\t R0:\t %d \n\t D0:\t %d\n\t Calculated R_0: %.4f\n",
           1.0 / params.tinc_inv, 1.0 / params.tinfmild_inv, 1.0 / params.tserint_inv, 1.0 / params.thosp2home_inv,
           1.0 / params.thome2hosp_inv, 1.0 / params.thosp2icu_inv, 1.0 / params.tinfasy_inv,
           1.0 / params.ticu2death_inv, params.cont_freq, params.alpha, params.beta, params.delta, params.rho,
           params.theta, (int)params.populations.get_total_t0(), (int)params.populations.get_exposed_t0(),
           (int)params.populations.get_carrier_t0(), (int)params.populations.get_infectious_t0(),
           (int)params.populations.get_hospitalized_t0(), (int)params.populations.get_icu_t0(),
           (int)params.populations.get_recovered_t0(), (int)params.populations.get_dead_t0(), params.base_reprod);
}

SecirParams::Populations::Populations()
{
    m_nb_total_t0 = 0;
    m_nb_exp_t0   = 0;
    m_nb_car_t0   = 0;
    m_nb_inf_t0   = 0;
    m_nb_hosp_t0  = 0;
    m_nb_icu_t0   = 0;
    m_nb_rec_t0   = 0;
    m_nb_dead_t0  = 0;
    m_nb_sus_t0   = 0;
}

void SecirParams::Populations::set_total_t0(double nb_total_t0)
{
    m_nb_total_t0 = nb_total_t0;
    SecirParams::Populations::set_suscetible_t0();
}

void SecirParams::Populations::set_exposed_t0(double nb_exp_t0)
{
    m_nb_exp_t0 = nb_exp_t0;
    SecirParams::Populations::set_suscetible_t0();
}

void SecirParams::Populations::set_carrier_t0(double nb_car_t0)
{
    m_nb_car_t0 = nb_car_t0;
    SecirParams::Populations::set_suscetible_t0();
}

void SecirParams::Populations::set_infectious_t0(double nb_inf_t0)
{
    m_nb_inf_t0 = nb_inf_t0;
    SecirParams::Populations::set_suscetible_t0();
}

void SecirParams::Populations::set_hospital_t0(double nb_hosp_t0)
{
    m_nb_hosp_t0 = nb_hosp_t0;
    SecirParams::Populations::set_suscetible_t0();
}

void SecirParams::Populations::set_icu_t0(double nb_icu_t0)
{
    m_nb_icu_t0 = nb_icu_t0;
    SecirParams::Populations::set_suscetible_t0();
}

void SecirParams::Populations::set_recovered_t0(double nb_rec_t0)
{
    m_nb_rec_t0 = nb_rec_t0;
    SecirParams::Populations::set_suscetible_t0();
}

void SecirParams::Populations::set_dead_t0(double nb_dead_t0)
{
    m_nb_dead_t0 = nb_dead_t0;
    SecirParams::Populations::set_suscetible_t0();
}

void SecirParams::Populations::set_suscetible_t0()
{
    m_nb_sus_t0 = m_nb_total_t0 - m_nb_exp_t0 - m_nb_car_t0 - m_nb_inf_t0 - m_nb_hosp_t0 - m_nb_icu_t0 - m_nb_rec_t0 -
                  m_nb_dead_t0;
}

double SecirParams::Populations::get_total_t0() const
{
    return m_nb_total_t0;
}

double SecirParams::Populations::get_exposed_t0() const
{
    return m_nb_exp_t0;
}

double SecirParams::Populations::get_carrier_t0() const
{
    return m_nb_car_t0;
}

double SecirParams::Populations::get_infectious_t0() const
{
    return m_nb_inf_t0;
}

double SecirParams::Populations::get_hospitalized_t0() const
{
    return m_nb_hosp_t0;
}

double SecirParams::Populations::get_icu_t0() const
{
    return m_nb_icu_t0;
}

double SecirParams::Populations::get_recovered_t0() const
{
    return m_nb_rec_t0;
}

double SecirParams::Populations::get_dead_t0() const
{
    return m_nb_dead_t0;
}

double SecirParams::Populations::get_suscetible_t0() const
{
    return m_nb_sus_t0;
}

SecirParams::SecirParams()
    : populations{}
{
}

SecirParams::SecirParams(double tinc, double tinfmild, double tserint, double thosp2home, double thome2hosp,
                         double thosp2icu, double ticu2home, double tinfasy, double ticu2death, double cont_freq_in,
                         double alpha_in, double beta_in, double delta_in, double rho_in, double theta_in)
    : populations{}
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

    // only for output purposes (not explicitly used as single parameter in SECIR)
    auto dummy_R3 = 0.5 / (tinfmild - tserint);
    base_reprod   = cont_freq * ((1 - rho) * tinfmild_inv + dummy_R3 * beta * (1 - alpha) + rho * thome2hosp) /
                  ((dummy_R3 * (1 - alpha) + alpha * tinfasy_inv) * (tinfmild_inv * (1 - rho) + rho * thome2hosp_inv)) *
                  populations.get_suscetible_t0() / populations.get_total_t0();
}

void secir_get_derivatives(const SecirParams& params, const std::vector<double>& y, double t, std::vector<double>& dydt)
{

    // 0: S,      1: E,     2: C,     3: I,     4: H,     5: U,     6: R,     7: D
    double cont_freq_eff = params.cont_freq * params.dampings.get_factor(t);
    double divN          = 1.0 / params.populations.get_total_t0();

    double dummy_S = cont_freq_eff * y[0] * divN * (y[2] + params.beta * y[3]);

    double dummy_R2 = 1.0 / (2 * (1.0 / params.tserint_inv) - (1.0 / params.tinfmild_inv)); // R2 = 1/(2SI-IP)
    double dummy_R3 = 0.5 / ((1.0 / params.tinfmild_inv) - (1.0 / params.tserint_inv)); // R3 = 1/(2(IP-SI))

    dydt[0] = -dummy_S; // -R1*(C+beta*I)*S/N0
    dydt[1] = dummy_S - dummy_R2 * y[1]; // R1*(C+beta*I)*S/N0-R2*E - R2*E
    dydt[2] = dummy_R2 * y[1] - ((1 - params.alpha) * dummy_R3 + params.alpha * params.tinfasy_inv) * y[2];
    dydt[3] = (1 - params.alpha) * dummy_R3 * y[2] -
              ((1 - params.rho) * params.tinfmild_inv + params.rho * params.thome2hosp_inv) * y[3];
    dydt[4] = params.rho * params.thome2hosp_inv * y[3] -
              ((1 - params.theta) * params.thosp2home_inv + params.theta * params.thosp2icu_inv) * y[4];
    dydt[5] = params.theta * params.thosp2icu_inv * y[4] -
              ((1 - params.delta) * params.ticu2home_inv + params.delta * params.ticu2death_inv) * y[5];
    dydt[6] = params.alpha * params.tinfasy_inv * y[2] + (1 - params.rho) * params.tinfmild_inv +
              (1 - params.theta) * params.thosp2home_inv * y[4] + (1 - params.delta) * params.ticu2home_inv * y[5];
    dydt[7] = params.delta * params.ticu2death_inv * y[5];
}

std::vector<double> simulate(double t0, double tmax, double dt, const SecirParams& params,
                             std::vector<std::vector<double>>& secir)
{
    size_t n_params = 8;
    secir           = std::vector<std::vector<double>>(1, std::vector<double>(n_params, 0.));

    //initial conditions
    secir[0][0] = params.populations.get_suscetible_t0();
    secir[0][1] = params.populations.get_exposed_t0();
    secir[0][2] = params.populations.get_carrier_t0();
    secir[0][3] = params.populations.get_infectious_t0();
    secir[0][4] = params.populations.get_hospitalized_t0();
    secir[0][5] = params.populations.get_icu_t0();
    secir[0][6] = params.populations.get_recovered_t0();
    secir[0][7] = params.populations.get_dead_t0();

    auto secir_fun = [&params](std::vector<double> const& y, const double t, std::vector<double>& dydt) {
        return secir_get_derivatives(params, y, t, dydt);
    };

#ifdef ARK_H
    double dtmin = 1e-3;
    double dtmax = 1.;
    RKIntegrator integrator(secir_fun, dtmin, dtmax);
    integrator.set_abs_tolerance(1e-1);
    integrator.set_rel_tolerance(1e-4);
#endif
#ifdef EULER_H
    double dtmin = dt;
    double dtmax = dt;
    EulerIntegrator integrator(secir_fun);
#endif

    return ode_integrate(t0, tmax, dt, integrator, secir);
}

} // namespace epi
