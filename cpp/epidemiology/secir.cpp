#include <epidemiology/secir.h>

//#include <epidemiology/euler.h>
#include <epidemiology/adapt_rk.h>

#include <cmath>
#include <cassert>
#include <cstdio>
#include <algorithm>

namespace epi
{

Damping::Damping(double day_in, double factor_in)
    : day(day_in)
    , factor(factor_in)
{
}

template <typename T, typename Pred>
typename std::vector<T>::iterator insert_sorted(std::vector<T>& vec, T const& item, Pred pred)
{
    return vec.insert(std::upper_bound(vec.begin(), vec.end(), item, pred), item);
}

void print_secir_params(const SecirParams& params)
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

SecirParams::SecirParams()
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

SecirParams::SecirParams(double tinc, double tinfmild, double base_reprod_in, double nb_total_t0_in,
                         double nb_exp_t0_in, double nb_inf_t0_in, double nb_rec_t0_in)
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

SecirParams::SecirParams(double tinc, double tinfmild, double tserint, double thosp2home, double thome2hosp,
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

void SecirParams::add_damping(const Damping& d)
{
    // make sure, the damping array is sorted
    insert_sorted(dampings, d, [](const Damping& d1, const Damping& d2) {
        return d1.day < d2.day;
    });
}

double SecirParams::get_damping_factor(double day) const
{
    // we assume, that the data_array is ordered in ascending order
    size_t ilow  = 0;
    size_t ihigh = dampings.size() - 1;

    // check extrapolation cases
    if (day < dampings[ilow].day) {
        return dampings[ilow].factor;
    }

    if (day >= dampings[ihigh].day) {
        return dampings[ihigh].factor;
    }

    // now do the search
    while (ilow < ihigh - 1) {
        size_t imid = (ilow + ihigh) / 2;
        if (dampings[ilow].day <= day && day < dampings[imid].day) {
            ihigh = imid;
        }
        else if (dampings[imid].day <= day && day < dampings[ihigh].day) {
            ilow = imid;
        }
        else {
            // this case can only occur, if
            // input data are not ordered
            return 1e16;
        }
    }

    assert(dampings[ilow].day <= day && day < dampings[ilow + 1].day);
    return dampings[ilow].factor;
}

void secir_get_derivatives(const SecirParams& params, const std::vector<double>& y, double t, std::vector<double>& dydt)
{

    // 0: S,      1: E,     2: C,     3: I,     4: H,     5: U,     6: R,     7: D
    double cont_freq_eff = params.cont_freq * params.get_damping_factor(t);
    double divN          = 1.0 / params.nb_total_t0;

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
                             std::vector<std::vector<double>>& seir)
{
    size_t n_params = params.model == 0 ? 4 : 8;
    seir            = std::vector<std::vector<double>>(1, std::vector<double>(n_params, 0.));

    if (params.model == 0) {
        //initial conditions
        seir[0][0] = params.nb_sus_t0;
        seir[0][1] = params.nb_exp_t0;
        seir[0][2] = params.nb_inf_t0;
        seir[0][3] = params.nb_rec_t0;
    }
    else {
        //initial conditions
        seir[0][0] = params.nb_sus_t0;
        seir[0][1] = params.nb_exp_t0;
        seir[0][2] = params.nb_car_t0;
        seir[0][3] = params.nb_inf_t0;
        seir[0][4] = params.nb_hosp_t0;
        seir[0][5] = params.nb_icu_t0;
        seir[0][6] = params.nb_rec_t0;
        seir[0][7] = params.nb_dead_t0;
    }

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
    EulerIntegrator<double> integrator(secir_fun);
#endif

    return ode_integrate(t0, tmax, dt, integrator, seir);
}

} // namespace epi
