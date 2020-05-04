#include <epidemiology/seir.h>

#include <epidemiology/euler.h>
// #include <epidemiology/adapt_rk.h>

#include <cmath>
#include <cassert>
#include <cstdio>
#include <algorithm>

namespace epi
{

void print_seir_params(const SeirParams& params)
{
    printf("\n SEIR model set.\n Parameters:\n\t Time incubation:\t %.4f \n\t Time infectious:\t %.4f \n\t contact "
           "rate:\t %.4f \n\t N:\t %d \n\t E0:\t %d \n\t I0:\t %d \n\t R0:\t %d\n",
           1.0 / params.tinc_inv, 1.0 / params.tinfmild_inv, params.cont_freq, (int)params.nb_total_t0,
           (int)params.nb_exp_t0, (int)params.nb_inf_t0, (int)params.nb_rec_t0);
}

SeirParams::SeirParams()
{
    // assume an incubation period of 5.2 days;
    // an infectious period of (nonhospitalized) people (after disease) of 6 days
    // and a basis reproduction number (R0) of 2.7
    tinc_inv     = 1.0 / 5.2; // 1.0/5.2 (in JS version)
    tinfmild_inv = 1.0 / 6.0; // 1.0/2.0 (in JS version)
    base_reprod  = 2.7; // 3.5 (in JS version)

    // Initial Population size
    nb_total_t0 = 10000;
    // Initial Number of exposed
    nb_exp_t0 = 50.0;
    // Intial Number of infected
    nb_inf_t0 = 10.0;
    // Initial Number of recovered
    nb_rec_t0 = 0.0;

    nb_sus_t0 = nb_total_t0 - nb_exp_t0 - nb_inf_t0 - nb_rec_t0;
    // List of damping initially empty
    dampings.push_back(Damping(0.0, 1.0));
}

SeirParams::SeirParams(double tinc, double tinfmild, double base_reprod_in, double nb_total_t0_in, double nb_exp_t0_in,
                       double nb_inf_t0_in, double nb_rec_t0_in)
{
    tinc_inv     = 1.0 / tinc;
    tinfmild_inv = 1.0 / tinfmild;
    base_reprod  = base_reprod_in; // only for output purposes

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

void SeirParams::add_damping(const Damping& d)
{
    // make sure, the damping array is sorted
    insert_sorted(dampings, d, [](const Damping& d1, const Damping& d2) {
        return d1.day < d2.day;
    });
}

double SeirParams::get_damping_factor(double day) const
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

void seir_get_derivatives(const SeirParams& params, const std::vector<double>& y, double t, std::vector<double>& dydt)
{

    // 0: S,      1: E,       2: I,     3: R
    double cont_freq_eff = params.cont_freq * params.get_damping_factor(t);
    double divN          = 1.0 / params.nb_total_t0;

    dydt[0] = -cont_freq_eff * y[0] * y[2] * divN;
    dydt[1] = cont_freq_eff * y[0] * y[2] * divN - params.tinc_inv * y[1];
    dydt[2] = params.tinc_inv * y[1] - params.tinfmild_inv * y[2];
    dydt[3] = params.tinfmild_inv * y[2];
}

std::vector<double> simulate(double t0, double tmax, double dt, const SeirParams& params,
                             std::vector<std::vector<double>>& seir)
{
    size_t n_params = 4;
    seir            = std::vector<std::vector<double>>(1, std::vector<double>(n_params, 0.));

    //initial conditions
    seir[0][0] = params.nb_sus_t0;
    seir[0][1] = params.nb_exp_t0;
    seir[0][2] = params.nb_inf_t0;
    seir[0][3] = params.nb_rec_t0;

    auto seir_fun = [&params](std::vector<double> const& y, const double t, std::vector<double>& dydt) {
        return seir_get_derivatives(params, y, t, dydt);
    };

#ifdef ARK_H
    double dtmin = 1e-3;
    double dtmax = 1.;
    RKIntegrator integrator(seir_fun, dtmin, dtmax);
    integrator.set_abs_tolerance(1e-1);
    integrator.set_rel_tolerance(1e-4);
#endif
#ifdef EULER_H
    double dtmin = dt;
    double dtmax = dt;
    EulerIntegrator integrator(seir_fun);
#endif

    return ode_integrate(t0, tmax, dt, integrator, seir);
}

} // namespace epi
