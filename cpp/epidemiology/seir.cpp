#include <epidemiology/seir.h>

#include <epidemiology/euler.h>
#include <epidemiology/adapt_rk.h>

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
           1.0 / params.times.get_incubation_inv(), 1.0 / params.times.get_infectious_inv(),
           params.times.get_cont_freq(), (int)params.populations.get_total_t0(),
           (int)params.populations.get_exposed_t0(), (int)params.populations.get_infectious_t0(),
           (int)params.populations.get_recovered_t0());
}

/**
 * @brief Initializes a time parameters' struct of the SEIR model
 */
SeirParams::StageTimes::StageTimes()
{
    // assume an incubation period of 5.2 days;
    // an infectious period of (nonhospitalized) people (after disease) of 6 days
    // and a basis reproduction number (R0) of 2.7
    m_tinc_inv     = 1.0 / 5.2; // 1.0/5.2 (in JS version)
    m_tinfmild_inv = 1.0 / 6.0; // 1.0/2.0 (in JS version)
    m_cont_freq    = 0.587; // taken from SECIR paper
    // base_reprod  = 2.7; // 3.5 (in JS version)
}

void SeirParams::StageTimes::set_cont_freq(double const& cont_freq)
{
    m_cont_freq = cont_freq;
}

void SeirParams::StageTimes::set_incubation(double const& tinc)
{
    m_tinc_inv = 1.0 / tinc;
}

void SeirParams::StageTimes::set_infectious(double const& tinfmild)
{
    m_tinfmild_inv = 1.0 / tinfmild;
}

double SeirParams::StageTimes::get_cont_freq() const
{
    return m_cont_freq;
}

double SeirParams::StageTimes::get_incubation_inv() const
{
    return m_tinc_inv;
}

double SeirParams::StageTimes::get_infectious_inv() const
{
    return m_tinfmild_inv;
}

SeirParams::Populations::Populations()
{
    m_nb_total_t0 = 0;
    m_nb_exp_t0   = 0;
    m_nb_inf_t0   = 0;
    m_nb_rec_t0   = 0;
    m_nb_sus_t0   = 0;
}

void SeirParams::Populations::set_total_t0(double nb_total_t0)
{
    m_nb_total_t0 = nb_total_t0;
    SeirParams::Populations::set_suscetible_t0();
}

void SeirParams::Populations::set_exposed_t0(double nb_exp_t0)
{
    m_nb_exp_t0 = nb_exp_t0;
    SeirParams::Populations::set_suscetible_t0();
}

void SeirParams::Populations::set_infectious_t0(double nb_inf_t0)
{
    m_nb_inf_t0 = nb_inf_t0;
    SeirParams::Populations::set_suscetible_t0();
}

void SeirParams::Populations::set_recovered_t0(double nb_rec_t0)
{
    m_nb_rec_t0 = nb_rec_t0;
    SeirParams::Populations::set_suscetible_t0();
}

void SeirParams::Populations::set_suscetible_t0()
{
    m_nb_sus_t0 = m_nb_total_t0 - m_nb_exp_t0 - m_nb_inf_t0 - m_nb_rec_t0;
}

double SeirParams::Populations::get_total_t0() const
{
    return m_nb_total_t0;
}

double SeirParams::Populations::get_exposed_t0() const
{
    return m_nb_exp_t0;
}

double SeirParams::Populations::get_infectious_t0() const
{
    return m_nb_inf_t0;
}

double SeirParams::Populations::get_recovered_t0() const
{
    return m_nb_rec_t0;
}

double SeirParams::Populations::get_suscetible_t0() const
{
    return m_nb_sus_t0;
}

SeirParams::SeirParams()
    : times{}
    , populations{}
{
}

void seir_get_derivatives(const SeirParams& params, const Eigen::VectorXd& y, double t, Eigen::VectorXd& dydt)
{
    // 0: S,      1: E,       2: I,     3: R
    double cont_freq_eff = params.times.get_cont_freq() * params.dampings.get_factor(t);
    double divN          = 1.0 / params.populations.get_total_t0();

    dydt[0] = -cont_freq_eff * y[0] * y[2] * divN;
    dydt[1] = cont_freq_eff * y[0] * y[2] * divN - params.times.get_incubation_inv() * y[1];
    dydt[2] = params.times.get_incubation_inv() * y[1] - params.times.get_infectious_inv() * y[2];
    dydt[3] = params.times.get_infectious_inv() * y[2];
}

std::vector<double> simulate(double t0, double tmax, double dt, const SeirParams& params,
                             std::vector<Eigen::VectorXd>& seir)
{
    size_t n_params = 4;
    seir            = std::vector<Eigen::VectorXd>(1, Eigen::VectorXd::Constant(n_params, 0));

    //initial conditions
    seir[0][0] = params.populations.get_suscetible_t0();
    seir[0][1] = params.populations.get_exposed_t0();
    seir[0][2] = params.populations.get_infectious_t0();
    seir[0][3] = params.populations.get_recovered_t0();

    auto seir_fun = [&params](Eigen::VectorXd const& y, const double t, Eigen::VectorXd& dydt) {
        return seir_get_derivatives(params, y, t, dydt);
    };

    // double dtmin = 1e-5;
    // double dtmax = 1.;
    // RKIntegrator integrator(seir_fun, dtmin, dtmax);
    // integrator.set_rel_tolerance(1e-6);
    // integrator.set_abs_tolerance(0);
    
    double dtmin = dt;
    double dtmax = dt;
    EulerIntegrator integrator(seir_fun);

    return ode_integrate(t0, tmax, dt, integrator, seir);
}

std::vector<double> simulate(double t0, double tmax, double dt, const std::vector<SeirParams>& group_params,
                             MigrationFunction migration_function, std::vector<Eigen::VectorXd>& result)
{
    auto num_groups     = group_params.size();
    auto num_vars       = 4;
    auto num_vars_total = num_vars * num_groups;
    auto num_steps      = static_cast<size_t>(std::ceil(tmax - t0)); //only report once each migration step
    auto migration      = Eigen::MatrixXd(num_groups, num_groups);
    result              = std::vector<Eigen::VectorXd>(1, Eigen::VectorXd::Constant(num_vars_total, 0));
    auto t              = std::vector<double>(1, 0);
    result.reserve(num_steps);
    t.reserve(num_steps);

    //for two groups: [S1, S2, E1, E2, I1, I2, R1, R2]
    for (size_t i = 0; i < num_groups; i++) {
        result[0][0 * num_groups + i] = group_params[i].populations.get_suscetible_t0();
        result[0][1 * num_groups + i]     = group_params[i].populations.get_exposed_t0();
        result[0][2 * num_groups + i]  = group_params[i].populations.get_infectious_t0();
        result[0][3 * num_groups + i]   = group_params[i].populations.get_recovered_t0();
    }

    while (t.back() < tmax) {
        //TODO: avoid copying the results of a single step/group
        //TODO: use vector slicing
        //TODO: variable migration interval size

        auto group_seir = Eigen::VectorXd::Zero(num_vars_total).eval();

        for (size_t i = 0; i < num_groups; i++) {
            auto seir = std::vector<Eigen::VectorXd>(1, Eigen::VectorXd::Constant(4, 0));            
            for (size_t j = 0; j < num_vars; j++) {
                seir[0][j] = result.back()[j * num_groups + i];
            }
            auto tmpT = simulate(t.back(), std::min(t.back() + 1, tmax), dt, group_params[i], seir);
            for (size_t j = 0; j < num_vars; j++) {
                group_seir[j * num_groups + i] = seir.back()[j];
            }
        }

        for (size_t i = 0; i < num_vars; i++) {
            migration_function(static_cast<SeirCompartment>(i), t.back(), migration);
            auto group_compartment = group_seir.segment(num_groups * i, num_groups);
            group_compartment      = migration * group_compartment;
        }

        result.push_back(group_seir);
        t.push_back(std::min(t.back() + 1, tmax));
    }

    return t;
}

} // namespace epi
