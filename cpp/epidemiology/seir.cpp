#include <epidemiology/seir.h>

#include <epidemiology/euler.h>
#include <epidemiology/adapt_rk.h>
#include <epidemiology/eigen_util.h>

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
           params.times.get_cont_freq(), (int)params.populations.get_total(),
           (int)params.populations.get({SeirCompartments::E}), (int)params.populations.get({SeirCompartments::I}),
           (int)params.populations.get({SeirCompartments::R}));
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

void seir_get_derivatives(const SeirParams& params, Eigen::Ref<const Eigen::VectorXd> y, double t, Eigen::Ref<Eigen::VectorXd> dydt)
{
    double cont_freq_eff = params.times.get_cont_freq() * params.dampings.get_factor(t);
    double divN          = 1.0 / params.populations.get_total();

    dydt[SeirCompartments::S] = -cont_freq_eff * y[SeirCompartments::S] * y[SeirCompartments::I] * divN;
    dydt[SeirCompartments::E] = cont_freq_eff * y[SeirCompartments::S] * y[SeirCompartments::I] * divN -
                                params.times.get_incubation_inv() * y[SeirCompartments::E];
    dydt[SeirCompartments::I] = params.times.get_incubation_inv() * y[SeirCompartments::E] -
                                params.times.get_infectious_inv() * y[SeirCompartments::I];
    dydt[SeirCompartments::R] = params.times.get_infectious_inv() * y[SeirCompartments::I];
}

std::vector<double> simulate(double t0, double tmax, double dt, const SeirParams& params,
                             std::vector<Eigen::VectorXd>& seir)
{
    SeirSimulation sim(params, t0, dt);
    sim.advance(tmax);
    auto& result = sim.get_result();
    std::vector<double> t(result.get_num_time_points());
    for (Eigen::Index i = 0; i < result.get_num_time_points(); i++)
    {
        t[i] = result.get_time(i);
    }
    std::transform(result.begin(), result.end(), std::back_inserter(seir), [](auto&& v_ref) { return v_ref.eval(); });
    return t;
}

SeirSimulation::SeirSimulation(const SeirParams& params, double t0, double dt_init)
    : m_integrator([params](auto&& y, auto&& t, auto&& dydt) { seir_get_derivatives(params, y, t, dydt); }, t0,
                   params.populations.get_compartments(), dt_init,
                   std::make_shared<EulerIntegratorCore>() /*std::make_shared<RkIntegratorCore>(1e-6, 1.)*/)
{
}

Eigen::Ref<Eigen::VectorXd> SeirSimulation::advance(double tmax)
{
    return m_integrator.advance(tmax);
}

} // namespace epi
