#include "epidemiology/secir/seir.h"
#include "epidemiology/math/euler.h"
#include "epidemiology/math/adapt_rk.h"
#include "epidemiology/utils/eigen_util.h"

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
           (int)params.populations.get({SeirCompartments::SeirE}), (int)params.populations.get({SeirCompartments::SeirI}),
           (int)params.populations.get({SeirCompartments::SeirR}));
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

void seir_get_derivatives(const SeirParams& params, Eigen::Ref<const Eigen::VectorXd> pop,
                          Eigen::Ref<const Eigen::VectorXd> y, double t, Eigen::Ref<Eigen::VectorXd> dydt)
{
    double cont_freq_eff = params.times.get_cont_freq() * params.dampings.get_factor(t);
    double divN          = 1.0 / params.populations.get_total();

    dydt[SeirCompartments::SeirS] = -cont_freq_eff * y[SeirCompartments::SeirS] * pop[SeirCompartments::SeirI] * divN;
    dydt[SeirCompartments::SeirE] = cont_freq_eff * y[SeirCompartments::SeirS] * pop[SeirCompartments::SeirI] * divN -
                                params.times.get_incubation_inv() * y[SeirCompartments::SeirE];
    dydt[SeirCompartments::SeirI] = params.times.get_incubation_inv() * y[SeirCompartments::SeirE] -
                                params.times.get_infectious_inv() * y[SeirCompartments::SeirI];
    dydt[SeirCompartments::SeirR] = params.times.get_infectious_inv() * y[SeirCompartments::SeirI];
}

TimeSeries<double> simulate(double t0, double tmax, double dt, const SeirParams& params)
{
    SeirSimulation sim(params, t0, dt);
    sim.advance(tmax);
    return sim.get_result();
}

SeirSimulation::SeirSimulation(const SeirParams& params, double t0, double dt_init)
    : m_params(params)
    , m_integrator(
          [params](auto&& y, auto&& t, auto&& dydt) {
              seir_get_derivatives(params, y, t, dydt);
          },
          t0, m_params.populations.get_compartments(), dt_init,
          std::make_shared<EulerIntegratorCore>() /*std::make_shared<RkIntegratorCore>(1e-6, 1.)*/)
{
}

Eigen::Ref<Eigen::VectorXd> SeirSimulation::advance(double tmax)
{
    return m_integrator.advance(tmax);
}

} // namespace epi
