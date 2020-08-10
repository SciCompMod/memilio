#include <epidemiology/secir.h>

// #include <epidemiology/euler.h>
#include <epidemiology/adapt_rk.h>
#include <epidemiology/eigen_util.h>
#include <epidemiology_io/secir_result_io.h>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <algorithm>

#include <iostream>
#include <fstream>

namespace epi
{

SecirParams::StageTimes::StageTimes()
    : m_tinc_inv{1.0}
    , m_tinfmild_inv{1.0}
    , m_tserint_inv{1.0}
    , m_thosp2home_inv{1.0}
    , m_thome2hosp_inv{1.0}
    , m_thosp2icu_inv{1.0}
    , m_ticu2home_inv{1.0}
    , m_tinfasy_inv{1.0}
    , m_ticu2death_inv{1.0}
{
}

void SecirParams::StageTimes::set_incubation(double const& tinc)
{
    m_tinc_inv = 1.0 / tinc;
}

void SecirParams::StageTimes::set_infectious_mild(double const& tinfmild)
{
    m_tinfmild_inv = 1.0 / tinfmild;
}

void SecirParams::StageTimes::set_serialinterval(double const& tserint)
{
    m_tserint_inv = 1.0 / tserint;
}

void SecirParams::StageTimes::set_hospitalized_to_home(double const& thosp2home)
{
    m_thosp2home_inv = 1.0 / thosp2home;
}

void SecirParams::StageTimes::set_home_to_hospitalized(double const& thome2hosp)
{
    m_thome2hosp_inv = 1.0 / thome2hosp;
}

void SecirParams::StageTimes::set_hospitalized_to_icu(double const& thosp2icu)
{
    m_thosp2icu_inv = 1.0 / thosp2icu;
}

void SecirParams::StageTimes::set_icu_to_home(double const& ticu2home)
{
    m_ticu2home_inv = 1.0 / ticu2home;
}

void SecirParams::StageTimes::set_infectious_asymp(double const& tinfasy)
{
    m_tinfasy_inv = 1.0 / tinfasy;
}

void SecirParams::StageTimes::set_icu_to_death(double const& ticu2death)
{
    m_ticu2death_inv = 1.0 / ticu2death;
}

const UncertainValue& SecirParams::StageTimes::get_incubation_inv() const
{
    return m_tinc_inv;
}

UncertainValue& SecirParams::StageTimes::get_incubation_inv()
{
    return m_tinc_inv;
}

const UncertainValue& SecirParams::StageTimes::get_infectious_mild_inv() const
{
    return m_tinfmild_inv;
}

UncertainValue& SecirParams::StageTimes::get_infectious_mild_inv()
{
    return m_tinfmild_inv;
}

const UncertainValue& SecirParams::StageTimes::get_serialinterval_inv() const
{
    return m_tserint_inv;
}

UncertainValue& SecirParams::StageTimes::get_serialinterval_inv()
{
    return m_tserint_inv;
}

const UncertainValue& SecirParams::StageTimes::get_hospitalized_to_home_inv() const
{
    return m_thosp2home_inv;
}

UncertainValue& SecirParams::StageTimes::get_hospitalized_to_home_inv()
{
    return m_thosp2home_inv;
}

const UncertainValue& SecirParams::StageTimes::get_home_to_hospitalized_inv() const
{
    return m_thome2hosp_inv;
}

UncertainValue& SecirParams::StageTimes::get_home_to_hospitalized_inv()
{
    return m_thome2hosp_inv;
}

const UncertainValue& SecirParams::StageTimes::get_hospitalized_to_icu_inv() const
{
    return m_thosp2icu_inv;
}

UncertainValue& SecirParams::StageTimes::get_hospitalized_to_icu_inv()
{
    return m_thosp2icu_inv;
}

const UncertainValue& SecirParams::StageTimes::get_icu_to_home_inv() const
{
    return m_ticu2home_inv;
}

UncertainValue& SecirParams::StageTimes::get_icu_to_home_inv()
{
    return m_ticu2home_inv;
}

const UncertainValue& SecirParams::StageTimes::get_infectious_asymp_inv() const
{
    return m_tinfasy_inv;
}

UncertainValue& SecirParams::StageTimes::get_infectious_asymp_inv()
{
    return m_tinfasy_inv;
}

const UncertainValue& SecirParams::StageTimes::get_icu_to_dead_inv() const
{
    return m_ticu2death_inv;
}

UncertainValue& SecirParams::StageTimes::get_icu_to_dead_inv()
{
    return m_ticu2death_inv;
}

SecirParams::Probabilities::Probabilities()
    : m_infprob{1}
    , m_asympinf{0}
    , m_risksymp{0}
    , m_hospinf{0}
    , m_icuhosp{0}
    , m_deathicu{0}
{
}

void SecirParams::Probabilities::set_infection_from_contact(double const& infprob)
{
    m_infprob = infprob;
}

void SecirParams::Probabilities::set_asymp_per_infectious(double const& asympinf)
{
    m_asympinf = asympinf;
}

void SecirParams::Probabilities::set_risk_from_symptomatic(double const& risksymp)
{
    m_risksymp = risksymp;
}

void SecirParams::Probabilities::set_hospitalized_per_infectious(double const& hospinf)
{
    m_hospinf = hospinf;
}

void SecirParams::Probabilities::set_icu_per_hospitalized(double const& icuhosp)
{
    m_icuhosp = icuhosp;
}

void SecirParams::Probabilities::set_dead_per_icu(double const& deathicu)
{
    m_deathicu = deathicu;
}

const UncertainValue& SecirParams::Probabilities::get_infection_from_contact() const
{
    return m_infprob;
}

UncertainValue& SecirParams::Probabilities::get_infection_from_contact()
{
    return m_infprob;
}

const UncertainValue& SecirParams::Probabilities::get_asymp_per_infectious() const
{
    return m_asympinf;
}

UncertainValue& SecirParams::Probabilities::get_asymp_per_infectious()
{
    return m_asympinf;
}

const UncertainValue& SecirParams::Probabilities::get_risk_from_symptomatic() const
{
    return m_risksymp;
}

UncertainValue& SecirParams::Probabilities::get_risk_from_symptomatic()
{
    return m_risksymp;
}

const UncertainValue& SecirParams::Probabilities::get_hospitalized_per_infectious() const
{
    return m_hospinf;
}

UncertainValue& SecirParams::Probabilities::get_hospitalized_per_infectious()
{
    return m_hospinf;
}

const UncertainValue& SecirParams::Probabilities::get_icu_per_hospitalized() const
{
    return m_icuhosp;
}

UncertainValue& SecirParams::Probabilities::get_icu_per_hospitalized()
{
    return m_icuhosp;
}

const UncertainValue& SecirParams::Probabilities::get_dead_per_icu() const
{
    return m_deathicu;
}

UncertainValue& SecirParams::Probabilities::get_dead_per_icu()
{
    return m_deathicu;
}

ContactFrequencyMatrix& SecirParams::get_cont_freq_matrix()
{
    return contact_patterns.get_cont_freq_mat();
}

ContactFrequencyMatrix const& SecirParams::get_cont_freq_matrix() const
{
    return contact_patterns.get_cont_freq_mat();
}

double get_reprod_rate(SecirParams const& params, double const t, std::vector<double> const& yt)
{
    if (params.size() == 1) {
        // (base_)reprod has to be computed time dependently !
        auto dummy_R3 =
            0.5 / (1.0 / params.times[0].get_infectious_mild_inv() - 1.0 / params.times[0].get_serialinterval_inv());

        double cont_freq_eff = params.get_cont_freq_matrix().get_cont_freq(0, 0) *
                               params.get_cont_freq_matrix().get_dampings(0, 0).get_factor(t);

        double nb_total = 0;
        for (size_t i = 0; i < yt.size(); i++) {
            nb_total += yt[0];
        }

        double reprod_rate =
            cont_freq_eff *
            ((1 - params.probabilities[0].get_hospitalized_per_infectious()) *
                 params.times[0].get_infectious_mild_inv() +
             dummy_R3 * params.probabilities[0].get_risk_from_symptomatic() *
                 (1 - params.probabilities[0].get_asymp_per_infectious()) +
             params.probabilities[0].get_hospitalized_per_infectious() *
                 params.times[0].get_home_to_hospitalized_inv()) /
            ((dummy_R3 * (1 - params.probabilities[0].get_asymp_per_infectious()) +
              params.probabilities[0].get_asymp_per_infectious() * params.times[0].get_infectious_asymp_inv()) *
             (params.times[0].get_infectious_mild_inv() *
                  (1 - params.probabilities[0].get_hospitalized_per_infectious()) +
              params.probabilities[0].get_hospitalized_per_infectious() *
                  params.times[0].get_home_to_hospitalized_inv())) *
            yt[0] / nb_total;

        return reprod_rate;
    }
    else {
        return -1;
    }
}

void secir_get_derivatives(SecirParams const& params, const Eigen::VectorXd& y, double t, Eigen::VectorXd& dydt)
{
    // alpha  // percentage of asymptomatic cases
    // beta // risk of infection from the infected symptomatic patients
    // rho   // hospitalized per infectious
    // theta // icu per hospitalized
    // delta  // deaths per ICUs
    // 0: S,      1: E,     2: C,     3: I,     4: H,     5: U,     6: R,     7: D
    size_t n_agegroups = params.size();

    for (size_t i = 0; i < n_agegroups; i++) {

        size_t Si = params.populations.get_flat_index({i, S});
        size_t Ei = params.populations.get_flat_index({i, E});
        size_t Ci = params.populations.get_flat_index({i, C});
        size_t Ii = params.populations.get_flat_index({i, I});
        size_t Hi = params.populations.get_flat_index({i, H});
        size_t Ui = params.populations.get_flat_index({i, U});
        size_t Ri = params.populations.get_flat_index({i, R});
        size_t Di = params.populations.get_flat_index({i, D});

        dydt[Si] = 0;
        dydt[Ei] = 0;
        for (size_t j = 0; j < n_agegroups; j++) {
            size_t Cj = params.populations.get_flat_index({j, C});
            size_t Ij = params.populations.get_flat_index({j, I});
            // effective contact rate by contact rate between groups i and j and damping j
            double cont_freq_eff = params.get_cont_freq_matrix().get_cont_freq(i, j) *
                                   params.get_cont_freq_matrix().get_dampings(i, j).get_factor(
                                       t); // get effective contact rate between i and j
            double divN    = 1.0 / params.populations.get_group_total(SecirCategory::AgeGroup, j); // precompute 1.0/Nj
            double dummy_S = y[Si] * cont_freq_eff * divN * params.probabilities[i].get_infection_from_contact() *
                             (y[Cj] + params.probabilities[j].get_risk_from_symptomatic() * y[Ij]);

            dydt[Si] -= dummy_S; // -R1*(C+beta*I)*S/N0
            dydt[Ei] += dummy_S; // R1*(C+beta*I)*S/N0-R2*E
        }

        double dummy_R2 = 1.0 / (2 * (1.0 / params.times[i].get_serialinterval_inv()) -
                                 (1.0 / params.times[i].get_incubation_inv())); // R2 = 1/(2SI-TINC)
        double dummy_R3 = 0.5 / ((1.0 / params.times[i].get_incubation_inv()) -
                                 (1.0 / params.times[i].get_serialinterval_inv())); // R3 = 1/(2(TINC-SI))

        dydt[Ei] -= dummy_R2 * y[Ei]; // only exchange of E and C done here
        dydt[Ci] = dummy_R2 * y[Ei] -
                   ((1 - params.probabilities[i].get_asymp_per_infectious()) * dummy_R3 +
                    params.probabilities[i].get_asymp_per_infectious() * params.times[i].get_infectious_asymp_inv()) *
                       y[Ci];
        dydt[Ii] = (1 - params.probabilities[i].get_asymp_per_infectious()) * dummy_R3 * y[Ci] -
                   ((1 - params.probabilities[i].get_hospitalized_per_infectious()) *
                        params.times[i].get_infectious_mild_inv() +
                    params.probabilities[i].get_hospitalized_per_infectious() *
                        params.times[i].get_home_to_hospitalized_inv()) *
                       y[Ii];
        dydt[Hi] =
            params.probabilities[i].get_hospitalized_per_infectious() * params.times[i].get_home_to_hospitalized_inv() *
                y[Ii] -
            ((1 - params.probabilities[i].get_icu_per_hospitalized()) * params.times[i].get_hospitalized_to_home_inv() +
             params.probabilities[i].get_icu_per_hospitalized() * params.times[i].get_hospitalized_to_icu_inv()) *
                y[Hi];
        dydt[Ui] =
            params.probabilities[i].get_icu_per_hospitalized() * params.times[i].get_hospitalized_to_icu_inv() * y[Hi] -
            ((1 - params.probabilities[i].get_dead_per_icu()) * params.times[i].get_icu_to_home_inv() +
             params.probabilities[i].get_dead_per_icu() * params.times[i].get_icu_to_dead_inv()) *
                y[Ui];
        dydt[Ri] =
            params.probabilities[i].get_asymp_per_infectious() * params.times[i].get_infectious_asymp_inv() * y[Ci] +
            (1 - params.probabilities[i].get_hospitalized_per_infectious()) *
                params.times[i].get_infectious_mild_inv() * y[Ii] +
            (1 - params.probabilities[i].get_icu_per_hospitalized()) * params.times[i].get_hospitalized_to_home_inv() *
                y[Hi] +
            (1 - params.probabilities[i].get_dead_per_icu()) * params.times[i].get_icu_to_home_inv() * y[Ui];
        dydt[Di] = params.probabilities[i].get_dead_per_icu() * params.times[i].get_icu_to_dead_inv() * y[Ui];
    }
}

std::vector<double> simulate(double t0, double tmax, double dt, SecirParams const& params,
                             std::vector<Eigen::VectorXd>& secir)
{
    SecirSimulation sim(params, t0, dt);
    sim.advance(tmax);
    secir = sim.get_y();

    return sim.get_t();
}

SecirSimulation::SecirSimulation(const SecirParams& params, double t0, double dt)
    : m_integratorCore(std::make_shared<RKIntegratorCore>(1e-3, 1.))
    , m_integrator(
          [params](auto&& y, auto&& t, auto&& dydt) {
              secir_get_derivatives(params, y, t, dydt);
          },
          t0, params.populations.get_compartments(), dt, m_integratorCore)
    , m_params(params)
{
    m_integratorCore->set_rel_tolerance(1e-4);
    m_integratorCore->set_abs_tolerance(1e-1);
}

Eigen::VectorXd& SecirSimulation::advance(double tmax)
{
    return m_integrator.advance(tmax);
}

} // namespace epi
