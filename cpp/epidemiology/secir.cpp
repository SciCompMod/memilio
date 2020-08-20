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
    : m_tinc{1.0}
    , m_tinfmild{1.0}
    , m_tserint{1.0}
    , m_thosp2home{1.0}
    , m_thome2hosp{1.0}
    , m_thosp2icu{1.0}
    , m_ticu2home{1.0}
    , m_tinfasy{1.0}
    , m_ticu2death{1.0}
{
}

void SecirParams::StageTimes::set_incubation(double tinc)
{
    m_tinc = tinc;
}

void SecirParams::StageTimes::set_incubation(ParameterDistribution const& tinc)
{
    m_tinc.set_distribution(tinc);
}

void SecirParams::StageTimes::set_infectious_mild(double tinfmild)
{
    m_tinfmild = tinfmild;
}

void SecirParams::StageTimes::set_infectious_mild(ParameterDistribution const& tinfmild)
{
    m_tinfmild.set_distribution(tinfmild);
}

void SecirParams::StageTimes::set_serialinterval(double tserint)
{
    m_tserint = tserint;
}

void SecirParams::StageTimes::set_serialinterval(ParameterDistribution const& tserint)
{
    m_tserint.set_distribution(tserint);
}

void SecirParams::StageTimes::set_hospitalized_to_home(double thosp2home)
{
    m_thosp2home = thosp2home;
}

void SecirParams::StageTimes::set_hospitalized_to_home(ParameterDistribution const& thosp2home)
{
    m_thosp2home.set_distribution(thosp2home);
}

void SecirParams::StageTimes::set_home_to_hospitalized(double thome2hosp)
{
    m_thome2hosp = thome2hosp;
}

void SecirParams::StageTimes::set_home_to_hospitalized(ParameterDistribution const& thome2hosp)
{
    m_thome2hosp.set_distribution(thome2hosp);
}

void SecirParams::StageTimes::set_hospitalized_to_icu(double thosp2icu)
{
    m_thosp2icu = thosp2icu;
}

void SecirParams::StageTimes::set_hospitalized_to_icu(ParameterDistribution const& thosp2icu)
{
    m_thosp2icu.set_distribution(thosp2icu);
}

void SecirParams::StageTimes::set_icu_to_home(double ticu2home)
{
    m_ticu2home = ticu2home;
}

void SecirParams::StageTimes::set_icu_to_home(ParameterDistribution const& ticu2home)
{
    m_ticu2home.set_distribution(ticu2home);
}

void SecirParams::StageTimes::set_infectious_asymp(double tinfasy)
{
    m_tinfasy = tinfasy;
}

void SecirParams::StageTimes::set_infectious_asymp(ParameterDistribution const& tinfasy)
{
    m_tinfasy.set_distribution(tinfasy);
}

void SecirParams::StageTimes::set_icu_to_death(double ticu2death)
{
    m_ticu2death = ticu2death;
}

void SecirParams::StageTimes::set_icu_to_death(ParameterDistribution const& ticu2death)
{
    m_ticu2death.set_distribution(ticu2death);
}

const UncertainValue& SecirParams::StageTimes::get_incubation() const
{
    return m_tinc;
}

UncertainValue& SecirParams::StageTimes::get_incubation()
{
    return m_tinc;
}

const UncertainValue& SecirParams::StageTimes::get_infectious_mild() const
{
    return m_tinfmild;
}

UncertainValue& SecirParams::StageTimes::get_infectious_mild()
{
    return m_tinfmild;
}

const UncertainValue& SecirParams::StageTimes::get_serialinterval() const
{
    return m_tserint;
}

UncertainValue& SecirParams::StageTimes::get_serialinterval()
{
    return m_tserint;
}

const UncertainValue& SecirParams::StageTimes::get_hospitalized_to_home() const
{
    return m_thosp2home;
}

UncertainValue& SecirParams::StageTimes::get_hospitalized_to_home()
{
    return m_thosp2home;
}

const UncertainValue& SecirParams::StageTimes::get_home_to_hospitalized() const
{
    return m_thome2hosp;
}

UncertainValue& SecirParams::StageTimes::get_home_to_hospitalized()
{
    return m_thome2hosp;
}

const UncertainValue& SecirParams::StageTimes::get_hospitalized_to_icu() const
{
    return m_thosp2icu;
}

UncertainValue& SecirParams::StageTimes::get_hospitalized_to_icu()
{
    return m_thosp2icu;
}

const UncertainValue& SecirParams::StageTimes::get_icu_to_home() const
{
    return m_ticu2home;
}

UncertainValue& SecirParams::StageTimes::get_icu_to_home()
{
    return m_ticu2home;
}

const UncertainValue& SecirParams::StageTimes::get_infectious_asymp() const
{
    return m_tinfasy;
}

UncertainValue& SecirParams::StageTimes::get_infectious_asymp()
{
    return m_tinfasy;
}

const UncertainValue& SecirParams::StageTimes::get_icu_to_dead() const
{
    return m_ticu2death;
}

UncertainValue& SecirParams::StageTimes::get_icu_to_dead()
{
    return m_ticu2death;
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

void SecirParams::Probabilities::set_infection_from_contact(double infprob)
{
    m_infprob = infprob;
}

void SecirParams::Probabilities::set_infection_from_contact(ParameterDistribution const& infprob)
{
    m_infprob.set_distribution(infprob);
}

void SecirParams::Probabilities::set_asymp_per_infectious(double asympinf)
{
    m_asympinf = asympinf;
}

void SecirParams::Probabilities::set_asymp_per_infectious(ParameterDistribution const& asympinf)
{
    m_asympinf.set_distribution(asympinf);
}

void SecirParams::Probabilities::set_risk_from_symptomatic(double risksymp)
{
    m_risksymp = risksymp;
}

void SecirParams::Probabilities::set_risk_from_symptomatic(ParameterDistribution const& risksymp)
{
    m_risksymp.set_distribution(risksymp);
}

void SecirParams::Probabilities::set_hospitalized_per_infectious(double hospinf)
{
    m_hospinf = hospinf;
}

void SecirParams::Probabilities::set_hospitalized_per_infectious(ParameterDistribution const& hospinf)
{
    m_hospinf.set_distribution(hospinf);
}

void SecirParams::Probabilities::set_icu_per_hospitalized(double icuhosp)
{
    m_icuhosp = icuhosp;
}

void SecirParams::Probabilities::set_icu_per_hospitalized(ParameterDistribution const& icuhosp)
{
    m_icuhosp.set_distribution(icuhosp);
}

void SecirParams::Probabilities::set_dead_per_icu(double deathicu)
{
    m_deathicu = deathicu;
}

void SecirParams::Probabilities::set_dead_per_icu(ParameterDistribution const& deathicu)
{
    m_deathicu.set_distribution(deathicu);
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

void SecirParams::set_contact_patterns(UncertainContactMatrix contact_patterns)
{
    m_contact_patterns = contact_patterns;
}

UncertainContactMatrix& SecirParams::get_contact_patterns()
{
    return m_contact_patterns;
}

UncertainContactMatrix const& SecirParams::get_contact_patterns() const
{
    return m_contact_patterns;
}

double get_reprod_rate(SecirParams const& params, double const t, std::vector<double> const& yt)
{
    if (params.get_num_groups() == 1) {
        // (base_)reprod has to be computed time dependently !
        auto dummy_R3 = 0.5 / (params.times[0].get_infectious_mild() - params.times[0].get_serialinterval());

        ContactFrequencyMatrix const& cont_freq_matrix = params.get_contact_patterns();

        double cont_freq_eff = cont_freq_matrix.get_cont_freq(0, 0) * cont_freq_matrix.get_dampings(0, 0).get_factor(t);

        double nb_total = 0;
        for (size_t i = 0; i < yt.size(); i++) {
            nb_total += yt[0];
        }

        double reprod_rate =
            cont_freq_eff *
            ((1 - params.probabilities[0].get_hospitalized_per_infectious()) / params.times[0].get_infectious_mild() +
             dummy_R3 * params.probabilities[0].get_risk_from_symptomatic() *
                 (1 - params.probabilities[0].get_asymp_per_infectious()) +
             params.probabilities[0].get_hospitalized_per_infectious() / params.times[0].get_home_to_hospitalized()) /
            ((dummy_R3 * (1 - params.probabilities[0].get_asymp_per_infectious()) +
              params.probabilities[0].get_asymp_per_infectious() / params.times[0].get_infectious_asymp()) *
             ((1 - params.probabilities[0].get_hospitalized_per_infectious()) / params.times[0].get_infectious_mild() +
              params.probabilities[0].get_hospitalized_per_infectious() / params.times[0].get_home_to_hospitalized())) *
            yt[0] / nb_total;

        return reprod_rate;
    }
    else {
        return -1;
    }
}

void secir_get_derivatives(SecirParams const& params, Eigen::Ref<const Eigen::VectorXd> y, double t,
                           Eigen::Ref<Eigen::VectorXd> dydt)
{
    // alpha  // percentage of asymptomatic cases
    // beta // risk of infection from the infected symptomatic patients
    // rho   // hospitalized per infectious
    // theta // icu per hospitalized
    // delta  // deaths per ICUs
    // 0: S,      1: E,     2: C,     3: I,     4: H,     5: U,     6: R,     7: D
    size_t n_agegroups = params.get_num_groups();

    ContactFrequencyMatrix const& cont_freq_matrix = params.get_contact_patterns();

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
            double cont_freq_eff =
                cont_freq_matrix.get_cont_freq(i, j) *
                cont_freq_matrix.get_dampings(i, j).get_factor(t); // get effective contact rate between i and j
            double divN    = 1.0 / params.populations.get_group_total(SecirCategory::AgeGroup, j); // precompute 1.0/Nj
            double dummy_S = y[Si] * cont_freq_eff * divN * params.probabilities[i].get_infection_from_contact() *
                             (y[Cj] + params.probabilities[j].get_risk_from_symptomatic() * y[Ij]);

            dydt[Si] -= dummy_S; // -R1*(C+beta*I)*S/N0
            dydt[Ei] += dummy_S; // R1*(C+beta*I)*S/N0-R2*E
        }

        double dummy_R2 = 1.0 / (2 * (params.times[i].get_serialinterval()) -
                                 (params.times[i].get_incubation())); // R2 = 1/(2SI-TINC)
        double dummy_R3 =
            0.5 / ((params.times[i].get_incubation()) - (params.times[i].get_serialinterval())); // R3 = 1/(2(TINC-SI))

        dydt[Ei] -= dummy_R2 * y[Ei]; // only exchange of E and C done here
        dydt[Ci] = dummy_R2 * y[Ei] -
                   ((1 - params.probabilities[i].get_asymp_per_infectious()) * dummy_R3 +
                    params.probabilities[i].get_asymp_per_infectious() / params.times[i].get_infectious_asymp()) *
                       y[Ci];
        dydt[Ii] =
            (1 - params.probabilities[i].get_asymp_per_infectious()) * dummy_R3 * y[Ci] -
            ((1 - params.probabilities[i].get_hospitalized_per_infectious()) / params.times[i].get_infectious_mild() +
             params.probabilities[i].get_hospitalized_per_infectious() / params.times[i].get_home_to_hospitalized()) *
                y[Ii];
        dydt[Hi] =
            params.probabilities[i].get_hospitalized_per_infectious() / params.times[i].get_home_to_hospitalized() *
                y[Ii] -
            ((1 - params.probabilities[i].get_icu_per_hospitalized()) / params.times[i].get_hospitalized_to_home() +
             params.probabilities[i].get_icu_per_hospitalized() / params.times[i].get_hospitalized_to_icu()) *
                y[Hi];
        dydt[Ui] =
            params.probabilities[i].get_icu_per_hospitalized() / params.times[i].get_hospitalized_to_icu() * y[Hi] -
            ((1 - params.probabilities[i].get_dead_per_icu()) / params.times[i].get_icu_to_home() +
             params.probabilities[i].get_dead_per_icu() / params.times[i].get_icu_to_dead()) *
                y[Ui];
        dydt[Ri] = params.probabilities[i].get_asymp_per_infectious() / params.times[i].get_infectious_asymp() * y[Ci] +
                   (1 - params.probabilities[i].get_hospitalized_per_infectious()) /
                       params.times[i].get_infectious_mild() * y[Ii] +
                   (1 - params.probabilities[i].get_icu_per_hospitalized()) /
                       params.times[i].get_hospitalized_to_home() * y[Hi] +
                   (1 - params.probabilities[i].get_dead_per_icu()) / params.times[i].get_icu_to_home() * y[Ui];
        dydt[Di] = params.probabilities[i].get_dead_per_icu() / params.times[i].get_icu_to_dead() * y[Ui];
    }
}

TimeSeries<double> simulate(double t0, double tmax, double dt, SecirParams const& params)
{
    SecirSimulation sim(params, t0, dt);
    sim.advance(tmax);
    return sim.get_result();
}

std::vector<double> simulate(double t0, double tmax, double dt, SecirParams const& params,
                             std::vector<Eigen::VectorXd>& secir)
{
    auto result = simulate(t0, tmax, dt, params);
    std::vector<double> t(result.get_times().begin(), result.get_times().end());
    std::transform(result.begin(), result.end(), std::back_inserter(secir), [](auto&& v_ref) { return v_ref.eval(); });

    return t;
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

Eigen::Ref<Eigen::VectorXd> SecirSimulation::advance(double tmax)
{
    return m_integrator.advance(tmax);
}

} // namespace epi
