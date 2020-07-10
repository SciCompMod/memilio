#include <epidemiology/secir.h>

// #include <epidemiology/euler.h>
#include <epidemiology/adapt_rk.h>
#include <epidemiology/eigen_util.h>
#include <epidemiology/save_result.h>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <algorithm>

#include <iostream>
#include <fstream>

namespace epi
{

void print_secir_params(ContactFrequencyMatrix const& cont_freq, std::vector<SecirParams> const& params)
{
    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs
    printf("\n SECIR model set. ");

    char output_file[35] = "output_secir_params_and_contacts";
    std::ofstream myfile(output_file);
    if (myfile.is_open()) {
        printf("Writing to file %s ", output_file);
        myfile << "Considering " << params.size() << " age groups.\n";
        for (size_t i = 0; i < params.size(); i++) {
            myfile << "\nParameters of Group " << i + 1 << "\n";
            myfile << "\t People at t0 \n";
            myfile << "\t\t Total: " << (int)params[i].populations.get_total_t0() << "\n";
            myfile << "\t\t Susceptible: " << (int)params[i].populations.get_suscetible_t0() << "\n";
            myfile << "\t\t Exposed: " << (int)params[i].populations.get_exposed_t0() << "\n";
            myfile << "\t\t Carrier: " << (int)params[i].populations.get_carrier_t0() << "\n";
            myfile << "\t\t Infectious: " << (int)params[i].populations.get_infectious_t0() << "\n";
            myfile << "\t\t Hospitalized: " << (int)params[i].populations.get_hospitalized_t0() << "\n";
            myfile << "\t\t ICU: " << (int)params[i].populations.get_icu_t0() << "\n";
            myfile << "\t\t Recovered: " << (int)params[i].populations.get_recovered_t0() << "\n";
            myfile << "\t\t Dead: " << (int)params[i].populations.get_dead_t0() << "\n";

            myfile << "\t Duration parameters \n";
            myfile << "\t\t Incubation time: \t" << 1.0 / params[i].times.get_incubation_inv() << "\n";
            myfile << "\t\t Infectious (mild) time: \t" << 1.0 / params[i].times.get_infectious_mild_inv() << "\n";
            myfile << "\t\t\t Serial Interval: \t" << 1.0 / params[i].times.get_serialinterval_inv() << "\n";
            myfile << "\t\t Hospitalized->Home time: \t" << 1.0 / params[i].times.get_hospitalized_to_home_inv()
                   << "\n";
            myfile << "\t\t Home->Hospitalized time: \t" << 1.0 / params[i].times.get_home_to_hospitalized_inv()
                   << "\n";
            myfile << "\t\t Infectious (asymp.) time: \t" << 1.0 / params[i].times.get_infectious_asymp_inv() << "\n";
            myfile << "\t\t Hospitalized->ICU time: \t" << 1.0 / params[i].times.get_hospitalized_to_icu_inv() << "\n";
            myfile << "\t\t ICU->Recovered time: \t" << 1.0 / params[i].times.get_icu_to_home_inv() << "\n";
            myfile << "\t\t ICU->Death time: \t" << 1.0 / params[i].times.get_icu_to_dead_inv() << "\n";

            myfile << "\t Probabilities \n";
            myfile << "\t\t Infect from contact: \t" << params[i].probabilities.get_infection_from_contact() << "\n";
            myfile << "\t\t Asymptomatic infections: \t" << params[i].probabilities.get_asymp_per_infectious() << "\n";
            myfile << "\t\t Risk of symptomatic contact: \t" << params[i].probabilities.get_risk_from_symptomatic()
                   << "\n";
            myfile << "\t\t Deaths per ICU care: \t" << params[i].probabilities.get_dead_per_icu() << "\n";
            myfile << "\t\t Hospitalized per Infection: \t" << params[i].probabilities.get_hospitalized_per_infectious()
                   << "\n";
            myfile << "\t\t ICU per Hospitalized: \t" << params[i].probabilities.get_icu_per_hospitalized() << "\n";
        }

        myfile << "\nContact frequency matrix \n\t";
        for (size_t i = 0; i < params.size(); i++) {
            myfile << "\t\t G" << i;
        }
        for (size_t i = 0; i < params.size(); i++) {
            myfile << "\n\t\t G" << i;
            for (size_t j = 0; j < params.size(); j++) {
                myfile << "\t" << cont_freq.get_cont_freq(static_cast<int>(i), static_cast<int>(j));
            }
        }

        myfile << "\n Dampings: \n\t";
        for (size_t i = 0; i < params.size(); i++) {
            myfile << "\n\t G" << i;
            for (size_t j = 0; j < params.size(); j++) {
                myfile << "\n\t\t G" << j;
                for (size_t k = 0;
                     k < cont_freq.get_dampings(static_cast<int>(i), static_cast<int>(j)).get_dampings_vector().size();
                     k++) {
                    myfile << "\t day: "
                           << cont_freq.get_dampings(static_cast<int>(i), static_cast<int>(j))
                                  .get_dampings_vector()
                                  .at(k)
                                  .day
                           << " fac: "
                           << cont_freq.get_dampings(static_cast<int>(i), static_cast<int>(j))
                                  .get_dampings_vector()
                                  .at(k)
                                  .factor;
                }
            }
        }

        myfile.close();
    }
}

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

double SecirParams::StageTimes::get_incubation_inv() const
{
    return m_tinc_inv;
}

double SecirParams::StageTimes::get_infectious_mild_inv() const
{
    return m_tinfmild_inv;
}

double SecirParams::StageTimes::get_serialinterval_inv() const
{
    return m_tserint_inv;
}

double SecirParams::StageTimes::get_hospitalized_to_home_inv() const
{
    return m_thosp2home_inv;
}

double SecirParams::StageTimes::get_home_to_hospitalized_inv() const
{
    return m_thome2hosp_inv;
}

double SecirParams::StageTimes::get_hospitalized_to_icu_inv() const
{
    return m_thosp2icu_inv;
}

double SecirParams::StageTimes::get_icu_to_home_inv() const
{
    return m_ticu2home_inv;
}

double SecirParams::StageTimes::get_infectious_asymp_inv() const
{
    return m_tinfasy_inv;
}

double SecirParams::StageTimes::get_icu_to_dead_inv() const
{
    return m_ticu2death_inv;
}

SecirParams::Populations::Populations()
    : m_nb_total_t0{0}
    , m_nb_exp_t0{0}
    , m_nb_car_t0{0}
    , m_nb_inf_t0{0}
    , m_nb_hosp_t0{0}
    , m_nb_icu_t0{0}
    , m_nb_rec_t0{0}
    , m_nb_dead_t0{0}
    , m_nb_sus_t0{0}
{
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

SecirParams::Probabilities::Probabilities()
    : m_infprob{1}
    , m_alpha{0}
    , m_beta{0}
    , m_rho{0}
    , m_theta{0}
    , m_delta{0}
{
}

void SecirParams::Probabilities::set_infection_from_contact(double const& infprob)
{
    m_infprob = infprob;
}

void SecirParams::Probabilities::set_asymp_per_infectious(double const& alpha)
{
    m_alpha = alpha;
}

void SecirParams::Probabilities::set_risk_from_symptomatic(double const& beta)
{
    m_beta = beta;
}

void SecirParams::Probabilities::set_hospitalized_per_infectious(double const& rho)
{
    m_rho = rho;
}

void SecirParams::Probabilities::set_icu_per_hospitalized(double const& theta)
{
    m_theta = theta;
}

void SecirParams::Probabilities::set_dead_per_icu(double const& delta)
{
    m_delta = delta;
}

double SecirParams::Probabilities::get_infection_from_contact() const
{
    return m_infprob;
}

double SecirParams::Probabilities::get_asymp_per_infectious() const
{
    return m_alpha;
}

double SecirParams::Probabilities::get_risk_from_symptomatic() const
{
    return m_beta;
}

double SecirParams::Probabilities::get_hospitalized_per_infectious() const
{
    return m_rho;
}

double SecirParams::Probabilities::get_icu_per_hospitalized() const
{
    return m_theta;
}

double SecirParams::Probabilities::get_dead_per_icu() const
{
    return m_delta;
}

ContactFrequencyMatrix::ContactFrequencyMatrix()
    : m_cont_freq{{1.0}}
    , m_dampings{{Dampings{}}}
{
}

ContactFrequencyMatrix::ContactFrequencyMatrix(size_t const nb_groups)
    : m_cont_freq{nb_groups, std::vector<double>(nb_groups, 0)}
    , m_dampings{nb_groups, std::vector<Dampings>(nb_groups, Dampings{})}
{
}

int ContactFrequencyMatrix::get_size() const
{
    return static_cast<int>(m_cont_freq.size());
}

void ContactFrequencyMatrix::set_cont_freq(double const cont_freq, int const self_group, int const contact_group)
{
    if (self_group <= contact_group) {
        m_cont_freq[self_group][contact_group] = cont_freq;
    }
    else {
        m_cont_freq[contact_group][self_group] = cont_freq;
    }
}

double ContactFrequencyMatrix::get_cont_freq(int self_group, int contact_group) const
{
    // prevent erroneous nonsymmetry
    return self_group <= contact_group ? m_cont_freq[self_group][contact_group]
                                       : m_cont_freq[contact_group][self_group];
}

void ContactFrequencyMatrix::set_dampings(Dampings const& damping, int self_group, int contact_group)
{
    if (self_group <= contact_group) {
        m_dampings[self_group][contact_group] = damping;
    }
    else {
        m_dampings[contact_group][self_group] = damping;
    }
}

const Dampings& ContactFrequencyMatrix::get_dampings(int self_group, int contact_group) const
{
    // prevent erroneous nonsymmetry
    return self_group <= contact_group ? m_dampings[self_group][contact_group] : m_dampings[contact_group][self_group];
}

void ContactFrequencyMatrix::add_damping(Damping const& damping, int self_group, int contact_group)
{
    if (self_group <= contact_group) {
        m_dampings[self_group][contact_group].add(damping);
    }
    else {
        m_dampings[contact_group][self_group].add(damping);
    }
}

double get_reprod_rate(ContactFrequencyMatrix const& cont_freq_matrix, std::vector<SecirParams> const& params,
                       double const t, std::vector<double> const& yt)
{
    if (params.size() == 1) {
        // (base_)reprod has to be computed time dependently !
        auto dummy_R3 =
            0.5 / (1.0 / params[0].times.get_infectious_mild_inv() - 1.0 / params[0].times.get_serialinterval_inv());

        double cont_freq_eff = cont_freq_matrix.get_cont_freq(0, 0) * cont_freq_matrix.get_dampings(0, 0).get_factor(t);

        double nb_total = 0;
        for (size_t i = 0; i < yt.size(); i++) {
            nb_total += yt[0];
        }

        double reprod_rate =
            cont_freq_eff *
            ((1 - params[0].probabilities.get_hospitalized_per_infectious()) *
                 params[0].times.get_infectious_mild_inv() +
             dummy_R3 * params[0].probabilities.get_risk_from_symptomatic() *
                 (1 - params[0].probabilities.get_asymp_per_infectious()) +
             params[0].probabilities.get_hospitalized_per_infectious() *
                 params[0].times.get_home_to_hospitalized_inv()) /
            ((dummy_R3 * (1 - params[0].probabilities.get_asymp_per_infectious()) +
              params[0].probabilities.get_asymp_per_infectious() * params[0].times.get_infectious_asymp_inv()) *
             (params[0].times.get_infectious_mild_inv() *
                  (1 - params[0].probabilities.get_hospitalized_per_infectious()) +
              params[0].probabilities.get_hospitalized_per_infectious() *
                  params[0].times.get_home_to_hospitalized_inv())) *
            yt[0] / nb_total;

        return reprod_rate;
    }
    else {
        return -1;
    }
}

void secir_get_derivatives(ContactFrequencyMatrix const& cont_freq_matrix, std::vector<SecirParams> const& params,
                           const Eigen::VectorXd& y, double t, Eigen::VectorXd& dydt)
{
    // alpha  // percentage of asymptomatic cases
    // beta // risk of infection from the infected symptomatic patients
    // rho   // hospitalized per infectious
    // theta // icu per hospitalized
    // delta  // deaths per ICUs
    // 0: S,      1: E,     2: C,     3: I,     4: H,     5: U,     6: R,     7: D
    size_t n_agegroups = params.size();

    for (int i = 0; i < n_agegroups; i++) {
        dydt[0 + 8 * i] = 0;
        dydt[1 + 8 * i] = 0;
        for (int j = 0; j < n_agegroups; j++) {
            // effective contact rate by contact rate between groups i and j and damping j
            double cont_freq_eff =
                cont_freq_matrix.get_cont_freq(i, j) *
                cont_freq_matrix.get_dampings(i, j).get_factor(t); // get effective contact rate between i and j
            double divN    = 1.0 / params[j].populations.get_total_t0(); // precompute 1.0/Nj
            double dummy_S = y[0 + 8 * i] * cont_freq_eff * divN *
                             params[i].probabilities.get_infection_from_contact() *
                             (y[2 + 8 * j] + params[j].probabilities.get_risk_from_symptomatic() * y[3 + 8 * j]);

            dydt[0 + 8 * i] -= dummy_S; // -R1*(C+beta*I)*S/N0
            dydt[1 + 8 * i] += dummy_S; // R1*(C+beta*I)*S/N0-R2*E
        }

        double dummy_R2 = 1.0 / (2 * (1.0 / params[i].times.get_serialinterval_inv()) -
                                 (1.0 / params[i].times.get_incubation_inv())); // R2 = 1/(2SI-TINC)
        double dummy_R3 = 0.5 / ((1.0 / params[i].times.get_incubation_inv()) -
                                 (1.0 / params[i].times.get_serialinterval_inv())); // R3 = 1/(2(TINC-SI))

        dydt[1 + 8 * i] -= dummy_R2 * y[1 + 8 * i]; // only exchange of E and C done here
        dydt[2 + 8 * i] =
            dummy_R2 * y[1 + 8 * i] -
            ((1 - params[i].probabilities.get_asymp_per_infectious()) * dummy_R3 +
             params[i].probabilities.get_asymp_per_infectious() * params[i].times.get_infectious_asymp_inv()) *
                y[2 + 8 * i];
        dydt[3 + 8 * i] = (1 - params[i].probabilities.get_asymp_per_infectious()) * dummy_R3 * y[2 + 8 * i] -
                          ((1 - params[i].probabilities.get_hospitalized_per_infectious()) *
                               params[i].times.get_infectious_mild_inv() +
                           params[i].probabilities.get_hospitalized_per_infectious() *
                               params[i].times.get_home_to_hospitalized_inv()) *
                              y[3 + 8 * i];
        dydt[4 + 8 * i] =
            params[i].probabilities.get_hospitalized_per_infectious() * params[i].times.get_home_to_hospitalized_inv() *
                y[3 + 8 * i] -
            ((1 - params[i].probabilities.get_icu_per_hospitalized()) * params[i].times.get_hospitalized_to_home_inv() +
             params[i].probabilities.get_icu_per_hospitalized() * params[i].times.get_hospitalized_to_icu_inv()) *
                y[4 + 8 * i];
        dydt[5 + 8 * i] = params[i].probabilities.get_icu_per_hospitalized() *
                              params[i].times.get_hospitalized_to_icu_inv() * y[4 + 8 * i] -
                          ((1 - params[i].probabilities.get_dead_per_icu()) * params[i].times.get_icu_to_home_inv() +
                           params[i].probabilities.get_dead_per_icu() * params[i].times.get_icu_to_dead_inv()) *
                              y[5 + 8 * i];
        dydt[6 + 8 * i] =
            params[i].probabilities.get_asymp_per_infectious() * params[i].times.get_infectious_asymp_inv() *
                y[2 + 8 * i] +
            (1 - params[i].probabilities.get_hospitalized_per_infectious()) *
                params[i].times.get_infectious_mild_inv() * y[3 + 8 * i] +
            (1 - params[i].probabilities.get_icu_per_hospitalized()) * params[i].times.get_hospitalized_to_home_inv() *
                y[4 + 8 * i] +
            (1 - params[i].probabilities.get_dead_per_icu()) * params[i].times.get_icu_to_home_inv() * y[5 + 8 * i];
        dydt[7 + 8 * i] =
            params[i].probabilities.get_dead_per_icu() * params[i].times.get_icu_to_dead_inv() * y[5 + 8 * i];
    }
}

namespace
{
    Eigen::VectorXd secir_get_initial_values(const SecirParams& params)
    {
        Eigen::VectorXd y(8);
        y[0] = params.populations.get_suscetible_t0();
        y[1] = params.populations.get_exposed_t0();
        y[2] = params.populations.get_carrier_t0();
        y[3] = params.populations.get_infectious_t0();
        y[4] = params.populations.get_hospitalized_t0();
        y[5] = params.populations.get_icu_t0();
        y[6] = params.populations.get_recovered_t0();
        y[7] = params.populations.get_dead_t0();
        return y;
    }

    Eigen::VectorXd secir_get_initial_values(const std::vector<SecirParams>& params)
    {
        Eigen::VectorXd y(params.size() * 8);
        for (size_t i = 0; i < params.size(); i++) {
            slice(y, Seq<Eigen::Index>{static_cast<Eigen::Index>(i) * 8, 8}) = secir_get_initial_values(params[i]);
        }
        return y;
    }
} // namespace

std::vector<double> simulate(double t0, double tmax, double dt, ContactFrequencyMatrix const& cont_freq_matrix,
                             std::vector<SecirParams> const& params, std::vector<Eigen::VectorXd>& secir)
{
    SecirSimulation sim(cont_freq_matrix, params, t0, dt);
    sim.advance(tmax);
    secir = sim.get_y();

    return sim.get_t();
}

SecirSimulation::SecirSimulation(const ContactFrequencyMatrix& cont_freq_matrix, const std::vector<SecirParams>& params,
                                 double t0, double dt_init)
    : m_integratorCore(std::make_shared<RKIntegratorCore>(1e-3, 1.))
    , m_integrator(
          [params, cont_freq_matrix](auto&& y, auto&& t, auto&& dydt) {
              secir_get_derivatives(cont_freq_matrix, params, y, t, dydt);
          },
          t0, secir_get_initial_values(params), dt_init, m_integratorCore)
{
    m_integratorCore->set_rel_tolerance(1e-4);
    m_integratorCore->set_abs_tolerance(1e-1);
    //print_secir_params(params, cont_freq_matrix);
}

Eigen::VectorXd& SecirSimulation::advance(double tmax)
{
    return m_integrator.advance(tmax);
}

} // namespace epi
