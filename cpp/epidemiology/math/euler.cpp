#include "epidemiology/math/euler.h"
#include "epidemiology/secir/secir.h"

namespace epi
{

bool EulerIntegratorCore::step(const DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, double& t, double& dt,
                               Eigen::Ref<Eigen::VectorXd> ytp1) const
{
    // we are misusing the next step y as temporary space to store the derivative
    f(yt, t, ytp1);
    ytp1 = yt + dt * ytp1;
    t += dt;
    return true;
}

ImplicitEulerIntegratorCore::ImplicitEulerIntegratorCore(SecirModel const& model)
    : m_model{model}
{
}

bool ImplicitEulerIntegratorCore::step(const DerivFunction& /*f*/, Eigen::Ref<const Eigen::VectorXd> yt, double& t,
                                       double& dt, Eigen::Ref<Eigen::VectorXd> ytp1) const
{
    // 0: S,      1: E,     2: C,     3: I,     4: H,     5: U,     6: R,     7: D

    const auto& params         = m_model.parameters;
    const auto& contact_matrix = params.get_contact_patterns().get_cont_freq_mat();

    auto yt_eval = yt.eval();

    Eigen::VectorXd yt_tilde;
    yt_tilde.resizeLike(yt);

    // go through the variables of the system: S, E, I, ....
    double divN = 1.0 / m_model.populations.get_total(); // precompute 1.0/N

    double cont_freq_eff = contact_matrix.get_matrix_at(t)(0, 0);

    double dummy_R2 =
        1.0 / (2 * (params.times[0].get_serialinterval()) - (params.times[0].get_incubation())); // R2 = 1/(2SI-TINC)

    double dummy_R3 =
        0.5 / ((params.times[0].get_incubation()) - (params.times[0].get_serialinterval())); // R3 = 1/(2(TINC-SI))

    double dummy_C =
        1. / (1. + dt * ((1. - params.probabilities[0].get_asymp_per_infectious()) * dummy_R3 +
                         params.probabilities[0].get_asymp_per_infectious() / params.times[0].get_infectious_asymp()));

    double dummy_I =
        1. /
        (1. +
         dt * ((1 - params.probabilities[0].get_hospitalized_per_infectious()) / params.times[0].get_infectious_mild() +
               params.probabilities[0].get_hospitalized_per_infectious() / params.times[0].get_home_to_hospitalized()));

    double dummy_H =
        1. /
        (1. +
         dt * ((1 - params.probabilities[0].get_icu_per_hospitalized()) / params.times[0].get_hospitalized_to_home() +
               params.probabilities[0].get_icu_per_hospitalized() / params.times[0].get_hospitalized_to_icu()));

    double dummy_U =
        1. / (1. + dt * ((1 - params.probabilities[0].get_dead_per_icu()) / params.times[0].get_icu_to_home() +
                         params.probabilities[0].get_dead_per_icu() / params.times[0].get_icu_to_dead()));

    // these temporary variables are used for the fix-point iteration
    double y_temp_C;
    double y_temp_I;
    yt_tilde[(size_t)InfectionState::Carrier] = yt_eval[(size_t)InfectionState::Carrier];
    yt_tilde[(size_t)InfectionState::Infected] = yt_eval[(size_t)InfectionState::Infected];

    // fix-point iteration since C and I are discretized explicitly (C^n, I^n), using this approach we obtain a close approximation to (C^{n+1}, I^{n+1})

    do {
        // S
        y_temp_C = yt_tilde[(size_t)InfectionState::Carrier];
        y_temp_I = yt_tilde[(size_t)InfectionState::Infected];

        double dummy_S = cont_freq_eff * divN * params.probabilities[0].get_infection_from_contact() *
                         (y_temp_C + params.probabilities[0].get_risk_from_symptomatic() * y_temp_I);

        yt_tilde[(size_t)InfectionState::Susceptible] = yt_eval[(size_t)InfectionState::Susceptible] / (1. + dt * dummy_S);

        // E (use new value for S)
        yt_tilde[(size_t)InfectionState::Exposed] =
            (yt_eval[(size_t)InfectionState::Exposed] + dt * dummy_S * yt_tilde[(size_t)InfectionState::Susceptible]) /
            (1. + dt * dummy_R2);

        // C (use new value for E)
        yt_tilde[(size_t)InfectionState::Carrier] =
            dummy_C * (yt_eval[(size_t)InfectionState::Carrier] + dt * dummy_R2 * yt_tilde[(size_t)InfectionState::Exposed]);

        // I (use new value for C)
        yt_tilde[(size_t)InfectionState::Infected] = dummy_I * (yt_eval[(size_t)InfectionState::Infected] +
                                                        dt * (1 - params.probabilities[0].get_asymp_per_infectious()) *
                                                            dummy_R3 * yt_tilde[(size_t)InfectionState::Carrier]);

    } while (std::fabs(yt_tilde[(size_t)InfectionState::Carrier] - y_temp_C) > m_abs_tol &&
             std::fabs(yt_tilde[(size_t)InfectionState::Infected] - y_temp_I) > m_abs_tol);

    // H (use new value for I, i.e. I^{n+1})
    yt_tilde[(size_t)InfectionState::Hospitalized] =
        dummy_H * (yt_eval[(size_t)InfectionState::Hospitalized] + dt * params.probabilities[0].get_hospitalized_per_infectious() /
                                                           params.times[0].get_home_to_hospitalized() *
                                                           yt_tilde[(size_t)InfectionState::Infected]);

    // U (use new value for H, i.e. H^{n+1})
    yt_tilde[(size_t)InfectionState::ICU] =
        dummy_U * (yt_eval[(size_t)InfectionState::ICU] + dt * params.probabilities[0].get_icu_per_hospitalized() /
                                                           params.times[0].get_hospitalized_to_icu() *
                                                           yt_tilde[(size_t)InfectionState::Hospitalized]);

    // R (use new values for C, I, H, U)
    yt_tilde[(size_t)InfectionState::Recovered] =
        yt_eval[(size_t)InfectionState::Recovered] +
        dt * (params.probabilities[0].get_asymp_per_infectious() / params.times[0].get_infectious_asymp() *
                  yt_tilde[(size_t)InfectionState::Carrier] +
              (1 - params.probabilities[0].get_hospitalized_per_infectious()) / params.times[0].get_infectious_mild() *
                  yt_tilde[(size_t)InfectionState::Infected] +
              (1 - params.probabilities[0].get_icu_per_hospitalized()) / params.times[0].get_hospitalized_to_home() *
                  yt_tilde[(size_t)InfectionState::Hospitalized] +
              (1 - params.probabilities[0].get_dead_per_icu()) / params.times[0].get_icu_to_home() *
                  yt_tilde[(size_t)InfectionState::ICU]);

    // D (use new value for U, i.e. U^{n+1})
    yt_tilde[(size_t)InfectionState::Dead] =
        yt_eval[(size_t)InfectionState::Dead] + dt * params.probabilities[0].get_dead_per_icu() /
                                                params.times[0].get_icu_to_dead() * yt_tilde[(size_t)InfectionState::ICU];

    ytp1 = yt_tilde;
    t += dt; // this is the new timestep t=(n+1)*\Delta t where ytp1 belongs to

    return true;
}

} // namespace epi
