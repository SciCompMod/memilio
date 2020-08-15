#include <epidemiology/parameter_studies/parameter_space.h>
#include <epidemiology/secir.h>
#include <epidemiology/parameter_studies/parameter_distributions.h>

namespace epi
{

ParameterSpace::ParameterSpace(SecirParams const& params)
    : m_params{params}
{
}

ParameterSpace::ParameterSpace(SecirParams const& params, double t0, double tmax, double dev_rel)
    : m_params{params}
{
    // ParameterSpace param_space{params};

    double min_val = 0.001;

    // populations
    for (size_t i = 0; i < m_params.size(); i++) {

        // variably sized groups
        // exposed
        double value_params = m_params.populations.get({i, SecirCompartments::E});
        m_params.populations.set({i, SecirCompartments::E},
                                 ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                                             (1 + dev_rel * 2.6) * value_params, value_params,
                                                             dev_rel * value_params));

        // carrier
        value_params = m_params.populations.get({i, SecirCompartments::C});
        m_params.populations.set({i, SecirCompartments::C},
                                 ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                                             (1 + dev_rel * 2.6) * value_params, value_params,
                                                             dev_rel * value_params));

        // infectious
        value_params = m_params.populations.get({i, SecirCompartments::I});
        m_params.populations.set({i, SecirCompartments::I},
                                 ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                                             (1 + dev_rel * 2.6) * value_params, value_params,
                                                             dev_rel * value_params));

        // hospitalized
        value_params = m_params.populations.get({i, SecirCompartments::H});
        m_params.populations.set({i, SecirCompartments::H},
                                 ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                                             (1 + dev_rel * 2.6) * value_params, value_params,
                                                             dev_rel * value_params));

        // icu
        value_params = m_params.populations.get({i, SecirCompartments::U});
        m_params.populations.set({i, SecirCompartments::U},
                                 ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                                             (1 + dev_rel * 2.6) * value_params, value_params,
                                                             dev_rel * value_params));

        // recovered
        value_params = m_params.populations.get({i, SecirCompartments::R});
        m_params.populations.set({i, SecirCompartments::R},
                                 ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                                             (1 + dev_rel * 2.6) * value_params, value_params,
                                                             dev_rel * value_params));
    }

    // times
    for (size_t i = 0; i < m_params.size(); i++) {
        // incubation time
        double value_params = m_params.times[i].get_incubation();
        m_params.times[i].set_incubation(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // serial interval
        value_params = m_params.times[i].get_serialinterval();
        m_params.times[i].set_serialinterval(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // infectious time (mild)
        value_params = m_params.times[i].get_infectious_mild();
        m_params.times[i].set_infectious_mild(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // infectious to recovered (hospitalized to home)
        value_params = m_params.times[i].get_hospitalized_to_home();
        m_params.times[i].set_hospitalized_to_home(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // infectious to hospitalized (home to hospitalized)
        value_params = m_params.times[i].get_home_to_hospitalized();
        m_params.times[i].set_home_to_hospitalized(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // infectious (asymptomatic)
        value_params = m_params.times[i].get_infectious_asymp();
        m_params.times[i].set_infectious_asymp(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // hospitalized to ICU
        value_params = m_params.times[i].get_hospitalized_to_icu();
        m_params.times[i].set_hospitalized_to_icu(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // ICU to recovered
        value_params = m_params.times[i].get_icu_to_home();
        m_params.times[i].set_icu_to_home(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // ICU to death
        value_params = m_params.times[i].get_icu_to_dead();
        m_params.times[i].set_icu_to_death(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));
    }

    // probabilities
    for (size_t i = 0; i < m_params.size(); i++) {
        // infection from contact
        double value_params = m_params.probabilities[i].get_infection_from_contact();
        m_params.probabilities[i].set_infection_from_contact(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // asymptomatic per infectious
        value_params = m_params.probabilities[i].get_asymp_per_infectious();
        m_params.probabilities[i].set_asymp_per_infectious(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // risk of infection from infectious
        value_params = m_params.probabilities[i].get_risk_from_symptomatic();
        m_params.probabilities[i].set_risk_from_symptomatic(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // deaths per icu treatments
        value_params = m_params.probabilities[i].get_dead_per_icu();
        m_params.probabilities[i].set_dead_per_icu(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // hospitalized per infections
        value_params = m_params.probabilities[i].get_hospitalized_per_infectious();
        m_params.probabilities[i].set_hospitalized_per_infectious(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // icu treatments per hospitalized
        value_params = m_params.probabilities[i].get_icu_per_hospitalized();
        m_params.probabilities[i].set_icu_per_hospitalized(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));
    }

    // maximum number of dampings; to avoid overfitting only allow one damping for every 10 days simulated
    // damping base values are between 0.1 and 1; diagonal values vary lie in the range of 0.6 to 1.4 times the base value
    // off diagonal values vary between 0.7 to 1.1 of the corresponding diagonal value (symmetrization is conducted)
    m_params.get_contact_patterns().set_dist_damp_nb(ParameterDistributionUniform(1, (tmax - t0) / 10));
    m_params.get_contact_patterns().set_dist_damp_days(ParameterDistributionUniform(t0, tmax));
    m_params.get_contact_patterns().set_dist_damp_diag_base(ParameterDistributionUniform(0.1, 1));
    m_params.get_contact_patterns().set_dist_damp_diag_rel(ParameterDistributionUniform(0.6, 1.4));
    m_params.get_contact_patterns().set_dist_damp_offdiag_rel(ParameterDistributionUniform(0.7, 1.1));
}

} // namespace epi