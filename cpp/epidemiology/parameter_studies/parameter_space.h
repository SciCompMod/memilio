#ifndef PARAMETER_SPACE_H
#define PARAMETER_SPACE_H

#include <assert.h>
#include <string>
#include <vector>
#include <random>
#include <memory>
#include <epidemiology/memory.h>
#include <epidemiology/secir.h>
#include <epidemiology/logging.h>
#include <epidemiology/parameter_studies/parameter_distributions.h>

namespace epi
{

/* The class ParameterSpace stores ranges of parameters
 * together with information on step sizes,
 * a start and end time as well as an initial time step.
 * The class provides an iterator that iterates over all
 * generated parameter combinations.
 *
 * Currently all parameters are of type double.
 */
class ParameterSpace
{
public:
    /* Constructor from given ContactFrequencyMatrix and SecirParams. Mainly used for testing.
     * @param[in] params SecirParams including ContactFrequencyMatrix for alle age groups
     * @param[in] t0 start time
     * @param[in] tmax end time
     * @param[in] dev_rel maximum relative deviation from particular value(s) given in params
     */
    ParameterSpace(SecirParams& params, double t0, double tmax, double dev_rel);

    const std::vector<double>& get_total() const
    {
        return m_total;
    }

    const std::vector<double>& get_hospitalized() const
    {
        return m_hospitalized;
    }

    const std::vector<double>& get_icu() const
    {
        return m_icu;
    }

    const std::vector<double>& get_dead() const
    {
        return m_dead;
    }

    const ParameterDistribution& get_dist_exposed(int group) const
    {
        return *m_exposed[group];
    }

    const ParameterDistribution& get_dist_infectious(int group) const
    {
        return *m_infectious[group];
    }

    const ParameterDistribution& get_dist_carrier(int group) const
    {
        return *m_carrier[group];
    }

    const ParameterDistribution& get_dist_recovered(int group) const
    {
        return *m_recovered[group];
    }

    const ParameterDistribution* get_dist_incubation(int group) const
    {
        return m_params.times[group].get_incubation().get_distribution().get();
    }

    const ParameterDistribution* get_dist_serial_int(int group) const
    {
        return m_params.times[group].get_serialinterval().get_distribution().get();
    }

    const ParameterDistribution* get_dist_inf_mild(int group) const
    {
        return m_params.times[group].get_infectious_mild().get_distribution().get();
    }

    const ParameterDistribution* get_dist_hosp_to_rec(int group) const
    {
        return m_params.times[group].get_hospitalized_to_home().get_distribution().get();
    }

    const ParameterDistribution* get_dist_inf_to_hosp(int group) const
    {
        return m_params.times[group].get_home_to_hospitalized().get_distribution().get();
    }

    const ParameterDistribution* get_dist_inf_asymp(int group) const
    {
        return m_params.times[group].get_infectious_asymp().get_distribution().get();
    }

    const ParameterDistribution* get_dist_hosp_to_icu(int group) const
    {
        return m_params.times[group].get_hospitalized_to_icu().get_distribution().get();
    }

    const ParameterDistribution* get_dist_icu_to_rec(int group) const
    {
        return m_params.times[group].get_icu_to_home().get_distribution().get();
    }

    const ParameterDistribution* get_dist_icu_to_death(int group) const
    {
        return m_params.times[group].get_icu_to_dead().get_distribution().get();
    }

    const ParameterDistribution* get_dist_inf_from_cont(int group) const
    {
        return m_params.probabilities[group].get_infection_from_contact().get_distribution().get();
    }

    const ParameterDistribution* get_dist_asymp_per_inf(int group) const
    {
        return m_params.probabilities[group].get_asymp_per_infectious().get_distribution().get();
    }

    const ParameterDistribution* get_dist_risk_from_symp(int group) const
    {
        return m_params.probabilities[group].get_risk_from_symptomatic().get_distribution().get();
    }

    const ParameterDistribution* get_dist_death_per_icu(int group) const
    {
        return m_params.probabilities[group].get_dead_per_icu().get_distribution().get();
    }

    const ParameterDistribution* get_dist_hosp_per_inf(int group) const
    {
        return m_params.probabilities[group].get_hospitalized_per_infectious().get_distribution().get();
    }

    const ParameterDistribution* get_dist_icu_per_hosp(int group) const
    {
        return m_params.probabilities[group].get_icu_per_hospitalized().get_distribution().get();
    }

    SecirParams get_secir_params_sample()
    {
        SecirParams secir_params_sample(m_params.size());
        for (size_t i = 0; i < m_params.size(); i++) {

            secir_params_sample.populations.set({i, SecirCompartments::H}, m_hospitalized[i]);
            secir_params_sample.populations.set({i, SecirCompartments::U}, m_icu[i]);
            secir_params_sample.populations.set({i, SecirCompartments::D}, m_dead[i]);

            secir_params_sample.populations.set({i, SecirCompartments::E}, m_exposed[i]->get_sample());
            secir_params_sample.populations.set({i, SecirCompartments::C}, m_carrier[i]->get_sample());
            secir_params_sample.populations.set({i, SecirCompartments::I}, m_infectious[i]->get_sample());
            secir_params_sample.populations.set({i, SecirCompartments::R}, m_recovered[i]->get_sample());
            secir_params_sample.populations.set_difference_from_group_total(
                {i, epi::SecirCompartments::S}, epi::SecirCategory::AgeGroup, i, m_total[i]);

            secir_params_sample.times[i].set_incubation(m_params.times[i].get_incubation().draw_sample());
            secir_params_sample.times[i].set_infectious_mild(m_params.times[i].get_infectious_mild().draw_sample());
            secir_params_sample.times[i].set_serialinterval(m_params.times[i].get_serialinterval().draw_sample());
            secir_params_sample.times[i].set_hospitalized_to_home(
                m_params.times[i].get_hospitalized_to_home().draw_sample()); // here: home=recovered
            secir_params_sample.times[i].set_home_to_hospitalized(
                m_params.times[i].get_home_to_hospitalized().draw_sample()); // here: home=infectious
            secir_params_sample.times[i].set_infectious_asymp(m_params.times[i].get_infectious_asymp().draw_sample());
            secir_params_sample.times[i].set_hospitalized_to_icu(
                m_params.times[i].get_hospitalized_to_icu().draw_sample());
            secir_params_sample.times[i].set_icu_to_death(m_params.times[i].get_icu_to_dead().draw_sample());
            secir_params_sample.times[i].set_icu_to_home(m_params.times[i].get_icu_to_home().draw_sample());

            secir_params_sample.probabilities[i].set_infection_from_contact(
                m_params.probabilities[i].get_infection_from_contact().draw_sample());
            secir_params_sample.probabilities[i].set_asymp_per_infectious(
                m_params.probabilities[i].get_asymp_per_infectious().draw_sample());
            secir_params_sample.probabilities[i].set_risk_from_symptomatic(
                m_params.probabilities[i].get_risk_from_symptomatic().draw_sample());
            secir_params_sample.probabilities[i].set_dead_per_icu(
                m_params.probabilities[i].get_dead_per_icu().draw_sample());
            secir_params_sample.probabilities[i].set_hospitalized_per_infectious(
                m_params.probabilities[i].get_hospitalized_per_infectious().draw_sample());
            secir_params_sample.probabilities[i].set_icu_per_hospitalized(
                m_params.probabilities[i].get_icu_per_hospitalized().draw_sample());
        }

        return secir_params_sample;
    }

private:
    SecirParams m_params;

    // populations
    std::vector<double> m_total;
    std::vector<std::unique_ptr<ParameterDistribution>> m_exposed;
    std::vector<std::unique_ptr<ParameterDistribution>> m_carrier;
    std::vector<std::unique_ptr<ParameterDistribution>> m_infectious;
    std::vector<double> m_hospitalized;
    std::vector<double> m_icu;
    std::vector<std::unique_ptr<ParameterDistribution>> m_recovered;
    std::vector<double> m_dead;
};

ParameterSpace::ParameterSpace(SecirParams& params, double t0, double tmax, double dev_rel)
    : m_params{params}
{
    double min_val = 0.001;

    // populations
    for (size_t i = 0; i < params.size(); i++) {

        // fixed size groups
        // total
        m_total.push_back(params.populations.get_group_total(SecirCategory::AgeGroup, i));
        m_hospitalized.push_back(params.populations.get({i, SecirCompartments::H}));
        m_icu.push_back(params.populations.get({i, SecirCompartments::U}));
        m_dead.push_back(params.populations.get({i, SecirCompartments::D}));

        // variably sized groups
        // exposed
        double value_params = params.populations.get({i, SecirCompartments::E});
        m_exposed.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // carrier
        value_params = params.populations.get({i, SecirCompartments::C});
        m_carrier.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // infectious
        value_params = params.populations.get({i, SecirCompartments::I});
        m_infectious.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));

        // recovered
        value_params = params.populations.get({i, SecirCompartments::R});
        m_recovered.push_back(std::make_unique<ParameterDistributionNormal>(
            std::max(min_val, (1 - dev_rel * 2.6) * value_params), (1 + dev_rel * 2.6) * value_params, value_params,
            dev_rel * value_params));
    }

    // times
    for (size_t i = 0; i < params.size(); i++) {
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
    for (size_t i = 0; i < params.size(); i++) {
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

#endif // PARAMETER_SPACE_H
