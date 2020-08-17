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
    /* Constructor for ParameterSpace from SecirParams; does not introduce new distributions 
     * nor change current ones
     * @param[in] params SecirParams including ContactFrequencyMatrix for alle age groups
     */
    ParameterSpace(SecirParams const& params);

    /* Constructor for ParameterSpace from SecirParams; introduces normal distributions
     * @param[in] params SecirParams including ContactFrequencyMatrix for alle age groups
     * @param[in] t0 start time
     * @param[in] tmax end time
     * @param[in] dev_rel maximum relative deviation from particular value(s) given in params
     */
    ParameterSpace(SecirParams const& params, double t0, double tmax, double dev_rel);

    double get_total(size_t group) const
    {
        return m_params.populations.get_group_total(AgeGroup, group);
    }

    double get_dead(size_t group) const
    {
        return m_params.populations.get({group, SecirCompartments::D});
    }

    const ParameterDistribution* get_distribution_exposed(size_t group) const
    {
        return m_params.populations.get({group, SecirCompartments::E}).get_distribution().get();
    }

    const ParameterDistribution* get_distribution_carrier(size_t group) const
    {
        return m_params.populations.get({group, SecirCompartments::C}).get_distribution().get();
    }

    const ParameterDistribution* get_distribution_infectious(size_t group) const
    {
        return m_params.populations.get({group, SecirCompartments::I}).get_distribution().get();
    }

    const ParameterDistribution* get_distribution_hospitalized(size_t group) const
    {
        return m_params.populations.get({group, SecirCompartments::H}).get_distribution().get();
    }

    const ParameterDistribution* get_distribution_icu(size_t group) const
    {
        return m_params.populations.get({group, SecirCompartments::U}).get_distribution().get();
    }

    const ParameterDistribution* get_distribution_recovered(size_t group) const
    {
        return m_params.populations.get({group, SecirCompartments::R}).get_distribution().get();
    }

    const ParameterDistribution* get_distribution_incubation(int group) const
    {
        return m_params.times[group].get_incubation().get_distribution().get();
    }

    const ParameterDistribution* get_distribution_serial_int(int group) const
    {
        return m_params.times[group].get_serialinterval().get_distribution().get();
    }

    const ParameterDistribution* get_distribution_inf_mild(int group) const
    {
        return m_params.times[group].get_infectious_mild().get_distribution().get();
    }

    const ParameterDistribution* get_distribution_hosp_to_rec(int group) const
    {
        return m_params.times[group].get_hospitalized_to_home().get_distribution().get();
    }

    const ParameterDistribution* get_distribution_inf_to_hosp(int group) const
    {
        return m_params.times[group].get_home_to_hospitalized().get_distribution().get();
    }

    const ParameterDistribution* get_distribution_inf_asymp(int group) const
    {
        return m_params.times[group].get_infectious_asymp().get_distribution().get();
    }

    const ParameterDistribution* get_distribution_hosp_to_icu(int group) const
    {
        return m_params.times[group].get_hospitalized_to_icu().get_distribution().get();
    }

    const ParameterDistribution* get_distribution_icu_to_rec(int group) const
    {
        return m_params.times[group].get_icu_to_home().get_distribution().get();
    }

    const ParameterDistribution* get_distribution_icu_to_death(int group) const
    {
        return m_params.times[group].get_icu_to_dead().get_distribution().get();
    }

    const ParameterDistribution* get_distribution_inf_from_cont(int group) const
    {
        return m_params.probabilities[group].get_infection_from_contact().get_distribution().get();
    }

    const ParameterDistribution* get_distribution_asymp_per_inf(int group) const
    {
        return m_params.probabilities[group].get_asymp_per_infectious().get_distribution().get();
    }

    const ParameterDistribution* get_distribution_risk_from_symp(int group) const
    {
        return m_params.probabilities[group].get_risk_from_symptomatic().get_distribution().get();
    }

    const ParameterDistribution* get_distribution_death_per_icu(int group) const
    {
        return m_params.probabilities[group].get_dead_per_icu().get_distribution().get();
    }

    const ParameterDistribution* get_distribution_hosp_per_inf(int group) const
    {
        return m_params.probabilities[group].get_hospitalized_per_infectious().get_distribution().get();
    }

    const ParameterDistribution* get_distribution_icu_per_hosp(int group) const
    {
        return m_params.probabilities[group].get_icu_per_hospitalized().get_distribution().get();
    }

    SecirParams& get_secir_params()
    {
        return m_params;
    }

    SecirParams const& get_secir_params() const
    {
        return m_params;
    }

    SecirParams draw_sample()
    {
        SecirParams secir_params_sample(m_params.size());
        for (size_t i = 0; i < m_params.size(); i++) {

            double group_total = m_params.populations.get_group_total(SecirCategory::AgeGroup, i);

            secir_params_sample.populations.set({i, SecirCompartments::E},
                                                m_params.populations.get({i, SecirCompartments::E}).draw_sample());
            secir_params_sample.populations.set({i, SecirCompartments::C},
                                                m_params.populations.get({i, SecirCompartments::C}).draw_sample());
            secir_params_sample.populations.set({i, SecirCompartments::I},
                                                m_params.populations.get({i, SecirCompartments::I}).draw_sample());
            secir_params_sample.populations.set({i, SecirCompartments::H},
                                                m_params.populations.get({i, SecirCompartments::H}).draw_sample());
            secir_params_sample.populations.set({i, SecirCompartments::U},
                                                m_params.populations.get({i, SecirCompartments::U}).draw_sample());
            secir_params_sample.populations.set({i, SecirCompartments::R},
                                                m_params.populations.get({i, SecirCompartments::R}).draw_sample());

            // no sampling for dead and total numbers
            secir_params_sample.populations.set({i, SecirCompartments::D},
                                                m_params.populations.get({i, SecirCompartments::D}));

            m_params.populations.set_difference_from_group_total({i, epi::SecirCompartments::S},
                                                                 epi::SecirCategory::AgeGroup, i, group_total);
            secir_params_sample.populations.set_difference_from_group_total(
                {i, epi::SecirCompartments::S}, epi::SecirCategory::AgeGroup, i,
                m_params.populations.get_group_total(SecirCategory::AgeGroup, i));

            secir_params_sample.times[i].set_incubation(m_params.times[i].get_incubation().draw_sample());
            double serint_dummy = m_params.times[i].get_serialinterval().draw_sample();
            if (serint_dummy > secir_params_sample.times[i].get_incubation() - 0.2) {
                serint_dummy = secir_params_sample.times[i].get_incubation() - 1;
                log_warning("To Do/Discussion: Redesign sample strategy for serial interval.");
            }
            else if (2 * serint_dummy < secir_params_sample.times[i].get_incubation()) {
                serint_dummy += 0.5;
                log_warning("To Do/Discussion: Redesign sample strategy for serial interval.");
            }
            secir_params_sample.times[i].set_serialinterval(serint_dummy);
            secir_params_sample.times[i].set_infectious_mild(m_params.times[i].get_infectious_mild().draw_sample());
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
};

} // namespace epi

#endif // PARAMETER_SPACE_H
