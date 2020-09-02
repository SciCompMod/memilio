#include <epidemiology/secir.h>
#include <epidemiology/parameter_studies/parameter_space.h>
#include <epidemiology/parameter_studies/parameter_distributions.h>

namespace epi
{

void set_params_distributions_normal(SecirParams& params, double t0, double tmax, double dev_rel)
{
    double min_val = 0.001;

    // populations
    for (size_t i = 0; i < params.get_num_groups(); i++) {

        // variably sized groups
        // exposed
        double value_params = params.populations.get({i, SecirCompartments::E});
        params.populations.set({i, SecirCompartments::E},
                               ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                                           (1 + dev_rel * 2.6) * value_params, value_params,
                                                           dev_rel * value_params));

        // carrier
        value_params = params.populations.get({i, SecirCompartments::C});
        params.populations.set({i, SecirCompartments::C},
                               ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                                           (1 + dev_rel * 2.6) * value_params, value_params,
                                                           dev_rel * value_params));

        // infectious
        value_params = params.populations.get({i, SecirCompartments::I});
        params.populations.set({i, SecirCompartments::I},
                               ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                                           (1 + dev_rel * 2.6) * value_params, value_params,
                                                           dev_rel * value_params));

        // hospitalized
        value_params = params.populations.get({i, SecirCompartments::H});
        params.populations.set({i, SecirCompartments::H},
                               ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                                           (1 + dev_rel * 2.6) * value_params, value_params,
                                                           dev_rel * value_params));

        // icu
        value_params = params.populations.get({i, SecirCompartments::U});
        params.populations.set({i, SecirCompartments::U},
                               ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                                           (1 + dev_rel * 2.6) * value_params, value_params,
                                                           dev_rel * value_params));

        // recovered
        value_params = params.populations.get({i, SecirCompartments::R});
        params.populations.set({i, SecirCompartments::R},
                               ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                                           (1 + dev_rel * 2.6) * value_params, value_params,
                                                           dev_rel * value_params));
    }

    // times
    for (size_t i = 0; i < params.get_num_groups(); i++) {
        // incubation time
        double value_params = params.times[i].get_incubation();
        params.times[i].set_incubation(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // serial interval
        value_params = params.times[i].get_serialinterval();
        params.times[i].set_serialinterval(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // infectious time (mild)
        value_params = params.times[i].get_infectious_mild();
        params.times[i].set_infectious_mild(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // infectious to recovered (hospitalized to home)
        value_params = params.times[i].get_hospitalized_to_home();
        params.times[i].set_hospitalized_to_home(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // infectious to hospitalized (home to hospitalized)
        value_params = params.times[i].get_home_to_hospitalized();
        params.times[i].set_home_to_hospitalized(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // infectious (asymptomatic)
        value_params = params.times[i].get_infectious_asymp();
        params.times[i].set_infectious_asymp(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // hospitalized to ICU
        value_params = params.times[i].get_hospitalized_to_icu();
        params.times[i].set_hospitalized_to_icu(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // ICU to recovered
        value_params = params.times[i].get_icu_to_home();
        params.times[i].set_icu_to_home(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // ICU to death
        value_params = params.times[i].get_icu_to_dead();
        params.times[i].set_icu_to_death(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));
    }

    // probabilities
    for (size_t i = 0; i < params.get_num_groups(); i++) {
        // infection from contact
        double value_params = params.probabilities[i].get_infection_from_contact();
        params.probabilities[i].set_infection_from_contact(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // asymptomatic per infectious
        value_params = params.probabilities[i].get_asymp_per_infectious();
        params.probabilities[i].set_asymp_per_infectious(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // risk of infection from infectious
        value_params = params.probabilities[i].get_risk_from_symptomatic();
        params.probabilities[i].set_risk_from_symptomatic(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // deaths per icu treatments
        value_params = params.probabilities[i].get_dead_per_icu();
        params.probabilities[i].set_dead_per_icu(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // hospitalized per infections
        value_params = params.probabilities[i].get_hospitalized_per_infectious();
        params.probabilities[i].set_hospitalized_per_infectious(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));

        // icu treatments per hospitalized
        value_params = params.probabilities[i].get_icu_per_hospitalized();
        params.probabilities[i].set_icu_per_hospitalized(
            ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * value_params),
                                        (1 + dev_rel * 2.6) * value_params, value_params, dev_rel * value_params));
    }

    // maximum number of dampings; to avoid overfitting only allow one damping for every 10 days simulated
    // damping base values are between 0.1 and 1; diagonal values vary lie in the range of 0.6 to 1.4 times the base value
    // off diagonal values vary between 0.7 to 1.1 of the corresponding diagonal value (symmetrization is conducted)
    params.get_contact_patterns().set_distribution_damp_nb(ParameterDistributionUniform(1, (tmax - t0) / 10));
    params.get_contact_patterns().set_distribution_damp_days(ParameterDistributionUniform(t0, tmax));
    params.get_contact_patterns().set_distribution_damp_diag_base(ParameterDistributionUniform(0.1, 1));
    params.get_contact_patterns().set_distribution_damp_diag_rel(ParameterDistributionUniform(0.6, 1.4));
    params.get_contact_patterns().set_distribution_damp_offdiag_rel(ParameterDistributionUniform(0.7, 1.1));
}

void draw_sample(SecirParams& params)
{
    for (size_t i = 0; i < params.get_num_groups(); i++) {

        double group_total = params.populations.get_group_total(SecirCategory::AgeGroup, i);

        params.populations.get({i, SecirCompartments::E}).draw_sample();
        params.populations.get({i, SecirCompartments::C}).draw_sample();
        params.populations.get({i, SecirCompartments::I}).draw_sample();
        params.populations.get({i, SecirCompartments::H}).draw_sample();
        params.populations.get({i, SecirCompartments::U}).draw_sample();
        params.populations.get({i, SecirCompartments::R}).draw_sample();

        // no sampling for dead and total numbers
        // [...]

        params.populations.set_difference_from_group_total({i, epi::SecirCompartments::S}, epi::SecirCategory::AgeGroup,
                                                           i, group_total);
        params.populations.set_difference_from_group_total(
            {i, epi::SecirCompartments::S}, epi::SecirCategory::AgeGroup, i,
            params.populations.get_group_total(SecirCategory::AgeGroup, i));

        params.times[i].get_incubation().draw_sample();
        params.times[i].get_serialinterval().draw_sample();
        params.times[i].get_infectious_mild().draw_sample();
        params.times[i].get_hospitalized_to_home().draw_sample(); // here: home=recovered
        params.times[i].get_home_to_hospitalized().draw_sample(); // here: home=infectious
        params.times[i].get_infectious_asymp().draw_sample();
        params.times[i].get_hospitalized_to_icu().draw_sample();
        params.times[i].get_icu_to_dead().draw_sample();
        params.times[i].get_icu_to_home().draw_sample();

        params.probabilities[i].get_infection_from_contact().draw_sample();
        params.probabilities[i].get_asymp_per_infectious().draw_sample();
        params.probabilities[i].get_risk_from_symptomatic().draw_sample();
        params.probabilities[i].get_dead_per_icu().draw_sample();
        params.probabilities[i].get_hospitalized_per_infectious().draw_sample();
        params.probabilities[i].get_icu_per_hospitalized().draw_sample();
    }

    params.get_contact_patterns().draw_sample();

    params.check_constraints();
}

} // namespace epi