#ifndef PARAMETER_SPACE_H
#define PARAMETER_SPACE_H

#include "epidemiology/utils/memory.h"
#include "epidemiology/utils/logging.h"
#include "epidemiology/utils/parameter_distributions.h"
#include "epidemiology/secir/secir.h"

#include <assert.h>
#include <string>
#include <vector>
#include <random>
#include <memory>

namespace epi
{
/* Sets alls SecirParams parameters normally distributed, 
*  using the current value and a given standard deviation
* @param[inout] params SecirParams including contact patterns for alle age groups
* @param[in] t0 start time
* @param[in] tmax end time
* @param[in] dev_rel maximum relative deviation from particular value(s) given in params
*/
template <class AgeGroup>
void set_params_distributions_normal(
    CompartmentalModel<Populations<AgeGroup, InfectionType>, SecirParams<(size_t)AgeGroup::Count>>& model, double t0,
    double tmax, double dev_rel)
{
    auto calc_distr = [dev_rel](double v, double min_val = 0.001){
        return ParameterDistributionNormal(std::max(min_val, (1 - dev_rel * 2.6) * v),
                                           (1 + dev_rel * 2.6) * v,
                                           v,
                                           dev_rel * v);
    };


    model.parameters.set_seasonality(  calc_distr(model.parameters.get_seasonality(), 0.0));
    model.parameters.set_icu_capacity( calc_distr(model.parameters.get_icu_capacity()));

    // populations
    for (size_t i = 0; i < model.parameters.get_num_groups(); i++) {
        for ( size_t j = 0; j < (size_t)InfectionType::Count; j++) {

            // don't touch S and D
            if ( j == (size_t)InfectionType::S || j == (size_t)InfectionType::D) {
                continue;
            }


            // variably sized groups
            model.populations.get(AgeGroup(i), InfectionType(j))
                .set_distribution(calc_distr(model.populations.get(AgeGroup(i), InfectionType(j))));
        }
    }

    // times
    for (size_t i = 0; i < model.parameters.get_num_groups(); i++) {

        model.parameters.times[i].set_incubation(           calc_distr(model.parameters.times[i].get_incubation()));
        model.parameters.times[i].set_serialinterval(       calc_distr(model.parameters.times[i].get_serialinterval()));
        model.parameters.times[i].set_infectious_mild(      calc_distr(model.parameters.times[i].get_infectious_mild()));
        model.parameters.times[i].set_hospitalized_to_home( calc_distr(model.parameters.times[i].get_hospitalized_to_home()));
        model.parameters.times[i].set_home_to_hospitalized( calc_distr(model.parameters.times[i].get_home_to_hospitalized()));
        model.parameters.times[i].set_infectious_asymp(     calc_distr(model.parameters.times[i].get_infectious_asymp()));
        model.parameters.times[i].set_hospitalized_to_icu(  calc_distr(model.parameters.times[i].get_hospitalized_to_icu()));
        model.parameters.times[i].set_icu_to_home(          calc_distr(model.parameters.times[i].get_icu_to_home()));
        model.parameters.times[i].set_icu_to_death(         calc_distr(model.parameters.times[i].get_icu_to_dead()));
    }

    // probabilities
    for (size_t i = 0; i < model.parameters.get_num_groups(); i++) {

        model.parameters.probabilities[i].set_infection_from_contact(      calc_distr(model.parameters.probabilities[i].get_infection_from_contact()));
        model.parameters.probabilities[i].set_carrier_infectability(       calc_distr(model.parameters.probabilities[i].get_carrier_infectability()));
        model.parameters.probabilities[i].set_asymp_per_infectious(        calc_distr(model.parameters.probabilities[i].get_asymp_per_infectious()));
        model.parameters.probabilities[i].set_risk_from_symptomatic(       calc_distr(model.parameters.probabilities[i].get_risk_from_symptomatic()));
        model.parameters.probabilities[i].set_dead_per_icu(                calc_distr(model.parameters.probabilities[i].get_dead_per_icu()));
        model.parameters.probabilities[i].set_hospitalized_per_infectious( calc_distr(model.parameters.probabilities[i].get_hospitalized_per_infectious()));
        model.parameters.probabilities[i].set_icu_per_hospitalized(        calc_distr(model.parameters.probabilities[i].get_icu_per_hospitalized()));
    }

    // maximum number of dampings; to avoid overfitting only allow one damping for every 10 days simulated
    // damping base values are between 0.1 and 1; diagonal values vary lie in the range of 0.6 to 1.4 times the base value
    // off diagonal values vary between 0.7 to 1.1 of the corresponding diagonal value (symmetrization is conducted)
    model.parameters.get_contact_patterns().set_distribution_damp_nb(ParameterDistributionUniform(1, (tmax - t0) / 10));
    model.parameters.get_contact_patterns().set_distribution_damp_days(ParameterDistributionUniform(t0, tmax));
    model.parameters.get_contact_patterns().set_distribution_damp_diag_base(ParameterDistributionUniform(0.1, 1));
    model.parameters.get_contact_patterns().set_distribution_damp_diag_rel(ParameterDistributionUniform(0.6, 1.4));
    model.parameters.get_contact_patterns().set_distribution_damp_offdiag_rel(ParameterDistributionUniform(0.7, 1.1));
}

/* Draws a sample from SecirParams parameter distributions and stores sample values
* as SecirParams parameter values (cf. UncertainValue and SecirParams classes)
* @param[inout] params SecirParams including contact patterns for alle age groups
*/
template <class AgeGroup>
void draw_sample(CompartmentalModel<Populations<AgeGroup, InfectionType>, SecirParams<(size_t)AgeGroup::Count>>& model)
{
    model.parameters.get_seasonality().draw_sample();
    model.parameters.get_icu_capacity().draw_sample();

    for (size_t i = 0; i < model.parameters.get_num_groups(); i++) {

        double group_total = model.populations.get_group_total(AgeGroup(i));

        model.populations.get(AgeGroup(i), InfectionType::E).draw_sample();
        model.populations.get(AgeGroup(i), InfectionType::C).draw_sample();
        model.populations.get(AgeGroup(i), InfectionType::I).draw_sample();
        model.populations.get(AgeGroup(i), InfectionType::H).draw_sample();
        model.populations.get(AgeGroup(i), InfectionType::U).draw_sample();
        model.populations.get(AgeGroup(i), InfectionType::R).draw_sample();

        // no sampling for dead and total numbers
        // [...]

        model.populations.set_difference_from_group_total(group_total, AgeGroup(i), AgeGroup(i), epi::InfectionType::S);
        model.populations.set_difference_from_group_total(model.populations.get_group_total(AgeGroup(i)), AgeGroup(i),
                                                          AgeGroup(i), epi::InfectionType::S);

        model.parameters.times[i].get_incubation().draw_sample();
        model.parameters.times[i].get_serialinterval().draw_sample();
        model.parameters.times[i].get_infectious_mild().draw_sample();
        model.parameters.times[i].get_hospitalized_to_home().draw_sample(); // here: home=recovered
        model.parameters.times[i].get_home_to_hospitalized().draw_sample(); // here: home=infectious
        model.parameters.times[i].get_infectious_asymp().draw_sample();
        model.parameters.times[i].get_hospitalized_to_icu().draw_sample();
        model.parameters.times[i].get_icu_to_dead().draw_sample();
        model.parameters.times[i].get_icu_to_home().draw_sample();

        model.parameters.probabilities[i].get_infection_from_contact().draw_sample();
        model.parameters.probabilities[i].get_asymp_per_infectious().draw_sample();
        model.parameters.probabilities[i].get_risk_from_symptomatic().draw_sample();
        model.parameters.probabilities[i].get_dead_per_icu().draw_sample();
        model.parameters.probabilities[i].get_hospitalized_per_infectious().draw_sample();
        model.parameters.probabilities[i].get_icu_per_hospitalized().draw_sample();
    }

    model.parameters.get_contact_patterns().draw_sample();

    model.apply_constraints();
}

} // namespace epi

#endif // PARAMETER_SPACE_H
