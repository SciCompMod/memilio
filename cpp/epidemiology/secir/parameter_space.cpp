#include "epidemiology/secir/parameter_space.h"
#include "epidemiology/utils/parameter_distributions.h"
#include "epidemiology/secir/secir.h"

namespace epi
{

void set_params_distributions_normal(
    SecirModel& model, double t0,
    double tmax, double dev_rel)
{
    auto set_distribution = [dev_rel](UncertainValue& v, double min_val = 0.001){
        v.set_distribution( ParameterDistributionNormal(std::max(min_val,
                                                       (1 - dev_rel * 2.6) * v),
                                                       (1 + dev_rel * 2.6) * v,
                                                       v,
                                                       dev_rel * v));
    };


    set_distribution(model.parameters.get_seasonality(), 0.0);
    set_distribution(model.parameters.get_icu_capacity());
    set_distribution(model.parameters.get_test_and_trace_capacity());

    // populations
    for (size_t i = 0; i < model.parameters.get_num_groups(); i++) {
        for ( size_t j = 0; j < (size_t)InfectionState::Count; j++) {

            // don't touch S and D
            if ( j == (size_t)InfectionState::Susceptible || j == (size_t)InfectionState::Dead) {
                continue;
            }


            // variably sized groups
            set_distribution(model.populations[{AgeGroup(i), InfectionState(j)}]);
        }
    }

    // times
    for (size_t i = 0; i < model.parameters.get_num_groups(); i++) {

        set_distribution(model.parameters.times[i].get_incubation());
        set_distribution(model.parameters.times[i].get_serialinterval());
        set_distribution(model.parameters.times[i].get_infectious_mild());
        set_distribution(model.parameters.times[i].get_hospitalized_to_home());
        set_distribution(model.parameters.times[i].get_home_to_hospitalized());
        set_distribution(model.parameters.times[i].get_infectious_asymp());
        set_distribution(model.parameters.times[i].get_hospitalized_to_icu());
        set_distribution(model.parameters.times[i].get_icu_to_home());
        set_distribution(model.parameters.times[i].get_icu_to_dead());
    }

    // probabilities
    for (size_t i = 0; i < model.parameters.get_num_groups(); i++) {

        set_distribution(model.parameters.probabilities[i].get_infection_from_contact());
        set_distribution(model.parameters.probabilities[i].get_carrier_infectability());
        set_distribution(model.parameters.probabilities[i].get_asymp_per_infectious());
        set_distribution(model.parameters.probabilities[i].get_risk_from_symptomatic());
        set_distribution(model.parameters.probabilities[i].get_test_and_trace_max_risk_from_symptomatic());
        set_distribution(model.parameters.probabilities[i].get_dead_per_icu());
        set_distribution(model.parameters.probabilities[i].get_hospitalized_per_infectious());
        set_distribution(model.parameters.probabilities[i].get_icu_per_hospitalized());
    }

    // maximum number of dampings; to avoid overfitting only allow one damping for every 10 days simulated
    // damping base values are between 0.1 and 1; diagonal values vary lie in the range of 0.6 to 1.4 times the base value
    // off diagonal values vary between 0.7 to 1.1 of the corresponding diagonal value (symmetrization is conducted)
    model.parameters.get_contact_patterns().set_distribution_damp_nb(ParameterDistributionUniform(1, (tmax - t0) / 10));
    model.parameters.get_contact_patterns().set_distribution_damp_days(ParameterDistributionUniform(t0, tmax));
    model.parameters.get_contact_patterns().set_distribution_damp_diag_base(ParameterDistributionUniform(0.0, 0.9));
    model.parameters.get_contact_patterns().set_distribution_damp_diag_rel(ParameterDistributionUniform(0.0, 0.4));
    model.parameters.get_contact_patterns().set_distribution_damp_offdiag_rel(ParameterDistributionUniform(0.0, 0.3));
}


void draw_sample_demographics(SecirModel& model)
{
    model.parameters.get_icu_capacity().draw_sample();
    model.parameters.get_test_and_trace_capacity().draw_sample();

    for (size_t i = 0; i < model.parameters.get_num_groups(); i++) {
        double group_total = model.populations.get_group_total(AgeGroup(i));

        model.populations[{AgeGroup(i), InfectionState::Exposed}].draw_sample();
        model.populations[{AgeGroup(i), InfectionState::Carrier}].draw_sample();
        model.populations[{AgeGroup(i), InfectionState::Infected}].draw_sample();
        model.populations[{AgeGroup(i), InfectionState::Hospitalized}].draw_sample();
        model.populations[{AgeGroup(i), InfectionState::ICU}].draw_sample();
        model.populations[{AgeGroup(i), InfectionState::Recovered}].draw_sample();

        // no sampling for dead and total numbers
        // [...]

        model.populations.set_difference_from_group_total<AgeGroup>({AgeGroup(i), epi::InfectionState::Susceptible}, group_total);
        model.populations.set_difference_from_group_total<AgeGroup>({AgeGroup(i), epi::InfectionState::Susceptible},
                                                                    model.populations.get_group_total(AgeGroup(i)));
    }
}

void draw_sample_infection(SecirModel& model)
{
    model.parameters.get_seasonality().draw_sample();

    //not age dependent
    model.parameters.times[0].get_incubation().draw_sample();
    model.parameters.times[0].get_serialinterval().draw_sample();
    model.parameters.times[0].get_infectious_mild().draw_sample();
    model.parameters.times[0].get_hospitalized_to_icu().draw_sample();
    model.parameters.probabilities[0].get_carrier_infectability().draw_sample();
    model.parameters.probabilities[0].get_risk_from_symptomatic().draw_sample();
    model.parameters.probabilities[0].get_test_and_trace_max_risk_from_symptomatic().draw_sample();

    for (size_t i = 0; i < model.parameters.get_num_groups(); i++) {
        //not age dependent
        model.parameters.times[i].get_incubation()          = model.parameters.times[0].get_incubation();
        model.parameters.times[i].get_serialinterval()      = model.parameters.times[0].get_serialinterval();
        model.parameters.times[i].get_infectious_mild()     = model.parameters.times[0].get_infectious_mild();
        model.parameters.times[i].get_hospitalized_to_icu() = model.parameters.times[0].get_hospitalized_to_icu();
        model.parameters.probabilities[i].get_carrier_infectability() =
            model.parameters.probabilities[0].get_carrier_infectability();
        model.parameters.probabilities[i].get_risk_from_symptomatic() =
            model.parameters.probabilities[0].get_risk_from_symptomatic();
        model.parameters.probabilities[i].get_test_and_trace_max_risk_from_symptomatic() =
            model.parameters.probabilities[0].get_test_and_trace_max_risk_from_symptomatic();

        //age dependent
        model.parameters.times[i].get_hospitalized_to_home().draw_sample(); // here: home=recovered
        model.parameters.times[i].get_home_to_hospitalized().draw_sample(); // here: home=infectious
        model.parameters.times[i].get_infectious_asymp().draw_sample();
        model.parameters.times[i].get_icu_to_dead().draw_sample();
        model.parameters.times[i].get_icu_to_home().draw_sample();

        model.parameters.probabilities[i].get_infection_from_contact().draw_sample();
        model.parameters.probabilities[i].get_asymp_per_infectious().draw_sample();
        model.parameters.probabilities[i].get_dead_per_icu().draw_sample();
        model.parameters.probabilities[i].get_hospitalized_per_infectious().draw_sample();
        model.parameters.probabilities[i].get_icu_per_hospitalized().draw_sample();
    }
}

void draw_sample(SecirModel& model)
{
    draw_sample_infection(model);
    draw_sample_demographics(model);
    model.parameters.get_contact_patterns().draw_sample();
    model.apply_constraints();
}

} // namespace epi
