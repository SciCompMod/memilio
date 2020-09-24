import epidemiology.secir as secir

def run_secir_simulation():
    params = secir.SecirParams()

    params.times[0].set_incubation(5.2)
    params.times[0].set_infectious_mild(6)
    params.times[0].set_serialinterval(4.2)
    params.times[0].set_hospitalized_to_home(12)
    params.times[0].set_home_to_hospitalized(5)
    params.times[0].set_hospitalized_to_icu(2)
    params.times[0].set_icu_to_home(8)
    params.times[0].set_icu_to_death(5)

    params.get_contact_patterns().get_cont_freq_mat().set_cont_freq(0.5, 0, 0)
    params.get_contact_patterns().get_cont_freq_mat().add_damping(secir.Damping(30, 0.3), 0, 0)

    params.populations.set([0, secir.SecirCompartments.E], 100)
    params.populations.set([0, secir.SecirCompartments.C], 50)
    params.populations.set([0, secir.SecirCompartments.I], 50)
    params.populations.set([0, secir.SecirCompartments.H], 20)
    params.populations.set([0, secir.SecirCompartments.U], 10)
    params.populations.set([0, secir.SecirCompartments.R], 10)
    params.populations.set([0, secir.SecirCompartments.D], 0)
    params.populations.set_difference_from_total([0, secir.SecirCompartments.S], 10000)

    params.probabilities[0].set_infection_from_contact(1.0)
    params.probabilities[0].set_asymp_per_infectious(0.09)
    params.probabilities[0].set_risk_from_symptomatic(0.25)
    params.probabilities[0].set_hospitalized_per_infectious(0.2)
    params.probabilities[0].set_icu_per_hospitalized(0.25)
    params.probabilities[0].set_dead_per_icu(0.3)

    params.apply_constraints()

    result = secir.simulate(0, 50, 0.1, params)
    print(result.get_last_value())

if __name__ == "__main__":
    run_secir_simulation()