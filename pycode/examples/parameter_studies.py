import epidemiology.secir as secir

def parameter_study():    
    #setup basic parameters
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

    #process the result of one run
    parameter_study.c = 0
    def handle_result(params, result, _):
        print("run {} with infection rate {:.2G}".format(parameter_study.c, params.probabilities[0].get_infection_from_contact().value))
        print("compartments at t = {}:".format(result.get_time(0)))
        print(result.get_value(0))
        print("compartments at t = {}:".format(result.get_last_time()))
        print(result.get_last_value())
        parameter_study.c += 1

    #study the effect of different infection rates
    infection_from_contact_uniform = secir.ParameterDistributionUniform(0.5, 1.0)
    params.probabilities[0].set_infection_from_contact(infection_from_contact_uniform)

    params.apply_constraints()
    
    t0 = 0
    tmax = 50
    study = secir.ParameterStudy(params, t0, tmax, num_runs = 3)
    study.run(handle_result)

if __name__ == "__main__":
    parameter_study()