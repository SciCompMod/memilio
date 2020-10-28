import epidemiology.secir as secir

def parameter_study():    
    #setup basic parameters
    model = secir.SecirModel1()

    model.parameters.times[0].set_incubation(5.2)
    model.parameters.times[0].set_infectious_mild(6)
    model.parameters.times[0].set_serialinterval(4.2)
    model.parameters.times[0].set_hospitalized_to_home(12)
    model.parameters.times[0].set_home_to_hospitalized(5)
    model.parameters.times[0].set_hospitalized_to_icu(2)
    model.parameters.times[0].set_icu_to_home(8)
    model.parameters.times[0].set_icu_to_death(5)

    model.parameters.get_contact_patterns().get_cont_freq_mat().set_cont_freq(0.5, 0, 0)
    model.parameters.get_contact_patterns().get_cont_freq_mat().add_damping(secir.Damping(30, 0.3), 0, 0)

    model.populations.set(100, secir.AgeGroup1.Group0, secir.InfectionType.E)
    model.populations.set( 50, secir.AgeGroup1.Group0, secir.InfectionType.C)
    model.populations.set( 50, secir.AgeGroup1.Group0, secir.InfectionType.I)
    model.populations.set( 20, secir.AgeGroup1.Group0, secir.InfectionType.H)
    model.populations.set( 10, secir.AgeGroup1.Group0, secir.InfectionType.U)
    model.populations.set( 10, secir.AgeGroup1.Group0, secir.InfectionType.R)
    model.populations.set(  0, secir.AgeGroup1.Group0, secir.InfectionType.D)
    model.populations.set_difference_from_total(10000, secir.AgeGroup1.Group0, secir.InfectionType.S)

    model.parameters.probabilities[0].set_infection_from_contact(1.0)
    model.parameters.probabilities[0].set_carrier_infectability(0.67)
    model.parameters.probabilities[0].set_asymp_per_infectious(0.09)
    model.parameters.probabilities[0].set_risk_from_symptomatic(0.25)
    model.parameters.probabilities[0].set_hospitalized_per_infectious(0.2)
    model.parameters.probabilities[0].set_icu_per_hospitalized(0.25)
    model.parameters.probabilities[0].set_dead_per_icu(0.3)

    #process the result of one run
    parameter_study.c = 0
    def handle_result(model, result, _):
        print("run {} with infection rate {:.2G}".format(parameter_study.c, model.parameters.probabilities[0].get_infection_from_contact().value))
        print("compartments at t = {}:".format(result.get_time(0)))
        print(result.get_value(0))
        print("compartments at t = {}:".format(result.get_last_time()))
        print(result.get_last_value())
        parameter_study.c += 1

    #study the effect of different infection rates
    infection_from_contact_uniform = secir.ParameterDistributionUniform(0.5, 1.0)
    model.parameters.probabilities[0].set_infection_from_contact(infection_from_contact_uniform)

    model.apply_constraints()
    
    t0 = 0
    tmax = 50
    study = secir.ParameterStudy1(model, t0, tmax, num_runs = 3)
    study.run(handle_result)

if __name__ == "__main__":
    parameter_study()