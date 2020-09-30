import epidemiology.secir as secir
import numpy as np

def parameter_study():    
    t0 = 0
    tmax = 50

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

    params.probabilities[0].set_infection_from_contact(1.0)
    params.probabilities[0].set_carrier_infectability(0.67)
    params.probabilities[0].set_asymp_per_infectious(0.09)
    params.probabilities[0].set_risk_from_symptomatic(0.25)
    params.probabilities[0].set_hospitalized_per_infectious(0.2)
    params.probabilities[0].set_icu_per_hospitalized(0.25)
    params.probabilities[0].set_dead_per_icu(0.3)

    #two regions with different populations and with some migration between them
    graph = secir.MigrationGraph()    
    params.populations.set([0, secir.SecirCompartments.E], 100)
    params.populations.set([0, secir.SecirCompartments.C], 50)
    params.populations.set([0, secir.SecirCompartments.I], 50)
    params.populations.set([0, secir.SecirCompartments.H], 20)
    params.populations.set([0, secir.SecirCompartments.U], 10)
    params.populations.set([0, secir.SecirCompartments.R], 10)
    params.populations.set([0, secir.SecirCompartments.D], 0)
    params.populations.set_difference_from_total([0, secir.SecirCompartments.S], 10000)
    params.apply_constraints()
    graph.add_node(params, t0)
    params.populations.set([0, secir.SecirCompartments.E], 0)
    params.populations.set([0, secir.SecirCompartments.C], 0)
    params.populations.set([0, secir.SecirCompartments.I], 0)
    params.populations.set([0, secir.SecirCompartments.H], 0)
    params.populations.set([0, secir.SecirCompartments.U], 0)
    params.populations.set([0, secir.SecirCompartments.R], 0)
    params.populations.set([0, secir.SecirCompartments.D], 0)
    params.populations.set_difference_from_total([0, secir.SecirCompartments.S], 2000)
    params.apply_constraints()
    graph.add_node(params, t0)
    migration_coefficients = 0.1 * np.ones(8)
    migration_params = secir.MigrationParams(migration_coefficients)
    graph.add_edge(0, 1, migration_params) #one coefficient per (age group x compartment)
    graph.add_edge(1, 0, migration_params) #directed graph -> add both directions so coefficients can be different
    
    #run simulation
    sim = secir.MigrationSimulation(graph, t0, dt = 0.5)
    sim.advance(tmax)

    #process results
    region0_result = sim.graph.get_node(0).get_result()
    region1_result = sim.graph.get_node(1).get_result()

if __name__ == "__main__":
    parameter_study()