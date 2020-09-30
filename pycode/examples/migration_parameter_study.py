import epidemiology.secir as secir
import numpy as np

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

    params.probabilities[0].set_infection_from_contact(1.0)
    params.probabilities[0].set_carrier_infectability(0.67)
    params.probabilities[0].set_asymp_per_infectious(0.09)
    params.probabilities[0].set_risk_from_symptomatic(0.25)
    params.probabilities[0].set_hospitalized_per_infectious(0.2)
    params.probabilities[0].set_icu_per_hospitalized(0.25)
    params.probabilities[0].set_dead_per_icu(0.3)

    #two regions with different populations and with some migration between them
    params_graph = secir.SecirParamsGraph()    
    params.populations.set([0, secir.SecirCompartments.E], 100)
    params.populations.set([0, secir.SecirCompartments.C], 50)
    params.populations.set([0, secir.SecirCompartments.I], 50)
    params.populations.set([0, secir.SecirCompartments.H], 20)
    params.populations.set([0, secir.SecirCompartments.U], 10)
    params.populations.set([0, secir.SecirCompartments.R], 10)
    params.populations.set([0, secir.SecirCompartments.D], 0)
    params.populations.set_difference_from_total([0, secir.SecirCompartments.S], 10000)
    params.apply_constraints()
    params_graph.add_node(params)
    params.populations.set([0, secir.SecirCompartments.E], 0)
    params.populations.set([0, secir.SecirCompartments.C], 0)
    params.populations.set([0, secir.SecirCompartments.I], 0)
    params.populations.set([0, secir.SecirCompartments.H], 0)
    params.populations.set([0, secir.SecirCompartments.U], 0)
    params.populations.set([0, secir.SecirCompartments.R], 0)
    params.populations.set([0, secir.SecirCompartments.D], 0)
    params.populations.set_difference_from_total([0, secir.SecirCompartments.S], 2000)
    params.apply_constraints()
    params_graph.add_node(params)
    migration_coefficients = 0.1 * np.ones(8)
    migration_params = secir.MigrationParams(migration_coefficients)
    params_graph.add_edge(0, 1, migration_params) #one coefficient per (age group x compartment)
    params_graph.add_edge(1, 0, migration_params) #directed graph -> add both directions so coefficients can be different

    #process the result of one run
    parameter_study.c = 0
    def handle_result(params, result, node_idx):
        if node_idx == 0:
            print("run {}".format(parameter_study.c))
        print("  node {}".format(node_idx))
        print("  initial carrier count {:.2G}.".format(params.populations.get([0, secir.SecirCompartments.C]).value))
        print("  compartments at t = {}:".format(result.get_time(0)))
        print("  ", result.get_value(0))
        print("  compartments at t = {}:".format(result.get_last_time()))
        print("  ", result.get_last_value())
        if node_idx == 1:
            parameter_study.c += 1

    #study with unknown number of undetected carriers
    carrier_distribution = secir.ParameterDistributionNormal(50, 2000, 200, 100)
    params_graph.get_node(0).populations.set([0, secir.SecirCompartments.C], carrier_distribution)
    
    t0 = 0
    tmax = 50
    study = secir.ParameterStudy(params_graph, t0, tmax, graph_dt = 1.0, num_runs = 3)
    study.run(handle_result)

if __name__ == "__main__":
    parameter_study()