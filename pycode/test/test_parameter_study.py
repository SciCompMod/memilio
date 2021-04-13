import unittest
import epidemiology.secir as secir
from epidemiology.secir import InfectionState as State
import numpy as np

class Test_ParameterStudy(unittest.TestCase):

    def _get_model(self):
        model = secir.SecirModel(1)

        model.parameters.times[0].set_incubation(5.2)
        model.parameters.times[0].set_infectious_mild(6)
        model.parameters.times[0].set_serialinterval(4.2)
        model.parameters.times[0].set_hospitalized_to_home(12)
        model.parameters.times[0].set_home_to_hospitalized(5)
        model.parameters.times[0].set_hospitalized_to_icu(2)
        model.parameters.times[0].set_icu_to_home(8)
        model.parameters.times[0].set_icu_to_death(5)

        model.parameters.get_contact_patterns().cont_freq_mat[0] = secir.ContactMatrix(np.r_[0.5])
        model.parameters.get_contact_patterns().cont_freq_mat.add_damping(secir.Damping(np.r_[0.7], 30.0))

        model.populations[secir.AgeGroup(0), secir.Index_InfectionState(State.Exposed)] = 100
        model.populations[secir.AgeGroup(0), secir.Index_InfectionState(State.Carrier)] = 50
        model.populations[secir.AgeGroup(0), secir.Index_InfectionState(State.Infected)] = 50
        model.populations[secir.AgeGroup(0), secir.Index_InfectionState(State.Hospitalized)] = 20
        model.populations[secir.AgeGroup(0), secir.Index_InfectionState(State.ICU)] = 10
        model.populations[secir.AgeGroup(0), secir.Index_InfectionState(State.Recovered)] = 10
        model.populations[secir.AgeGroup(0), secir.Index_InfectionState(State.Dead)] = 0
        model.populations.set_difference_from_total((secir.AgeGroup(0), secir.Index_InfectionState(State.Susceptible)), 10000)

        model.parameters.probabilities[0].set_infection_from_contact(1.0)
        model.parameters.probabilities[0].set_asymp_per_infectious(0.09)
        model.parameters.probabilities[0].set_risk_from_symptomatic(0.25)
        model.parameters.probabilities[0].set_hospitalized_per_infectious(0.2)
        model.parameters.probabilities[0].set_icu_per_hospitalized(0.25)
        model.parameters.probabilities[0].set_dead_per_icu(0.3)

        model.apply_constraints()
        return model


    def test_run(self):
        model = self._get_model()

        t0 = 1
        tmax = 10
        num_runs = 3
        study = secir.ParameterStudy(model, t0, tmax, num_runs)

        self.assertEqual(study.t0, t0)
        self.assertEqual(study.tmax, tmax)
        self.assertEqual(study.num_runs, num_runs)

        #mock callback
        def handle_result_func(graph):
            handle_result_func.c += 1
            self.assertAlmostEqual(graph.get_node(0).property.result.get_time(0), t0)
            self.assertAlmostEqual(graph.get_node(0).property.result.get_last_time(), tmax)

        handle_result_func.c = 0
        result = study.run(handle_result_func)

        self.assertEqual(handle_result_func.c, num_runs)
        self.assertEqual(len(result), num_runs)
        
        handle_result_func.c = 0
        result = study.run_single(handle_result_func)

        self.assertEqual(handle_result_func.c, num_runs)
        self.assertEqual(len(result), num_runs)

    def test_graph(self):        
        model = self._get_model()
        graph = secir.SecirModelGraph()
        graph.add_node(0, model)
        graph.add_node(1, model)
        graph.add_edge(0, 1, 0.01 * np.ones(8))
        graph.add_edge(1, 0, 0.01 * np.ones(8))

        study = secir.ParameterStudy(graph, t0 = 1, tmax = 10, dt = 0.5, num_runs = 3)

        self.assertEqual(study.secir_model_graph.num_nodes, 2)
        self.assertEqual(study.secir_model_graph.num_edges, 2)

if __name__ == '__main__':
    unittest.main()
