import unittest
import epidemiology.secir as secir
import numpy as np

class Test_ParameterStudy(unittest.TestCase):
    def _get_params(self):
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
        return params

    def test_run(self):
        params = self._get_params()

        t0 = 1
        tmax = 10
        num_runs = 3
        study = secir.ParameterStudy(params, t0, tmax, num_runs)

        self.assertEqual(study.t0, t0)
        self.assertEqual(study.tmax, tmax)
        self.assertEqual(study.num_runs, num_runs)

        #mock callback
        def handle_result_func(params, result, node_idx):
            handle_result_func.c += 1
            self.assertAlmostEqual(result.get_time(0), t0)
            self.assertAlmostEqual(result.get_last_time(), tmax)
            self.assertEqual(node_idx, 0)

        handle_result_func.c = 0
        result = study.run(handle_result_func)

        self.assertEqual(handle_result_func.c, num_runs)
        self.assertEqual(len(result), num_runs)
        
        handle_result_func.c = 0
        result = study.run_single(handle_result_func)

        self.assertEqual(handle_result_func.c, num_runs)
        self.assertEqual(len(result), num_runs)

    def test_graph(self):        
        params = self._get_params()
        graph = secir.SecirParamsGraph()
        graph.add_node(params)
        graph.add_node(params)
        graph.add_edge(0, 1, 0.01 * np.ones(8))
        graph.add_edge(1, 0, 0.01 * np.ones(8))

        t0 = 1
        tmax = 10
        num_runs = 3
        study = secir.ParameterStudy(graph, t0, tmax, num_runs)

        self.assertEqual(study.secir_params_graph.num_nodes, 2)
        self.assertEqual(study.secir_params_graph.num_edges, 2)

if __name__ == '__main__':
    unittest.main()