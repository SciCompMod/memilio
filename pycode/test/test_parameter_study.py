import unittest
import epidemiology.secir as secir
from epidemiology.secir import InfectionState as State
import numpy as np

class Test_ParameterStudy(unittest.TestCase):

    def _get_model(self):
        model = secir.SecirModel(1)

        A0 = secir.AgeGroup(0)

        model.parameters.IncubationTime[A0] = 5.2
        model.parameters.InfectiousTimeMild[A0] = 6
        model.parameters.SerialInterval[A0] = 4.2
        model.parameters.HospitalizedToHomeTime[A0] = 12
        model.parameters.HomeToHospitalizedTime[A0] = 5
        model.parameters.HospitalizedToICUTime[A0] = 2
        model.parameters.ICUToHomeTime[A0] = 8
        model.parameters.ICUToDeathTime[A0] = 5

        model.parameters.ContactPatterns.cont_freq_mat[0] = secir.ContactMatrix(np.r_[0.5])
        model.parameters.ContactPatterns.cont_freq_mat.add_damping(secir.Damping(np.r_[0.7], 30.0))

        model.populations[A0, secir.Index_InfectionState(State.Exposed)] = 100
        model.populations[A0, secir.Index_InfectionState(State.Carrier)] = 50
        model.populations[A0, secir.Index_InfectionState(State.Infected)] = 50
        model.populations[A0, secir.Index_InfectionState(State.Hospitalized)] = 20
        model.populations[A0, secir.Index_InfectionState(State.ICU)] = 10
        model.populations[A0, secir.Index_InfectionState(State.Recovered)] = 10
        model.populations[A0, secir.Index_InfectionState(State.Dead)] = 0
        model.populations.set_difference_from_total((A0, secir.Index_InfectionState(State.Susceptible)), 10000)

        model.parameters.InfectionProbabilityFromContact[A0] = 1.0
        model.parameters.AsymptoticCasesPerInfectious[A0] = 0.09
        model.parameters.RiskOfInfectionFromSympomatic[A0] = 0.25
        model.parameters.HospitalizedCasesPerInfectious[A0] = 0.2
        model.parameters.ICUCasesPerHospitalized[A0] = 0.25
        model.parameters.DeathsPerHospitalized[A0] = 0.3

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

        #run as graph
        def handle_result_func(graph):
            handle_result_func.c += 1
            handle_result_func.results.append(graph)
            self.assertAlmostEqual(graph.get_node(0).property.result.get_time(0), t0)
            self.assertAlmostEqual(graph.get_node(0).property.result.get_last_time(), tmax)

        handle_result_func.c = 0
        handle_result_func.results = []
        study.run(handle_result_func)

        self.assertEqual(handle_result_func.c, num_runs)
        
        #run as single node
        def handle_single_result_func(sim):
            handle_single_result_func.c += 1
            handle_single_result_func.results.append(sim)
            self.assertAlmostEqual(sim.result.get_time(0), t0)
            self.assertAlmostEqual(sim.result.get_last_time(), tmax)

        handle_single_result_func.c = 0
        handle_single_result_func.results = []
        study.run_single(handle_single_result_func)

        self.assertEqual(handle_single_result_func.c, num_runs)

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
