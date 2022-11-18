#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors:
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################
import unittest

import numpy as np

import memilio.simulation as mio
import memilio.simulation.secir as secir


class Test_ParameterStudy(unittest.TestCase):

    def _get_model(self):
        model = secir.Model(1)

        A0 = secir.AgeGroup(0)

        model.parameters.IncubationTime[A0] = 5.2
        model.parameters.TimeInfectedSymptoms[A0] = 6
        model.parameters.SerialInterval[A0] = 4.2
        model.parameters.TimeInfectedSevere[A0] = 12
        model.parameters.TimeInfectedCritical[A0] = 8

        model.parameters.ContactPatterns.cont_freq_mat[0] = mio.ContactMatrix(np.r_[
                                                                              0.5])
        model.parameters.ContactPatterns.cont_freq_mat.add_damping(
            mio.Damping(np.r_[0.7], 30.0))

        model.populations[A0, secir.InfectionState.Exposed] = 100
        model.populations[A0, secir.InfectionState.InfectedNoSymptoms] = 50
        model.populations[A0, secir.InfectionState.InfectedSymptoms] = 50
        model.populations[A0, secir.InfectionState.InfectedSevere] = 20
        model.populations[A0, secir.InfectionState.InfectedCritical] = 10
        model.populations[A0, secir.InfectionState.Recovered] = 10
        model.populations[A0, secir.InfectionState.Dead] = 0
        model.populations.set_difference_from_total(
            (A0, secir.InfectionState.Susceptible), 10000)

        model.parameters.TransmissionProbabilityOnContact[A0] = 1.0
        model.parameters.RecoveredPerInfectedNoSymptoms[A0] = 0.09
        model.parameters.RiskOfInfectionFromSymptomatic[A0] = 0.25
        model.parameters.SeverePerInfectedSymptoms[A0] = 0.2
        model.parameters.CriticalPerSevere[A0] = 0.25
        model.parameters.DeathsPerCritical[A0] = 0.3

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

        # run as graph
        def handle_result_func(graph):
            handle_result_func.c += 1
            handle_result_func.results.append(graph)
            self.assertAlmostEqual(graph.get_node(
                0).property.result.get_time(0), t0)
            self.assertAlmostEqual(graph.get_node(
                0).property.result.get_last_time(), tmax)

        handle_result_func.c = 0
        handle_result_func.results = []
        study.run(handle_result_func)

        self.assertEqual(handle_result_func.c, num_runs)

        # run as single node
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
        graph = secir.ModelGraph()
        graph.add_node(0, model)
        graph.add_node(1, model)
        graph.add_edge(0, 1, 0.01 * np.ones(8))
        graph.add_edge(1, 0, 0.01 * np.ones(8))

        study = secir.ParameterStudy(graph, t0=1, tmax=10, dt=0.5, num_runs=3)

        self.assertEqual(study.secir_model_graph.num_nodes, 2)
        self.assertEqual(study.secir_model_graph.num_edges, 2)


if __name__ == '__main__':
    unittest.main()
