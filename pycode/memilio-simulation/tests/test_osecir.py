#############################################################################
# Copyright (C) 2020-2026 MEmilio
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

import os
import numpy as np
import pandas as pd

from memilio.simulation import AgeGroup, ContactMatrix, Damping, UncertainContactMatrix
from memilio.simulation.osecir import Index_InfectionState
from memilio.simulation.osecir import InfectionState as State
from memilio.simulation.osecir import Model, Simulation, simulate, simulate_flows
from memilio.simulation.osecir import ParameterStudy, GraphParameterStudy, ModelGraph


class Test_osecir_integration(unittest.TestCase):
    """ """

    here = os.path.dirname(os.path.abspath(__file__))

    def setUp(self):
        """ """

        model = Model(1)

        self.t0 = 0
        self.tmax = 50
        self.dt = 0.1

        cont_freq = 10
        nb_total_t0, nb_exp_t0, nb_inf_t0, nb_car_t0, nb_hosp_t0, nb_icu_t0, nb_rec_t0, nb_dead_t0 = 10000, 100, 50, 50, 20, 10, 10, 0

        A0 = AgeGroup(0)

        model.parameters.StartDay = 60
        model.parameters.Seasonality.value = 0.2

        model.parameters.TimeExposed[A0] = 3.2
        model.parameters.TimeInfectedNoSymptoms[A0] = 2.
        model.parameters.TimeInfectedSymptoms[A0] = 5.8
        model.parameters.TimeInfectedSevere[A0] = 9.5
        model.parameters.TimeInfectedCritical[A0] = 7.1

        model.populations[A0, State.Exposed] = nb_exp_t0
        model.populations[A0, State.InfectedNoSymptoms] = nb_car_t0
        model.populations[A0, State.InfectedNoSymptomsConfirmed] = 0
        model.populations[A0, State.InfectedSymptoms] = nb_inf_t0
        model.populations[A0, State.InfectedSymptomsConfirmed] = 0
        model.populations[A0, State.InfectedSevere] = nb_hosp_t0
        model.populations[A0, State.InfectedCritical] = nb_icu_t0
        model.populations[A0, State.Recovered] = nb_rec_t0
        model.populations[A0, State.Dead] = nb_dead_t0
        model.populations.set_difference_from_total(
            (A0, State.Susceptible), nb_total_t0)

        model.parameters.TransmissionProbabilityOnContact[A0] = 0.05
        model.parameters.RelativeTransmissionNoSymptoms[A0] = 0.7
        model.parameters.RecoveredPerInfectedNoSymptoms[A0] = 0.09
        model.parameters.RiskOfInfectionFromSymptomatic[A0] = 0.25
        model.parameters.MaxRiskOfInfectionFromSymptomatic[A0] = 0.45
        model.parameters.TestAndTraceCapacity.value = 35
        model.parameters.SeverePerInfectedSymptoms[A0] = 0.2
        model.parameters.CriticalPerSevere[A0] = 0.25
        model.parameters.DeathsPerCritical[A0] = 0.3

        contacts = ContactMatrix(np.r_[cont_freq])
        contacts.add_damping(
            Damping(coeffs=np.r_[0.7], t=30.0, level=0, type=0))
        model.parameters.ContactPatterns.cont_freq_mat[0] = contacts

        model.apply_constraints()

        self.model = model

    def test_simulate_simple(self):
        """ """
        result = simulate(t0=0., tmax=100., dt=0.1, model=self.model)
        self.assertAlmostEqual(result.get_time(0), 0.)
        self.assertAlmostEqual(result.get_time(1), 0.1)
        self.assertAlmostEqual(result.get_last_time(), 100.)

    def test_flow_simulation_simple(self):
        """ """
        flow_sim_results = simulate_flows(
            t0=0., tmax=100., dt=0.1, model=self.model)
        flows = flow_sim_results[1]
        self.assertEqual(flows.get_time(0), 0.)
        self.assertEqual(flows.get_last_time(), 100.)
        self.assertEqual(len(flows.get_last_value()), 15)

        compartments = flow_sim_results[0]
        self.assertEqual(compartments.get_time(0), 0.)
        self.assertEqual(compartments.get_last_time(), 100.)
        self.assertEqual(len(compartments.get_last_value()), 10)

    def test_simulation_simple(self):
        """ """
        sim = Simulation(self.model, t0=0., dt=0.1)
        sim.advance(tmax=100.)
        result = sim.result
        self.assertAlmostEqual(result.get_time(0), 0.)
        self.assertAlmostEqual(result.get_time(1), 0.1)
        self.assertAlmostEqual(result.get_last_time(), 100.)

    def test_compare_with_cpp(self):
        """Tests the correctness of the python bindings. The results of a simulation
        in python get compared to the results of a cpp simulation. Cpp simulation
        results contained in the file ode-secihurd-compare.csv.
        If cpp model changes this test needs to be adjusted accordingly.


        """
        refData = pd.read_csv(
            os.path.join(self.here + '/data/ode-secihurd-compare.csv'),
            sep=r'(?<!#)\s+', engine='python')
        refData.columns = pd.Series(refData.columns.str.replace(
            r"#\s", "", regex=True))

        result = simulate(t0=self.t0, tmax=self.tmax,
                          dt=self.dt, model=self.model)

        # compare num elements
        for index_timestep, timestep in refData.iterrows():
            # compare num elements
            t = float(timestep.at['t'])
            self.assertAlmostEqual(
                t, result.get_time(index_timestep),
                delta=1e-9)

            for index_compartment in range(0, 10):
                self.assertAlmostEqual(
                    timestep[index_compartment + 1],
                    result[index_timestep][index_compartment], delta=1e-9)

    def test_study(self):
        """ Runs a parameterstudy with a single model, to check that it is possible """

        t0 = 1
        tmax = 10
        dt = 0.5
        num_runs = 3

        study = ParameterStudy(self.model, t0, tmax, dt, num_runs)

        results = study.run()

        self.assertEqual(len(results), num_runs)
        self.assertEqual(results[0].result.get_last_time(), tmax)

    def test_study_graph(self):
        """ Runs a parameterstudy with a graph, to check that it is possible """

        t0 = 1
        tmax = 10
        dt = 0.5
        num_runs = 3

        graph = ModelGraph()
        graph.add_node(0, self.model)

        study = GraphParameterStudy(graph, t0, tmax, dt, num_runs)

        results = study.run()

        self.assertEqual(len(results), num_runs)
        self.assertEqual(results[0].get_node(
            0).property.result.get_last_time(), tmax)


if __name__ == '__main__':
    unittest.main()
