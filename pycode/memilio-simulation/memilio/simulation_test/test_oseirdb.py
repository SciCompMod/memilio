#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Henrik Zunker
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

from memilio.simulation import AgeGroup, Damping
from memilio.simulation.oseirdb import (InfectionState as State, Model,
                                        simulate, simulate_flows)


class Test_oseirdb_integration(unittest.TestCase):
    """ """

    def setUp(self):
        """ """

        model = Model(1)
        A0 = AgeGroup(0)

        self.t0 = 0.0
        self.tmax = 50.0
        self.dt = 0.5
        total_population = 10000

        model.populations[A0, State.Exposed] = 1000
        model.populations[A0, State.Infected] = 100
        model.populations[A0, State.Recovered] = 100
        model.populations[A0, State.Dead] = 25
        model.populations[A0, State.Buried] = 20
        model.populations.set_difference_from_total(
            (A0, State.Susceptible), total_population)

        model.parameters.TransmissionProbabilityOnContact[A0] = 0.1
        model.parameters.TransmissionProbabilityFromDead[A0] = 0.01
        model.parameters.ProbabilityToRecover[A0] = 0.6
        model.parameters.TimeExposed[A0] = 5.2
        model.parameters.TimeInfected[A0] = 2.0
        model.parameters.TimeToBurial[A0] = 2.0

        model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
            (1, 1)) * 2.0
        model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.zeros(
            (1, 1))
        model.parameters.ContactPatterns.cont_freq_mat.add_damping(
            Damping(coeffs=np.r_[0.5], t=10.0, level=0, type=0))

        model.check_constraints()

        self.model = model

    def test_simulate_simple(self):
        """ """
        result = simulate(t0=0.0, tmax=10.0, dt=0.5, model=self.model)
        self.assertAlmostEqual(result.get_time(0), 0.0)
        self.assertAlmostEqual(result.get_time(1), 0.5)
        self.assertAlmostEqual(result.get_last_time(), 10.0)
        self.assertAlmostEqual(
            result.get_last_value().sum(),
            result.get_value(0).sum(),
            delta=1e-6)

    def test_flow_simulation_simple(self):
        """ """
        flow_sim_results = simulate_flows(
            t0=0.0, tmax=5.0, dt=0.5, model=self.model)
        flows = flow_sim_results[1]
        self.assertEqual(flows.get_time(0), 0.0)
        self.assertEqual(flows.get_last_time(), 5.0)
        # S->E, E->I, I->R, I->D, D->B
        self.assertEqual(len(flows.get_last_value()), 5)

        compartments = flow_sim_results[0]
        self.assertEqual(compartments.get_time(0), 0.0)
        self.assertEqual(compartments.get_last_time(), 5.0)
        self.assertEqual(len(compartments.get_last_value()), 6)

    def test_check_constraints_parameters(self):
        """ """

        model = Model(1)
        A0 = AgeGroup(0)

        model.parameters.TimeExposed[A0] = 5.2
        model.parameters.TimeInfected[A0] = 6.0
        model.parameters.TimeToBurial[A0] = 2.0
        model.parameters.TransmissionProbabilityOnContact[A0] = 0.04
        model.parameters.TransmissionProbabilityFromDead[A0] = 0.03
        model.parameters.ProbabilityToRecover[A0] = 0.4
        self.assertEqual(model.parameters.check_constraints(), 0)

        model.parameters.TimeExposed[A0] = -1.0
        self.assertEqual(model.parameters.check_constraints(), 1)

        model.parameters.TimeExposed[A0] = 5.2
        model.parameters.TimeInfected[A0] = 0
        self.assertEqual(model.parameters.check_constraints(), 1)

        model.parameters.TimeInfected[A0] = 6.0
        model.parameters.TimeToBurial[A0] = 0
        self.assertEqual(model.parameters.check_constraints(), 1)

        model.parameters.TimeToBurial[A0] = 2.0
        self.assertEqual(model.parameters.check_constraints(), 0)
        model.parameters.TransmissionProbabilityOnContact[A0] = -1.0
        self.assertEqual(model.parameters.check_constraints(), 1)

        model.parameters.TransmissionProbabilityOnContact[A0] = 0.1
        self.assertEqual(model.parameters.check_constraints(), 0)
        model.parameters.TransmissionProbabilityFromDead[A0] = -1.0
        self.assertEqual(model.parameters.check_constraints(), 1)

        model.parameters.TransmissionProbabilityFromDead[A0] = 0.01
        self.assertEqual(model.parameters.check_constraints(), 0)
        model.parameters.ProbabilityToRecover[A0] = -0.1
        self.assertEqual(model.parameters.check_constraints(), 1)


if __name__ == '__main__':
    unittest.main()
