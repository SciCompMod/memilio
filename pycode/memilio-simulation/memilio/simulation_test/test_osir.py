#############################################################################
# Copyright (C) 2020-2025 MEmilio
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

from memilio.simulation import AgeGroup, Damping
from memilio.simulation.osir import Index_InfectionState
from memilio.simulation.osir import InfectionState as State
from memilio.simulation.osir import Model, simulate


class Test_osir_integration(unittest.TestCase):
    """ """

    def setUp(self):
        """ """

        model = Model(1)
        A0 = AgeGroup(0)

        model.parameters.TimeInfected[A0] = 6.
        model.parameters.TransmissionProbabilityOnContact[A0] = 1.

        model.populations[A0, State.Susceptible] = 4800
        model.populations[A0, State.Infected] = 50
        model.populations[A0, State.Recovered] = 50

        model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
            (1, 1))
        model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.zeros(
            (1, 1))
        model.parameters.ContactPatterns.cont_freq_mat.add_damping(
            Damping(coeffs=np.r_[0.9], t=30.0, level=0, type=0))

        model.check_constraints()

        self.model = model

    def test_simulate_simple(self):
        """ """
        result = simulate(t0=0., tmax=100., dt=0.1, model=self.model)
        self.assertAlmostEqual(result.get_time(0), 0.)
        self.assertAlmostEqual(result.get_time(1), 0.1)
        self.assertAlmostEqual(result.get_last_time(), 100.)

    def test_check_constraints_parameters(self):
        """ """

        model = Model(1)
        A0 = AgeGroup(0)

        model.parameters.TimeInfected[A0] = 6.
        model.parameters.TransmissionProbabilityOnContact[A0] = 1.

        model.parameters.TimeInfected[A0] = 6.
        model.parameters.TransmissionProbabilityOnContact[A0] = 1.
        self.assertEqual(model.parameters.check_constraints(), 0)

        model.parameters.TimeInfected[A0] = 0
        self.assertEqual(model.parameters.check_constraints(), 1)

        model.parameters.TimeInfected[A0] = 6.
        model.parameters.TransmissionProbabilityOnContact[A0] = -1.
        self.assertEqual(model.parameters.check_constraints(), 1)


if __name__ == '__main__':
    unittest.main()
