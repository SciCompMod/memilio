#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Maximilian Betz
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
from memilio.simulation.osir import Model, simulate, InfectionState


class Test_Integrators(unittest.TestCase):
    """ """

    def test_euler_step(self):
        """ """

        # Define a simple deriv function for step
        def deriv_function(y, t, dydt):
            """

            :param y: 
            :param t: 
            :param dydt: 

            """
            dydt[:] = -y

        integrator = mio.EulerIntegratorCore()
        yt = np.array([1.0, 2.0])
        t = 0.0
        dt = 0.1
        ytp1 = np.zeros_like(yt)

        result = integrator.step({deriv_function}, yt, t, dt, ytp1)
        self.assertTrue(result)
        self.assertTrue((ytp1 == [0.9, 1.8]).all())

    def test_model_integration(self):
        """ """

        model = Model(1)
        A0 = mio.AgeGroup(0)

        # Compartment transition duration
        model.parameters.TimeInfected[A0] = 6.

        # Compartment transition propabilities
        model.parameters.TransmissionProbabilityOnContact[A0] = 1.

        # Initial number of people in each compartment
        model.populations[A0, InfectionState.Infected] = 50
        model.populations[A0, InfectionState.Recovered] = 10
        model.populations.set_difference_from_total(
            (A0, InfectionState.Susceptible), 8000)

        integrator = mio.RKIntegratorCore()
        result1 = simulate(0, 5, 1, model, integrator)

        dt_max = 0.5
        integrator.dt_max = dt_max
        result2 = simulate(0, 5, 0.5, model, integrator)

        self.assertTrue((dt_max >= np.diff(result2.as_ndarray()[0, :])).all())
        self.assertFalse((result1.get_last_value() ==
                         result2.get_last_value()).all())


if __name__ == '__main__':
    unittest.main()
