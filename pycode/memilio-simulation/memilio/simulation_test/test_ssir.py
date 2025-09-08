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

from memilio.simulation import Damping, LogLevel, set_log_level
from memilio.simulation.ssir import InfectionState as State
from memilio.simulation.ssir import (
    Model,
    simulate_stochastic,
    interpolate_simulation_result,
)


class Test_ssir_integration(unittest.TestCase):
    """Integration tests for the SDE SIR Python bindings."""

    set_log_level(LogLevel.Off)

    def setUp(self):
        # Initialize model similar to the example
        self.population = 10000
        self.days = 5.0
        self.dt = 0.1

        model = Model()

        # Set parameters
        model.parameters.TimeInfected.value = 10.0
        model.parameters.TransmissionProbabilityOnContact.value = 1.0

        # Initial compartment counts
        model.populations[State.Infected] = 100
        model.populations[State.Recovered] = 1000
        model.populations.set_difference_from_total(
            State.Susceptible, self.population)

        # Contacts and damping
        model.parameters.ContactPatterns.baseline = np.ones((1, 1)) * 2.7
        model.parameters.ContactPatterns.minimum = np.zeros((1, 1))
        model.parameters.ContactPatterns.add_damping(
            Damping(coeffs=np.r_[0.6], t=12.5, level=0, type=0)
        )

        model.check_constraints()
        self.model = model

    def test_simulate_simple(self):
        result = simulate_stochastic(0.0, self.days, self.dt, self.model)
        self.assertAlmostEqual(result.get_time(0), 0.0)
        self.assertAlmostEqual(result.get_time(1), self.dt)
        self.assertAlmostEqual(result.get_last_time(), self.days)
        self.assertEqual(len(result.get_last_value()), 3)

    def test_check_constraints_parameters(self):
        model = Model()
        # valid
        model.parameters.TimeInfected.value = 6.0
        model.parameters.TransmissionProbabilityOnContact.value = 1.0
        self.assertEqual(model.parameters.check_constraints(), 0)
        # invalid: time <= 0.1
        model.parameters.TimeInfected.value = 0.0
        self.assertEqual(model.parameters.check_constraints(), 1)
        # invalid: transmission < 0
        model.parameters.TimeInfected.value = 6.0
        model.parameters.TransmissionProbabilityOnContact.value = -1.0
        self.assertEqual(model.parameters.check_constraints(), 1)
        # invalid: transmission > 1
        model.parameters.TransmissionProbabilityOnContact.value = 2.0
        self.assertEqual(model.parameters.check_constraints(), 1)

    def test_apply_constraints_parameters(self):
        model = Model()
        # set invalid values
        model.parameters.TimeInfected.value = 0.0
        model.parameters.TransmissionProbabilityOnContact.value = -0.5
        corrected = model.parameters.apply_constraints()
        self.assertTrue(corrected)
        # after apply_constraints(), values should be set according to constraints
        self.assertAlmostEqual(
            model.parameters.TimeInfected.value, 0.1, places=6)
        self.assertAlmostEqual(
            model.parameters.TransmissionProbabilityOnContact.value, 0.0,
            places=6)

    def test_population_conservation(self):
        result = simulate_stochastic(0.0, self.days, self.dt, self.model)
        total0 = self.population
        atol = 1e-10 * total0
        for i in range(int(self.days / self.dt) + 1):
            vals = result[i]
            self.assertEqual(len(vals), 3)
            # Non-negativity
            self.assertTrue(np.all(np.array(vals) >= -1e-10))
            # Conservation of total population
            self.assertAlmostEqual(np.sum(vals), total0, delta=atol)

    def test_interpolate_times(self):
        result = simulate_stochastic(0.0, self.days, self.dt, self.model)
        times = [0.0, 1.0, 2.0, 3.0, self.days]
        interp = interpolate_simulation_result(result, times)
        self.assertEqual(interp.get_time(0), times[0])
        self.assertEqual(interp.get_last_time(), times[-1])
        self.assertEqual(len(interp.get_value(0)), 3)
        self.assertEqual(len(interp.get_value(len(times) - 1)), 3)


if __name__ == '__main__':
    unittest.main()
