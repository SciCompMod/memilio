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

from memilio.simulation.omseirs4 import InfectionState as State
from memilio.simulation.omseirs4 import Model, Simulation, simulate
import memilio.simulation as mio


class Test_omseirs4_integration(unittest.TestCase):
    """Python binding tests for ODE MSEIRS4 model."""

    def setUp(self):
        # create a simple single-population model
        model = Model()

        # parameters
        model.parameters.BaseTransmissionRate.value = 0.2
        model.parameters.SeasonalAmplitude.value = 0.0
        model.parameters.SeasonalPhase.value = 0.0
        model.parameters.NaturalBirthDeathRate.value = 0.0
        model.parameters.LossMaternalImmunityRate.value = 1.0 / 90.0
        model.parameters.ProgressionRate.value = 1.0 / 7.0
        model.parameters.RecoveryRate.value = 1.0 / 14.0
        model.parameters.ImmunityWaningRate.value = 0.0
        # Beta2/3/4Factor keep defaults

        # initialize populations
        N = 100000.0
        assigned = 0.0
        model.populations[State.E1] = 20.0
        model.populations[State.I1] = 10.0
        assigned += 30.0
        model.populations[State.R1] = 100.0
        model.populations[State.R2] = 50.0
        assigned += 150.0
        model.populations[State.S2] = 500.0
        model.populations[State.S3] = 200.0
        model.populations[State.S4] = 100.0
        assigned += 800.0
        model.populations[State.MaternalImmune] = 50.0
        assigned += 50.0
        model.populations[State.S1] = max(0.0, N - assigned)

        model.check_constraints()

        self.model = model

    def test_simulate_simple(self):
        integrator = mio.EulerIntegratorCore()
        result = simulate(t0=0.0, tmax=50.0, dt=0.5,
                          model=self.model, integrator=integrator)
        self.assertAlmostEqual(result.get_time(0), 0.0)
        self.assertAlmostEqual(result.get_time(1), 0.5)
        self.assertAlmostEqual(result.get_last_time(), 50.0)
        # vector length equals number of infection states
        self.assertEqual(len(result.get_last_value()),
                         17)  # 17 infection states

    def test_simulation_simple(self):
        sim = Simulation(self.model, t0=0.0, dt=1.0)
        sim.advance(tmax=10.0)
        res = sim.result
        self.assertAlmostEqual(res.get_time(0), 0.0)
        self.assertAlmostEqual(res.get_last_time(), 10.0)
        self.assertEqual(len(res.get_last_value()), 17)  # 17 infection states

    def test_check_and_apply_constraints_parameters(self):
        model = Model()
        p = model.parameters

        # default valid
        self.assertEqual(p.check_constraints(), False)

        p.SeasonalAmplitude.value = -0.1
        self.assertEqual(p.check_constraints(), True)
        p.SeasonalAmplitude.value = 0.5
        self.assertEqual(p.check_constraints(), False)
        p.SeasonalAmplitude.value = 1.2
        self.assertEqual(p.check_constraints(), True)

        p.BaseTransmissionRate.value = -0.3
        p.ProgressionRate.value = -1.0
        p.RecoveryRate.value = -1.0
        p.ImmunityWaningRate.value = -1.0
        p.NaturalBirthDeathRate.value = -1.0
        p.LossMaternalImmunityRate.value = -1.0
        p.Beta2Factor.value = -0.1
        p.Beta3Factor.value = -0.2
        p.Beta4Factor.value = -0.3
        p.SeasonalAmplitude.value = 1.5
        changed = p.apply_constraints()
        self.assertTrue(changed)
        self.assertAlmostEqual(p.BaseTransmissionRate.value, 0.0)
        self.assertAlmostEqual(p.ProgressionRate.value, 0.0)
        self.assertAlmostEqual(p.RecoveryRate.value, 0.0)
        self.assertAlmostEqual(p.ImmunityWaningRate.value, 0.0)
        self.assertAlmostEqual(p.NaturalBirthDeathRate.value, 0.0)
        self.assertAlmostEqual(p.LossMaternalImmunityRate.value, 0.0)
        self.assertAlmostEqual(p.Beta2Factor.value, 0.0)
        self.assertAlmostEqual(p.Beta3Factor.value, 0.0)
        self.assertAlmostEqual(p.Beta4Factor.value, 0.0)
        self.assertAlmostEqual(p.SeasonalAmplitude.value, 1.0)


if __name__ == '__main__':
    unittest.main()
