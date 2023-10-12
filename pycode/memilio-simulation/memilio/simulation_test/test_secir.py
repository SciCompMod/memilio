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

from memilio.simulation import ContactMatrix, Damping, UncertainContactMatrix
from memilio.simulation.secir import AgeGroup, Index_InfectionState
from memilio.simulation.secir import InfectionState as State
from memilio.simulation.secir import Model, Simulation, simulate


class Test_secir_integration(unittest.TestCase):

    def setUp(self):

        model = Model(1)

        A0 = AgeGroup(0)

        model.parameters.IncubationTime[A0] = 5.2
        model.parameters.TimeInfectedSymptoms[A0] = 6.
        model.parameters.SerialInterval[A0] = 4.2
        model.parameters.TimeInfectedSevere[A0] = 12.
        model.parameters.TimeInfectedCritical[A0] = 8.

        model.parameters.TransmissionProbabilityOnContact[A0] = 1.0
        model.parameters.RecoveredPerInfectedNoSymptoms[A0] = 0.09
        model.parameters.RiskOfInfectionFromSymptomatic[A0] = 0.25
        model.parameters.SeverePerInfectedSymptoms[A0] = 0.2
        model.parameters.CriticalPerSevere[A0] = 0.25
        model.parameters.DeathsPerCritical[A0] = 0.3

        model.populations[A0, State.Susceptible] = 7600
        model.populations[A0, State.Exposed] = 100
        model.populations[A0, State.InfectedNoSymptoms] = 50
        model.populations[A0, State.InfectedSymptoms] = 50
        model.populations[A0, State.InfectedSevere] = 20
        model.populations[A0, State.InfectedCritical] = 10
        model.populations[A0, State.Recovered] = 10
        model.populations[A0, State.Dead] = 0

        contacts = ContactMatrix(np.r_[0.5])
        contacts.add_damping(
            Damping(coeffs=np.r_[0.0], t=0.0, level=0, type=0))
        model.parameters.ContactPatterns.cont_freq_mat[0] = contacts

        model.apply_constraints()

        self.model = model

    def test_simulate_simple(self):
        result = simulate(t0=0., tmax=100., dt=0.1, model=self.model)
        self.assertAlmostEqual(result.get_time(0), 0.)
        self.assertAlmostEqual(result.get_time(1), 0.1)
        self.assertAlmostEqual(result.get_last_time(), 100.)

    def test_simulation_simple(self):
        sim = Simulation(self.model, t0=0., dt=0.1)
        sim.advance(tmax=100.)
        result = sim.result
        self.assertAlmostEqual(result.get_time(0), 0.)
        self.assertAlmostEqual(result.get_time(1), 0.1)
        self.assertAlmostEqual(result.get_last_time(), 100.)


if __name__ == '__main__':
    unittest.main()
