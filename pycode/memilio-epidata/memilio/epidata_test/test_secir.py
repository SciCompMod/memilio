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

from epidemiology.secir import (UncertainContactMatrix, ContactMatrix, Damping, SecirModel,
                                simulate, AgeGroup, Index_InfectionState, SecirSimulation)
from epidemiology.secir import InfectionState as State
import numpy as np

class Test_secir_integration(unittest.TestCase):

    def setUp(self):

        model = SecirModel(1)

        A0 = AgeGroup(0)

        model.parameters.IncubationTime[A0] = 5.2  # R_2^(-1)+R_3^(-1)
        model.parameters.InfectiousTimeMild[A0] =  6.  # 4-14  (=R4^(-1))
        model.parameters.SerialInterval[A0] = 4.2   # 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        model.parameters.HospitalizedToHomeTime[A0] = 12.  # 7-16 (=R5^(-1))
        model.parameters.HomeToHospitalizedTime[A0] = 5.  # 2.5-7 (=R6^(-1))
        model.parameters.HospitalizedToICUTime[A0] = 2.  # 1-3.5 (=R7^(-1))
        model.parameters.ICUToHomeTime[A0] = 8.  # 5-16 (=R8^(-1))
        model.parameters.ICUToDeathTime[A0] = 5.  # 3.5-7 (=R5^(-1))

        model.parameters.InfectionProbabilityFromContact[A0] = 1.0
        model.parameters.AsymptoticCasesPerInfectious[A0] = 0.09  # 0.01-0.16
        model.parameters.RiskOfInfectionFromSympomatic[A0] = 0.25  # 0.05-0.5
        model.parameters.HospitalizedCasesPerInfectious[A0] = 0.2  # 0.1-0.35
        model.parameters.ICUCasesPerHospitalized[A0] = 0.25  # 0.15-0.4
        model.parameters.DeathsPerHospitalized[A0] = 0.3  # 0.15-0.77

        model.populations[A0, Index_InfectionState(State.Susceptible)] = 7600
        model.populations[A0, Index_InfectionState(State.Exposed)] = 100
        model.populations[A0, Index_InfectionState(State.Carrier)] = 50
        model.populations[A0, Index_InfectionState(State.Infected)] = 50
        model.populations[A0, Index_InfectionState(State.Hospitalized)] = 20
        model.populations[A0, Index_InfectionState(State.ICU)] = 10
        model.populations[A0, Index_InfectionState(State.Recovered)] = 10
        model.populations[A0, Index_InfectionState(State.Dead)] = 0

        contacts = ContactMatrix(np.r_[0.5])
        contacts.add_damping(Damping(coeffs = np.r_[0.0], t = 0.0, level = 0, type = 0))
        model.parameters.ContactPatterns.cont_freq_mat[0] = contacts

        model.apply_constraints()

        self.model = model

    def test_simulate_simple(self):
      result = simulate(t0=0., tmax=100., dt=0.1, model=self.model)
      self.assertAlmostEqual(result.get_time(0), 0.)
      self.assertAlmostEqual(result.get_time(1), 0.1)
      self.assertAlmostEqual(result.get_last_time(), 100.)

    def test_simulation_simple(self):
        sim = SecirSimulation(self.model, t0 = 0., dt = 0.1)
        sim.advance(tmax = 100.)
        result = sim.result
        self.assertAlmostEqual(result.get_time(0), 0.)
        self.assertAlmostEqual(result.get_time(1), 0.1)
        self.assertAlmostEqual(result.get_last_time(), 100.)

if __name__ == '__main__':
    unittest.main()
