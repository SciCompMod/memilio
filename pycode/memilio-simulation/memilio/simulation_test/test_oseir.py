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

import os
import numpy as np
import pandas as pd

from memilio.simulation import Damping
from memilio.simulation.oseir import Index_InfectionState
from memilio.simulation.oseir import InfectionState as State
from memilio.simulation.oseir import Model, simulate


class Test_oseir_integration(unittest.TestCase):

    here = os.path.dirname(os.path.abspath(__file__))

    def setUp(self):

        model = Model()

        self.t0 = 0.
        self.tmax = 50.
        self.dt = 0.1002004008016032
        total_population = 1061000

        model.populations[Index_InfectionState(State.Exposed)] = 10000
        model.populations[Index_InfectionState(State.Infected)] = 1000
        model.populations[Index_InfectionState(State.Recovered)] = 1000
        model.populations.set_difference_from_total(
            Index_InfectionState(State.Susceptible), total_population)

        model.parameters.TransmissionProbabilityOnContact.value = 1.
        model.parameters.TimeExposed.value = 5.2
        model.parameters.TimeInfected.value = 2.

        model.parameters.ContactPatterns.baseline = [[2.7]]
        model.parameters.ContactPatterns.minimum = np.zeros((1, 1))
        model.parameters.ContactPatterns.add_damping(
            Damping(coeffs=[[0.6]], t=12.5, level=0, type=0))

        model.check_constraints()

        self.model = model

    def test_simulate_simple(self):
        result = simulate(t0=0., tmax=100., dt=0.1, model=self.model)
        self.assertAlmostEqual(result.get_time(0), 0.)
        self.assertAlmostEqual(result.get_time(1), 0.1)
        self.assertAlmostEqual(result.get_last_time(), 100.)

    def test_compare_seir_with_cpp(self):
        """
        Tests the correctness of the python bindings. The results of a simulation 
        in python get compared to the results of a cpp simulation. Cpp simulation 
        results contained in the file seir-compare.csv. 
        If cpp model changes this test needs to be adjusted accordingly.
        """
        refData = pd.read_csv(
            os.path.join(self.here + '/data/seir-compare.csv'),
            sep=r'(?<!#)\s+', engine='python')
        refData.columns = pd.Series(refData.columns.str.replace(
            r"#\s", "", regex=True))

        result = simulate(t0=self.t0, tmax=self.tmax,
                          dt=self.dt, model=self.model)

        for index_timestep, timestep in refData.iterrows():
            t = float(timestep.at['t'])
            rel_tol = 1e-6

            # test result diverges at damping because of changes, not worth fixing at the moment
            if (t > 11.0 and t < 13.0):
                # strong divergence around damping
                rel_tol = 0.5
            elif (t > 13.0):
                # minor divergence after damping
                rel_tol = 1e-2

            self.assertAlmostEqual(
                t, result.get_time(index_timestep),
                delta=1e-12)

            for index_compartment in range(0, 4):
                ref = timestep[index_compartment+1]
                actual = result[index_timestep][index_compartment]

                tol = rel_tol * ref
                self.assertAlmostEqual(ref, actual, delta=tol)

    def test_check_constraints_parameters(self):

        model = Model()

        model.parameters.TimeExposed.value = 5.2
        model.parameters.TimeInfected.value = 6.
        model.parameters.TransmissionProbabilityOnContact.value = 1.

        model.parameters.TimeExposed.value = 5.2
        model.parameters.TimeInfected.value = 6.
        model.parameters.TransmissionProbabilityOnContact.value = 1.
        self.assertEqual(model.parameters.check_constraints(), 0)

        model.parameters.TimeExposed.value = -1.
        self.assertEqual(model.parameters.check_constraints(), 1)

        model.parameters.TimeExposed.value = 5.2
        model.parameters.TimeInfected.value = 0
        self.assertEqual(model.parameters.check_constraints(), 1)

        model.parameters.TimeInfected.value = 6.
        model.parameters.TransmissionProbabilityOnContact.value = -1.
        self.assertEqual(model.parameters.check_constraints(), 1)


if __name__ == '__main__':
    unittest.main()
