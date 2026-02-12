#############################################################################
# Copyright (C) 2020-2026 MEmilio
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

import os
import numpy as np
import pandas as pd

from memilio.simulation import ContactMatrix, Damping, UncertainContactMatrix, AgeGroup
from memilio.simulation.osecir import Index_InfectionState
from memilio.simulation.osecir import InfectionState as State
from memilio.simulation.osecir import Model, Simulation, simulate


class Test_osecir_integration(unittest.TestCase):
    """ """

    here = os.path.dirname(os.path.abspath(__file__))

    def setUp(self):
        """ """

        self.t0 = 0
        self.tmax = 50
        self.dt = 0.1

        cont_freq = 10
        nb_total_t0, nb_exp_t0, nb_inf_t0, nb_car_t0, nb_hosp_t0, nb_icu_t0, nb_rec_t0, nb_dead_t0 = 10000, 100, 50, 50, 20, 10, 10, 0

        self.nb_comp = 10
        self.nb_groups = 3
        fact = 1.0/self.nb_groups

        model = Model(self.nb_groups)

        model.parameters.StartDay = 60
        model.parameters.Seasonality.value = 0.2
        model.parameters.TestAndTraceCapacity.value = 35

        for i in range(0, self.nb_groups):
            Ai = AgeGroup(i)

            model.parameters.TimeExposed[Ai] = 3.2
            model.parameters.TimeInfectedNoSymptoms[Ai] = 2.
            model.parameters.TimeInfectedSymptoms[Ai] = 5.8
            model.parameters.TimeInfectedSevere[Ai] = 9.5
            model.parameters.TimeInfectedCritical[Ai] = 7.1

            model.populations[Ai, State.Exposed] = fact * nb_exp_t0
            model.populations[Ai, State.InfectedNoSymptoms] = fact * nb_car_t0
            model.populations[Ai, State.InfectedNoSymptomsConfirmed] = 0
            model.populations[Ai, State.InfectedSymptoms] = fact * nb_inf_t0
            model.populations[Ai, State.InfectedSymptomsConfirmed] = 0
            model.populations[Ai, State.InfectedSevere] = fact * nb_hosp_t0
            model.populations[Ai, State.InfectedCritical] = fact * nb_icu_t0
            model.populations[Ai, State.Recovered] = fact * nb_rec_t0
            model.populations[Ai, State.Dead] = fact * nb_dead_t0
            model.populations.set_difference_from_group_total_AgeGroup(
                (Ai, State.Susceptible), fact * nb_total_t0)

            model.parameters.TransmissionProbabilityOnContact[Ai] = 0.05
            model.parameters.RelativeTransmissionNoSymptoms[Ai] = 0.7
            model.parameters.RecoveredPerInfectedNoSymptoms[Ai] = 0.09
            model.parameters.RiskOfInfectionFromSymptomatic[Ai] = 0.25
            model.parameters.MaxRiskOfInfectionFromSymptomatic[Ai] = 0.45
            model.parameters.SeverePerInfectedSymptoms[Ai] = 0.2
            model.parameters.CriticalPerSevere[Ai] = 0.25
            model.parameters.DeathsPerCritical[Ai] = 0.3

        model.apply_constraints()

        contacts = ContactMatrix(
            np.full((self.nb_groups, self.nb_groups), fact * cont_freq))
        contacts.add_damping(
            Damping(
                coeffs=np.full(
                    (self.nb_groups, self.nb_groups),
                    0.7),
                t=30.0, level=0, type=0))
        model.parameters.ContactPatterns.cont_freq_mat[0] = contacts

        self.model = model

    def test_simulate_simple(self):
        """ """
        result = simulate(t0=0., tmax=100., dt=0.1, model=self.model)
        self.assertAlmostEqual(result.get_time(0), 0.)
        self.assertAlmostEqual(result.get_time(1), 0.1)
        self.assertAlmostEqual(result.get_last_time(), 100.)

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
        results contained in the file ode-secihurd-ageres-compare.csv.
        If cpp model changes this test needs to be adjusted accordingly.


        """
        refData = pd.read_csv(
            os.path.join(self.here + '/data/ode-secihurd-ageres-compare.csv'),
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
                delta=1e-10)

            for index_compartment in range(0, self.nb_comp):
                dummy = 0
                for index_agegroup in range(0, self.nb_groups):
                    dummy += result[index_timestep][
                        index_compartment + self.nb_comp * index_agegroup]

                self.assertAlmostEqual(
                    timestep[index_compartment + 1],
                    dummy, delta=1e-10)


if __name__ == '__main__':
    unittest.main()
