#############################################################################
# Copyright (C) 2020-2026 MEmilio
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

from memilio.simulation import AgeGroup, Region
from memilio.simulation.oseir_metapop import (
    FlowSimulation,
    InfectionState as State,
    Model,
    simulate_flows,
)


class Test_oseir_metapop_integration(unittest.TestCase):
    """Python binding tests for ODE SEIR metapopulation model."""

    def setUp(self):
        model = Model(3, 1)
        A0 = AgeGroup(0)

        for region in range(3):
            model.populations[Region(region), A0, State.Susceptible] = 1000.0

        model.populations[Region(0), A0, State.Susceptible] = 990.0
        model.populations[Region(0), A0, State.Infected] = 10.0

        model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
            (1, 1)) * 2.7
        model.set_commuting_strengths_identity()

        self.model = model

    def test_population_after_commuting_is_bound(self):
        A0 = AgeGroup(0)
        population_after_commuting = self.model.parameters.PopulationAfterCommuting

        self.assertAlmostEqual(
            population_after_commuting[Region(0), A0].value, 1000.0)
        self.assertAlmostEqual(
            population_after_commuting.get_group_total_Region(Region(0)),
            1000.0)

    def test_simulate_flows_simple(self):
        flow_sim_results = simulate_flows(
            t0=0.0, tmax=1.0, dt=0.1, model=self.model)
        compartments = flow_sim_results[0]
        flows = flow_sim_results[1]

        self.assertAlmostEqual(compartments.get_time(0), 0.0)
        self.assertAlmostEqual(compartments.get_last_time(), 1.0)
        self.assertEqual(len(compartments.get_last_value()), 12)
        self.assertEqual(len(flows.get_last_value()), 9)

    def test_flow_simulation_simple(self):
        sim = FlowSimulation(self.model, t0=0.0, dt=0.1)
        sim.advance(tmax=1.0)

        self.assertAlmostEqual(sim.result.get_last_time(), 1.0)
        self.assertEqual(len(sim.result.get_last_value()), 12)
        self.assertEqual(len(sim.flows.get_last_value()), 9)


if __name__ == '__main__':
    unittest.main()
