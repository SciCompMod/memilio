#############################################################################
# Copyright (C) 2020-2024 MEmilio
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

import memilio.simulation as mio
import memilio.simulation.secir as secir


class Test_Movement(unittest.TestCase):
    def test_params(self):
        coeffs = mio.MovementCoefficientGroup(1, 10)
        coeffs[0] = mio.MovementCoefficients(np.ones(10))
        coeffs[0].add_damping(mio.MovementDamping(0.5 * np.ones(10), t=1.0))
        params = mio.MovementParameters(coeffs)
        self.assertTrue(
            (params.coefficients.get_matrix_at(0) == np.ones(10)).all())
        self.assertTrue((params.coefficients.get_matrix_at(2)
                        == 0.5 * np.ones(10)).all())

    def test_params_graph(self):
        graph = secir.ModelGraph()
        graph.add_node(0, secir.Model(1))
        graph.add_node(1, secir.Model(1))
        graph.add_edge(0, 1, np.ones(10))
        self.assertEqual(graph.num_nodes, 2)
        self.assertEqual(graph.num_edges, 1)
        self.assertEqual(graph.get_num_out_edges(0), 1)
        self.assertEqual(graph.get_num_out_edges(1), 0)

    def test_sim_graph(self):
        graph = secir.MovementGraph()
        graph.add_node(0, secir.Model(1), 0, 0.1)
        graph.add_node(1, secir.Model(1), 0)
        graph.add_edge(0, 1, np.ones(10))
        self.assertEqual(graph.num_nodes, 2)
        self.assertEqual(graph.num_edges, 1)
        self.assertEqual(graph.get_num_out_edges(0), 1)
        self.assertEqual(graph.get_num_out_edges(1), 0)

    def test_movement_sim(self):
        graph = secir.MovementGraph()
        graph.add_node(0, secir.Model(1), 0, 0.1)
        graph.add_node(1, secir.Model(1), 0)
        graph.add_edge(0, 1, np.ones(10))

        sim = secir.MovementSimulation(graph, t0=0.0)
        sim.advance(2)

        # integration does adaptive time steps so exact count is unknown
        self.assertGreaterEqual(sim.graph.get_node(
            0).property.result.get_num_time_points(), 3)


if __name__ == '__main__':
    unittest.main()
