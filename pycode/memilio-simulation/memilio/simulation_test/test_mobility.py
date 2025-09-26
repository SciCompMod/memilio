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
import tempfile

import numpy as np

import memilio.simulation as mio
import memilio.simulation.osecir as osecir


class Test_Mobility(unittest.TestCase):
    """ """

    def test_params(self):
        """ """
        coeffs = mio.MobilityCoefficientGroup(1, 10)
        coeffs[0] = mio.MobilityCoefficients(np.ones(10))
        coeffs[0].add_damping(mio.MobilityDamping(0.5 * np.ones(10), t=1.0))
        params = mio.MobilityParameters(coeffs)
        self.assertTrue(
            (params.coefficients.get_matrix_at(0) == np.ones(10)).all())
        self.assertTrue((params.coefficients.get_matrix_at(2)
                        == 0.5 * np.ones(10)).all())

    def test_params_graph(self):
        """ """
        graph = osecir.ModelGraph()
        graph.add_node(0, osecir.Model(1))
        graph.add_node(1, osecir.Model(1))
        graph.add_edge(0, 1, np.ones(10))
        self.assertEqual(graph.num_nodes, 2)
        self.assertEqual(graph.num_edges, 1)
        self.assertEqual(graph.get_num_out_edges(0), 1)
        self.assertEqual(graph.get_num_out_edges(1), 0)

    def test_sim_graph(self):
        """ """
        graph = osecir.MobilityGraph()
        graph.add_node(0, osecir.Model(1), 0, 0.1)
        graph.add_node(1, osecir.Model(1), 0)
        graph.add_edge(0, 1, np.ones(10))
        self.assertEqual(graph.num_nodes, 2)
        self.assertEqual(graph.num_edges, 1)
        self.assertEqual(graph.get_num_out_edges(0), 1)
        self.assertEqual(graph.get_num_out_edges(1), 0)

    def test_mobility_sim(self):
        """ """
        graph = osecir.MobilityGraph()
        graph.add_node(0, osecir.Model(1), 0, 0.1)
        graph.add_node(1, osecir.Model(1), 0)
        graph.add_edge(0, 1, np.ones(10))

        sim = osecir.MobilitySimulation(graph, t0=0.0)
        sim.advance(2)

        # integration does adaptive time steps so exact count is unknown
        self.assertGreaterEqual(sim.graph.get_node(
            0).property.result.get_num_time_points(), 3)

    def test_write_read_graph_simple(self):
        # build a simple model graph
        model = osecir.Model(1)
        model.parameters.TestAndTraceCapacity.value = 42
        model.apply_constraints()

        graph = osecir.ModelGraph()
        graph.add_node(0, model)
        graph.add_node(1, model)

        num_compartments = 10
        graph.add_edge(0, 1, 0.1 * np.ones(num_compartments))
        graph.add_edge(1, 0,  0.1 * np.ones(num_compartments))

        with tempfile.TemporaryDirectory() as tmpdir:
            # save graph
            osecir.write_graph(graph, tmpdir)
            # read graph back
            g_read = osecir.read_graph(tmpdir)

        # basic structure
        self.assertEqual(graph.num_nodes, g_read.num_nodes)
        self.assertEqual(graph.num_edges, g_read.num_edges)

        # check one parameter
        self.assertEqual(
            g_read.get_node(0).property.parameters.TestAndTraceCapacity.value,
            42,
        )


if __name__ == '__main__':
    unittest.main()
