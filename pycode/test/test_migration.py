import unittest
import epidemiology.secir as secir
import numpy as np

class Test_Migration(unittest.TestCase):
    def test_params_graph(self):
        graph = secir.SecirParamsGraph()
        graph.add_node(0, secir.SecirParams())
        graph.add_node(1, secir.SecirParams())
        graph.add_edge(0, 1, np.ones(8))
        self.assertEqual(graph.num_nodes, 2)
        self.assertEqual(graph.num_edges, 1)
        self.assertEqual(graph.get_num_out_edges(0), 1)
        self.assertEqual(graph.get_num_out_edges(1), 0)

    def test_sim_graph(self):
        graph = secir.MigrationGraph()
        graph.add_node(0, secir.SecirParams(), 0, 0.1)
        graph.add_node(1, secir.SecirParams(), 0)
        graph.add_edge(0, 1, np.ones(8))
        self.assertEqual(graph.num_nodes, 2)
        self.assertEqual(graph.num_edges, 1)
        self.assertEqual(graph.get_num_out_edges(0), 1)
        self.assertEqual(graph.get_num_out_edges(1), 0)

    def test_migration_sim(self):        
        graph = secir.MigrationGraph()
        graph.add_node(0, secir.SecirParams(), 0, 0.1)
        graph.add_node(1, secir.SecirParams(), 0)
        graph.add_edge(0, 1, np.ones(8))

        sim = secir.MigrationSimulation(graph, t0 = 0.0)
        sim.advance(2)

        #integration does adaptive time steps so exact count is unknown
        self.assertGreaterEqual(sim.graph.get_node(0).property.result.get_num_time_points(), 3)

if __name__ == '__main__':
    unittest.main()
