import unittest
import epidemiology.secir as secir
import numpy as np

class Test_AnalyzeResult(unittest.TestCase):
    def test_interpolate_time_series(self):
        ts = secir.TimeSeries(1)
        ts.add_time_point(0.0, np.r_[0.0])
        ts.add_time_point(0.5, np.r_[1.0])
        ts.add_time_point(1.5, np.r_[2.0])
        interpolated = secir.interpolate_simulation_result(ts)
        self.assertEqual(interpolated.get_num_time_points(), 3)
        self.assertEqual(interpolated.get_time(0), 0.0)
        self.assertEqual(interpolated.get_time(1), 1.0)
        self.assertEqual(interpolated.get_time(2), 2.0)

    def test_ensemble(self):
        ts = secir.TimeSeries(1)
        ts.add_time_point(0.0, np.r_[0.0])
        ts.add_time_point(0.5, np.r_[1.0])
        ts.add_time_point(1.5, np.r_[2.0])
        interpolated = secir.interpolate_ensemble_results([ts, ts])
        self.assertEqual(interpolated[1].get_num_time_points(), 3)
        self.assertEqual(interpolated[1].get_time(0), 0.0)
        self.assertEqual(interpolated[1].get_time(1), 1.0)
        self.assertEqual(interpolated[1].get_time(2), 2.0)
    
    def test_ensemble_graph(self):
        model = secir.SecirModel1()
        graph = secir.MigrationGraph1()
        graph.add_node(model)
        sim = secir.MigrationSimulation1(graph, t0 = 0, dt = 0.5)
        sim.advance(2)
        interpolated = secir.interpolate_ensemble_results([sim.graph, sim.graph])
        self.assertEqual(interpolated[0][0].get_time(0), 0.0)
        self.assertEqual(interpolated[0][0].get_time(1), 1.0)
        self.assertEqual(interpolated[0][0].get_time(2), 2.0)

    def test_mean(self):
        ts1 = secir.TimeSeries(1)
        ts1.add_time_point(0.0, np.r_[0.0])
        ts1.add_time_point(1.0, np.r_[1.0])
        ts2 = secir.TimeSeries(1)
        ts2.add_time_point(0.0, np.r_[1.0])
        ts2.add_time_point(1.0, np.r_[2.0])
        mean = secir.ensemble_mean([[ts1], [ts2]])
        self.assertEqual(mean[0][0], np.r_[0.5])
        self.assertEqual(mean[0][1], np.r_[1.5])

    def test_percentile(self):
        ts1 = secir.TimeSeries(1)
        ts1.add_time_point(0.0, np.r_[0.0])
        ts1.add_time_point(1.0, np.r_[1.0])
        ts2 = secir.TimeSeries(1)
        ts2.add_time_point(0.0, np.r_[1.0])
        ts2.add_time_point(1.0, np.r_[0.0])
        percentile = secir.ensemble_percentile([[ts1], [ts2]], 0.5)
        self.assertEqual(percentile[0][0], np.r_[1.0])
        self.assertEqual(percentile[0][1], np.r_[1.0])
        
if __name__ == '__main__':
    unittest.main()