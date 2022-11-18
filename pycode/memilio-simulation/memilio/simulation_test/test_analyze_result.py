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

import memilio.simulation as mio
import memilio.simulation.secir as secir


class Test_AnalyzeResult(unittest.TestCase):
    def test_interpolate_time_series(self):
        ts = mio.TimeSeries(1)
        ts.add_time_point(0.0, np.r_[0.0])
        ts.add_time_point(0.5, np.r_[1.0])
        ts.add_time_point(1.5, np.r_[2.0])
        interpolated = secir.interpolate_simulation_result(ts)
        self.assertEqual(interpolated.get_num_time_points(), 2)
        self.assertEqual(interpolated.get_time(0), 0.0)
        self.assertEqual(interpolated.get_time(1), 1.0)

    def test_ensemble(self):
        ts = mio.TimeSeries(1)
        ts.add_time_point(0.0, np.r_[0.0])
        ts.add_time_point(0.5, np.r_[1.0])
        ts.add_time_point(1.5, np.r_[2.0])
        interpolated = secir.interpolate_ensemble_results([ts, ts])
        self.assertEqual(interpolated[1].get_num_time_points(), 2)
        self.assertEqual(interpolated[1].get_time(0), 0.0)
        self.assertEqual(interpolated[1].get_time(1), 1.0)

    def test_ensemble_graph(self):
        model = secir.Model(1)
        graph = secir.ModelGraph()
        graph.add_node(0, model)
        graph.add_node(1, model)
        graph.add_edge(0, 1, 0.01 * np.ones(8))
        graph.add_edge(1, 0, 0.01 * np.ones(8))

        study = secir.ParameterStudy(graph, t0=0, tmax=2, dt=0.5, num_runs=3)
        r = study.run()
        interpolated = secir.interpolate_ensemble_results(r)
        self.assertEqual(interpolated[0][0].get_time(0), 0.0)
        self.assertEqual(interpolated[0][0].get_time(1), 1.0)

    def test_mean(self):
        ts1 = mio.TimeSeries(1)
        ts1.add_time_point(0.0, np.r_[0.0])
        ts1.add_time_point(1.0, np.r_[1.0])
        ts2 = mio.TimeSeries(1)
        ts2.add_time_point(0.0, np.r_[1.0])
        ts2.add_time_point(1.0, np.r_[2.0])
        mean = secir.ensemble_mean([[ts1], [ts2]])
        self.assertEqual(mean[0][0], np.r_[0.5])
        self.assertEqual(mean[0][1], np.r_[1.5])

    def test_percentile(self):
        ts1 = mio.TimeSeries(1)
        ts1.add_time_point(0.0, np.r_[0.0])
        ts1.add_time_point(1.0, np.r_[1.0])
        ts2 = mio.TimeSeries(1)
        ts2.add_time_point(0.0, np.r_[1.0])
        ts2.add_time_point(1.0, np.r_[0.0])
        percentile = secir.ensemble_percentile([[ts1], [ts2]], 0.5)
        self.assertEqual(percentile[0][0], np.r_[1.0])
        self.assertEqual(percentile[0][1], np.r_[1.0])


if __name__ == '__main__':
    unittest.main()
