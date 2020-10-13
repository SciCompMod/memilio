import unittest
from numpy.testing import assert_array_equal
import epidemiology.secir as secir
import numpy as np

class Test_TimeSeries(unittest.TestCase):
    def test_add_time_point(self):
        ts = secir.TimeSeries(1)
        ts.add_time_point(2, np.r_[1])
        self.assertEqual(ts.get_num_time_points(), 1)
        self.assertEqual(ts.get_time(0), 2)
        self.assertEqual(ts.get_value(0), np.r_[1])

        ts.add_time_point(3.5)
        self.assertEqual(ts.get_num_time_points(), 2)
        self.assertEqual(ts.get_last_time(), 3.5)
    
    def test_set_value(self):
        ts = secir.TimeSeries(1)
        ts.add_time_point(0.0, np.r_[1.0])
        ts.get_value(0)[:] = np.r_[2.0]
        self.assertEqual(ts.get_value(0), np.r_[2.0])

    def test_ndarray(self):
        ts = secir.TimeSeries(2)
        ts.add_time_point()
        ts.add_time_point()
        arr = ts.as_ndarray()
        self.assertEqual(arr.shape, (3, 2))
        arr[:, 1] = np.r_[1.0, 1.1, 1.2]
        assert_array_equal(ts.get_last_value(), np.r_[1.1, 1.2])
        assert_array_equal(ts.get_last_time(), 1.0)

if __name__ == '__main__':
    unittest.main()