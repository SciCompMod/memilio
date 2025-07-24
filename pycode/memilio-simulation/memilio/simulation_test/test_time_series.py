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
import os
import tempfile

import numpy as np
from numpy.testing import assert_array_equal

import memilio.simulation as mio


class Test_TimeSeries(unittest.TestCase):
    """ """

    def test_add_time_point(self):
        """ """
        ts = mio.TimeSeries(1)
        ts.add_time_point(2, np.r_[1])
        self.assertEqual(ts.get_num_time_points(), 1)
        self.assertEqual(ts.get_time(0), 2)
        self.assertEqual(ts.get_value(0), np.r_[1])

        ts.add_time_point(3.5)
        self.assertEqual(ts.get_num_time_points(), 2)
        self.assertEqual(ts.get_last_time(), 3.5)

    def test_set_value(self):
        """ """
        ts = mio.TimeSeries(1)
        ts.add_time_point(0.0, np.r_[1.0])
        ts.get_value(0)[:] = np.r_[2.0]
        self.assertEqual(ts.get_value(0), np.r_[2.0])

    def test_ndarray(self):
        """ """
        ts = mio.TimeSeries(2)
        ts.add_time_point()
        ts.add_time_point()
        arr = ts.as_ndarray()
        self.assertEqual(arr.shape, (3, 2))
        arr[:, 1] = np.r_[1.0, 1.1, 1.2]
        assert_array_equal(ts.get_last_value(), np.r_[1.1, 1.2])
        assert_array_equal(ts.get_last_time(), 1.0)

    def test_print_table(self):
        """ """
        ts = mio.TimeSeries(1)
        ts.add_time_point(2, np.r_[1])
        ts.add_time_point(3.5, np.r_[2])
        output = ts.print_table(True, ["a", "b"], 2, 2)
        self.assertEqual(
            output, '\nTime a \n2.00 1.00\n3.50 2.00\n')

    def test_print_table_with_separator(self):
        """Test print_table with custom separator"""
        ts = mio.TimeSeries(1)
        ts.add_time_point(2, np.r_[1])
        ts.add_time_point(3.5, np.r_[2])
        output = ts.print_table(True, ["a"], 4, 1, ',', "# ")
        self.assertEqual(
            output, '# Time,a   \n 2.0, 1.0\n 3.5, 2.0\n')

    def test_print_table_console_output(self):
        """Test the new print_table function that prints directly to console"""
        ts = mio.TimeSeries(1)
        ts.add_time_point(2, np.r_[1])
        ts.add_time_point(3.5, np.r_[2])

        # Test that print_table() without return_string parameter returns None
        # (meaning it prints directly to console)
        result = ts.print_table()
        self.assertIsNone(result)

        # Test that print_table() with custom parameters still prints to console
        result = ts.print_table(["Column1"], 10, 2)
        self.assertIsNone(result)

    def test_export_csv(self):
        """Test export_csv functionality"""

        ts = mio.TimeSeries(2)
        ts.add_time_point(1.0, np.r_[10.0, 20.0])
        ts.add_time_point(2.5, np.r_[15.0, 25.0])

        # Create a temporary file for testing
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as temp:
            temp_filename = temp.name

        # Export to the temporary CSV file
        ts.export_csv(temp_filename, ["Col1", "Col2"], ',', 2)

        # Read and verify the file contents
        with open(temp_filename) as f:
            content = f.read()

        expected = "Time,Col1,Col2\n1.00,10.00,20.00\n2.50,15.00,25.00\n"
        self.assertEqual(content, expected)

        # Clean up the temporary file
        if os.path.exists(temp_filename):
            os.remove(temp_filename)


if __name__ == '__main__':
    unittest.main()
