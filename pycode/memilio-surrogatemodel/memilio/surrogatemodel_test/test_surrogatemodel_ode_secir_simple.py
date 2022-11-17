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
from pyfakefs import fake_filesystem_unittest
import os

from memilio.surrogatemodel.ode_secir_simple.data_generation import generate_data, run_secir_simulation


class TestSurrogatemodelOdeSecirSimple(fake_filesystem_unittest.TestCase):

    path = '/home/'

    def setUp(self):
        self.setUpPyfakefs()

    def test_simulation_run(self):

        days_1 = 10
        days_2 = 30
        days_3 = 50

        simulation_1 = run_secir_simulation(days_1)
        simulation_2 = run_secir_simulation(days_2)
        simulation_3 = run_secir_simulation(days_3)

        self.assertEqual(len(simulation_1), days_1+1)
        self.assertEqual(len(simulation_2), days_2+1)
        self.assertEqual(len(simulation_3), days_3+1)

    def test_data_generation_runs(self):

        input_width_1 = 1
        input_width_2 = 3

        label_width_1 = 1
        label_width_2 = 10

        num_runs_1 = 1
        num_runs_2 = 5

        data_1 = generate_data(num_runs_1, self.path, input_width_1,
                               label_width_1, save_data=False)
        self.assertEqual(len(data_1['inputs']), num_runs_1)
        self.assertEqual(len(data_1['inputs'][0]), input_width_1)
        self.assertEqual(len(data_1['inputs'][0][0]), 8)
        self.assertEqual(len(data_1['labels']), num_runs_1)
        self.assertEqual(len(data_1['labels'][0]), label_width_1)
        self.assertEqual(len(data_1['labels'][0][0]), 8)

        data_2 = generate_data(num_runs_2, self.path, input_width_2,
                               label_width_2, save_data=False)
        self.assertEqual(len(data_2['inputs']), num_runs_2)
        self.assertEqual(len(data_2['inputs'][0]), input_width_2)
        self.assertEqual(len(data_2['inputs'][0][0]), 8)
        self.assertEqual(len(data_2['labels']), num_runs_2)
        self.assertEqual(len(data_2['labels'][0]), label_width_2)
        self.assertEqual(len(data_2['labels'][0][0]), 8)

    def test_data_generation_save(self):

        input_width = 5
        label_width = 10
        num_runs = 15

        generate_data(num_runs, self.path, input_width,
                      label_width)
        self.assertEqual(len(os.listdir(self.path)), 1)

        self.assertEqual(os.listdir(self.path),
                         ['data_secir_simple.pickle'])


if __name__ == '__main__':
    unittest.main()
