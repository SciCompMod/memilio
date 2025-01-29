#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
from pyfakefs import fake_filesystem_unittest

from memilio.surrogatemodel.GNN import (
    data_generation_nodamp, data_generation_withdamp)
import memilio.simulation as mio
import memilio.simulation.osecir as osecir

from unittest.mock import patch
import os
import unittest

import numpy as np

import tensorflow as tf
tf.get_logger().setLevel('ERROR')


class TestSurrogatemodelGNN(fake_filesystem_unittest.TestCase):

    path = '/home/'

    num_groups = 6
    model = osecir.Model(num_groups)
    graph = osecir.ModelGraph()
    graph.add_node(0, model)
    graph.add_node(1, model)
    mobility_coefficients = 0.01 * np.ones(model.populations.numel())
    for i in range(num_groups):
        flat_index = model.populations.get_flat_index(
            osecir.MultiIndex_PopulationsArray(mio.AgeGroup(i), osecir.InfectionState.Dead))
        mobility_coefficients[flat_index] = 0
    graph.add_edge(0, 1, mobility_coefficients)
    graph.add_edge(1, 0, mobility_coefficients)

    def setUp(self):
        self.setUpPyfakefs()

#### test simulation no damp ####

    @patch('memilio.surrogatemodel.GNN.data_generation_nodamp.getBaselineMatrix',
           return_value=0.6 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.GNN.data_generation_nodamp.make_graph',
           return_value=graph)
    @patch('memilio.epidata.transformMobilityData.updateMobility2022')
    @patch('memilio.surrogatemodel.GNN.data_generation_nodamp.transform_mobility_directory',
           autospec=True)
    def test_simulation_run_nodamp(self, mock_transform_mobility,
                                   mock_update_mobility, mock_baseline, mock_graph):

        mock_transform_mobility.side_effect = lambda: None
        mock_update_mobility.side_effect = lambda directory, mobility_file: None

        days_1 = 10
        days_2 = 30
        days_3 = 50

        population = [[5256.0, 10551, 32368.5, 43637.833333333336, 22874.066666666666, 8473.6],
                      [4000, 8000, 40000, 28000, 15000, 6000]]

        # Call the actual function being tested
        simulation_1 = data_generation_nodamp.run_secir_groups_simulation(
            days_1, population)
        simulation_2 = data_generation_nodamp.run_secir_groups_simulation(
            days_2, population)
        simulation_3 = data_generation_nodamp.run_secir_groups_simulation(
            days_3, population)

        # Ensure that the results are of the correct length.
        self.assertEqual(len(simulation_1[0]), days_1 + 1)
        self.assertEqual(len(simulation_2[0]), days_2 + 1)
        self.assertEqual(len(simulation_3[0]), days_3 + 1)


#### test data genertion no damp ####


    @patch('memilio.surrogatemodel.GNN.data_generation_nodamp.getBaselineMatrix',
           return_value=0.6 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.GNN.data_generation_nodamp.make_graph',
           return_value=graph)
    @patch('memilio.surrogatemodel.GNN.data_generation_nodamp.get_population',
           return_value=np.random.randint(0, 700001, size=(400, 6)))
    @patch('memilio.epidata.transformMobilityData.updateMobility2022')
    @patch('memilio.surrogatemodel.GNN.data_generation_nodamp.transform_mobility_directory',
           autospec=True)
    def test_data_generation_runs_nodamp(
        self, mock_transform_mobility, mock_update_mobility,
            mock_baseline, mock_graph, mock_population):

        mock_transform_mobility.side_effect = lambda: None
        mock_update_mobility.side_effect = lambda directory, mobility_file: None

        input_width_1 = 1
        input_width_2 = 5

        label_width_1 = 1
        label_width_2 = 10

        num_runs_1 = 1
        num_runs_2 = 2

        # test data generation without dampings
        data_1 = data_generation_nodamp.generate_data(
            num_runs_1, self.path, input_width_1, label_width_1,
            save_data=False)
        self.assertEqual(len(data_1['inputs']), num_runs_1)
        self.assertEqual(len(data_1['inputs'][0]), input_width_1)
        self.assertEqual(len(data_1['inputs'][0][0]), 48)
        self.assertEqual(len(data_1['labels']), num_runs_1)
        self.assertEqual(len(data_1['labels'][0]), label_width_1)
        self.assertEqual(len(data_1['labels'][0][0]), 48)

        data_2 = data_generation_nodamp.generate_data(
            num_runs_2, self.path, input_width_2, label_width_2,
            save_data=False)
        self.assertEqual(len(data_2['inputs']), num_runs_2)
        self.assertEqual(len(data_2['inputs'][0]), input_width_2)
        self.assertEqual(len(data_2['inputs'][0][0]), 48)
        self.assertEqual(len(data_2['labels']), num_runs_2)
        self.assertEqual(len(data_2['labels'][0]), label_width_2)
        self.assertEqual(len(data_2['labels'][0][0]), 48)

    @patch('memilio.surrogatemodel.GNN.data_generation_nodamp.getBaselineMatrix',
           return_value=0.6 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.GNN.data_generation_nodamp.get_population',
           return_value=np.random.randint(0, 700001, size=(400, 6)))
    # create mock graph
    @patch('memilio.surrogatemodel.GNN.data_generation_nodamp.make_graph',
           return_value=graph)
    @patch('memilio.surrogatemodel.GNN.data_generation_nodamp.scale_data',
           autospec=True)
    @patch('memilio.epidata.transformMobilityData.updateMobility2022')
    @patch('memilio.surrogatemodel.GNN.data_generation_nodamp.transform_mobility_directory',
           autospec=True)
    def test_data_generation_save_nodamp(
            self, mock_transform_mobility, mock_update_mobility, mock_scale_data,
            mock_population, mock_baseline, mock_graph):

        mock_transform_mobility.side_effect = lambda: None
        mock_update_mobility.side_effect = lambda directory, mobility_file: None

        # Mock the return value of scale_data with dummy inputs and labels
        # Assuming scaled inputs and labels are 4D arrays based on the original function's output
        mock_scale_data.return_value = (np.random.rand(
            10, 8, 10, 6), np.random.rand(10, 8, 10, 6))

        input_width = 2
        label_width = 3
        num_runs = 1

        # Call the function being tested
        data_generation_nodamp.generate_data(num_runs, self.path, input_width,
                                             label_width)

        # Check the number of generated files
        self.assertEqual(len(os.listdir(self.path)), 1)

        # Check the contents of the directory
        self.assertEqual(os.listdir(self.path),
                         ['data_secir_age_groups.pickle'])

    @patch('memilio.surrogatemodel.GNN.data_generation_withdamp.getMinimumMatrix',
           return_value=0 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.GNN.data_generation_withdamp.getBaselineMatrix',
           return_value=0.6 * np.ones((6, 6)))
    # create mock graph
    @patch('memilio.surrogatemodel.GNN.data_generation_withdamp.make_graph',
           return_value=graph)
    @patch('memilio.epidata.transformMobilityData.updateMobility2022')
    @patch('memilio.surrogatemodel.GNN.data_generation_withdamp.transform_mobility_directory',
           autospec=True)
    def test_simulation_run_withdamp(self, mock_transform_mobility, mock_update_mobility,
                                     mock_minimum, mock_baseline, mock_graph):

        mock_transform_mobility.side_effect = lambda: None
        mock_update_mobility.side_effect = lambda directory, mobility_file: None

        days_1 = 10
        days_2 = 30
        days_3 = 50

        dampings1 = [5]
        dampings2 = [6, 15]
        dampings3 = [8, 18, 35]

        population = [[5256.0, 10551, 32368.5,
                      43637.833333333336, 22874.066666666666, 8473.6],
                      [4000, 8000, 40000,
                      28000, 15000, 6000]]

        dataset_entry1, damped_matrices1, num_damp1, damping_coeff1 = data_generation_withdamp.run_secir_groups_simulation(
            days_1, population, dampings1)
        dataset_entry2, damped_matrices2, num_damp2, damping_coeff2 = data_generation_withdamp.run_secir_groups_simulation(
            days_2, population, dampings2)
        dataset_entry3, damped_matrices3, num_damp3, damping_coeff3 = data_generation_withdamp.run_secir_groups_simulation(
            days_3, population, dampings3)

        # result length
        self.assertEqual(len(dataset_entry1[0]), days_1+1)
        self.assertEqual(len(dataset_entry2[0]), days_2+1)
        self.assertEqual(len(dataset_entry3[0]), days_3+1)

        baseline = data_generation_withdamp.getBaselineMatrix()
        # damping factor
        self.assertEqual(damped_matrices1[0].all(),
                         (baseline * damping_coeff1[0]).all())
        self.assertEqual(
            damped_matrices2[1].all(),
            (baseline * damping_coeff2[1]).all())
        self.assertEqual(
            damped_matrices3[2].all(),
            (baseline * damping_coeff3[2]).all())

        # number of dampings length
        self.assertEqual(len(damping_coeff1), len(dampings1))
        self.assertEqual(len(damping_coeff2), len(dampings2))
        self.assertEqual(len(damping_coeff3), len(dampings3))

# test data generation with dampings

    @patch('memilio.surrogatemodel.GNN.data_generation_withdamp.getMinimumMatrix',
           return_value=0 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.GNN.data_generation_withdamp.getBaselineMatrix',
           return_value=0.6 * np.ones((6, 6)))
    # create mock graph
    @patch('memilio.surrogatemodel.GNN.data_generation_withdamp.make_graph',
           return_value=graph)
    @patch('memilio.surrogatemodel.GNN.data_generation_withdamp.get_population',
           return_value=np.random.randint(0, 700001, size=(400, 6)))
    # mock transform directory function
    @patch('memilio.epidata.transformMobilityData.updateMobility2022')
    @patch('memilio.surrogatemodel.GNN.data_generation_withdamp.transform_mobility_directory', autospec=True)
    def test_data_generation_runs_withdamp(
            self, mock_transform_mobility, mock_update_mobility,
            mock_minimum, mock_baseline, mock_graph, mock_population):

        mock_transform_mobility.side_effect = lambda: None
        mock_update_mobility.side_effect = lambda directory, mobility_file: None

        input_width_1 = 1
        input_width_2 = 5

        label_width_1 = 10
        label_width_2 = 30

        num_runs_1 = 1
        num_runs_2 = 2

        damping1 = 1
        damping2 = 2

        data_1 = data_generation_withdamp.generate_data(
            num_runs_1, self.path, input_width_1, label_width_1,
            damping1, save_data=False)
        self.assertEqual(len(data_1['inputs']), num_runs_1)
        self.assertEqual(len(data_1['inputs'][0]), input_width_1)
        self.assertEqual(len(data_1['inputs'][0][0]), 48)
        self.assertEqual(len(data_1['labels']), num_runs_1)
        self.assertEqual(len(data_1['labels'][0]), label_width_1)
        self.assertEqual(len(data_1['labels'][0][0]), 48)

        data_2 = data_generation_withdamp.generate_data(
            num_runs_2, self.path, input_width_2, label_width_2,
            damping2, save_data=False)
        self.assertEqual(len(data_2['inputs']), num_runs_2)
        self.assertEqual(len(data_2['inputs'][0]), input_width_2)
        self.assertEqual(len(data_2['inputs'][0][0]), 48)
        self.assertEqual(len(data_2['labels']), num_runs_2)
        self.assertEqual(len(data_2['labels'][0]), label_width_2)
        self.assertEqual(len(data_2['labels'][0][0]), 48)

 # test saving for model with dampings
    @patch('memilio.surrogatemodel.GNN.data_generation_withdamp.getMinimumMatrix',
           return_value=0 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.GNN.data_generation_withdamp.getBaselineMatrix',
           return_value=0.6 * np.ones((6, 6)))
    # create mock graph
    @patch('memilio.surrogatemodel.GNN.data_generation_withdamp.make_graph',
           return_value=graph)
    @patch('memilio.surrogatemodel.GNN.data_generation_withdamp.get_population',
           return_value=np.random.randint(0, 700001, size=(400, 6)))
    @patch('memilio.surrogatemodel.GNN.data_generation_withdamp.scale_data',
           autospec=True)
    @patch('memilio.epidata.transformMobilityData.updateMobility2022')
    @patch('memilio.surrogatemodel.GNN.data_generation_withdamp.transform_mobility_directory',
           autospec=True)
    def test_data_generation_save_withdamp(
            self, mock_transform_mobility, mock_update_mobility, mock_scale_data,
            mock_population, mock_baseline, mcok_minimum, mock_graph):

        mock_transform_mobility.side_effect = lambda: None
        mock_update_mobility.side_effect = lambda directory, mobility_file: None

        mock_scale_data.return_value = (np.random.rand(
            10, 8, 10, 6), np.random.rand(10, 8, 10, 6))

        input_width = 5
        label_width = 20
        num_runs = 2
        num_dampings = 2
        data = data_generation_withdamp.generate_data(num_runs, self.path, input_width,
                                                      label_width, num_dampings, save_data=True)

        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(os.listdir(self.path),
                         ['data_secir_age_groups.pickle'])


if __name__ == '__main__':
    unittest.main()
