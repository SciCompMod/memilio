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

from memilio.surrogatemodel.ode_secir_simple import data_generation
from memilio.surrogatemodel.ode_secir_simple import different_networks
from memilio.surrogatemodel.ode_secir_simple import model


class TestSurrogatemodelOdeSecirSimple(fake_filesystem_unittest.TestCase):

    path = '/home/'

    def setUp(self):
        self.setUpPyfakefs()

    def test_simulation_run(self):

        days_1 = 10
        days_2 = 30
        days_3 = 50

        simulation_1 = data_generation.run_secir_simulation(days_1)
        simulation_2 = data_generation.run_secir_simulation(days_2)
        simulation_3 = data_generation.run_secir_simulation(days_3)

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

        data_1 = data_generation.generate_data(
            num_runs_1, self.path, input_width_1, label_width_1,
            save_data=False)
        self.assertEqual(len(data_1['inputs']), num_runs_1)
        self.assertEqual(len(data_1['inputs'][0]), input_width_1)
        self.assertEqual(len(data_1['inputs'][0][0]), 8)
        self.assertEqual(len(data_1['labels']), num_runs_1)
        self.assertEqual(len(data_1['labels'][0]), label_width_1)
        self.assertEqual(len(data_1['labels'][0][0]), 8)

        data_2 = data_generation.generate_data(
            num_runs_2, self.path, input_width_2, label_width_2,
            save_data=False)
        self.assertEqual(len(data_2['inputs']), num_runs_2)
        self.assertEqual(len(data_2['inputs'][0]), input_width_2)
        self.assertEqual(len(data_2['inputs'][0][0]), 8)
        self.assertEqual(len(data_2['labels']), num_runs_2)
        self.assertEqual(len(data_2['labels'][0]), label_width_2)
        self.assertEqual(len(data_2['labels'][0][0]), 8)

    def test_data_generation_save(self):

        input_width = 2
        label_width = 3
        num_runs = 1

        data_generation.generate_data(num_runs, self.path, input_width,
                                      label_width)
        self.assertEqual(len(os.listdir(self.path)), 1)

        self.assertEqual(os.listdir(self.path),
                         ['data_secir_simple.pickle'])

    def test_network_fit(self):

        max_epochs = 5

        # models with single output
        model_single = different_networks.single_output()
        model_dense_multi_input = different_networks.multilayer_multi_input()
        model_lstm_single = different_networks.lstm_network_multi_input()

        # no existing dataset
        with self.assertRaises(FileNotFoundError) as error:
            model.network_fit(
                self.path, model=model_single, max_epochs=max_epochs)
        error_message = "[Errno 2] No such file or directory in the fake filesystem: '" + \
            self.path + "data_secir_simple.pickle'"
        self.assertEqual(str(error.exception), error_message)

        # generate dataset with single output
        input_width = 5
        label_width = 1
        num_runs = 15
        data_generation.generate_data(num_runs, self.path, input_width,
                                      label_width)

        # test different models
        single_output = model.network_fit(
            self.path, model=model_single, max_epochs=max_epochs, plot=False)
        self.assertEqual(
            len(single_output.history['val_loss']), max_epochs)
        dense_output = model.network_fit(
            self.path, model=model_dense_multi_input, max_epochs=max_epochs,
            plot=False)
        self.assertEqual(
            len(dense_output.history['val_loss']), max_epochs)
        lstm_single_output = model.network_fit(
            self.path, model=model_lstm_single, max_epochs=max_epochs, plot=False)
        self.assertEqual(
            len(lstm_single_output.history['val_loss']), max_epochs)

        # generate dataset with multiple output
        input_width = 5
        label_width = 10
        num_runs = 25
        data_generation.generate_data(num_runs, self.path, input_width,
                                      label_width)

        # models with multiple outputs
        model_cnn = different_networks.cnn_multi_output(label_width)
        model_lstm = different_networks.lstm_multi_output(label_width)

        cnn_output = model.network_fit(
            self.path, model=model_cnn, max_epochs=max_epochs, plot=False)
        self.assertEqual(
            cnn_output.model.output_shape[1], label_width)
        self.assertEqual(
            len(cnn_output.history['val_loss']), max_epochs)

        lstm_output = model.network_fit(
            self.path, model=model_lstm, max_epochs=max_epochs, plot=False)
        self.assertEqual(
            lstm_output.model.output_shape[1], label_width)
        self.assertEqual(
            len(lstm_output.history['val_loss']), max_epochs)


if __name__ == '__main__':
    unittest.main()
