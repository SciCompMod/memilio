#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
#
# Authors: Agatha Schmidt, Henrik Zunker
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
import tensorflow as tf


def mlp_multi_input_single_output():
    """! Simple MLP Network which takes the compartments for multiple time steps as input and returns the 8 compartments for one single time step.

    Reshaping adds an extra dimension to the output, so the shape of the output is 1x8. This makes the shape comparable to that of the multi-output models.
    """
    model = tf.keras.Sequential([
        tf.keras.layers.Flatten(),
        tf.keras.layers.Dense(units=32, activation='relu'),
        tf.keras.layers.Dense(units=32, activation='relu'),
        tf.keras.layers.Dense(units=8),
        tf.keras.layers.Reshape([1, -1]), ])
    return model


def lstm_network_multi_input_single_output():
    """! LSTM Network which uses multiple time steps as input and returns the 8 compartments for one single time step in the future.

    Input and output have shape [number of expert model simulations, time points in simulation, number of individuals in infection states].
    """
    model = tf.keras.models.Sequential([
        tf.keras.layers.LSTM(32, return_sequences=True),
        tf.keras.layers.Dense(units=8)
    ])
    return model


def cnn_multi_input_multi_output(label_width):
    """! CNN Network which uses multiple time steps as input and returns the 8 compartments for multiple time step in the future.

    Input and output have shape [batch, time, features].

    @param label_width Number of time steps in the output.
    """
    CONV_WIDTH = 3
    num_outputs = 8
    model = tf.keras.Sequential([
        # Lambda layer transforms input to shape [batch, CONV_WIDTH, features]
        tf.keras.layers.Lambda(lambda x: x[:, -CONV_WIDTH:, :]),
        tf.keras.layers.Conv1D(256, activation='relu',
                               kernel_size=(CONV_WIDTH)),
        tf.keras.layers.Dense(label_width*num_outputs,
                              kernel_initializer=tf.initializers.zeros()),
        tf.keras.layers.Reshape([label_width, num_outputs])
    ])
    return model


def lstm_multi_input_multi_output(label_width):
    """! LSTM Network which uses multiple time steps as input and returns the 8 compartments for one single time step in the future.

    Input and output have shape [batch, time, features].

    @param label_width Number of time steps in the output.
    """
    num_outputs = 8
    model = tf.keras.Sequential([
        tf.keras.layers.LSTM(32, return_sequences=False),
        tf.keras.layers.Dense(label_width*num_outputs,
                              kernel_initializer=tf.initializers.zeros()),
        tf.keras.layers.Reshape([label_width, num_outputs])])
    return model
