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


def mlp_multi_input_single_output(num_outputs=8):
    """! Simple MLP Network which takes the compartments for multiple time steps as input and returns the 8 compartments for one single time step.

    Reshaping adds an extra dimension to the output, so the shape of the output is 1x8. This makes the shape comparable to that of the multi-output models.
    @param num_outputs [Default: 8] Number of compartments. Default value is reached when aggregating the confirmed compartments.
    """
    model = tf.keras.Sequential([
        tf.keras.layers.Flatten(),
        tf.keras.layers.Dense(units=32, activation='relu'),
        tf.keras.layers.Dense(units=32, activation='relu'),
        tf.keras.layers.Dense(units=num_outputs),
        tf.keras.layers.Reshape([1, -1]), ])
    return model


def lstm_network_multi_input_single_output(num_outputs=8):
    """! LSTM Network which uses multiple time steps as input and returns the 8 compartments for one single time step in the future.

    Input and output have shape [number of expert model simulations, time points in simulation, number of individuals in infection states].

    @param num_outputs [Default: 8] Number of compartments. Default value is reached when aggregating the confirmed compartments.
    """
    model = tf.keras.models.Sequential([
        tf.keras.layers.LSTM(32, return_sequences=True),
        tf.keras.layers.Dense(units=num_outputs)
    ])
    return model


def cnn_multi_input_multi_output(label_width, conv_size=3, num_outputs=8):
    """! CNN Network which uses multiple time steps as input and returns the 8 compartments for multiple time step in the future.

    Input and output have shape [number of expert model simulations, time points in simulation, number of individuals in infection states].
    The parameter conv_size describes the kernel_size of the 1d Conv layer.
    We also use the parameter in combination with a lambda layer to transform the input to shape [batch, CONV_WIDTH, features].  


    @param label_width Number of time steps in the output.
    @param conv_size [Default: 3] Convolution kernel width which is 3 per default.
    @param num_outputs [Default: 8] Number of compartments. Default value is reached when aggregating the confirmed compartments.
    """

    model = tf.keras.Sequential([
        tf.keras.layers.Lambda(lambda x: x[:, -conv_size:, :]),
        tf.keras.layers.Conv1D(256, activation='relu',
                               kernel_size=(conv_size)),
        tf.keras.layers.Dense(label_width*num_outputs,
                              kernel_initializer=tf.initializers.zeros()),
        tf.keras.layers.Reshape([label_width, num_outputs])
    ])
    return model


def lstm_multi_input_multi_output(label_width, num_outputs=8):
    """! LSTM Network which uses multiple time steps as input and returns the 8 compartments for one single time step in the future.

    Input and output have shape [number of expert model simulations, time points in simulation, number of individuals in infection states].

    @param label_width Number of time steps in the output.
    @param num_outputs [Default: 8] Number of compartments. Default value is reached when aggregating the confirmed compartments.
    """
    model = tf.keras.Sequential([
        tf.keras.layers.LSTM(32, return_sequences=False),
        tf.keras.layers.Dense(label_width*num_outputs,
                              kernel_initializer=tf.initializers.zeros()),
        tf.keras.layers.Reshape([label_width, num_outputs])])
    return model
