#############################################################################
# Copyright (C) 2020-2025 MEmilio
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


def mlp_multi_input_single_output(num_outputs=8, num_hidden_layers=3, num_neurons_per_layer = 32, activation = 'relu'):
    """ Simple MLP Network which takes the compartments for multiple time steps as input and returns the 8 compartments for one single time step.

    Reshaping adds an extra dimension to the output, so the shape of the output is 1x8. This makes the shape comparable to that of the multi-output models.

    :param num_outputs: Default: 8 Number of compartments. Default value is reached when aggregating the confirmed compartments.
    :param num_hidden_layers: Number of hidden dense layers in the MLP architecture 
    :param num_neurons_per_layer: Number of neurons per hidden layer 
    :param activation: name of the used activation function
    """
    #model = tf.keras.Sequential([
    #    tf.keras.layers.Flatten(),
    #    tf.keras.layers.Dense(units=32, activation='relu'),
    #    tf.keras.layers.Dense(units=32, activation='relu'),
    #    tf.keras.layers.Dense(units=num_outputs),
    #    tf.keras.layers.Reshape([1, -1]), ])
    #return model

    # Catching unallowed inputs
    if num_outputs < 1:
        raise ValueError("Output dimension must be at least 1, here %d"%(num_outputs))
    if num_hidden_layers < 0: 
        raise ValueError("Number of layers must be at least 0, here %d"%(num_hidden_layers))
    if num_neurons_per_layer < 1: 
        raise ValueError("Number of neurons per layer must be at least 1, here %d"%(num_neurons_per_layer))

    # Setting up basic model
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Flatten())

    # Adding new hidden layers 
    for _ in range(num_hidden_layers):
        model.add(tf.keras.layers.Dense(units = num_neurons_per_layer, activation = activation))
    
    # Adding output layer and reshaping output 
    model.add(tf.keras.layers.Dense(units = num_outputs))
    model.add(tf.keras.layers.Reshape([1,-1]))

    return model 

def mlp_multi_input_multi_output(label_width, num_outputs=8, num_hidden_layers=3, num_neurons_per_layer = 32, activation = 'relu'):
    """ Simple MLP Network which takes the compartments for multiple time steps as input and returns the 8 compartments for multiple time step.

    Reshaping adds an extra dimension to the output, so the shape of the output is 1x8. This makes the shape comparable to that of the multi-output models.

    :param label_width: Number of time steps in the output.
    :param num_outputs: Default: 8 Number of compartments. Default value is reached when aggregating the confirmed compartments.
    :param num_hidden_layers: Number of hidden dense layers in the MLP architecture 
    :param num_neurons_per_layer: Number of neurons per hidden layer 
    :param activation: name of the used activation function
    """
    # Catching unallowed inputs
    if label_width < 1: 
        raise ValueError("label width has to be a positive integer")
    if num_outputs < 1:
        raise ValueError("Output dimension must be at least 1, here %d"%(num_outputs))
    if num_hidden_layers < 0: 
        raise ValueError("Number of layers must be at least 0, here %d"%(num_hidden_layers))
    if num_neurons_per_layer < 1: 
        raise ValueError("Number of neurons per layer must be at least 1, here %d"%(num_neurons_per_layer))

    # Setting up basic model
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Flatten())

    # Adding new hidden layers 
    for _ in range(num_hidden_layers):
        model.add(tf.keras.layers.Dense(units = num_neurons_per_layer, activation = activation))
    
    # Adding output layer and reshaping output 
    model.add(tf.keras.layers.Dense(label_width*num_outputs))
    model.add(tf.keras.layers.Reshape([label_width, num_outputs]))

    return model 

def lstm_network_multi_input_single_output(num_outputs=8, internal_dimension = 32,  num_hidden_layers = 1, num_neurons_per_layer = 32, activation = "relu"):
    """ LSTM Network which uses multiple time steps as input and returns the 8 compartments for one single time step in the future.

    Input and output have shape [number of expert model simulations, time points in simulation, number of individuals in infection states].

    :param num_outputs: Default: 8 Number of compartments. Default value is reached when aggregating the confirmed compartments.
    :param internal_dimension: Output dimension of the LSTM-layer. 
    :param num_hidden_layers: Number of hidden layers in the dense network. 
    :param num_neurons_per_layer: Number of neurons per hidden layer.
    :param activation: Name of the used activation function. 
   
    """
    # Catching unallowed inputs
    if num_outputs < 1:
        raise ValueError("Output dimension must be at least 1, here %d"%(num_outputs))
    if internal_dimension < 1: 
        raise ValueError("Internal dimension must be at least 1, here %d"%(internal_dimension))
   
    # Defining final model starting with LSTM-layer 
    model = tf.keras.models.Sequential([
        tf.keras.layers.LSTM(internal_dimension, return_sequences=True)
        ])
    
    # Adding hidden layers in the dense network 
    for _ in range(num_hidden_layers):
        model.add(tf.keras.layers.Dense(units = num_neurons_per_layer, activation = activation))

    # Adding output layer and reshaping 
    model.add(tf.keras.layers.Dense(units=num_outputs))
    model.add(tf.keras.layers.Reshape([1,-1]))
    return model


def cnn_multi_input_multi_output(label_width, conv_size=3, num_outputs=8, num_filters = 256, num_hidden_layers = 1, num_neurons_per_layer = 256, activation = "relu"):
    """ CNN Network which uses multiple time steps as input and returns the 8 compartments for multiple time step in the future.

    Input and output have shape [number of expert model simulations, time points in simulation, number of individuals in infection states].
    The parameter conv_size describes the kernel_size of the 1d Conv layer.
    We also use the parameter in combination with a lambda layer to transform the input to shape [batch, CONV_WIDTH, features].

    :param label_width: Number of time steps in the output.
    :param conv_size: Default: 3 Convolution kernel width which is 3 per default.
    :param num_outputs: Default: 8 Number of compartments. Default value is reached when aggregating the confirmed compartments.
    :param num_filters: Number of different filters used in the Conv1D-Layer
    :param num_hidden_layers: Number of layers in the dense network following the convolution layer 
    :param num_neurons_per_layer: Number of neurons in each of the hidden layers (except the output layer)
    :param activation: activation function used in the hidden MLP-layers. 
    """

    #model = tf.keras.Sequential([
    #    tf.keras.layers.Lambda(lambda x: x[:, -conv_size:, :]),
    #    tf.keras.layers.Conv1D(256, activation='relu',
    #                           kernel_size=(conv_size)),
    #    tf.keras.layers.Dense(label_width*num_outputs,
    #                          kernel_initializer=tf.initializers.zeros()),
    #    tf.keras.layers.Reshape([label_width, num_outputs])
    #])
    #return model

    # Catching unallowed inputs
    if label_width < 1: 
        raise ValueError("label width has to be a positive integer")
    if conv_size < 2: 
        raise ValueError("Size of the convolution kernel has to be larger than 1, here %d"%(conv_size))
    if num_outputs < 1:
        raise ValueError("Output dimension must be at least 1, here %d"%(num_outputs))
    if num_hidden_layers < 0: 
        raise ValueError("Number of hidden layers must be at least 0, here %d"%(num_hidden_layers))
    if num_neurons_per_layer < 1: 
        raise ValueError("Number of neurons per layer must be at least 1, here %d"%(num_neurons_per_layer))

    # Defining convolutional layer 
    model = tf.keras.Sequential([
        tf.keras.layers.Lambda(lambda x: x[:, -conv_size:, :]),
        tf.keras.layers.Conv1D(num_filters, activation='relu',
                               kernel_size=(conv_size))
    ])
    # Constructing the hidden layers 
    for _ in range(num_hidden_layers):
        model.add(tf.keras.layers.Dense(units= num_neurons_per_layer, activation = activation))

    # Adding the output layer 
    model.add(tf.keras.layers.Dense(label_width*num_outputs,
                              kernel_initializer=tf.initializers.zeros()))
    # Adding reshaping layer 
    model.add(tf.keras.layers.Reshape([label_width, num_outputs]))
        
    return model


def lstm_multi_input_multi_output(label_width, num_outputs=8, internal_dimension = 32, num_hidden_layers = 1, num_neurons_per_layer = 32, activation = "relu"):
    """ LSTM Network which uses multiple time steps as input and returns the 8 compartments for one single time step in the future.

    Input and output have shape [number of expert model simulations, time points in simulation, number of individuals in infection states].

    :param label_width: Number of time steps in the output.
    :param num_outputs: Default: 8 Number of compartments. Default value is reached when aggregating the confirmed compartments.
    :param internal_dimension: Output dimension of the LSTM-layer. 
    :param num_hidden_layers: Number of hidden layers in the dense network 
    :param num_neurons_per_layer: Number of neurons per hidden layer 
    :param activation: Name of the used activation function
    """
    # Catching unallowed inputs
    if label_width < 1: 
        raise ValueError("label width has to be a positive integer")
    if num_outputs < 1:
        raise ValueError("Output dimension must be at least 1, here %d"%(num_outputs))
    if internal_dimension < 1: 
        raise ValueError("Internal dimension must be at least 1, here %d"%(internal_dimension))
    if num_hidden_layers < 0: 
        raise ValueError("Number of hidden layers must be at least 0, here %d"%(num_hidden_layers))
    if num_neurons_per_layer < 1: 
        raise ValueError("Number of neurons per layer must be at least 1, here %d"%(num_neurons_per_layer))

    # Defining the model with one LSTM layer 
    model = tf.keras.Sequential([
        tf.keras.layers.LSTM(internal_dimension, return_sequences=False)
        ])
    
    # Constructing the hidden layers 
    for _ in range(num_hidden_layers):
        model.add(tf.keras.layers.Dense(units= num_neurons_per_layer, activation = activation))

    # Adding output and reshape layer 
    model.add(tf.keras.layers.Dense(label_width*num_outputs,
                              kernel_initializer=tf.initializers.zeros()))
    model.add(tf.keras.layers.Reshape([label_width, num_outputs]))
    return model

