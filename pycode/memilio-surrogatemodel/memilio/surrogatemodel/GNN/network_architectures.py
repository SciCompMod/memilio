#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Agatha Schmidt, Henrik Zunker, Manuel Heger
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
import numpy as np
import tensorflow as tf

import spektral.layers as spektral_layers
import spektral.utils.convolution as spektral_convolution


# Function to generate a model class dynamically

def generate_model_class(name, layer_types, num_repeat, num_output, transform=None):
    '''Generates a custom Keras model class with specified layers.

    :param name: Name of the generated class.
    :param layer_types: List of layer types (classes or callables) to include in the model.
    :param num_repeat: List of integers specifying how many times to repeat each layer type.
    :param num_output: Number of output units in the final layer.
    :param transform: Optional function to transform the adjacency matrix before passing it to layers.
    :return: A dynamically created Keras model class.'''

    def __init__(self):
        if len(layer_types) != len(num_repeat):
            raise ValueError(
                "layer_types and num_repeat must have the same length.")
        if any(n < 1 for n in num_repeat):
            raise ValueError("All values in num_repeat must be at least 1.")

        super(type(self), self).__init__()
        self.layer_seq = []
        for i, layer_type in enumerate(layer_types):
            for _ in range(num_repeat[i]):
                layer = layer_type() if callable(layer_type) else layer_type
                self.layer_seq.append(layer)
        # Final output layer
        self.output_layer = tf.keras.layers.Dense(
            num_output, activation="relu")

    def call(self, inputs):
        # Input consists of a tuple (x, a) where x is the node features and a is the adjacency matrix
        x, a = inputs

        # Apply transformation on adjacency matrix if provided
        if not tf.is_symbolic_tensor(a):
            if transform is not None:
                a = transform(a)

        # Pass through the layers
        for layer in self.layer_seq:
            if type(layer).__module__.startswith("spektral.layers"):
                x = layer([x, a])  # Pass both `x` and `a` to the layer
            else:
                x = layer(x)  # Pass only `x` to the layer

        output = self.output_layer(x)
        return output

    # Define the methods
    class_dict = {
        '__init__': __init__,
        'call': call
    }
    return type(name, (tf.keras.Model,), class_dict)


def get_model(layer_type, num_layers, num_channels, activation, num_output=1):
    """Generates a GNN model based on the specified parameters. 
    :param layer_type: Name of GNN layer to use (possible  'ARMAConv', 'GCSConv', 'GATConv', 
            'GCNConv', 'APPNPConv'), provided as a string.
    :param num_layers: Number of hidden layers in the model.
    :param num_channels: Number of channels (units) in each hidden layer.
    :param activation: Activation function to use in the hidden layers (e.g., 'relu', 'elu', 'tanh', 'sigmoid').
    :param num_output: Number of output units in the final layer.
    :return: A Keras model instance with the specified architecture.
    """
    if layer_type not in ["ARMAConv", "GCSConv", "GATConv", "GCNConv", "APPNPConv"]:
        raise ValueError(
            f"Unsupported layer_type: {layer_type}. Supported types are 'ARMAConv', 'GCSConv', 'GATConv', 'GCNConv', 'APPNPConv'.")
    if num_layers < 1:
        raise ValueError("num_layers must be at least 1.")
    if num_channels < 1:
        raise ValueError("num_channels must be at least 1.")
    if not isinstance(activation, str):
        raise ValueError(
            "activation must be a string representing the activation function.")
    if num_output < 1:
        raise ValueError("num_output must be at least 1.")

    # Define the layer based on the specified type
    if layer_type == "ARMAConv":
        layer_name = spektral_layers.ARMAConv

        def transform(a):
            a = np.array(a)
            a = spektral_convolution.rescale_laplacian(
                spektral_convolution.normalized_laplacian(a))
            return tf.convert_to_tensor(a, dtype=tf.float32)
    elif layer_type == "GCSConv":
        layer_name = spektral_layers.GCSConv

        def transform(a):
            a = a.numpy()
            a = spektral_convolution.normalized_adjacency(a)
            return tf.convert_to_tensor(a, dtype=tf.float32)
    elif layer_type == "GATConv":
        layer_name = spektral_layers.GATConv

        def transform(a):
            a = a.numpy()
            a = spektral_convolution.normalized_adjacency(a)
            return tf.convert_to_tensor(a, dtype=tf.float32)
    elif layer_type == "GCNConv":
        layer_name = spektral_layers.GCNConv

        def transform(a):
            a = a.numpy()
            a = spektral_convolution.gcn_filter(a)
            return tf.convert_to_tensor(a, dtype=tf.float32)
    elif layer_type == "APPNPConv":
        layer_name = spektral_layers.APPNPConv

        def transform(a):
            a = a.numpy()
            a = spektral_convolution.gcn_filter(a)
            return tf.convert_to_tensor(a, dtype=tf.float32)
    # Generate the input for generate_model_class
    layer_types = [lambda: layer_name(
        num_channels, activation=activation, kernel_initializer=tf.keras.initializers.GlorotUniform())]
    num_repeat = [num_layers]
    model_class = generate_model_class(
        "CustomModel", layer_types, num_repeat, num_output, transform=transform)
    return model_class()


if __name__ == "__main__":
    layer_types = [
        # Dense layer (only uses x)
        lambda: tf.keras.layers.Dense(10, activation="relu"),
        # Dense layer (only uses x)
        lambda: tf.keras.layers.Dense(20, activation="relu"),
        lambda: tf.keras.layers.Dense(30, activation="relu")
    ]
    num_repeat = [2, 3, 1]

    testclass = generate_model_class(
        "testclass", layer_types, num_repeat, num_output=2)
    model = testclass()

    # Example inputs
    x = tf.random.normal([10, 20])  # Node features
    a = tf.random.normal([10, 10])  # Adjacency matrix
    labels = tf.random.normal([10, 2])  # Example labels

    model.compile(optimizer="adam", loss="mse")
    model.fit([x, a], labels, epochs=5)

    model2 = get_model("ARMAConv", 3, 16, "relu", num_output=2)

    model2.compile(optimizer="adam", loss="mse")
    model2.fit([x, a], labels, epochs=5)
