import os
import pickle

import pandas as pd
import numpy as np
from sklearn.preprocessing import FunctionTransformer
import tensorflow as tf
from tensorflow.keras.layers import Dense
from tensorflow.keras.losses import MeanAbsolutePercentageError
from tensorflow.keras.metrics import mean_absolute_percentage_error
from tensorflow.keras.models import Model
import tensorflow.keras.initializers as initializers

from tensorflow.keras.optimizers import Adam, Nadam, RMSprop, SGD, Adagrad
import spektral.layers as spektral_layers
import spektral.utils.convolution as spektral_convolution


from sklearn.model_selection import KFold

# from spektral.data import Dataset, DisjointLoader, Graph, Loader, BatchLoader, MixedLoader
# from spektral.layers import GCSConv, GlobalAvgPool, GlobalAttentionPool, ARMAConv, AGNNConv, APPNPConv, CrystalConv, GATConv, GINConv, XENetConv, GCNConv, GCSConv
# from spektral.transforms.normalize_adj import NormalizeAdj
# from spektral.utils.convolution import gcn_filter, normalized_laplacian, rescale_laplacian, normalized_adjacency


# Function to generate a model class dynamically

def generate_model_class(name, layer_types, num_repeat, num_output, transform=None):
    def __init__(self):
        super(type(self), self).__init__()
        self.layer_seq = []
        for i, layer_type in enumerate(layer_types):
            for _ in range(num_repeat[i]):
                layer = layer_type() if callable(layer_type) else layer_type
                self.layer_seq.append(layer)
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

    layer_types = [lambda: layer_name(num_channels, activation=activation)]
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
