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
import spektral.layers

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


if __name__ == "__main__":
    layer_types = [
        # Dense layer (only uses x)
        lambda: tf.keras.layers.Dense(10, activation="relu"),
        # Dense layer (only uses x)
        lambda: tf.keras.layers.Dense(20, activation="relu"),
        lambda: tf.keras.layers.Dense(30, activation="relu")
    ]
    num_repeat = [2, 3, 1]

    Karlheinz = generate_model_class(
        "Rosmarie", layer_types, num_repeat, num_output=2)
    model = Karlheinz()

    # Example inputs
    x = tf.random.normal([10, 20])  # Node features
    a = tf.random.normal([10, 10])  # Adjacency matrix
    labels = tf.random.normal([10, 2])  # Example labels

    model.compile(optimizer="adam", loss="mse")
    model.fit([x, a], labels, epochs=5)

    # Print trainable variables
    print(model.trainable_variables)
