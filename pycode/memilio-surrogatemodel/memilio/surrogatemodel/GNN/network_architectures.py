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
"""
Network architecture module for GNN-based surrogate models.

This file provides functionality to generate Graph Neural Network
architectures with various layer types, configurations, and preprocessing
transformations.
"""

import numpy as np
import tensorflow as tf
from spektral import layers as spektral_layers
from spektral.utils import convolution as spektral_convolution


# Supported GNN layer types
SUPPORTED_LAYER_TYPES = [
    "ARMAConv",
    "GCSConv",
    "GATConv",
    "GCNConv",
    "APPNPConv"
]


def _rescale_laplacian(laplacian):
    """Rescales a Laplacian matrix while remaining compatible with newer SciPy."""
    laplacian = laplacian.toarray() if hasattr(
        laplacian, "toarray") else np.asarray(laplacian)
    try:
        return spektral_convolution.rescale_laplacian(laplacian)
    except TypeError as err:
        if "eigvals" not in str(err):
            raise

    lmax = np.linalg.eigvalsh(laplacian)[-1]
    if lmax <= 0:
        lmax = 2.0
    identity = np.eye(laplacian.shape[0], dtype=laplacian.dtype)
    return (2.0 / lmax) * laplacian - identity


def _apply_graph_transform(adjacency, single_transform):
    """Applies a transform function to adjacency tensors with rank 2 or 3."""
    adj_array = adjacency.numpy()

    if adj_array.ndim == 2:
        transformed = np.asarray(single_transform(adj_array), dtype=np.float32)
        return tf.convert_to_tensor(transformed, dtype=tf.float32)

    if adj_array.ndim == 3:
        transformed = [
            np.asarray(single_transform(graph_adj), dtype=np.float32)
            for graph_adj in adj_array
        ]
        return tf.convert_to_tensor(np.stack(transformed), dtype=tf.float32)

    raise ValueError(
        f"Adjacency tensor must be rank 2 or 3, got rank {adj_array.ndim}."
    )


def generate_model_class(
        name,
        layer_types,
        num_repeat,
        num_output,
        transform=None):
    """Dynamically generates a custom Keras GNN model class with specified layer configuration.

    This function creates a new Keras Model class with a configurable sequence of layers.
    Each layer type can be repeated multiple times, allowing for flexible architecture design.
    The generated model supports both graph-based layers (Spektral) and standard layers.

    :param name: Name for the generated model class.
    :param layer_types: List of layer constructors.
                       Each element should be a callable that instantiates a layer.
    :param num_repeat: List of integers specifying repetition count for each layer type.
                      Must have same length as layer_types.
    :param num_output: Number of output units in the final dense layer.
    :param transform: Optional function to preprocess the adjacency matrix before passing
                     it to graph layers. Should accept a tensor and return a tensor.
    :returns: Created Keras Model class.
    :raises ValueError: If layer_types and num_repeat have different lengths, or if
                       any num_repeat value is less than 1.

    """

    def __init__(self):
        """Initializes the model with the specified layer sequence."""
        if len(layer_types) != len(num_repeat):
            raise ValueError(
                f"layer_types and num_repeat must have the same length. "
                f"Got {len(layer_types)} and {len(num_repeat)}."
            )
        if any(n < 1 for n in num_repeat):
            raise ValueError(
                "All values in num_repeat must be at least 1."
            )

        super(type(self), self).__init__()

        # Build sequence of hidden layers
        self.layer_seq = []
        for layer_idx, layer_type in enumerate(layer_types):
            for _ in range(num_repeat[layer_idx]):
                # Instantiate layer from callable
                layer = layer_type() if callable(layer_type) else layer_type
                self.layer_seq.append(layer)

        # Final output layer with ReLU activation
        self.output_layer = tf.keras.layers.Dense(
            num_output, activation="relu"
        )

    def call(self, inputs, mask=None):
        """Forward pass through the model.

        :param inputs: Tuple of (node_features, adjacency_matrix) where:
                      - node_features: [batch_size, num_nodes, num_features]
                      - adjacency_matrix: [num_nodes, num_nodes] or [batch_size, num_nodes, num_nodes]
        :param mask: Optional node mask from data loader. Can be None or a list [node_mask, None].
        :returns: Model output tensor with shape [batch_size, num_nodes, num_output].

        """
        # Unpack inputs
        x, a = inputs

        # Extract and prepare node mask for masking operations
        node_mask = None
        if mask is not None:
            # Spektral typically provides masks as [node_mask, None] for [x, a]
            node_mask = mask[0] if isinstance(
                mask, (tuple, list)) else mask

            if node_mask is not None:
                # Ensure mask has shape [batch_size, num_nodes, 1]
                if tf.rank(node_mask) == 2:
                    node_mask = tf.expand_dims(node_mask, axis=-1)

        # Create default mask if none provided
        if node_mask is None:
            # x has shape [batch_size, num_nodes, features]
            x_shape = tf.shape(x)
            node_mask = tf.ones([x_shape[0], x_shape[1], 1], dtype=tf.float32)

        # Apply adjacency matrix transformation if provided
        if not tf.is_symbolic_tensor(a):
            if transform is not None:
                a = transform(a)

        # Forward pass through layer sequence
        for layer in self.layer_seq:
            # Check if layer is a Spektral graph layer
            if type(layer).__module__.startswith("spektral.layers"):
                # Graph layers need both features and adjacency matrix
                x = layer([x, a], mask=[node_mask, None])
            else:
                # Standard layers only need features
                x = layer(x)

        # Apply final output layer
        output = self.output_layer(x)
        return output

    # Create class dictionary with methods
    class_dict = {
        '__init__': __init__,
        'call': call
    }

    return type(name, (tf.keras.Model,), class_dict)


def _get_layer_config(layer_type):
    """Returns layer class and transformation function for a given layer type.

    Each GNN layer type requires specific preprocessing of the adjacency matrix.
    This function encapsulates the layer-specific configuration.

    :param layer_type: String identifier for the GNN layer type.
    :returns: Tuple of (layer_class, transform_function).
    :raises ValueError: If layer_type is not supported.

    """
    if layer_type == "ARMAConv":
        layer_class = spektral_layers.ARMAConv

        def transform(adjacency):
            """Applies rescaled Laplacian transformation for ARMA convolution."""

            def single_transform(adj_array):
                laplacian = spektral_convolution.normalized_laplacian(
                    adj_array)
                return _rescale_laplacian(laplacian)

            return _apply_graph_transform(adjacency, single_transform)

        return layer_class, transform

    elif layer_type == "GCSConv":
        layer_class = spektral_layers.GCSConv

        def transform(adjacency):
            """Applies normalized adjacency for GCS convolution."""
            return _apply_graph_transform(
                adjacency,
                spektral_convolution.normalized_adjacency
            )

        return layer_class, transform

    elif layer_type == "GATConv":
        layer_class = spektral_layers.GATConv

        def transform(adjacency):
            """Applies normalized adjacency for GAT convolution."""
            return _apply_graph_transform(
                adjacency,
                spektral_convolution.normalized_adjacency
            )

        return layer_class, transform

    elif layer_type == "GCNConv":
        layer_class = spektral_layers.GCNConv

        def transform(adjacency):
            """Applies GCN filter for GCN convolution."""
            return _apply_graph_transform(
                adjacency,
                spektral_convolution.gcn_filter
            )

        return layer_class, transform

    elif layer_type == "APPNPConv":
        layer_class = spektral_layers.APPNPConv

        def transform(adjacency):
            """Applies GCN filter for APPNP convolution."""
            return _apply_graph_transform(
                adjacency,
                spektral_convolution.gcn_filter
            )

        return layer_class, transform

    else:
        raise ValueError(
            f"Unsupported layer_type: '{layer_type}'. "
            f"Supported types are: {', '.join(SUPPORTED_LAYER_TYPES)}"
        )


def get_model(
        layer_type,
        num_layers,
        num_channels,
        activation,
        num_output=1):
    """Generates a GNN model instance with specified architecture.

    Creates a Graph Neural Network model with a layer type repeated
    multiple times. The model includes appropriate preprocessing for the adjacency
    matrix based on the layer type.

    :param layer_type: Type of GNN layer to use. Must be one of: 'ARMAConv', 'GCSConv',
                      'GATConv', 'GCNConv', 'APPNPConv'.
    :param num_layers: Number of hidden GNN layers to stack.
    :param num_channels: Number of channels (units/features) in each hidden layer.
    :param activation: Activation function for hidden layers (e.g., 'relu', 'elu', 'tanh', 'sigmoid').
    :param num_output: Number of output units in the final dense layer (default: 1).
    :returns: Instantiated Keras Model ready for training.
    :raises ValueError: If parameters are invalid or layer_type is not supported.

    """
    # Validate inputs
    if layer_type not in SUPPORTED_LAYER_TYPES:
        raise ValueError(
            f"Unsupported layer_type: '{layer_type}'. "
            f"Supported types are: {', '.join(SUPPORTED_LAYER_TYPES)}"
        )

    if num_layers < 1:
        raise ValueError(
            f"num_layers must be at least 1, got {num_layers}."
        )

    if num_channels < 1:
        raise ValueError(
            f"num_channels must be at least 1, got {num_channels}."
        )

    if not isinstance(activation, str):
        raise ValueError(
            f"activation must be a string, got {type(activation).__name__}."
        )

    if num_output < 1:
        raise ValueError(
            f"num_output must be at least 1, got {num_output}."
        )

    # Get layer configuration for the specified type
    layer_class, transform_fn = _get_layer_config(layer_type)

    # Create layer constructor with specified parameters
    def create_layer():
        return layer_class(
            num_channels,
            activation=activation,
            kernel_initializer=tf.keras.initializers.GlorotUniform()
        )

    # Configure model structure
    layer_types = [create_layer]
    num_repeat = [num_layers]

    # Generate model class
    model_class = generate_model_class(
        name="CustomGNNModel",
        layer_types=layer_types,
        num_repeat=num_repeat,
        num_output=num_output,
        transform=transform_fn
    )

    return model_class()


def main():
    """Main function demonstrating network architecture usage.
    """
    print("=" * 70)
    print("GNN Network Architecture - Examples")
    print("=" * 70)

    # Example 1: Custom model with mixed layer types
    print("\nExample 1: Custom model with mixed Dense layers")
    layer_types = [
        lambda: tf.keras.layers.Dense(10, activation="relu"),
        lambda: tf.keras.layers.Dense(20, activation="relu"),
        lambda: tf.keras.layers.Dense(30, activation="relu")
    ]
    num_repeat = [2, 3, 1]  # 2 layers of 10 units, 3 of 20, 1 of 30

    custom_model_class = generate_model_class(
        name="CustomMixedModel",
        layer_types=layer_types,
        num_repeat=num_repeat,
        num_output=2
    )
    model1 = custom_model_class()

    # Example inputs for dense layers (no graph structure)
    batch_size = 8
    num_nodes = 20
    num_features = 5

    # Node features: [batch_size, num_nodes, features]
    x1 = tf.random.normal([batch_size, num_nodes, num_features])
    # Adjacency matrix per sample (not used by dense layers but required for input signature)
    a1_single = tf.random.normal([num_nodes, num_nodes])
    a1 = tf.tile(tf.expand_dims(a1_single, axis=0), [batch_size, 1, 1])
    # Labels: [batch_size, num_nodes, outputs]
    labels1 = tf.random.normal([batch_size, num_nodes, 2])

    print("Compiling and training custom model...")
    # Build the model by calling it once
    _ = model1([x1, a1])
    model1.compile(optimizer="adam", loss="mse")
    model1.fit([x1, a1], labels1, epochs=5, verbose=0)
    print("Custom model trained successfully")

    # Example 2: Pre-configured GNN model
    print("\nExample 2: Pre-configured ARMA GNN model")
    model2 = get_model(
        layer_type="ARMAConv",
        num_layers=3,
        num_channels=16,
        activation="relu",
        num_output=2
    )

    # Example inputs for GNN (proper graph structure)
    # Node features: [batch_size, num_nodes, features]
    x2 = tf.random.normal([batch_size, num_nodes, num_features])
    # Adjacency matrix per sample
    a2_single = tf.eye(
        num_nodes) + 0.1 * tf.random.normal([num_nodes, num_nodes])
    a2 = tf.tile(tf.expand_dims(a2_single, axis=0), [batch_size, 1, 1])
    # Labels: [batch_size, num_nodes, outputs]
    labels2 = tf.random.normal([batch_size, num_nodes, 2])

    print("Compiling and training ARMA model...")
    # Build the model by calling it once
    _ = model2([x2, a2])
    model2.compile(optimizer="adam", loss="mse")
    model2.fit([x2, a2], labels2, epochs=5, verbose=0)
    print("ARMA model trained successfully")

    print("\n" + "=" * 70)
    print("Examples completed successfully!")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
