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
Grid search module for GNN-based surrogate model hyperparameter optimization.

This module provides functionality to perform systematic hyperparameter search
over GNN architectures, training configurations, and model parameters to identify
optimal configurations for epidemiological surrogate models.
"""

from pathlib import Path

import pandas as pd
import tensorflow as tf
from spektral.data import MixedLoader
from tensorflow.keras.losses import MeanAbsolutePercentageError
from tensorflow.keras.metrics import MeanAbsolutePercentageError as MeanAbsolutePercentageErrorMetric
from tensorflow.keras.optimizers import Adam

from memilio.surrogatemodel.GNN.evaluate_and_train import (
    create_dataset,
    train_and_evaluate
)
from memilio.surrogatemodel.GNN.network_architectures import get_model


# Default hyperparameter grid
DEFAULT_LAYER_TYPES = [
    "ARMAConv",
    "GCSConv",
    "GATConv",
    "GCNConv",
    "APPNPConv"
]

DEFAULT_NUM_LAYERS = [2, 3, 4, 5, 6, 7]
DEFAULT_NUM_CHANNELS = [2, 3, 4, 5, 6, 7]
DEFAULT_ACTIVATION_FUNCTIONS = ["elu", "relu", "tanh", "sigmoid"]

# Default training configuration
DEFAULT_BATCH_SIZE = 32
DEFAULT_MAX_EPOCHS = 100
DEFAULT_ES_PATIENCE = 30


def generate_parameter_grid(
        layer_types=None,
        num_layers_options=None,
        num_channels_options=None,
        activation_functions=None):
    """Generates a grid of model parameter combinations for hyperparameter search.

    :param layer_types: List of GNN layer types to test (default: DEFAULT_LAYER_TYPES).
    :param num_layers_options: List of layer counts to test (default: DEFAULT_NUM_LAYERS).
    :param num_channels_options: List of channel counts to test (default: DEFAULT_NUM_CHANNELS).
    :param activation_functions: List of activation functions to test (default: DEFAULT_ACTIVATION_FUNCTIONS).
    :returns: List of tuples, each containing (layer_type, num_layers, num_channels, activation).

    """
    if layer_types is None:
        layer_types = DEFAULT_LAYER_TYPES
    if num_layers_options is None:
        num_layers_options = DEFAULT_NUM_LAYERS
    if num_channels_options is None:
        num_channels_options = DEFAULT_NUM_CHANNELS
    if activation_functions is None:
        activation_functions = DEFAULT_ACTIVATION_FUNCTIONS

    parameter_grid = [
        (layer_type, num_layers, num_channels, activation)
        for layer_type in layer_types
        for num_layers in num_layers_options
        for num_channels in num_channels_options
        for activation in activation_functions
    ]

    return parameter_grid


def perform_grid_search(
        data,
        parameter_grid,
        save_dir,
        batch_size=DEFAULT_BATCH_SIZE,
        max_epochs=DEFAULT_MAX_EPOCHS,
        es_patience=DEFAULT_ES_PATIENCE,
        learning_rate=0.001):
    """Performs systematic grid search over GNN hyperparameters.

    Trains and evaluates models for each parameter combination in the grid,
    tracking performance metrics and saving results.
    Results are stored in a CSV file for later analysis.

    :param data: Spektral dataset containing training, validation, and test samples.
    :param parameter_grid: List of tuples with (layer_type, num_layers, num_channels, activation) combinations.
    :param save_dir: Directory to save results.
    :param batch_size: Batch size for training (default: 32).
    :param max_epochs: Maximum number of training epochs per configuration (default: 100).
    :param es_patience: Early stopping in epochs (default: 30).
    :param learning_rate: Learning rate for Adam optimizer (default: 0.001).
    :returns: DataFrame containing all grid search results.
    :raises ValueError: If data is empty or parameter_grid is invalid.

    """
    if not data or len(data) == 0:
        raise ValueError("Dataset must contain at least one sample.")

    if not parameter_grid or len(parameter_grid) == 0:
        raise ValueError(
            "Parameter grid must contain at least one configuration.")

    # Convert save_dir to Path if it's a string
    save_dir = Path(save_dir) / "saves"

    # Determine output dimension from data
    output_dim = data[0].y.shape[-1]

    # Initialize loss function and optimizer
    loss_function = MeanAbsolutePercentageError()

    # Initialize results DataFrame
    results_df = pd.DataFrame(
        columns=[
            'model',
            'optimizer',
            'number_of_hidden_layers',
            'number_of_channels',
            'activation',
            'mean_train_loss',
            'mean_validation_loss',
            'mean_test_loss',
            'mean_test_loss_orig',
            'training_time',
            'train_losses',
            'val_losses'
        ]
    )

    save_dir.mkdir(parents=True, exist_ok=True)
    results_file = save_dir / 'grid_search_results.csv'

    print(f"\n{'=' * 70}")
    print("GNN Grid Search - Hyperparameter Optimization")
    print(f"{'=' * 70}")
    print(f"Total configurations to test: {len(parameter_grid)}")
    print(f"Results will be saved to: {results_file}")
    print(f"{'=' * 70}\n")

    # Iterate through all parameter combinations
    for idx, (layer_type, num_layers, num_channels, activation) in enumerate(
            parameter_grid, 1):
        print(f"\n[{idx}/{len(parameter_grid)}] Training configuration:")
        print(f"  Layer type: {layer_type}")
        print(f"  Number of layers: {num_layers}")
        print(f"  Number of channels: {num_channels}")
        print(f"  Activation function: {activation}")

        try:
            # Create model instance
            model = get_model(
                layer_type=layer_type,
                num_layers=num_layers,
                num_channels=num_channels,
                activation=activation,
                num_output=output_dim
            )

            # Initialize optimizer for this configuration
            optimizer = Adam(learning_rate=learning_rate)

            # Build model by passing a sample batch through it
            build_loader = MixedLoader(
                data, batch_size=batch_size, epochs=1, shuffle=False)
            build_inputs, _ = next(build_loader)
            model(build_inputs)

            # Initialize optimizer variables
            optimizer.build(model.trainable_variables)

            model.compile(
                optimizer=optimizer,
                loss=MeanAbsolutePercentageError(),
                metrics=[MeanAbsolutePercentageErrorMetric()]
            )

            # Train and evaluate model
            results = train_and_evaluate(
                data=data,
                batch_size=batch_size,
                epochs=max_epochs,
                model=model,
                loss_fn=loss_function,
                optimizer=optimizer,
                es_patience=es_patience,
                save_name="",  # Dont save individual models during grid search
                save_dir=None  # Dont save individual results
            )

            # Store results
            results_df.loc[len(results_df)] = [
                layer_type,
                optimizer.__class__.__name__,
                num_layers,
                num_channels,
                activation,
                results["mean_train_loss"],
                results["mean_val_loss"],
                results["mean_test_loss"],
                results["mean_test_loss_orig"],
                results["training_time"],
                results["train_losses"],
                results["val_losses"]
            ]

            # Save intermediate results after each configuration
            results_df.to_csv(results_file, index=False)
            print(
                f"Configuration complete. Results saved to {results_file}")

        except Exception as e:
            print(f"Error training configuration: {e}")
            # Continue with next configuration rather than failing entire search

        finally:
            # Clear TensorFlow session to free memory
            tf.keras.backend.clear_session()

    print(f"\n{'=' * 70}")
    print("Grid Search Complete!")
    print(f"{'=' * 70}")
    print(f"Total configurations tested: {len(results_df)}")
    print(f"Results saved to: {results_file}")

    # Print best configuration
    if len(results_df) > 0:
        best_idx = results_df['mean_validation_loss'].idxmin()
        best_config = results_df.loc[best_idx]
        print(f"\nBest Configuration:")
        print(f"  Model: {best_config['model']}")
        print(f"  Layers: {best_config['number_of_hidden_layers']}")
        print(f"  Channels: {best_config['number_of_channels']}")
        print(f"  Activation: {best_config['activation']}")
        print(f"  Validation Loss: {best_config['mean_validation_loss']:.4f}")
        print(f"  Test Loss: {best_config['mean_test_loss']:.4f}")

    print(f"{'=' * 70}\n")

    return results_df


def main():
    """Main function demonstrating grid search usage.
    """
    # Dataset path
    dataset_path = Path.cwd() / "generated_datasets" / \
        "GNN_data_30days_3dampings_classic5.pickle"
    mobility_dir = Path.cwd() / "data" / "Germany" / "mobility"

    # Output directory for results
    output_dir = Path.cwd() / "grid_search_results"

    print("=" * 70)
    print("GNN Grid Search - Configuration")
    print("=" * 70)
    print(f"Dataset: {dataset_path}")
    print(f"Mobility data: {mobility_dir}")
    print(f"Output directory: {output_dir}")
    print("=" * 70)

    # Verify files exist
    if not dataset_path.exists():
        raise FileNotFoundError(f"Dataset not found: {dataset_path}")
    if not mobility_dir.exists():
        raise FileNotFoundError(
            f"Mobility directory not found: {mobility_dir}")

    # Load dataset
    print("\nLoading dataset...")
    data = create_dataset(str(dataset_path), str(mobility_dir))
    print(f"Dataset loaded: {len(data)} samples")

    # Generate parameter grid
    # For demonstration, use a smaller grid. Remove restrictions for full search.
    parameter_grid = generate_parameter_grid(
        layer_types=["ARMAConv", "GCNConv"],
        num_layers_options=[3, 5],
        num_channels_options=[4, 6],
        activation_functions=["elu", "relu"]
    )

    print(
        f"Generated parameter grid with {len(parameter_grid)} configurations")

    # Perform grid search
    results = perform_grid_search(
        data=data,
        parameter_grid=parameter_grid,
        save_dir=str(output_dir),
        batch_size=32,
        max_epochs=100,
        es_patience=30,
        learning_rate=0.001
    )

    print(f"\nGrid search complete. Results shape: {results.shape}")


if __name__ == "__main__":
    main()
