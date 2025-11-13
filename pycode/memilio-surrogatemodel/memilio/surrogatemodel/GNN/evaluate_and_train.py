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
Training and evaluation module for GNN-based surrogate models.

This module loads GNN training data, trains/evaluates surrogate models, and saves weights plus metrics.
"""

import pickle
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import spektral


from tensorflow.keras.optimizers import Adam
from tensorflow.keras.losses import MeanAbsolutePercentageError
import tensorflow.keras.initializers as initializers

import tensorflow as tf
import memilio.surrogatemodel.GNN.network_architectures as network_architectures
from memilio.surrogatemodel.utils.helper_functions import (calc_split_index)

from spektral.data import MixedLoader
from spektral.layers import ARMAConv
from spektral.utils.convolution import normalized_laplacian, rescale_laplacian


def _iterate_batches(loader):
    """Yield all batches from a Spektral loader for one epoch."""
    for _ in range(loader.steps_per_epoch):
        yield next(loader)


class StaticGraphDataset(spektral.data.Dataset):
    """Spektral dataset wrapper for samples that share a single adjacency matrix."""

    def __init__(
            self,
            node_features,
            node_labels,
            adjacency):
        self._node_features = node_features
        self._node_labels = node_labels
        self._adjacency = adjacency.astype(np.float32)
        super().__init__()
        # This must be set AFTER calling super().__init__()
        self.a = self._adjacency

    def read(self):
        """Create one Graph object per sample.

        Note: For MixedLoader, individual Graph objects should not have an 
        adjacency matrix (a). The adjacency matrix is stored at the dataset level.
        """
        return [
            spektral.data.Graph(
                x=feature.astype(np.float32),
                y=label.astype(np.float32),
            )
            for feature, label in zip(self._node_features, self._node_labels)
        ]


@dataclass
class TrainingSummary:
    """Container for aggregated training and evaluation metrics."""

    model_name: str
    mean_train_loss: float
    mean_val_loss: float
    mean_test_loss: float
    mean_test_loss_orig: float
    training_time: float
    train_losses: List[List[float]]
    val_losses: List[List[float]]


def load_gnn_dataset(
        dataset_path,
        mobility_dir,
        number_of_nodes=400,
        mobility_filename="commuter_mobility_2022.txt"):
    """Load serialized samples and mobility data into a Spektral dataset.

    Args:
        dataset_path: Pickle file containing the generated training samples.
        mobility_dir: Directory containing the commuter mobility file.
        number_of_nodes: Number of spatial nodes to retain from the mobility data.
        mobility_filename: Mobility file name to use.

    Returns:
        Spektral dataset with one `Graph` per sample sharing a common adjacency matrix.
    """
    dataset_path = Path(dataset_path)
    mobility_dir = Path(mobility_dir)

    if not dataset_path.exists():
        raise FileNotFoundError(f"Dataset not found: {dataset_path}")

    mobility_path = mobility_dir / mobility_filename
    if not mobility_path.exists():
        raise FileNotFoundError(f"Mobility file not found: {mobility_path}")

    with dataset_path.open("rb") as fp:
        data = pickle.load(fp)

    if "inputs" not in data or "labels" not in data:
        raise KeyError(
            f"Dataset at {dataset_path} must contain 'inputs' and 'labels' keys."
        )

    inputs = np.asarray(data["inputs"])
    labels = np.asarray(data["labels"])
    if inputs.shape[0] == 0:
        raise ValueError(
            "Loaded dataset is empty; expected at least one sample.")

    # Flatten temporal dimensions into feature vectors per node.
    num_samples, input_width, num_nodes, num_features = inputs.shape
    _, label_width, _, label_features = labels.shape

    if num_nodes != number_of_nodes:
        raise ValueError(
            f"Number of nodes in dataset ({num_nodes}) does not match expected "
            f"value ({number_of_nodes}).")

    node_features = inputs.transpose(0, 2, 1, 3).reshape(
        num_samples, number_of_nodes, input_width * num_features)
    node_labels = labels.transpose(0, 2, 1, 3).reshape(
        num_samples, number_of_nodes, label_width * label_features)

    commuter_data = pd.read_csv(mobility_path, sep=" ", header=None)
    adjacency_matrix = commuter_data.iloc[
        :number_of_nodes, :number_of_nodes
    ].to_numpy()
    adjacency_matrix = (adjacency_matrix > 0).astype(np.float32)
    adjacency_matrix = np.maximum(adjacency_matrix, adjacency_matrix.T)

    return StaticGraphDataset(node_features, node_labels, adjacency_matrix)


def create_dataset(path_cases, path_mobility, number_of_nodes=400):
    """Compatibility wrapper around `load_gnn_dataset`."""
    return load_gnn_dataset(
        Path(path_cases),
        Path(path_mobility),
        number_of_nodes=number_of_nodes)


def _train_step_impl(model, optimizer, loss_fn, inputs, target):
    """Perform one optimization step."""
    with tf.GradientTape() as tape:
        predictions = model(inputs, training=True)
        loss = loss_fn(target, predictions)
        if model.losses:
            loss += tf.add_n(model.losses)

    gradients = tape.gradient(loss, model.trainable_variables)
    optimizer.apply_gradients(zip(gradients, model.trainable_variables))

    metric = tf.reduce_mean(loss_fn(target, predictions))
    return loss, metric


def train_step(*args, **kwargs):
    """Wrapper that accepts both old and new call signatures."""
    if len(args) == 5 and not kwargs:
        first, second, third, fourth, fifth = args
        if isinstance(first, (tuple, list)):
            inputs, target = first, second
            loss_fn = third
            model = fourth
            optimizer = fifth
        else:
            model, optimizer, loss_fn, inputs, target = args
    else:
        raise TypeError(
            "train_step expects either (inputs, target, loss_fn, model, optimizer) "
            "or (model, optimizer, loss_fn, inputs, target).")

    return _train_step_impl(model, optimizer, loss_fn, inputs, target)


def evaluate(loader, model, loss_fn, retransform=False):
    """Evaluate the model on the dataset provided by loader."""
    total_loss = 0.0
    total_metric = 0.0
    total_samples = 0

    for inputs, target in _iterate_batches(loader):
        predictions = model(inputs, training=False)

        target_tensor = tf.convert_to_tensor(target, dtype=tf.float32)
        prediction_tensor = tf.cast(predictions, tf.float32)
        if retransform:
            target_tensor = tf.math.expm1(target_tensor)
            prediction_tensor = tf.math.expm1(prediction_tensor)

        batch_losses = loss_fn(target_tensor, prediction_tensor)
        batch_loss = tf.reduce_mean(batch_losses)
        batch_metric = tf.reduce_mean(batch_losses)

        batch_size_tensor = tf.cast(tf.shape(target_tensor)[0], tf.float32)
        batch_size = float(batch_size_tensor.numpy())
        total_loss += float(batch_loss.numpy()) * batch_size
        total_metric += float(batch_metric.numpy()) * batch_size
        total_samples += int(round(batch_size))

    if total_samples == 0:
        return 0.0, 0.0

    mean_loss = total_loss / total_samples
    mean_metric = total_metric / total_samples
    return mean_loss, mean_metric


def train_and_evaluate(
        data, batch_size, epochs, model, loss_fn, optimizer, es_patience,
        save_dir=None, save_name="model"):
    """Train the provided GNN model."""
    dataset_size = len(data)
    if dataset_size == 0:
        raise ValueError("Dataset must contain at least one sample.")

    n_train, n_valid, n_test = calc_split_index(
        dataset_size, split_train=0.7, split_valid=0.2, split_test=0.1)
    if n_train == 0 or n_valid == 0 or n_test == 0:
        raise ValueError(
            "Dataset split produced empty partitions. Provide a larger dataset "
            "or adjust the split configuration.")

    train_data = data[:n_train]
    valid_data = data[n_train:n_train + n_valid]
    test_data = data[n_train + n_valid:]

    # Build the model by passing a single batch through it.
    build_loader = MixedLoader(train_data, batch_size=min(
        batch_size, max(1, n_train)), epochs=1, shuffle=False)
    build_inputs, _ = next(build_loader)
    model(build_inputs)

    def _make_loader(dataset, *, batch_size, shuffle=False):
        return MixedLoader(
            dataset, batch_size=batch_size, epochs=1, shuffle=shuffle)

    best_val_loss = np.inf
    best_weights = model.get_weights()
    patience_counter = es_patience

    epoch_train_losses: List[float] = []
    epoch_val_losses: List[float] = []

    start_time = time.perf_counter()

    for epoch in range(1, epochs + 1):
        train_loader = _make_loader(
            train_data, batch_size=batch_size, shuffle=True)
        batch_losses = []
        for inputs, target in _iterate_batches(train_loader):
            loss, _ = train_step(model, optimizer, loss_fn, inputs, target)
            batch_losses.append(float(loss.numpy()))

        epoch_train_loss = float(np.mean(batch_losses)
                                 ) if batch_losses else 0.0
        epoch_train_losses.append(epoch_train_loss)

        val_loader = _make_loader(valid_data, batch_size=min(
            batch_size, max(1, n_valid)), shuffle=False)
        val_loss, _ = evaluate(val_loader, model, loss_fn)
        epoch_val_losses.append(val_loss)

        print(
            f"Epoch {epoch:02d} | train_loss={epoch_train_loss:.4f} "
            f"| val_loss={val_loss:.4f}"
        )

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_weights = model.get_weights()
            patience_counter = es_patience
            print(f"  ↳ New best validation loss: {best_val_loss:.4f}")
        else:
            patience_counter -= 1
            if patience_counter == 0:
                print("Early stopping triggered.")
                break

    elapsed = time.perf_counter() - start_time

    # Restore best weights and evaluate on the test set.
    model.set_weights(best_weights)
    test_loader = _make_loader(
        test_data, batch_size=min(batch_size, max(1, n_test)), shuffle=False)
    test_loss, _ = evaluate(test_loader, model, loss_fn)

    test_loader_retransform = _make_loader(
        test_data, batch_size=min(batch_size, max(1, n_test)), shuffle=False)
    test_loss_orig, _ = evaluate(
        test_loader_retransform, model, loss_fn, retransform=True)

    print(f"Test loss (log space): {test_loss:.4f}")
    print(f"Test loss (original scale): {test_loss_orig:.4f}")
    print(f"Training runtime: {elapsed:.2f}s ({elapsed / 60:.2f} min)")

    summary = TrainingSummary(
        model_name=save_name, mean_train_loss=float(
            np.min(epoch_train_losses))
        if epoch_train_losses else float("nan"),
        mean_val_loss=float(np.min(epoch_val_losses))
        if epoch_val_losses else float("nan"), mean_test_loss=float(test_loss),
        mean_test_loss_orig=float(test_loss_orig),
        training_time=elapsed / 60, train_losses=[epoch_train_losses],
        val_losses=[epoch_val_losses])

    if save_dir:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)

        metrics_df = pd.DataFrame(columns=[
            "train_loss", "val_loss", "test_loss",
            "test_loss_orig", "training_time",
            "loss_history", "val_loss_history"])
        metrics_df.loc[len(metrics_df.index)] = [
            summary.mean_train_loss,
            summary.mean_val_loss,
            summary.mean_test_loss,
            summary.mean_test_loss_orig,
            summary.training_time,
            summary.train_losses,
            summary.val_losses
        ]

        weights_dir = save_dir / "saved_weights"
        weights_dir.mkdir(parents=True, exist_ok=True)

        weights_filename = save_name if save_name.endswith(
            ".pickle") else f"{save_name}.pickle"
        weights_path = weights_dir / weights_filename
        with weights_path.open("wb") as fp:
            pickle.dump(best_weights, fp)

        results_dir = save_dir / "model_evaluations_paper"
        results_dir.mkdir(parents=True, exist_ok=True)
        results_path = results_dir / \
            weights_filename.replace(".pickle", ".csv")
        metrics_df.to_csv(results_path, index=False)
        print(f"Saved weights to {weights_path}")
        print(f"Saved evaluation metrics to {results_path}")

    return asdict(summary)


if __name__ == "__main__":
    start_hyper = time.perf_counter()
    epochs = 10
    batch_size = 2
    es_patience = 10
    optimizer = Adam(learning_rate=0.001)
    loss_fn = MeanAbsolutePercentageError()

    repo_root = Path(__file__).resolve().parents[4]
    artifacts_root = repo_root / "artifacts"

    dataset_path = artifacts_root / \
        "generated_datasets" / "GNN_data_30days_3dampings_classic5.pickle"

    mobility_dir = repo_root / "data" / "Germany" / "mobility"
    data = load_gnn_dataset(dataset_path, mobility_dir)

    # Define the model architecture
    def transform_a(adjacency_matrix):
        a = adjacency_matrix.numpy()
        a = rescale_laplacian(normalized_laplacian(a))
        return tf.convert_to_tensor(a, dtype=tf.float32)

    layer_types = [
        # Dense layer (only uses x)
        lambda: ARMAConv(512, activation='elu',
                         kernel_initializer=initializers.GlorotUniform(seed=None))
    ]
    num_repeat = [7]

    model_class = network_architectures.generate_model_class(
        "ARMA", layer_types, num_repeat, num_output=1440, transform=transform_a)

    model = model_class()

    save_name = 'GNN_30days'  # name for model
    save_dir = artifacts_root / "model_results"

    train_and_evaluate(
        data,
        batch_size,
        epochs,
        model,
        loss_fn,
        optimizer,
        es_patience,
        save_dir=save_dir,
        save_name=save_name)

    elapsed_hyper = time.perf_counter() - start_hyper
    print(
        "Time for hyperparameter testing: {:.4f} minutes".format(
            elapsed_hyper / 60))
    print(
        "Time for hyperparameter testing: {:.4f} hours".format(
            elapsed_hyper / 60 / 60))
