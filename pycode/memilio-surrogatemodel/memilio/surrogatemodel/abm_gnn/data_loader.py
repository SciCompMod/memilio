#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Henrik Zunker
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
Data loader for ABM-generated GNN training data.

Reads the output of the C++ graph_abm_data_generation executable
and formats it for GNN surrogate model training.

Output format is compatible with the existing GNN training pipeline
"""

import glob
import os
import pickle
from typing import Dict, List, Optional, Tuple

import numpy as np


# Constants matching the C++ data generator
NUM_AGE_GROUPS = 6
NUM_INFECTION_STATES = 8  # S, E, INS, ISy, ISev, ICrit, R, D
NUM_FEATURES_PER_NODE = NUM_AGE_GROUPS * NUM_INFECTION_STATES  # 48

# Calibration parameter names must match the order in graph_abm_data_generation.cpp
CALIB_PARAM_NAMES = [
    "transmission_rate_scale",
    "contact_rate_work",
    "contact_rate_school",
    "contact_rate_social",
    "contact_rate_shop",
    "initial_infected_fraction",
    "initial_recovered_fraction",
    "commuter_fraction",
]
NUM_CALIB_PARAMS = len(CALIB_PARAM_NAMES)


def load_single_run(run_dir: str) -> Optional[Dict]:
    """Loads data from a single ABM simulation run.

    :param run_dir: Path to a run directory (e.g., 'abm_gnn_data/run_0000/').
    :returns: Dictionary with keys:
        - 'timeseries': np.ndarray of shape [num_timesteps, num_nodes, 48]
        - 'adjacency': np.ndarray of shape [num_nodes, num_nodes]
        - 'dampings': list of (day, severity) tuples
        - 'metadata': dict of run configuration
        Returns None if the run directory is incomplete.
    """
    if not os.path.isdir(run_dir):
        return None

    # Load adjacency matrix
    adj_path = os.path.join(run_dir, "adjacency.csv")
    if not os.path.exists(adj_path):
        return None
    adjacency = np.loadtxt(adj_path, delimiter=",")

    num_nodes = adjacency.shape[0]

    # Load per-node time series
    node_timeseries = []
    for node_idx in range(num_nodes):
        ts_path = os.path.join(
            run_dir, f"timeseries_node_{node_idx:03d}.csv")
        if not os.path.exists(ts_path):
            return None
        # Skip header row, first column is time
        data = np.loadtxt(ts_path, delimiter=",", skiprows=1)
        times = data[:, 0]
        values = data[:, 1:]  # shape: [num_timesteps, 48]
        node_timeseries.append(values)

    # Stack: [num_timesteps, num_nodes, 48]
    timeseries = np.stack(node_timeseries, axis=1)

    # Load dampings
    dampings = []
    damp_path = os.path.join(run_dir, "dampings.csv")
    if os.path.exists(damp_path):
        damp_data = np.loadtxt(damp_path, delimiter=",", skiprows=1)
        if damp_data.ndim == 1 and damp_data.size > 0:
            damp_data = damp_data.reshape(1, -1)
        if damp_data.size > 0:
            dampings = [(int(row[0]), float(row[1])) for row in damp_data]

    # Load metadata
    metadata = {}
    meta_path = os.path.join(run_dir, "metadata.txt")
    if os.path.exists(meta_path):
        with open(meta_path) as f:
            for line in f:
                line = line.strip()
                if "=" in line:
                    key, val = line.split("=", 1)
                    metadata[key] = val

    # Load sampled calibration parameters (the GNN input vector)
    parameters = None
    param_path = os.path.join(run_dir, "parameters.csv")
    if os.path.exists(param_path):
        param_dict = {}
        with open(param_path) as f:
            next(f)  # skip header
            for line in f:
                line = line.strip()
                if not line:
                    continue
                name, val = line.split(",", 1)
                param_dict[name] = float(val)
        parameters = np.array(
            [param_dict.get(n, np.nan) for n in CALIB_PARAM_NAMES],
            dtype=np.float32,
        )

    return {
        "timeseries": timeseries,
        "times": times,
        "adjacency": adjacency,
        "dampings": dampings,
        "parameters": parameters,
        "metadata": metadata,
    }


def load_dataset(
    data_dir: str,
    input_width: int = 5,
    label_width: int = 25,
    max_runs: Optional[int] = None,
    transform: bool = True,
) -> Dict[str, np.ndarray]:
    """Loads all ABM simulation runs and formats for GNN training.

    Splits each run's time series into input (first input_width days)
    and label (next label_width days) windows.

    :param data_dir: Root directory containing run_XXXX subdirectories.
    :param input_width: Number of input time steps (default: 5).
    :param label_width: Number of prediction time steps (default: 25).
    :param max_runs: Maximum number of runs to load (None = all).
    :param transform: Whether to apply log(1+x) transformation.
    :returns: Dictionary with keys:
        - 'inputs': [num_runs, input_width, num_nodes, 48]
        - 'labels': [num_runs, label_width, num_nodes, 48]
        - 'parameters': [num_runs, NUM_CALIB_PARAMS] sampled calibration vectors
        - 'adjacency': [num_nodes, num_nodes] (from first run, assumed same)
        - 'dampings': list of damping info per run
    """
    # Find all run directories
    run_dirs = sorted(glob.glob(os.path.join(data_dir, "run_*")))

    if max_runs is not None:
        run_dirs = run_dirs[:max_runs]

    if not run_dirs:
        raise FileNotFoundError(
            f"No run directories found in {data_dir}")

    total_width = input_width + label_width
    inputs_list = []
    labels_list = []
    dampings_list = []
    parameters_list = []
    adjacency = None

    loaded = 0
    skipped = 0

    for run_dir in run_dirs:
        run_data = load_single_run(run_dir)
        if run_data is None:
            skipped += 1
            continue

        ts = run_data["timeseries"]  # [T, N, 48]
        num_timesteps = ts.shape[0]

        if adjacency is None:
            adjacency = run_data["adjacency"]

        # We need to sample daily values from the (potentially sub-daily) ABM output.
        # The ABM logs at exchange intervals (e.g., every 12h).
        # Resample to daily by taking one sample per day.
        times = run_data["times"]
        daily_indices = _resample_to_daily(times, total_width)

        if daily_indices is None or len(daily_indices) < total_width:
            skipped += 1
            continue

        ts_daily = ts[daily_indices]  # [total_width, N, 48]

        inputs_list.append(ts_daily[:input_width])
        labels_list.append(ts_daily[input_width:total_width])
        dampings_list.append(run_data["dampings"])
        parameters_list.append(
            run_data["parameters"]
            if run_data["parameters"] is not None
            else np.full(NUM_CALIB_PARAMS, np.nan, dtype=np.float32)
        )
        loaded += 1

    if loaded == 0:
        raise RuntimeError(
            f"No valid runs loaded from {data_dir} "
            f"(skipped {skipped})")

    inputs = np.array(inputs_list)            # [R, input_width, N, 48]
    labels = np.array(labels_list)            # [R, label_width, N, 48]
    parameters = np.array(parameters_list)    # [R, NUM_CALIB_PARAMS]

    if transform:
        inputs = np.log1p(inputs)
        labels = np.log1p(labels)

    print(f"Loaded {loaded} runs ({skipped} skipped)")
    print(f"  Inputs shape: {inputs.shape}")
    print(f"  Labels shape: {labels.shape}")
    print(f"  Parameters shape: {parameters.shape}")
    print(f"  Adjacency shape: {adjacency.shape}")

    return {
        "inputs": inputs,
        "labels": labels,
        "parameters": parameters,
        "adjacency": adjacency,
        "dampings": dampings_list,
    }


def _resample_to_daily(times: np.ndarray, min_days: int) -> Optional[np.ndarray]:
    """Selects one index per integer day from the time array.

    The ABM may output at sub-daily intervals (e.g., every 12 hours).
    This selects the index closest to each integer day value.

    :param times: Array of time points (in days, float).
    :param min_days: Minimum number of daily samples required.
    :returns: Array of indices, one per day, or None if insufficient data.
    """
    if len(times) == 0:
        return None

    max_day = int(np.floor(times[-1]))

    if max_day < min_days:
        return None

    daily_indices = []
    for day in range(min_days):
        # Find index closest to this day
        idx = np.argmin(np.abs(times - day))
        daily_indices.append(idx)

    return np.array(daily_indices)


def get_binary_adjacency(adjacency: np.ndarray) -> np.ndarray:
    """Converts weighted adjacency matrix to binary (0/1).

    :param adjacency: Weighted adjacency matrix.
    :returns: Binary adjacency matrix.
    """
    return (adjacency > 0).astype(np.float32)


def save_as_pickle(data: Dict, output_path: str, filename: str):
    """Saves the loaded dataset as a pickle file compatible with the
    existing GNN training pipeline.

    :param data: Dataset dictionary from load_dataset().
    :param output_path: Output directory.
    :param filename: Output filename.
    """
    os.makedirs(output_path, exist_ok=True)

    # Reshape to match existing pipeline format:
    # inputs: [num_runs, num_nodes, features, input_width]
    # labels: [num_runs, num_nodes, features, label_width]
    inputs_reshaped = data["inputs"].transpose(0, 2, 3, 1)
    labels_reshaped = data["labels"].transpose(0, 2, 3, 1)

    save_data = {
        "inputs": inputs_reshaped,
        "labels": labels_reshaped,
        "parameters": data.get("parameters"),
        "param_names": CALIB_PARAM_NAMES,
        "adjacency": get_binary_adjacency(data["adjacency"]),
        "dampings": data["dampings"],
    }

    filepath = os.path.join(output_path, filename)
    with open(filepath, "wb") as f:
        pickle.dump(save_data, f)

    print(f"Dataset saved to {filepath}")


def main():
    """Example usage: load ABM data and prepare for GNN training."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Load ABM simulation data for GNN surrogate training.")
    parser.add_argument("--data_dir", type=str, required=True,
                        help="Directory containing run_XXXX subdirectories")
    parser.add_argument("--input_width", type=int, default=5,
                        help="Number of input time steps (default: 5)")
    parser.add_argument("--label_width", type=int, default=25,
                        help="Number of prediction time steps (default: 25)")
    parser.add_argument("--max_runs", type=int, default=None,
                        help="Max runs to load (default: all)")
    parser.add_argument("--output_dir", type=str, default=None,
                        help="Save as pickle to this directory")
    parser.add_argument("--no_transform", action="store_true",
                        help="Disable log(1+x) transformation")
    args = parser.parse_args()

    data = load_dataset(
        data_dir=args.data_dir,
        input_width=args.input_width,
        label_width=args.label_width,
        max_runs=args.max_runs,
        transform=not args.no_transform,
    )

    if args.output_dir:
        label_w = args.label_width
        n_runs = data["inputs"].shape[0]
        filename = f"ABM_GNN_data_{label_w}days_{n_runs}runs.pickle"
        save_as_pickle(data, args.output_dir, filename)


if __name__ == "__main__":
    main()
