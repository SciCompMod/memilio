import os
from typing import Any
import numpy as np
import pickle

from prior import PriorScaler

from bayesflow.simulation import GenerativeModel


def configure_input(forward_dict: dict[str, Any], prior_scaler: PriorScaler) -> dict[str, Any]:
    """
    Function to configure the simulated quantities (i.e., simulator outputs)
    into a neural network-friendly (BayesFlow) format.
    """

    # Prepare placeholder dict
    out_dict = {}

    # Remove a batch if it contains negative values
    data = forward_dict["sim_data"]
    idx_keep = np.all((data >= 0), axis=(1, 2))
    if not np.all(idx_keep):
        print("Invalid value encountered...removing from batch")

    # Convert data to logscale
    logdata = np.log1p(data[idx_keep]).astype(np.float32)

    # Extract prior draws and z-standardize with previously computed means
    params = forward_dict["prior_draws"][idx_keep].astype(np.float32)
    params = prior_scaler.transform(params)

    # Remove a batch if it contains nan, inf or -inf
    idx_keep = np.all(np.isfinite(logdata), axis=(1, 2))
    if not np.all(idx_keep):
        print("Invalid value encountered...removing from batch")

    # Add to keys
    out_dict["summary_conditions"] = logdata[idx_keep]
    out_dict["parameters"] = params[idx_keep]

    return out_dict


def generate_offline_data(output_folder_path: os.PathLike, generative_model: GenerativeModel, batch_size: int) -> dict[str, Any]:
    offline_data_file_path = os.path.join(
        output_folder_path, "offline_data.pkl")
    offline_data = None
    if os.path.isfile(offline_data_file_path):
        with open(offline_data_file_path, 'rb') as file:
            offline_data = pickle.load(file)
    else:
        offline_data = generative_model(batch_size)
        with open(offline_data_file_path, 'wb') as file:
            pickle.dump(offline_data, file)
    return offline_data
