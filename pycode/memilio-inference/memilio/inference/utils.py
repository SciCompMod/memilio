import os
from typing import Any, Callable
import datetime
import numpy as np
import pandas as pd
import pickle
import logging

from memilio.inference.prior import PriorScaler

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
    # Logger init
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    offline_data_file_path = os.path.join(
        output_folder_path, "offline_data.pkl")
    offline_data = None
    if os.path.isfile(offline_data_file_path):
        logger.info(f"Training data loaded from {offline_data_file_path}")
        with open(offline_data_file_path, 'rb') as file:
            offline_data = pickle.load(file)
    else:
        logger.info(
            f"Generate Training data and save to {offline_data_file_path}")
        offline_data = generative_model(batch_size)
        with open(offline_data_file_path, 'wb') as file:
            pickle.dump(offline_data, file)
    return offline_data


def load_data_rki_sir(date_data_begin: datetime.date, T: int, data_path: str) -> np.ndarray:
    """Helper function to load cumulative cases and transform them to new cases."""

    # Use right corona data based on the model (either reporting or reference date)
    confirmed_cases_json = data_path
    confirmed_cases = pd.read_json(confirmed_cases_json)
    confirmed_cases = confirmed_cases.set_index('Date')

    date_data_end = date_data_begin + datetime.timedelta(T)
    cases_obs = np.array(
        confirmed_cases.loc[date_data_begin:date_data_end]
    ).flatten()
    new_cases_obs = np.diff(cases_obs)
    return new_cases_obs


def load_data_synthetic(simulator_fun: Callable[[list[float]], np.ndarray], params_synthetic_data: list[float]) -> np.ndarray:
    """Helper function to generate new cases from ."""
    new_cases_obs = simulator_fun(
        params=params_synthetic_data).flatten()
    return new_cases_obs


def start_training(trainer, epochs, batch_size, offline_data=None, iterations_per_epoch=None, **kwargs):

    # either use offline_data for offline training or iterations_per_epoch for online

    epochs -= trainer.loss_history.latest

    # check if epochs were already reached with checkpoints
    if epochs <= 0:
        return trainer.loss_history.get_plottable()

    # Train
    if offline_data is not None:
        history = trainer.train_offline(
            offline_data, epochs=epochs, batch_size=batch_size, **kwargs)
    elif iterations_per_epoch is not None:
        history = trainer.train_online(
            epochs=epochs, iterations_per_epoch=500, **kwargs)
    else:
        raise Exception(
            "No offline_data for offline training or iterations_per_epoch for online training were defined.")

    return history
