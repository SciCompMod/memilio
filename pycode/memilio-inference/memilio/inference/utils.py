import os
from typing import Any, Callable
import datetime
import numpy as np
import pandas as pd
import pickle
import logging
import tqdm
import sys
from contextlib import nullcontext
from multiprocessing import Pool
from functools import partial

from memilio.inference.prior import PriorScaler

from bayesflow.simulation import GenerativeModel, TwoLevelGenerativeModel


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


class RedirectStdout:
    def __init__(self, log_file):
        self.log_file_handle = open(log_file, "a")
        self.original_stdout = sys.stdout
        self.original_stderr = sys.stderr

    def __enter__(self):
        class Logger:
            def __init__(self, file_handle, terminal_stdout, terminal_stderr):
                self.terminal_stdout = terminal_stdout
                self.terminal_stderr = terminal_stderr
                self.file = file_handle

            def write(self, message):
                # Write to both file and terminal
                if self.terminal_stdout:  # For stdout messages
                    self.terminal_stdout.write(message)
                if self.terminal_stderr:  # For stderr messages (like tqdm)
                    self.terminal_stderr.write(message)
                self.file.write(message)

            def flush(self):
                # Ensure both file and terminal are flushed
                if self.terminal_stdout:
                    self.terminal_stdout.flush()
                if self.terminal_stderr:
                    self.terminal_stderr.flush()
                self.file.flush()

        logger = Logger(self.log_file_handle,
                        self.original_stdout, self.original_stderr)
        sys.stdout = logger
        sys.stderr = logger  # Redirect tqdm output (stderr) to the same logger

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stdout = self.original_stdout
        sys.stderr = self.original_stderr
        self.log_file_handle.close()


def generate_offline_data(output_folder_path: os.PathLike, generative_model: GenerativeModel | TwoLevelGenerativeModel, batch_size: int, **kwargs) -> dict[str, Any]:
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
        offline_data = generative_model(
            batch_size, **kwargs.pop("generative_model_args", {}))

        logger.info(
            f"Generate Training data finished")
        with open(offline_data_file_path, 'wb') as file:
            pickle.dump(offline_data, file)
        logger.info(
            f"Save Training data finished")
    return offline_data


def generate_offline_data_splitwise(output_folder_path: os.PathLike, generative_model: GenerativeModel | TwoLevelGenerativeModel, batch_size: int, step_size: int, **kwargs) -> dict[str, Any]:
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
        n_steps = batch_size // step_size
        last_step = batch_size % step_size
        offline_data = {}

        pbar = tqdm.tqdm("Generate Trainings data",
                         total=n_steps + 1 if last_step else n_steps)
        for _ in range(n_steps):
            step_data = generative_model(
                step_size, **kwargs.get("generative_model_args", {}))
            for key, value in step_data.items():
                if key not in offline_data.keys():
                    offline_data[key] = value
                elif value is None:
                    continue
                else:
                    offline_data[key] = np.append(
                        offline_data[key], value, axis=0)

            pbar.update()
        if last_step:
            step_data = generative_model(
                last_step, **kwargs.get("generative_model_args", {}))
            for key, value in step_data.items():
                if key not in offline_data.keys():
                    offline_data[key] = value
                elif value is None:
                    continue
                else:
                    offline_data[key] = np.append(
                        offline_data[key], value, axis=0)
            pbar.update()
        pbar.close()

        logger.info(
            f"Generate Training data finished")
        with open(offline_data_file_path, 'wb') as file:
            pickle.dump(offline_data, file)
        logger.info(
            f"Save Training data finished")
    return offline_data


def parallelize_simulation(num_threads=4, use_batchable_context=False, use_non_batchable_context=False):
    """
    Decorator to optionally parallelize the simulator function based on flags.

    From 7:20 -> 2:50 for 4 threads. 
    """
    def decorator(simulator):
        def wrapper(params, *args, **kwargs):
            # Extract batch size
            # Always assume first dimension is batch dimension
            # Handle cases with multiple inputs to simulator
            if isinstance(params, tuple) or isinstance(params, list):
                batch_size = params[0].shape[0]
                non_batched_params = [[params[i][b]
                                       for i in range(len(params))] for b in range(batch_size)]
            # Handle all other cases or fail gently
            else:
                # expand dimension by one to handle both cases in the same way
                batch_size = params.shape[0]
                non_batched_params = params

            # No context type
            if (
                not use_batchable_context
                and not use_non_batchable_context
            ):
                batch_args = [
                    (non_batched_params[b], *args)
                    for b in range(batch_size)
                ]

            # Only batchable context
            elif not use_non_batchable_context:
                batch_args = [
                    (non_batched_params[b],
                     args[0][b], *args[1:])
                    for b in range(batch_size)
                ]
            # Only non-batchable context
            elif not use_batchable_context:
                batch_args = [
                    (non_batched_params[b],
                     args[0], *args[1:])
                    for b in range(batch_size)
                ]

            # Both batchable and non_batchable context
            else:
                batch_args = [
                    (non_batched_params[b], args[0][b],
                     args[1], *args[2:])
                    for b in range(batch_size)
                ]

            # Parallel execution using Pool
            with Pool(num_threads) as pool:
                results = pool.map(
                    partial(process_batch, simulator=simulator, **kwargs), batch_args)

            return np.array(results)

        return wrapper
    return decorator


def process_batch(batch_args, simulator, **kwargs):
    # Unpack batch_args and call the simulator
    return simulator(*batch_args, **kwargs)


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


def start_training(trainer, epochs, batch_size, offline_data=None, iterations_per_epoch=None, save_logs=True, output_folder_path="", **kwargs):

    # either use offline_data for offline training or iterations_per_epoch for online
    epochs -= trainer.loss_history.latest

    # check if epochs were already reached with checkpoints
    if epochs <= 0:
        return trainer.loss_history.get_plottable()

    # Context to redirect output into a file
    with RedirectStdout(os.path.join(
            output_folder_path, "training_output.log")) if save_logs else nullcontext():

        # Train
        if offline_data is not None:
            # history = trainer.train_offline(
            #     offline_data, epochs=epochs, batch_size=batch_size, **kwargs)
            history = trainer.train_offline(
                offline_data, epochs=epochs, batch_size=batch_size, use_autograph=True, **kwargs)  # For Debug
        elif iterations_per_epoch is not None:
            history = trainer.train_online(
                epochs=epochs, iterations_per_epoch=500, **kwargs)
        else:
            raise Exception(
                "No offline_data for offline training or iterations_per_epoch for online training were defined.")

    return history
