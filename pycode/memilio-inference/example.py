import datetime
from functools import partial
import os
from typing import Callable, Any, Self
from inspect import isclass
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt

from plots import *
from prior import ModelPriorBuilder, PriorScaler
from sir import *
from config import InferenceConfig

import bayesflow.diagnostics as diag
from bayesflow.amortizers import AmortizedPosterior
from bayesflow.networks import InvertibleNetwork, SequenceNetwork
from bayesflow.simulation import GenerativeModel, Simulator, Prior
from bayesflow.trainers import Trainer

RNG = np.random.default_rng(2023)
FILE_PATH = os.path.dirname(os.path.abspath(__file__))


def load_data_rki(T: int, data_path: str) -> np.ndarray:
    """Helper function to load cumulative cases and transform them to new cases."""

    # Use right corona data based on the model (either reporting or reference date)
    confirmed_cases_json = data_path
    confirmed_cases = pd.read_json(confirmed_cases_json)
    confirmed_cases = confirmed_cases.set_index('Date')

    date_data_begin = datetime.date(2020, 3, 1)
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


def run_inference(output_folder_path: os.PathLike) -> None:

    fixed_parameters = {ParameterNamesSir.I0.value: 200}

    config = InferenceConfig(T=81, N=83e6, intervention_model=False,
                             observation_model=False)

    # create dir if it does not exist
    if not os.path.exists(output_folder_path):
        os.path.makedirs(output_folder_path)

    # Create Priors
    prior_builder = ModelPriorBuilder(SIRStrategy)
    prior_builder.add_base()
    if config["intervention_model"]:
        prior_builder.add_intervention()
    if config["observation_model"]:
        prior_builder.add_observation()
    prior = prior_builder.set_fixed_parameters(fixed_parameters).build()
    prior_scaler = PriorScaler().fit(prior)
    # stationary_SIR(prior(1)["prior_draws"][0], 80000, 81, True, True,
    #                param_names=prior.param_names, fixed_params=fixed_parameters) # Debug

    # Create GenerativeModel
    simulator_function = partial(
        stationary_SIR, N=config["N"], T=config["T"], intervention_model=config["intervention_model"],
        observation_model=config["observation_model"], param_names=prior.param_names, fixed_params=fixed_parameters)
    simulator = Simulator(simulator_fun=simulator_function)
    generative_model = GenerativeModel(
        prior, simulator, name="sir_covid_simulator")

    offline_data = generate_offline_data(
        output_folder_path, generative_model, 100)

    # Create Trainer
    # Should not be lower than number of parameters
    summary_net = SequenceNetwork(summary_dim=4)
    # could try setting coupling_design=spline for more superior performance on lower-dimensional problems.
    inference_net = InvertibleNetwork(num_params=len(
        prior.param_names), num_coupling_layers=6)
    amortizer = AmortizedPosterior(
        inference_net, summary_net, name="covid_amortizer")
    trainer = Trainer(amortizer=amortizer, generative_model=generative_model,
                      configurator=partial(
                          configure_input, prior_scaler=prior_scaler), memory=True,
                      checkpoint_path=os.path.join(output_folder_path, "checkpoint"))

    # Train
    history = trainer.train_offline(
        offline_data, epochs=5, batch_size=320, validation_sims=200)

    f = diag.plot_losses(history["train_losses"],
                         history["val_losses"], moving_average=True)
    plt.savefig(os.path.join(output_folder_path, "losses.png"))

    # config["obs_data"] = load_data_rki(T, os.path.join(os.path.dirname(os.path.abspath(
    #     __file__)), "../../data/pydata/Germany/cases_infected_repdate.json"))
    params_synthetic_data = [0.8, 3]
    config["obs_data"] = load_data_synthetic(
        simulator_function, params_synthetic_data)
    # Format data into a 3D array of shape (1, n_time_steps, 1) and perform log transform
    obs_data = np.log1p(config["obs_data"])[
        np.newaxis, :, np.newaxis].astype(np.float32)
    # Obtain 500 posterior draws given real data
    post_samples = amortizer.sample({"summary_conditions": obs_data}, 5000)
    # Undo standardization to get parameters on their original (unstandardized) scales
    post_samples = prior_scaler.inverse_transform(post_samples)

    f = diag.plot_posterior_2d(post_samples, prior=prior)
    plt.savefig(os.path.join(output_folder_path, "posterior_2d.png"))

    f = plot_ppc(config, post_samples, simulator_function)
    plt.savefig(os.path.join(output_folder_path, "ppc.png"))


if __name__ == "__main__":
    run_inference(output_folder_path=os.path.join(FILE_PATH, "output"))
