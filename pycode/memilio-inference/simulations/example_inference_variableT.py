import datetime
from functools import partial
import os
from typing import Callable, Any
import numpy as np
import pandas as pd

import memilio.simulation as mio
from memilio.inference.plotting import Plotting
from memilio.inference.prior import (ModelPriorBuilder, PriorScaler, UnboundParameter, LambdaParameter, DelayParameter,
                                     InterventionChangePointParameter, WeeklyModulationAmplitudeParameter,
                                     ScaleMultiplicativeReportingNoiseParameter, ModelStrategy)
from memilio.inference.sir import SIRStrategy, simulator_SIR, DEFAULT_PRIORS
from memilio.inference.config import InferenceConfig, TrainerParameters
from memilio.inference.utils import generate_offline_data

from bayesflow.amortizers import AmortizedPosterior
from bayesflow.networks import InvertibleNetwork, SequenceNetwork
from bayesflow.simulation import GenerativeModel, Simulator, ContextGenerator
from bayesflow.trainers import Trainer

RNG = np.random.default_rng(2023)
FILE_PATH = os.path.dirname(os.path.abspath(__file__))
mio.set_log_level(mio.LogLevel.Warning)

VARIABLE_T_PRIORS = {
    "LAMBDA_1": {
        "type": "prior",
        "parameter_class": LambdaParameter,
        "distribution": np.random.lognormal,
        "params": {"mean": np.log(0.6), "sigma": 0.5},
        "name": r"$\lambda_1$",
        "description": "infectionrate after first intervention"
    },
    "LAMBDA_2": {
        "type": "prior",
        "parameter_class": LambdaParameter,
        "distribution": np.random.lognormal,
        "params": {"mean": np.log(0.6), "sigma": 0.5},
        "name": r"$\lambda_2$",
        "description": "infectionrate after second intervention"
    },
    "LAMBDA_3": {
        "type": "prior",
        "parameter_class": LambdaParameter,
        "distribution": np.random.lognormal,
        "params": {"mean": np.log(0.6), "sigma": 0.5},
        "name": r"$\lambda_3$",
        "description": "infectionrate after third intervention"
    },
    "LAMBDA_4": {
        "type": "prior",
        "parameter_class": LambdaParameter,
        "distribution": np.random.lognormal,
        "params": {"mean": np.log(0.6), "sigma": 0.5},
        "name": r"$\lambda_4$",
        "description": "infectionrate after fourth intervention"
    },
    "T1": {
        "type": "prior",
        "parameter_class": InterventionChangePointParameter,
        "distribution": np.random.uniform,
        "params": {"low": 0, "high": 180},
        "name": r"$t_1$",
        "description": "first InterventionChangePoint"
    },
    "T2": {
        "type": "prior",
        "parameter_class": InterventionChangePointParameter,
        "distribution": np.random.uniform,
        "params": {"low": 0, "high": 180},
        "name": r"$t_2$",
        "description": "second InterventionChangePoint"
    },
    "T3": {
        "type": "prior",
        "parameter_class": InterventionChangePointParameter,
        "distribution": np.random.uniform,
        "params": {"low": 0, "high": 180},
        "name": r"$t_3$",
        "description": "third InterventionChangePoint"
    },
    "T4": {
        "type": "prior",
        "parameter_class": InterventionChangePointParameter,
        "distribution": np.random.uniform,
        "params": {"low": 0, "high": 180},
        "name": r"$t_4$",
        "description": "fourth InterventionChangePoint"
    }
}


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

    # Make inference network aware of varying numbers of trials
    # We create a vector of shape (batch_size, 1) by repeating the sqrt(num_obs)
    vec_num_obs = forward_dict["sim_non_batchable_context"] * \
        np.ones((data.shape[0], 1))
    out_dict["direct_conditions"] = np.sqrt(vec_num_obs).astype(np.float32)

    # Add to keys
    out_dict["summary_conditions"] = logdata[idx_keep]
    out_dict["parameters"] = params[idx_keep]

    return out_dict


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


def build_amortizer(trainer_parameters, prior):
    # Create Trainer
    # Should not be lower than number of parameters
    summary_net = SequenceNetwork(summary_dim=trainer_parameters.summary_dim)
    # could try setting coupling_design=spline for more superior performance on lower-dimensional problems.
    inference_net = InvertibleNetwork(num_params=len(
        prior.param_names), num_coupling_layers=trainer_parameters.num_coupling_layers, coupling_design=trainer_parameters.coupling_design)
    amortizer = AmortizedPosterior(
        inference_net, summary_net, name="covid_amortizer")

    return amortizer


def run_inference(output_folder_path: os.PathLike, config: InferenceConfig) -> None:

    # create dir if it does not exist
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    # save current config
    # config.save(output_folder_path)

    # Change default values of prior
    for key, value in VARIABLE_T_PRIORS.items():
        DEFAULT_PRIORS[key] = value

    # Create Priors
    prior_builder = ModelPriorBuilder(SIRStrategy)
    prior_builder.add_base()
    if config["intervention_model"]:
        prior_builder.add_intervention()
    if config["observation_model"]:
        prior_builder.add_observation()
    prior = prior_builder.build()
    prior_scaler = PriorScaler().fit(prior)
    # simulator_SIR(prior(1)["prior_draws"][0], 80000, 81, True, True,
    #                param_names=prior.param_names, fixed_params=fixed_parameters) # Debug

    def random_num_days(min_days=70, max_days=140):
        """Draws a random number of observations for all simulations in a batch."""
        return RNG.integers(low=min_days, high=max_days + 1)

    context_gen = ContextGenerator(
        non_batchable_context_fun=random_num_days
    )

    # Create GenerativeModel
    simulator_function = partial(
        simulator_SIR, N=config["N"], intervention_model=config["intervention_model"],
        observation_model=config["observation_model"], param_names=prior.param_names)
    simulator = Simulator(simulator_fun=simulator_function,
                          context_generator=context_gen)
    generative_model = GenerativeModel(
        prior, simulator, name="sir_covid_simulator")

    # offline_data = generate_offline_data(
    #     output_folder_path, generative_model, 100000)

    # Create Trainer
    amortizer = build_amortizer(config.trainer_parameters, prior)
    trainer = Trainer(amortizer=amortizer, generative_model=generative_model,
                      configurator=partial(
                          configure_input, prior_scaler=prior_scaler), memory=True,
                      checkpoint_path=os.path.join(output_folder_path, "checkpoint"))

    # # Train
    # history = trainer.train_offline(
    #     offline_data, epochs=config.trainer_parameters.epochs, batch_size=320, validation_sims=2000)
    # history = trainer.train_online(
    #     epochs=config.trainer_parameters.epochs, iterations_per_epoch=500, batch_size=32, validation_sims=2000)

    config["obs_data"] = load_data_rki(81, os.path.join(os.path.dirname(os.path.abspath(
        __file__)), "../../data/pydata/Germany/cases_infected_repdate.json"))

    history = trainer.loss_history
    Plotting(output_folder_path).plot_all(history, config, prior, prior_scaler,
                                          simulator_function, generative_model, trainer, amortizer)


if __name__ == "__main__":
    trainer_parameters = TrainerParameters(epochs=20,
                                           summary_dim=16, num_coupling_layers=8)
    config = InferenceConfig(T=None, N=83e6, intervention_model=True,
                             observation_model=True, trainer_parameters=trainer_parameters)

    run_inference(output_folder_path=os.path.join(
        FILE_PATH, "output/variableT"), config=config)
