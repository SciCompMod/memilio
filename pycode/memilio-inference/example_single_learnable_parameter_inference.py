import datetime
from functools import partial
import os
from typing import Callable, Any
import numpy as np
import pandas as pd

from plotting import Plotting
from prior import ModelPriorBuilder, PriorScaler
from sir import ParameterNamesSir, SIRStrategy, simulator_SIR
from config import InferenceConfig
from utils import generate_offline_data, configure_input

from bayesflow.amortizers import AmortizedPosterior
from bayesflow.networks import InvertibleNetwork, SequenceNetwork
from bayesflow.simulation import GenerativeModel, Simulator
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


def run_single_learnable_inference(output_folder_path: os.PathLike, config: InferenceConfig) -> None:
    fixed_parameter_values = {ParameterNamesSir.LAMBDA_0.value: 0.7,
                              ParameterNamesSir.MU.value: 4,
                              ParameterNamesSir.I0.value: 200, }

    learnable_parameters_list = [ParameterNamesSir.I0.value]
    fixed_parameters_list = {k: fixed_parameter_values[k]
                             for k in fixed_parameter_values.keys() - list(learnable_parameters_list)}

    # create dir if it does not exist
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    # Create Priors
    prior_builder = ModelPriorBuilder(SIRStrategy)
    prior_builder.add_base()
    if config["intervention_model"]:
        prior_builder.add_intervention()
    if config["observation_model"]:
        prior_builder.add_observation()
    prior_builder.add_dummy()
    prior = prior_builder.set_fixed_parameters(
        fixed_parameters_list).build()
    prior_scaler = PriorScaler().fit(prior)

    # Create GenerativeModel
    simulator_function = partial(
        simulator_SIR, N=config["N"], T=config["T"], intervention_model=config["intervention_model"],
        observation_model=config["observation_model"], param_names=prior.param_names, fixed_params=fixed_parameters_list)
    simulator = Simulator(simulator_fun=simulator_function)
    generative_model = GenerativeModel(
        prior, simulator, name="sir_covid_simulator")

    offline_data = generate_offline_data(
        output_folder_path, generative_model, 100000)

    # Create Trainer
    # Should not be lower than number of parameters
    summary_net = SequenceNetwork(summary_dim=5)
    # could try setting coupling_design=spline for more superior performance on lower-dimensional problems.
    inference_net = InvertibleNetwork(num_params=len(
        prior.param_names), num_coupling_layers=8)
    amortizer = AmortizedPosterior(
        inference_net, summary_net, name="covid_amortizer")
    trainer = Trainer(amortizer=amortizer, generative_model=generative_model,
                      configurator=partial(
                          configure_input, prior_scaler=prior_scaler), memory=True,
                      checkpoint_path=os.path.join(output_folder_path, "checkpoint"))

    # Train
    # history = trainer.train_offline(
    #     offline_data, epochs=10, batch_size=320, validation_sims=200)

    # synthetic_case_parameter_values = [
    #     fixed_parameter_values[learnable_parameter] for learnable_parameter in learnable_parameters_list]
    # synthetic_case_parameter_values.append(1.0)

    # config["obs_data"] = load_data_synthetic(
    #     simulator_function, synthetic_case_parameter_values)

    # Plotting(output_folder_path).plot_all(history, config, prior, prior_scaler,
    #                                       simulator_function, generative_model, trainer, amortizer)


if __name__ == "__main__":
    config = InferenceConfig(T=81, N=83e6, intervention_model=False,
                             observation_model=False)

    run_single_learnable_inference(output_folder_path=os.path.join(
        FILE_PATH, "output/learnable_I_0"), config=config)
