import itertools
import os
from typing import Callable
from functools import partial
import re

from memilio.inference.utils import generate_offline_data, configure_input
from memilio.inference.sir import SIRStrategy, simulator_SIR, DEFAULT_PRIORS
from memilio.inference.plotting import Plotting
from memilio.inference.prior import ModelPriorBuilder, PriorScaler
from memilio.inference.config import InferenceConfig

import bayesflow.diagnostics as diag
from bayesflow.amortizers import AmortizedPosterior
from bayesflow.networks import InvertibleNetwork, SequenceNetwork
from bayesflow.simulation import GenerativeModel, Simulator
from bayesflow.trainers import Trainer
import numpy as np

RNG = np.random.default_rng(2023)
FILE_PATH = os.path.dirname(os.path.abspath(__file__))


def load_data_synthetic(simulator_fun: Callable[[list[float]], np.ndarray], params_synthetic_data: list[float]) -> np.ndarray:
    """Helper function to generate new cases from ."""
    new_cases_obs = simulator_fun(
        params=params_synthetic_data).flatten()
    return new_cases_obs


def single_inference_run(single_output_folder_path, config, fixed_parameter_values, learnable_parameters_list, add_dummy):

    fixed_parameters_list = {k: fixed_parameter_values[k]
                             for k in fixed_parameter_values.keys() - list(learnable_parameters_list)}

    # Create Priors
    prior_builder = ModelPriorBuilder(SIRStrategy)
    prior_builder.add_base()
    if config["intervention_model"]:
        prior_builder.add_intervention()
    if config["observation_model"]:
        prior_builder.add_observation()
    if add_dummy:
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
        single_output_folder_path, generative_model, 100000)

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
                      checkpoint_path=os.path.join(single_output_folder_path, "checkpoint"))

    # Train
    history = trainer.train_offline(
        offline_data, epochs=10, batch_size=320, validation_sims=200)  # big batches run into performance problems when plotting

    synthetic_case_parameter_values = [
        fixed_parameter_values[learnable_parameter] for learnable_parameter in learnable_parameters_list]
    if add_dummy:
        synthetic_case_parameter_values.append(1.0)

    config["obs_data"] = load_data_synthetic(
        simulator_function, synthetic_case_parameter_values)

    Plotting(single_output_folder_path).plot_all(history, config, prior, prior_scaler,
                                                 simulator_function, generative_model, trainer, amortizer)


def run_pairwise_parameter_inference(output_folder_path: os.PathLike, config: InferenceConfig, fixed_parameter_values: dict[str, float], add_dummy: bool = False):

    # create dir if it does not exist
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    with open(os.path.join(output_folder_path, "config.txt"), 'w') as file:
        print(fixed_parameter_values, file=file)
        print(config, file=file)

    for learnable_parameter_pair in itertools.combinations(
            fixed_parameter_values, 2):

        single_output_folder_path = os.path.join(
            output_folder_path, "output_" + re.sub(r'\\|\$', '', '_'.join(learnable_parameter_pair)))

        if add_dummy:
            single_output_folder_path += "_dummy"

        # create dir if it does not exist
        if not os.path.exists(single_output_folder_path):
            os.makedirs(single_output_folder_path)

        single_inference_run(single_output_folder_path, config,
                             fixed_parameter_values, learnable_parameter_pair, add_dummy)


if __name__ == "__main__":
    fixed_parameter_values = {DEFAULT_PRIORS["LAMBDA_0"]["name"]: 0.7,
                              DEFAULT_PRIORS["MU"]["name"]: 4,
                              DEFAULT_PRIORS["I0"]["name"]: 200, }

    config = InferenceConfig(T=81, N=83e6, intervention_model=False,
                             observation_model=False)

    run_pairwise_parameter_inference(
        output_folder_path=os.path.join(FILE_PATH, "output_pairwise_dummy"),
        config=config, fixed_parameter_values=fixed_parameter_values, add_dummy=True)
