import itertools
import os
from typing import Callable
from functools import partial
import re

from utils import generate_offline_data, configure_input
from sir import ParameterNamesSir, SIRStrategy, stationary_SIR
from plots import *
from prior import ModelPriorBuilder, PriorScaler
from config import InferenceConfig

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


def plot_all(single_output_folder_path, history, config, prior, prior_scaler, simulator_function, generative_model, trainer, amortizer):
    f = prior.plot_prior2d()
    plt.savefig(os.path.join(single_output_folder_path, "prior2D.png"))

    f = trainer.diagnose_latent2d()
    plt.savefig(os.path.join(single_output_folder_path, "latent2d.png"))

    f = diag.plot_losses(history["train_losses"],
                         history["val_losses"], moving_average=True)
    plt.savefig(os.path.join(single_output_folder_path, "losses.png"))

    f = trainer.diagnose_sbc_histograms()
    plt.savefig(os.path.join(single_output_folder_path, "sbc_histograms.png"))

    # Generate some validation data
    validation_sims = trainer.configurator(generative_model(batch_size=300))

    # Generate posterior draws for all simulations
    post_samples = amortizer.sample(validation_sims, n_samples=100)

    # Create ECDF plot
    f = diag.plot_sbc_ecdf(
        post_samples, validation_sims["parameters"], param_names=prior.param_names)
    plt.savefig(os.path.join(single_output_folder_path, "sbc_ecdf.png"))

    f = diag.plot_sbc_ecdf(
        post_samples, validation_sims["parameters"], stacked=True, difference=True, legend_fontsize=12, fig_size=(6, 5)
    )
    plt.savefig(os.path.join(
        single_output_folder_path, "sbc_ecdf_stacked.png"))

    f = diag.plot_sbc_histograms(
        post_samples, validation_sims["parameters"], param_names=prior.param_names)
    plt.savefig(os.path.join(
        single_output_folder_path, "sbc_ecdf_histograms.png"))

    # Format data into a 3D array of shape (1, n_time_steps, 1) and perform log transform
    obs_data = np.log1p(config["obs_data"])[
        np.newaxis, :, np.newaxis].astype(np.float32)
    # Obtain 500 posterior draws given real data
    post_samples = amortizer.sample({"summary_conditions": obs_data}, 5000)
    # Undo standardization to get parameters on their original (unstandardized) scales
    post_samples = prior_scaler.inverse_transform(post_samples)

    f = diag.plot_posterior_2d(post_samples, param_names=prior.param_names)
    plt.savefig(os.path.join(single_output_folder_path, "posterior_2d.png"))

    f = diag.plot_posterior_2d(post_samples, prior=prior)
    plt.savefig(os.path.join(single_output_folder_path,
                "posterior_2d_with_prior.png"))

    f = plot_ppc(config, post_samples, simulator_function)
    plt.savefig(os.path.join(single_output_folder_path, "ppc.png"))


def single_inference_run(single_output_folder_path, config, fixed_parameter_values, learnable_parameters_list):

    fixed_parameters_list = {k: fixed_parameter_values[k]
                             for k in fixed_parameter_values.keys() - list(learnable_parameters_list)}

    # Create Priors
    prior_builder = ModelPriorBuilder(SIRStrategy)
    prior_builder.add_base()
    if config["intervention_model"]:
        prior_builder.add_intervention()
    if config["observation_model"]:
        prior_builder.add_observation()
    prior = prior_builder.set_fixed_parameters(
        fixed_parameters_list).build()
    prior_scaler = PriorScaler().fit(prior)

    # Create GenerativeModel
    simulator_function = partial(
        stationary_SIR, N=config["N"], T=config["T"], intervention_model=config["intervention_model"],
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
        offline_data, epochs=20, batch_size=3200, validation_sims=200)

    config["obs_data"] = load_data_synthetic(
        simulator_function, [fixed_parameter_values[learnable_parameter] for learnable_parameter in learnable_parameters_list])

    plot_all(single_output_folder_path, history, config, prior, prior_scaler,
             simulator_function, generative_model, trainer, amortizer)


def run_pairwise_parameter_inference(output_folder_path, config, fixed_parameter_values):

    # create dir if it does not exist
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    with open(os.path.join(output_folder_path, "config.txt"), 'w') as file:
        print(fixed_parameter_values, file=file)
        print(config, file=file)

    # for fixed_parameter_pair in list(map(dict, itertools.combinations(
    #         fixed_parameter_values.items(), 2))):
    for learnable_parameter_pair in itertools.combinations(
            fixed_parameter_values, 2):

        single_output_folder_path = os.path.join(
            output_folder_path, "output_" + re.sub(r'\\|\$', '', '_'.join(learnable_parameter_pair)))
        # create dir if it does not exist
        if not os.path.exists(single_output_folder_path):
            os.makedirs(single_output_folder_path)

        single_inference_run(single_output_folder_path, config,
                             fixed_parameter_values, learnable_parameter_pair)


if __name__ == "__main__":
    fixed_parameter_values = {ParameterNamesSir.LAMBDA_0.value: 0.7,
                              ParameterNamesSir.MU.value: 4,
                              ParameterNamesSir.I0.value: 200, }

    config = InferenceConfig(T=81, N=83e6, intervention_model=False,
                             observation_model=False)

    run_pairwise_parameter_inference(
        output_folder_path=os.path.join(FILE_PATH, "output_pairwise"),
        config=config, fixed_parameter_values=fixed_parameter_values)
