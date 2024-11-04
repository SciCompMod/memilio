import datetime
from functools import partial
import os
from typing import Callable, Any
import numpy as np

from memilio.inference.plotting import Plotting
from memilio.inference.prior import ModelPriorBuilder, PriorScaler
from memilio.inference.sir_dampingIntervals import ParameterNamesSir, SIRStrategy, simulator_SIR
from memilio.inference.config import InferenceConfig, TrainerParameters
from memilio.inference.utils import generate_offline_data, configure_input, load_data_rki_sir, start_training

from bayesflow.amortizers import AmortizedPosterior
from bayesflow.networks import InvertibleNetwork, SequenceNetwork
from bayesflow.simulation import GenerativeModel, Simulator, ContextGenerator
from bayesflow.trainers import Trainer

RNG = np.random.default_rng(2023)
FILE_PATH = os.path.dirname(os.path.abspath(__file__))


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
    config.save(output_folder_path)

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

    # Create GenerativeModel
    simulator_function = partial(
        simulator_SIR, T=config["T"], N=config["N"], intervention_model=config["intervention_model"],
        observation_model=config["observation_model"], param_names=prior.param_names)
    simulator = Simulator(simulator_fun=simulator_function)
    generative_model = GenerativeModel(
        prior, simulator, name="sir_covid_simulator")

    offline_data = generate_offline_data(
        output_folder_path, generative_model, 1000000)

    # Create Trainer
    amortizer = build_amortizer(config.trainer_parameters, prior)
    trainer = Trainer(amortizer=amortizer, generative_model=generative_model,
                      configurator=partial(
                          configure_input, prior_scaler=prior_scaler), memory=True,
                      checkpoint_path=os.path.join(output_folder_path, "checkpoint"))

    # Train
    history = start_training(trainer, epochs=config.trainer_parameters.epochs,
                             batch_size=320, offline_data=offline_data, early_stopping=True, validation_sims=2000)
    # history = start_training(trainer, epochs=config.trainer_parameters.epochs,
    #                          batch_size=320, iterations_per_epoch=500, early_stopping=True, validation_sims=2000)

    config["obs_data"] = load_data_rki_sir(datetime.date(2020, 3, 1), config.T, os.path.join(os.path.dirname(os.path.abspath(
        __file__)), "../../data/pydata/Germany/cases_infected_repdate.json"))

    Plotting(output_folder_path).plot_all(history, config, prior, prior_scaler,
                                          simulator_function, generative_model, trainer, amortizer)


if __name__ == "__main__":
    trainer_parameters = TrainerParameters(epochs=2,
                                           summary_dim=16, num_coupling_layers=8)
    config = InferenceConfig(T=81, N=83e6, intervention_model=True,
                             observation_model=True, trainer_parameters=trainer_parameters)

    run_inference(output_folder_path=os.path.join(
        FILE_PATH, "output/dynamicDamping4"), config=config)
