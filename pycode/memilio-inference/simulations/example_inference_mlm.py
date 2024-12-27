import datetime
from functools import partial
import os
import numpy as np
from typing import Any, Callable
import pandas as pd
import tensorflow as tf

from memilio.inference.plotting import MLMPlotting
from memilio.inference.prior import ModelPriorBuilder, TwoLevelPriorScaler
from memilio.inference.sir_mlm import HyperparameterNamesSir, ParameterNamesSir, SIRStrategy, simulator_SIR, get_population, create_prior
from memilio.inference.config import InferenceConfig, TrainerParameters
from memilio.inference.utils import generate_offline_data, generate_offline_data_splitwise, start_training, RedirectStdout, parallelize_simulation
from memilio.inference.networks import TwoLevelSequenceNetwork

from bayesflow.amortizers import TwoLevelAmortizedPosterior, AmortizedPosterior
from bayesflow.networks import InvertibleNetwork, SequenceNetwork, DeepSet, HierarchicalNetwork
from bayesflow.simulation import TwoLevelGenerativeModel, Simulator
from bayesflow.trainers import Trainer

RNG = np.random.default_rng(2023)
FILE_PATH = os.path.dirname(os.path.abspath(__file__))

# Check intra-op threads
intra_threads = tf.config.threading.get_intra_op_parallelism_threads()
print("Default intra-op threads:", intra_threads)

# Check inter-op threads
inter_threads = tf.config.threading.get_inter_op_parallelism_threads()
print("Default inter-op threads:", inter_threads)

# In TensorFlow 2.x, use the same configuration when starting a session:
tf.config.threading.set_intra_op_parallelism_threads(16)
tf.config.threading.set_inter_op_parallelism_threads(16)


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


def load_data_synthetic(simulator_fun: Callable[..., np.ndarray], **kwargs: Any) -> np.ndarray:
    """Helper function to generate new cases from ."""
    new_cases_obs = simulator_fun(**kwargs)
    return new_cases_obs


def configure_input(forward_dict: dict[str, Any], prior_scaler: TwoLevelPriorScaler, direct_conditions: dict[str, Any] = {}) -> dict[str, Any]:
    """
    Function to configure the simulated quantities (i.e., simulator outputs)
    into a neural network-friendly (BayesFlow) format.
    """

    # Prepare placeholder dict
    out_dict = {}

    # Remove a batch if it contains negative values
    data = forward_dict["sim_data"]
    idx_keep = np.all((data >= 0), axis=(1, 2, 3))
    if not np.all(idx_keep):
        print("Invalid value encountered...removing from batch")

    # Convert data to logscale
    logdata = np.log1p(data[idx_keep]).astype(np.float32)

    # Extract prior draws and z-standardize with previously computed means
    params = {}
    params["hyper_parameters"] = forward_dict["hyper_prior_draws"][idx_keep]
    params["local_parameters"] = forward_dict["local_prior_draws"][idx_keep]
    params["shared_parameters"] = forward_dict["shared_prior_draws"][idx_keep]
    params = prior_scaler.transform(params)

    # Remove a batch if it contains nan, inf or -inf
    idx_keep = np.all(np.isfinite(logdata), axis=(1, 2, 3))
    if not np.all(idx_keep):
        print("Invalid value encountered...removing from batch")

    # Add to keys
    out_dict["summary_conditions"] = logdata[idx_keep]
    out_dict["hyper_parameters"] = params["hyper_parameters"][idx_keep].astype(
        np.float32)
    out_dict["local_parameters"] = params["local_parameters"][idx_keep].transpose((0, 2, 1)).astype(
        np.float32)
    out_dict["shared_parameters"] = params["shared_parameters"][idx_keep].astype(
        np.float32)

    out_dict["direct_local_conditions"] = np.repeat([direct_conditions.get(
        "direct_local_conditions", {})], out_dict["summary_conditions"].shape[0], 0)
    out_dict["direct_global_conditions"] = np.repeat([direct_conditions.get(
        "direct_global_conditions", {})], out_dict["summary_conditions"].shape[0], 0)

    return out_dict


def build_amortizer(trainer_parameters, prior):
    # # Should not be lower than number of parameters
    # summary_net = SequenceNetwork(summary_dim=trainer_parameters.summary_dim)

    # could try setting coupling_design=spline for more superior performance on lower-dimensional problems.
    # infere local parameters on local data
    local_inference_net = InvertibleNetwork(
        num_params=2, num_coupling_layers=trainer_parameters.num_coupling_layers, coupling_design=trainer_parameters.coupling_design)
    local_amortizer = AmortizedPosterior(
        local_inference_net, name="local_amortizer")

    # could try setting coupling_design=spline for more superior performance on lower-dimensional problems.
    # infere hyper and shared parameters on global data
    global_inference_net = InvertibleNetwork(
        num_params=7, num_coupling_layers=trainer_parameters.num_coupling_layers, coupling_design=trainer_parameters.coupling_design)
    global_amortizer = AmortizedPosterior(
        global_inference_net, name="global_amortizer")

    summary_net = HierarchicalNetwork(
        [TwoLevelSequenceNetwork(summary_dim=10), DeepSet(summary_dim=20)])

    amortizer = TwoLevelAmortizedPosterior(
        local_amortizer, global_amortizer, summary_net, name="covid_amortizer")

    return amortizer


def run_inference(output_folder_path: os.PathLike, config: InferenceConfig) -> None:

    # create dir if it does not exist
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    # save current config
    config.save(output_folder_path)

    # Define synthetic case
    # n_regions = 4
    # N = [20000000, 22000000, 25000000, 20000000]
    # mobility_params = [[0] * n_regions] * n_regions
    # mobility_condition = [[0, 0]] * n_regions

    # Define synthetic with commuting
    n_regions = 4
    N = [20000000, 22000000, 25000000, 20000000]
    mobility_params = [[0, 20000, 15000, 18000],
                       [20000, 0, 10000, 14000],
                       [25000, 22000, 0, 28000],
                       [5000, 8000, 1000, 0]]
    mobility_condition = np.transpose(
        [np.sum(mobility_params, axis=0), np.sum(mobility_params, axis=1)])

    # ([[lambd0], [I0]], [mu, f_i, phi_i, D_i, scale_I])
    combined_params = ([[0.7, 0.75, 0.6, 0.8], [500, 2000, 500, 100]], [
                       6, 0.7, 2.5, 7, 2])
    direct_conditions = {"direct_local_conditions": np.concatenate([np.log1p(np.array(N)[:, np.newaxis]), np.log1p(mobility_condition)], axis=1).astype(np.float32),
                         "direct_global_conditions": np.log1p(np.sum(N, keepdims=True)).astype(np.float32)}

    # Define germany states
    # n_regions = 16
    # path_population_data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
    #                                     "..", "data", "pydata", "Germany",
    #                                     "county_current_population.json")
    # N = get_population(path_population_data)

    # prior
    prior = create_prior()
    prior_args = {"local_args": {"n_regions": n_regions}}
    # prior(10, **prior_args)
    prior_scaler = TwoLevelPriorScaler().fit(
        prior, **prior_args)

    # Create GenerativeModel
    simulator_function = partial(
        simulator_SIR, n_regions=n_regions, N=N, T=config[
            "T"], mobility_params=mobility_params, intervention_model=config["intervention_model"],
        observation_model=config["observation_model"], param_names=[])

    batch_simulator_fun = parallelize_simulation(
        num_threads=16)(simulator=simulator_function)
    simulator = Simulator(batch_simulator_fun=batch_simulator_fun)
    generative_model = TwoLevelGenerativeModel(
        prior, simulator, skip_test=True, name="sir_covid_simulator")

    # generative_model(1, prior_args=prior_args)
    # offline_data = generate_offline_data_splitwise(
    #     output_folder_path, generative_model, 1000000, 100000, generative_model_args={"prior_args": prior_args})
    offline_data = generate_offline_data(
        output_folder_path, generative_model, 1000000, generative_model_args={"prior_args": prior_args})

    # offline_data = generate_offline_data_splitwise(
    #     output_folder_path, generative_model, 1000, 100, generative_model_args={"prior_args": prior_args})

    # Create Trainer
    amortizer = build_amortizer(config.trainer_parameters, prior)
    trainer = Trainer(amortizer=amortizer, generative_model=generative_model,
                      configurator=partial(
                          configure_input, prior_scaler=prior_scaler, direct_conditions=direct_conditions), memory=True,
                      checkpoint_path=os.path.join(output_folder_path, "checkpoint"), skip_checks=True)

    # Train **kwargs.pop("val_model_args", {})
    history = start_training(trainer, epochs=config.trainer_parameters.epochs,
                             batch_size=25000, offline_data=offline_data, early_stopping=True, validation_sims=2000, val_model_args={"prior_args": prior_args}, save_logs=False, output_folder_path=output_folder_path)
    # history = start_training(trainer, epochs=config.trainer_parameters.epochs,
    #                          batch_size=3200, iterations_per_epoch=500, early_stopping=True, validation_sims=5000)

    # obs_data: (n_groups, n_time_steps, n_time_series)
    # config["obs_data"] = load_data_rki_sir(datetime.date(2020, 3, 1), config.T, os.path.join(os.path.dirname(os.path.abspath(
    #     __file__)), "../../../data/pydata/Germany/cases_infected_repdate.json"))
    config["obs_data"] = load_data_synthetic(
        simulator_function, combined_params=combined_params)

    config["obs_input"] = {
        "summary_conditions": np.log1p(config["obs_data"])[np.newaxis].astype(np.float32),
        "direct_local_conditions": [direct_conditions["direct_local_conditions"]],
        "direct_global_conditions": [direct_conditions["direct_global_conditions"]]}

    # history = trainer.loss_history
    MLMPlotting(output_folder_path).plot_all(history, config, prior, prior_scaler,
                                             simulator_function, generative_model, trainer, prior_args=prior_args)


if __name__ == "__main__":
    trainer_parameters = TrainerParameters(epochs=20, coupling_design="spline",
                                           summary_dim=16, num_coupling_layers=10)
    config = InferenceConfig(T=81, N=83e6, intervention_model=False,
                             observation_model=True, trainer_parameters=trainer_parameters)

    run_inference(output_folder_path=os.path.join(
        FILE_PATH, "output_mlm/mlm_resimulation_mobility_spline"), config=config)
