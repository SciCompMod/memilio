import datetime
from functools import partial
import os
import numpy as np
from typing import Any, Callable
import pandas as pd
import tensorflow as tf
import yaml

from memilio.inference.plotting import MLMPlotting
from memilio.inference.prior import ModelPriorBuilder, TwoLevelPriorScaler
from memilio.inference.sir_mlm import DEFAULT_PRIORS, simulator_SIR, get_population, create_prior, build_mlm_amortizer
from memilio.inference.utils import generate_offline_data, generate_offline_data_splitwise, start_training, RedirectStdout, parallelize_simulation

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

# import torch
# from torch_geometric.nn import GCNConv
# import torch.nn as nn
# class MobilityGNN(nn.Module):
#     def __init__(self, node_feature_dim, hidden_dim, output_dim):
#         super().__init__()
#         self.conv1 = GCNConv(node_feature_dim, hidden_dim)
#         self.conv2 = GCNConv(hidden_dim, output_dim)

#     def forward(self, node_features, edge_index, edge_weight=None):
#         x = self.conv1(node_features, edge_index, edge_weight)
#         x = torch.relu(x)
#         x = self.conv2(x, edge_index, edge_weight)
#         return x  # Node embeddings


# def prepare_graph_data(mobility_matrix, population):
#     """
#     Converts the mobility matrix and population data into graph representation.
#     Args:
#         mobility_matrix: [n_regions, n_regions] Mobility matrix.
#         population: [n_regions] Population data for each region.
#     Returns:
#         node_features: [n_regions, node_feature_dim] Node feature matrix.
#         edge_index: [2, num_edges] Edge indices for the graph.
#         edge_weight: [num_edges] Weights corresponding to edges (mobility values).
#     """
#     from torch_geometric.utils import dense_to_sparse

#     # Node features: Add population as a feature
#     node_features = torch.tensor(
#         population, dtype=torch.float32).unsqueeze(-1)  # [n_regions, 1]

#     # Edge indices and weights from mobility matrix
#     mobility_tensor = torch.tensor(mobility_matrix, dtype=torch.float32)
#     edge_index, edge_weight = dense_to_sparse(
#         mobility_tensor)  # Sparse representation

#     return node_features, edge_index, edge_weight


# gnn_model = MobilityGNN(node_feature_dim=1, hidden_dim=32, output_dim=16)
# gnn_model.eval()


def configure_direct_conditions(population, mobility_params, mobility_condition_strategy):

    if mobility_condition_strategy == "accumulated_flows":
        mobility_condition = np.transpose(
            [np.sum(mobility_params, axis=0), np.sum(mobility_params, axis=1)])
        return {"direct_local_conditions": np.concatenate([np.log1p(np.array(population)[:, np.newaxis]), np.log1p(mobility_condition)], axis=1).astype(np.float32),
                "direct_global_conditions": np.log1p(np.sum(population, keepdims=True)).astype(np.float32)}
    elif mobility_condition_strategy == "flow_distributions":
        mobility_condition = np.transpose(
            [np.sum(mobility_params, axis=0), np.sum(mobility_params, axis=1)])
        return {"direct_local_conditions": np.concatenate([np.log1p(np.array(population)[:, np.newaxis]), np.log1p(mobility_condition)], axis=1).astype(np.float32),
                "direct_global_conditions": np.log1p(np.sum(population, keepdims=True)).astype(np.float32)}
    # elif mobility_condition_strategy == "gnn":
    #     # Process the mobility matrix and population with GNN (if necessary)
    #     node_features, edge_index, edge_weight = prepare_graph_data(
    #         mobility_params, population)

    #     # Run the GNN model to get node embeddings
    #     with torch.no_grad():
    #         node_embeddings = gnn_model(node_features, edge_index, edge_weight)
    #     return node_embeddings


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


def run_inference(output_folder_path: str, settings: dict) -> None:

    # create dir if it does not exist
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    # Define synthetic case
    settings["task"]["population"] = [20000000, 22000000, 25000000, 20000000]
    settings["task"]["mobility_params"] = [[0, 20000, 15000, 18000],
                                           [20000, 0, 10000, 14000],
                                           [25000, 22000, 0, 28000],
                                           [5000, 8000, 1000, 0]]

    # ([[lambd0], [I0]], [mu, f_i, phi_i, D_i, scale_I])
    combined_params = ([[0.7, 0.75, 0.6, 0.8], [500, 2000, 500, 100]], [
        6, 0.7, 2.5, 7, 2])
    settings["task"]["direct_conditions"] = configure_direct_conditions(
        settings["task"]["population"], settings["task"]["mobility_params"], settings["task"]["mobility_condition_strategy"])

    # save current settings
    with open(os.path.join(output_folder_path, "settings_mlm.yaml"), "w") as f:
        yaml.dump(settings, f)

    # Define germany states
    # n_regions = 16
    # path_population_data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
    #                                     "..", "data", "pydata", "Germany",
    #                                     "county_current_population.json")
    # population = get_population(path_population_data)

    # prior
    prior = create_prior()
    prior_args = {"local_args": {"n_regions": settings["task"]["n_regions"]}}
    # prior(10, **prior_args)
    prior_scaler = TwoLevelPriorScaler().fit(
        prior, **prior_args)

    # Create GenerativeModel
    simulator_function = partial(
        simulator_SIR, n_regions=settings["task"]["n_regions"], population=settings["task"]["population"], T=settings["task"][
            "T"], mobility_params=settings["task"]["mobility_params"], intervention_model=settings["task"]["intervention_model"],
        observation_model=settings["task"]["observation_model"], param_names=[])

    batch_simulator_fun = parallelize_simulation(
        num_threads=16)(simulator=simulator_function)
    simulator = Simulator(batch_simulator_fun=batch_simulator_fun)
    generative_model = TwoLevelGenerativeModel(
        prior, simulator, skip_test=True, name="sir_covid_simulator")

    # Create Trainer
    amortizer = build_mlm_amortizer(settings, local_params=2, global_params=7)
    trainer = Trainer(amortizer=amortizer, generative_model=generative_model,
                      configurator=partial(
                          configure_input, prior_scaler=prior_scaler, direct_conditions=settings["task"]["direct_conditions"]),
                      checkpoint_path=os.path.join(output_folder_path, "checkpoint"), default_lr=settings["training"]["trainer"]["default_lr"],
                      skip_checks=settings["training"]["trainer"]["skip_checks"], memory=settings["training"]["trainer"]["memory"])

    # generative_model(1, prior_args=prior_args)
    # offline_data = generate_offline_data_splitwise(
    #     output_folder_path, generative_model, 1000000, 100000, generative_model_args={"prior_args": prior_args})
    # offline_data = generate_offline_data(
    #     output_folder_path, generative_model, 1000000, generative_model_args={"prior_args": prior_args})
    offline_data = generate_offline_data(
        output_folder_path, generative_model, settings["task"]["num_samples"], generative_model_args={"prior_args": prior_args})

    # offline_data = generate_offline_data_splitwise(
    #     output_folder_path, generative_model, 1000, 100, generative_model_args={"prior_args": prior_args})

    # Train **kwargs.pop("val_model_args", {})
    # history = start_training(trainer, offline_data=offline_data, epochs=settings["training"]["epochs"],
    #                         batch_size=25000, early_stopping=True, validation_sims=2000, val_model_args={"prior_args": prior_args}, save_logs=False, output_folder_path=output_folder_path)
    history = start_training(trainer, offline_data=offline_data, epochs=settings["training"]["epochs"],
                             batch_size=settings["training"]["batch_size"], early_stopping=settings["training"][
        "early_stopping"], validation_sims=settings["training"]["epochs"], val_model_args={"prior_args": prior_args},
        save_logs=False, output_folder_path=output_folder_path)
    # history = start_training(trainer, epochs=settings["training"]["batch_size"],
    #                          batch_size=3200, iterations_per_epoch=500, early_stopping=True, validation_sims=5000)

    # obs_data: (n_groups, n_time_steps, n_time_series)
    # settings["obs_data"] = load_data_rki_sir(datetime.date(2020, 3, 1), settings["task"]["T"], os.path.join(os.path.dirname(os.path.abspath(
    #     __file__)), "../../../data/pydata/Germany/cases_infected_repdate.json"))
    settings["obs_data"] = load_data_synthetic(
        simulator_function, combined_params=combined_params)

    settings["obs_input"] = {
        "summary_conditions": np.log1p(settings["obs_data"])[np.newaxis].astype(np.float32),
        "direct_local_conditions": [settings["task"]["direct_conditions"]["direct_local_conditions"]],
        "direct_global_conditions": [settings["task"]["direct_conditions"]["direct_global_conditions"]]}

    # history = trainer.loss_history
    MLMPlotting(output_folder_path).plot_all(history, settings, prior, prior_scaler,
                                             simulator_function, generative_model, trainer, prior_args=prior_args)


if __name__ == "__main__":
    with open(os.path.join("settings_mlm.yaml")) as f:
        settings = yaml.safe_load(f)

    run_inference(output_folder_path=os.path.join(
        FILE_PATH, settings["task"]["output_path"], settings["task"]["name"]), settings=settings)
