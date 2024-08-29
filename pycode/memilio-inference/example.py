import datetime
from functools import partial

import numpy as np
import pandas as pd
from scipy import stats
import pickle
import matplotlib.pyplot as plt

from plots import *
from prior import ModelPriorBuilder
from sir import *

import bayesflow.diagnostics as diag
from bayesflow.amortizers import AmortizedPosterior
from bayesflow.networks import InvertibleNetwork, SequenceNetwork
from bayesflow.simulation import GenerativeModel, Simulator
from bayesflow.trainers import Trainer

RNG = np.random.default_rng(2023)


def load_data(T):
    """Helper function to load cumulative cases and transform them to new cases."""

    # Use right corona data based on the model (either reporting or reference date)
    confirmed_cases_json = "/home/betz_mx/Project_Memilio/experimental-master/data/pydata/Germany/cases_infected_repdate.json"
    confirmed_cases = pd.read_json(confirmed_cases_json)
    confirmed_cases = confirmed_cases.set_index('Date')

    date_data_begin = datetime.date(2020, 3, 1)
    date_data_end = date_data_begin + datetime.timedelta(T)

    cases_obs = np.array(
        confirmed_cases.loc[date_data_begin:date_data_end]
    ).flatten()
    new_cases_obs = np.diff(cases_obs)
    return new_cases_obs


if __name__ == "__main__":

    # Building a SEIR model
    sir_builder = ModelPriorBuilder(SIRStrategy())
    sir_builder.add_base().add_intervention().add_observation()

    fixed_parameters = {ParameterNamesSir.I0.value: 40,
                        ParameterNamesSir.MU.value: 0.5}
    prior = sir_builder.set_fixed_parameters(fixed_parameters).build()
    prior_means, prior_stds = prior.estimate_means_and_stds()
    for idx, prior_std in enumerate(prior_stds[0]):
        if prior_std == 0:
            prior_stds[0, idx] += 10 ** -14

    stationary_SIR(prior(1)["prior_draws"][0], 80000, 81, True,
                   param_names=prior.param_names, fixed_params=fixed_parameters)

    def configure_input(forward_dict):
        """Function to configure the simulated quantities (i.e., simulator outputs)
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
        params = (params - prior_means) / prior_stds

        # Remove a batch if it contains nan, inf or -inf
        idx_keep = np.all(np.isfinite(logdata), axis=(1, 2))
        if not np.all(idx_keep):
            print("Invalid value encountered...removing from batch")

        # Add to keys
        out_dict["summary_conditions"] = logdata[idx_keep]
        out_dict["parameters"] = params[idx_keep]

        return out_dict

    T = 81
    config = {"T": T, "N": 83e6, "observation_model": True,
              "obs_data": load_data(T)}

    simulator = Simulator(simulator_fun=partial(
        stationary_SIR, T=config["T"], N=config["N"], observation_model=config["observation_model"], param_names=prior.param_names, fixed_params=fixed_parameters))
    model = GenerativeModel(prior, simulator, name="sir_covid_simulator")

    # params = model_prior()
    # params = np.array([ 9.48431235e+00,  1.55999015e+01,  1.96929563e+01,  6.69060727e+01,
    #     2.68912110e-01,  6.37664326e-01,  5.61281235e-01,  6.49926284e-03,
    #     7.13581645e-03,  7.99039524e+00,  7.93024437e-01, -2.67922548e+00,
    #     7.68616098e+00,  8.62828212e+01,  3.21992246e+00])
    # params = np.array([8.18585338e+00, 1.37430922e+01, 2.35805176e+01, 6.68899973e+01,
    #                 7.78653693e-01, 3.21214143e-01, 8.36588495e-01, 1.20224499e-02, 1.32240528e-02,
    #                 9.53819259e+00, 4.83260432e-01, 1.65663618e+00, 9.15314481e+00, 1.60834143e+02, 1.03021437e+01])
    # result = stationary_SIR(params, config["N"], config["T"], observation_model=False)

    # f = plt.plot(result)
    # plt.savefig("results.png")

    offline_data = model(1000)
    with open('offline_obs_data.pkl', 'wb') as file:
        pickle.dump(offline_data, file)

    # with open('offline_obs_data.pkl', 'rb') as file:
    #     offline_data = pickle.load(file)

    # Should not be lower than number of parameters
    summary_net = SequenceNetwork(summary_dim=15)
    # could try setting coupling_design=spline for more superior performance on lower-dimensional problems.
    inference_net = InvertibleNetwork(num_params=len(
        prior.param_names), num_coupling_layers=6)
    amortizer = AmortizedPosterior(
        inference_net, summary_net, name="covid_amortizer")
    trainer = Trainer(amortizer=amortizer, generative_model=model,
                      configurator=configure_input, memory=True)

    # Train
    history = trainer.train_offline(
        offline_data, epochs=5, batch_size=320, validation_sims=200)

    f = diag.plot_losses(history["train_losses"],
                         history["val_losses"], moving_average=True)
    plt.savefig("losses.png")

    # Format data into a 3D array of shape (1, n_time_steps, 1) and perform log transform
    obs_data = np.log1p(config["obs_data"])[
        np.newaxis, :, np.newaxis].astype(np.float32)
    # Obtain 500 posterior draws given real data
    post_samples = amortizer.sample({"summary_conditions": obs_data}, 5000)
    # Undo standardization to get parameters on their original (unstandardized) scales
    post_samples = prior_means + post_samples * prior_stds

    f = diag.plot_posterior_2d(post_samples, prior=prior)
    plt.savefig("posterior_2d.png")

    f = plot_ppc(config, post_samples, stationary_SIR)
    plt.savefig("ppc.png")
