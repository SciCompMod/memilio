import os
from typing import Self
import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bayesflow.simulation import Prior, TwoLevelPrior
import bayesflow.diagnostics as diag
import tensorflow as tf
import seaborn as sns


class MLMPlotting():
    fig = plt.figure()

    def __init__(self: Self, output_folder_path: os.PathLike = None, overwrite_figures: bool = True):
        self.output_folder_path = output_folder_path
        self.overwrite_figures = overwrite_figures

    @staticmethod
    def figure_wrapper(func):
        def wrapper(self: Self, *args, **kwargs):
            plt.clf()
            plt.close(self.fig)
            self.fig, name = func(self, *args, **kwargs)
            file_path = os.path.join(
                self.output_folder_path, name)
            if self.output_folder_path is None:
                return
            if not self.overwrite_figures and os.path.isfile(file_path):
                print(f"Figure {name} already exists at {file_path}.")
                return
            plt.savefig(file_path)
        return wrapper

    def plot_ppc(self: Self, config, post_samples, ode_model, logscale=True, color="#8f2727", dummy=True, figsize=(12, 6), font_size=18):

        @self.figure_wrapper
        def helper_plot_ppc(self, key, config, sims, obs_data):
            """
            Helper function to perform some plotting of the posterior predictive.
            """
            # Plot settings
            plt.rcParams["font.size"] = font_size

            f, ax = plt.subplots(1, 1, figsize=figsize)

            # Compute quantiles for each t = 1,...,T
            qs_50 = np.quantile(sims, q=[0.25, 0.75], axis=0)
            qs_90 = np.quantile(sims, q=[0.05, 0.95], axis=0)
            qs_95 = np.quantile(sims, q=[0.025, 0.975], axis=0)

            # Plot median predictions and observed data
            ax.plot(np.median(sims, axis=0),
                    label="Median predicted cases", color=color)
            ax.plot(obs_data, marker="o", label="Reported cases",
                    color="black", linestyle="dashed", alpha=0.8)

            # Add compatibility intervals (also called credible intervals)
            ax.fill_between(range(config["T"]), qs_50[0],
                            qs_50[1], color=color, alpha=0.5, label="50% CI")
            ax.fill_between(range(config["T"]), qs_90[0],
                            qs_90[1], color=color, alpha=0.3, label="90% CI")
            ax.fill_between(range(config["T"]), qs_95[0],
                            qs_95[1], color=color, alpha=0.1, label="95% CI")

            # Grid and schmuck
            ax.grid(color="grey", linestyle="-", linewidth=0.25, alpha=0.5)
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)
            ax.set_xlabel("Days since pandemic onset")
            ax.set_ylabel("Number of cases")
            ax.minorticks_off()
            if logscale:
                ax.set_yscale('symlog', linthresh=1)
            ax.legend(fontsize=font_size)
            return f, "ppc_" + key + ".png"

        # Re-simulations
        sims = []
        for params in zip(post_samples["local_parameters"], post_samples["shared_parameters"]):
            # Remove parameters < 0
            if any((arr < 0).any() for arr in params):
                continue
            # Note - simulator returns 2D arrays of shape (T, 1), so we remove trailing dim
            sims.append(ode_model(params)[:, :, 0])
        sims = np.array(sims)

        for region_idx in range(sims.shape[1]):
            helper_plot_ppc(
                self, "region" + str(region_idx), config, sims[:, region_idx, :], config["obs_data"][region_idx, :, 0])

        helper_plot_ppc(self, "global", config, np.sum(sims, axis=1),
                        np.sum(config["obs_data"][:, :, 0], axis=0))
        return

    def plot_prior2d(self: Self, prior: TwoLevelPrior, n_samples=2000, **kwargs):
        """Creates pair-plots for a given joint prior.

        Parameters
        ----------
        prior_samples       : callable
            The prior object which takes a single integer argument and generates random draws.
        param_names : list of str or None, optional, default None
            An optional list of strings which
        n_samples   : int, optional, default: 1000
            The number of random draws from the joint prior
        height      : float, optional, default: 2.5
            The height of the pair plot
        color       : str, optional, default : '#8f2727'
            The color of the plot
        **kwargs    : dict, optional
            Additional keyword arguments passed to the sns.PairGrid constructor

        Returns
        -------
        f : plt.Figure - the figure instance for optional saving
        """

        @self.figure_wrapper
        def helper_plot_prior2d(self, key, prior_samples, param_names=None, height=2.5, color="#8f2727", **kwargs):

            # Get latent dimensionality and prepare titles
            dim = prior_samples.shape[-1]

            # Convert samples to a pandas data frame
            if param_names is None:
                titles = [f"Prior Param. {i}" for i in range(1, dim + 1)]
            else:
                titles = [f"Prior {p}" for p in param_names]
            data_to_plot = pd.DataFrame(prior_samples, columns=titles)

            # Generate plots
            g = sns.PairGrid(data_to_plot, height=height, **kwargs)
            g.map_diag(sns.histplot, fill=True,
                       color=color, alpha=0.9, kde=True)

            # Kernel density estimation (KDE) may not always be possible
            # (e.g. with parameters whose correlation is close to 1 or -1).
            # In this scenario, a scatter-plot is generated instead.
            try:
                g.map_lower(sns.kdeplot, fill=True, color=color, alpha=0.9)
            except Exception as e:
                logging.warning("KDE failed due to the following exception:\n" +
                                repr(e) + "\nSubstituting scatter plot.")
                g.map_lower(sns.scatterplot, alpha=0.6,
                            s=40, edgecolor="k", color=color)
            g.map_upper(sns.scatterplot, alpha=0.6,
                        s=40, edgecolor="k", color=color)

            # Add grids
            for i in range(dim):
                for j in range(dim):
                    g.axes[i, j].grid(alpha=0.5)
            g.tight_layout()

            return g.figure, "prior2d_" + key + ".png"

        # Generate prior draws
        prior_out = prior(n_samples, **kwargs.get("prior_args", {}))
        prior_samples_combined = {}
        prior_samples_combined["hyper_parameters"] = prior_out["hyper_parameters"]
        prior_samples_combined["local_parameters"] = prior_out["local_parameters"]
        prior_samples_combined["shared_parameters"] = prior_out["shared_parameters"]

        for key, prior_samples in prior_samples_combined.items():
            if len(prior_samples.shape) == 3:  # If n_groups exists
                # Combine values across groups: reshape to (batch * n_groups, values)
                prior_samples = prior_samples.transpose(
                    0, 2, 1).reshape(-1, prior_samples.shape[1])

            helper_plot_prior2d(self, key, prior_samples)

        return

    def diagnose_latent2d(self: Self, validation_sims, trainer, **kwargs):

        @self.figure_wrapper
        def helper_diagnose_latent2d(self, key, z, **kwargs):
            f = diag.plot_latent_space_2d(z, **kwargs.pop("plot_args", {}))
            return f, "latent2d_" + key + ".png"

        # shoudl this be per region or just for local and global
        local_out, global_out = trainer.amortizer(validation_sims)
        z = {}
        z["local"] = tf.reshape(local_out[0], (-1, *local_out[0].shape[2:]))
        z["global"] = global_out[0]

        for key, value in z.items():
            helper_diagnose_latent2d(self, key, value, **kwargs)
        return

    @ figure_wrapper
    def plot_losses(self: Self, history):
        f = diag.plot_losses(history["train_losses"],
                             history["val_losses"], moving_average=True)
        return f, "losses.png"

    def diagnose_sbc_histograms(self: Self, validation_sims, trainer, n_samples=None, **kwargs):
        @ self.figure_wrapper
        def helper_diagnose_sbc_histograms(self, key, post_samples, prior_samples, **kwargs):
            # Check for prior names and override keyword if available
            plot_kwargs = kwargs.pop("plot_args", {})
            f = diag.plot_sbc_histograms(
                post_samples, prior_samples, **plot_kwargs)
            return f, "sbc_histogram_" + key + ".png"

        # Heuristically determine the number of posterior samples
        if n_samples is None:
            n_samples = int(
                np.ceil(validation_sims["shared_parameters"].shape[0] / 20))

        # Do inference
        # Generate posterior draws for all simulations
        post_samples = {"global_samples": [],
                        "local_samples": []}
        for i in range(validation_sims["shared_parameters"].shape[0]):
            input_dict = {}
            for key in validation_sims.keys():
                input_dict[key] = validation_sims[key][i:i+1].copy()
            post_sample = trainer.amortizer.sample(
                input_dict, n_samples=n_samples)
            post_samples["global_samples"].append(
                post_sample["global_samples"])
            post_samples["local_samples"].append(
                post_sample["local_samples"])

        sbc_input = {}
        local_shape = np.shape(post_samples["local_samples"])
        sbc_input["local_parameters"] = (
            np.transpose(post_samples["local_samples"], (0, 2, 3, 1)).reshape((local_shape[0], local_shape[2], -1)), np.reshape(validation_sims["local_parameters"], (local_shape[0], -1)))
        sbc_input["global_parameters"] = (
            post_samples["global_samples"], np.concatenate((validation_sims["hyper_parameters"], validation_sims["shared_parameters"]), axis=1))
        for key, (post_sample, prior_sample) in sbc_input.items():
            helper_diagnose_sbc_histograms(
                self, key, np.array(post_sample), np.array(prior_sample), **kwargs)
        return

    def plot_sbc_ecdf(self: Self, post_samples, validation_sims, prior: Prior):

        @ self.figure_wrapper
        def helper_plot_sbc_ecdf(self, key, post_sample, prior_sample):
            f = diag.plot_sbc_ecdf(
                post_sample, prior_sample)
            return f, "sbc_ecdf_" + key + ".png"

        sbc_input = {}
        local_shape = np.shape(post_samples["local_samples"])
        sbc_input["local_parameters"] = (
            np.transpose(post_samples["local_samples"], (0, 2, 3, 1)).reshape((local_shape[0], local_shape[2], -1)), np.reshape(validation_sims["local_parameters"], (local_shape[0], -1)))
        sbc_input["global_parameters"] = (
            post_samples["global_samples"], np.concatenate((validation_sims["hyper_parameters"], validation_sims["shared_parameters"]), axis=1))
        for key, (post_sample, prior_sample) in sbc_input.items():
            helper_plot_sbc_ecdf(self, key, np.array(
                post_sample), np.array(prior_sample))
        return

    def plot_sbc_ecdf_stacked(self: Self, post_samples, validation_sims):

        @ self.figure_wrapper
        def helper_plot_sbc_ecdf_stacked(self, key, post_sample, prior_sample):
            f = diag.plot_sbc_ecdf(
                post_sample, prior_sample, stacked=True, difference=True, legend_fontsize=12, fig_size=(6, 5))
            return f, "sbc_ecdf_stacked_" + key + ".png"

        sbc_input = {}
        local_shape = np.shape(post_samples["local_samples"])
        sbc_input["local_parameters"] = (
            np.transpose(post_samples["local_samples"], (0, 2, 3, 1)).reshape((local_shape[0], local_shape[2], -1)), np.reshape(validation_sims["local_parameters"], (local_shape[0], -1)))
        sbc_input["global_parameters"] = (
            post_samples["global_samples"], np.concatenate((validation_sims["hyper_parameters"], validation_sims["shared_parameters"]), axis=1))
        for key, (post_sample, prior_sample) in sbc_input.items():
            helper_plot_sbc_ecdf_stacked(self, key, np.array(
                post_sample), np.array(prior_sample))
        return

    def plot_sbc_histograms(self: Self, post_samples, validation_sims, prior: Prior):

        @ self.figure_wrapper
        def helper_plot_sbc_histograms(self, key, post_sample, prior_sample):
            f = diag.plot_sbc_histograms(
                post_sample, prior_sample)
            return f, "sbc_ecdf_histograms_" + key + ".png"

        sbc_input = {}
        local_shape = np.shape(post_samples["local_samples"])
        sbc_input["local_parameters"] = (
            np.transpose(post_samples["local_samples"], (0, 2, 3, 1)).reshape((local_shape[0], local_shape[2], -1)), np.reshape(validation_sims["local_parameters"], (local_shape[0], -1)))
        sbc_input["global_parameters"] = (
            post_samples["global_samples"], np.concatenate((validation_sims["hyper_parameters"], validation_sims["shared_parameters"]), axis=1))
        for key, (post_sample, prior_sample) in sbc_input.items():
            helper_plot_sbc_histograms(self, key, np.array(
                post_sample), np.array(prior_sample))
        return

    def plot_recovery(self: Self,  post_samples, validation_sims, prior):

        @ self.figure_wrapper
        def helper_plot_recovery(self, key, post_sample, prior_sample):
            f = diag.plot_recovery(
                post_sample, prior_sample)
            return f, "recovery_" + key + ".png"

        recovery_input = {}
        local_shape = np.shape(post_samples["local_samples"])
        recovery_input["local_parameters"] = (
            np.transpose(post_samples["local_samples"], (0, 2, 3, 1)).reshape((local_shape[0], local_shape[2], -1)), np.reshape(validation_sims["local_parameters"], (local_shape[0], -1)))
        recovery_input["global_parameters"] = (
            post_samples["global_samples"], np.concatenate((validation_sims["hyper_parameters"], validation_sims["shared_parameters"]), axis=1))
        for key, (post_sample, prior_sample) in recovery_input.items():
            helper_plot_recovery(self, key, np.array(
                post_sample), np.array(prior_sample))
        return

    def plot_z_score_contraction(self: Self, post_samples, validation_sims, prior):

        @ self.figure_wrapper
        def helper_plot_z_score_contraction(self, key, post_sample, prior_sample):
            f = diag.plot_z_score_contraction(
                post_sample, prior_sample)
            return f, "z_score_contraction_" + key + ".png"

        zscore_input = {}
        local_shape = np.shape(post_samples["local_samples"])
        zscore_input["local_parameters"] = (
            np.transpose(post_samples["local_samples"], (0, 2, 3, 1)).reshape((local_shape[0], local_shape[2], -1)), np.reshape(validation_sims["local_parameters"], (local_shape[0], -1)))
        zscore_input["global_parameters"] = (
            post_samples["global_samples"], np.concatenate((validation_sims["hyper_parameters"], validation_sims["shared_parameters"]), axis=1))
        for key, (post_sample, prior_sample) in zscore_input.items():
            helper_plot_z_score_contraction(
                self, key, np.array(post_sample), np.array(prior_sample))

    def plot_posterior_2d(self: Self, post_samples, prior: Prior, with_prior: bool = False, **kwargs):

        @ self.figure_wrapper
        def helper_plot_posterior_2d(self, key, single_post_samples, single_prior_draws):
            f = diag.plot_posterior_2d(
                single_post_samples, prior_draws=single_prior_draws)
            return f, "posterior_2d_" + key + ".png"

        # Generate prior draws
        prior_out = {}
        if with_prior:
            prior_out = prior(
                post_samples["hyper_parameters"].shape[0], **kwargs.get("prior_args", {}))
        post_samples_combined = {}
        post_samples_combined["hyper_parameters"] = post_samples["hyper_parameters"]
        post_samples_combined["local_parameters"] = post_samples["local_parameters"]
        post_samples_combined["shared_parameters"] = post_samples["shared_parameters"]

        for key, single_post_samples in post_samples_combined.items():
            if len(single_post_samples.shape) == 3:  # If n_groups exists
                # Plot values for each group:
                for i in range(single_post_samples.shape[2]):
                    if with_prior:
                        helper_plot_posterior_2d(
                            self, "with_prior_" + key + "_region" + str(i), single_post_samples[:, :, i], prior_out[key][:, :, i])
                    else:
                        helper_plot_posterior_2d(
                            self, key + "_region" + str(i), single_post_samples[:, :, i], {})
            else:
                helper_plot_posterior_2d(
                    self, key, single_post_samples, prior_out.get(key, {}))

        return

    def plot_all(self: Self, history, config, prior, prior_scaler, simulator_function, generative_model, trainer, **kwargs):
        self.plot_prior2d(prior, **kwargs)

        self.plot_losses(history)

        # Generate some validation data
        batch_size = 16000  # batch_size*50
        validation_sims = trainer.configurator(
            generative_model(batch_size=batch_size, prior_args=kwargs.get("prior_args", {})))

        self.diagnose_latent2d(validation_sims, trainer)

        self.diagnose_sbc_histograms(validation_sims, trainer)

        # Generate some validation data
        batch_size = 3000
        validation_sims = trainer.configurator(
            generative_model(batch_size=batch_size, prior_args=kwargs.get("prior_args", {})))

        # Generate posterior draws for all simulations
        post_samples = {"global_samples": [],
                        "local_samples": []}
        for i in range(batch_size):
            input_dict = {}
            for key in validation_sims.keys():
                input_dict[key] = validation_sims[key][i:i+1].copy()
            post_sample = trainer.amortizer.sample(input_dict, n_samples=100)
            post_samples["global_samples"].append(
                post_sample["global_samples"])
            post_samples["local_samples"].append(
                post_sample["local_samples"])

        # Create ECDF plot
        self.plot_sbc_ecdf(
            post_samples, validation_sims, prior)

        self.plot_sbc_ecdf_stacked(
            post_samples, validation_sims,)

        self.plot_sbc_histograms(
            post_samples, validation_sims, prior)

        # Generate posterior draws for all simulations
        post_samples = {"global_samples": [],
                        "local_samples": []}
        for i in range(batch_size):
            input_dict = {}
            for key in validation_sims.keys():
                input_dict[key] = validation_sims[key][i:i+1].copy()
            post_sample = trainer.amortizer.sample(input_dict, n_samples=1000)
            post_samples["global_samples"].append(
                post_sample["global_samples"])
            post_samples["local_samples"].append(
                post_sample["local_samples"])
        self.plot_recovery(
            post_samples, validation_sims, prior)

        self.plot_z_score_contraction(
            post_samples, validation_sims, prior)

        # Local samples will hold an array-like of shape (num_regions, num_samples, num_local)
        # and global samples will hold an array-like of shape (num_samples, num_hyper + num_shared)
        post_samples = trainer.amortizer.sample(config["obs_input"], 5000)

        # Update the dictionary
        local_samples = post_samples["local_samples"]
        global_samples = post_samples["global_samples"]

        # Overwrite with new keys
        num_hyper = 2
        post_samples = {
            # Renamed local_samples
            "local_parameters": local_samples.transpose((1, 2, 0)),
            # First num_hyper columns
            "hyper_parameters": global_samples[:, :num_hyper],
            # Remaining columns
            "shared_parameters": global_samples[:, num_hyper:],
        }
        # Undo standardization to get parameters on their original (unstandardized) scales
        post_samples = prior_scaler.inverse_transform(post_samples)

        self.plot_posterior_2d(post_samples, prior)

        self.plot_posterior_2d(post_samples, prior, with_prior=True, **kwargs)

        self.plot_ppc(config, post_samples, simulator_function)


class Plotting():
    fig = plt.figure()

    def __init__(self: Self, output_folder_path: os.PathLike = None, overwrite_figures: bool = True):
        self.output_folder_path = output_folder_path
        self.overwrite_figures = overwrite_figures

    @ staticmethod
    def figure_wrapper(func):
        def wrapper(self: Self, *args, **kwargs):
            plt.clf()
            plt.close(self.fig)
            self.fig, name = func(self, *args, **kwargs)
            file_path = os.path.join(
                self.output_folder_path, name)
            if self.output_folder_path is None:
                return
            if not self.overwrite_figures and os.path.isfile(file_path):
                print(f"Figure {name} already exists at {file_path}.")
                return
            plt.savefig(file_path)
        return wrapper

    @ figure_wrapper
    def plot_ppc(self: Self, config, post_samples, ode_model, logscale=True, color="#8f2727", dummy=True, figsize=(12, 6), font_size=18):
        """
        Helper function to perform some plotting of the posterior predictive.
        """
        # Plot settings
        plt.rcParams["font.size"] = font_size

        # Remove parameters < 0
        samples = post_samples[np.sum(post_samples < 0, axis=1) == 0]

        f, ax = plt.subplots(1, 1, figsize=figsize)

        # Re-simulations
        sims = []
        for i in range(samples.shape[0]):
            # Note - simulator returns 2D arrays of shape (T, 1), so we remove trailing dim
            sims.append(ode_model(samples[i])[:, 0])
        sims = np.array(sims)

        # Compute quantiles for each t = 1,...,T
        qs_50 = np.quantile(sims, q=[0.25, 0.75], axis=0)
        qs_90 = np.quantile(sims, q=[0.05, 0.95], axis=0)
        qs_95 = np.quantile(sims, q=[0.025, 0.975], axis=0)

        # Plot median predictions and observed data
        ax.plot(np.median(sims, axis=0),
                label="Median predicted cases", color=color)
        ax.plot(config["obs_data"], marker="o", label="Reported cases",
                color="black", linestyle="dashed", alpha=0.8)

        # Add compatibility intervals (also called credible intervals)
        ax.fill_between(range(config["T"]), qs_50[0],
                        qs_50[1], color=color, alpha=0.5, label="50% CI")
        ax.fill_between(range(config["T"]), qs_90[0],
                        qs_90[1], color=color, alpha=0.3, label="90% CI")
        ax.fill_between(range(config["T"]), qs_95[0],
                        qs_95[1], color=color, alpha=0.1, label="95% CI")

        # Grid and schmuck
        ax.grid(color="grey", linestyle="-", linewidth=0.25, alpha=0.5)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.set_xlabel("Days since pandemic onset")
        ax.set_ylabel("Number of cases")
        ax.minorticks_off()
        if logscale:
            ax.set_yscale("log")
        ax.legend(fontsize=font_size)
        return f, "ppc.png"

    @ figure_wrapper
    def plot_prior2d(self: Self, prior: Prior):
        f = prior.plot_prior2d()
        return f, "prior2D.png"

    @ figure_wrapper
    def diagnose_latent2d(self: Self, trainer):
        f = trainer.diagnose_latent2d()
        return f, "latent2d.png"

    @ figure_wrapper
    def plot_losses(self: Self, history):
        f = diag.plot_losses(history["train_losses"],
                             history["val_losses"], moving_average=True)
        return f, "losses.png"

    @ figure_wrapper
    def diagnose_sbc_histograms(self: Self, trainer):
        f = trainer.diagnose_sbc_histograms()
        return f, "sbc_histograms.png"

    @ figure_wrapper
    def plot_sbc_ecdf(self: Self, post_samples, validation_sims, prior: Prior):
        f = diag.plot_sbc_ecdf(
            post_samples, validation_sims["parameters"], param_names=prior.param_names)
        return f, "sbc_ecdf.png"

    @ figure_wrapper
    def plot_sbc_ecdf_stacked(self: Self, post_samples, validation_sims):
        f = diag.plot_sbc_ecdf(
            post_samples, validation_sims["parameters"], stacked=True, difference=True, legend_fontsize=12, fig_size=(6, 5)
        )
        return f, "sbc_ecdf_stacked.png"

    @ figure_wrapper
    def plot_sbc_histograms(self: Self, post_samples, validation_sims, prior: Prior):
        f = diag.plot_sbc_histograms(
            post_samples, validation_sims["parameters"], param_names=prior.param_names)
        return f, "sbc_ecdf_histograms.png"

    @ figure_wrapper
    def plot_recovery(self: Self,  post_samples, validation_sims, prior):
        f = diag.plot_recovery(
            post_samples, validation_sims["parameters"], param_names=prior.param_names)
        return f, "recovery.png"

    @ figure_wrapper
    def plot_z_score_contraction(self: Self, post_samples, validation_sims, prior):
        f = diag.plot_z_score_contraction(
            post_samples, validation_sims["parameters"], param_names=prior.param_names)
        return f, "z_score_contraction.png"

    @ figure_wrapper
    def plot_posterior_2d(self: Self, post_samples, prior: Prior, with_prior: bool = False):
        if with_prior:
            f = diag.plot_posterior_2d(post_samples, prior=prior)
            return f, "posterior_2d_with_prior.png"
        else:
            f = diag.plot_posterior_2d(
                post_samples, param_names=prior.param_names)
            return f, "posterior_2d.png"

    def plot_all(self: Self, history, config, prior, prior_scaler, simulator_function, generative_model, trainer):
        # self.plot_prior2d(prior)

        # self.diagnose_latent2d(trainer)

        self.plot_losses(history)

        # self.diagnose_sbc_histograms(trainer)

        # # Generate some validation data
        # validation_sims = trainer.configurator(
        #     generative_model(batch_size=3000))

        # # Generate posterior draws for all simulations
        # post_samples = trainer.amortizer.sample(validation_sims, n_samples=100)

        # # Create ECDF plot
        # self.plot_sbc_ecdf(
        #     post_samples, validation_sims, prior)

        # self.plot_sbc_ecdf_stacked(
        #     post_samples, validation_sims,)

        # self.plot_sbc_histograms(
        #     post_samples, validation_sims, prior)

        # post_samples = trainer.amortizer.sample(validation_sims, n_samples=1000)
        # self.plot_recovery(
        #     post_samples, validation_sims, prior)

        # self.plot_z_score_contraction(
        #     post_samples, validation_sims, prior)

        # # Format data into a 3D array of shape (1, n_time_steps, 1) and perform log transform
        # obs_data = np.log1p(config["obs_data"])[
        #     np.newaxis, :, np.newaxis].astype(np.float32)

        # # vec_num_obs = config["obs_data"].shape[0] * \
        # #     np.ones((obs_data.shape[0], 1))
        # # # Obtain 500 posterior draws given real data
        # # post_samples = trainer.amortizer.sample({"summary_conditions": obs_data, "direct_conditions": np.sqrt(
        # #     vec_num_obs).astype(np.float32)}, 5000)

        # post_samples = trainer.amortizer.sample({"summary_conditions": obs_data}, 5000)
        # # Undo standardization to get parameters on their original (unstandardized) scales
        # post_samples = prior_scaler.inverse_transform(post_samples)

        # self.plot_posterior_2d(post_samples, prior)

        # self.plot_posterior_2d(post_samples, prior, with_prior=True)

        # self.plot_ppc(config, post_samples, simulator_function)
