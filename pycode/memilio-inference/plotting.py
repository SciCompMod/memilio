import os
from typing import Self

import matplotlib.pyplot as plt
import numpy as np

from bayesflow.simulation import Prior
import bayesflow.diagnostics as diag


class Plotting():
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

    @figure_wrapper
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

    @figure_wrapper
    def plot_prior2d(self: Self, prior: Prior):
        f = prior.plot_prior2d()
        return f, "prior2D.png"

    @figure_wrapper
    def diagnose_latent2d(self: Self, trainer):
        f = trainer.diagnose_latent2d()
        return f, "latent2d.png"

    @figure_wrapper
    def plot_losses(self: Self, history):
        f = diag.plot_losses(history["train_losses"],
                             history["val_losses"], moving_average=True)
        return f, "losses.png"

    @figure_wrapper
    def diagnose_sbc_histograms(self: Self, trainer):
        f = trainer.diagnose_sbc_histograms()
        return f, "sbc_histograms.png"

    @figure_wrapper
    def plot_sbc_ecdf(self: Self, post_samples, validation_sims, prior: Prior):
        f = diag.plot_sbc_ecdf(
            post_samples, validation_sims["parameters"], param_names=prior.param_names)
        return f, "sbc_ecdf.png"

    @figure_wrapper
    def plot_sbc_ecdf_stacked(self: Self, post_samples, validation_sims):
        f = diag.plot_sbc_ecdf(
            post_samples, validation_sims["parameters"], stacked=True, difference=True, legend_fontsize=12, fig_size=(6, 5)
        )
        return f, "sbc_ecdf_stacked.png"

    @figure_wrapper
    def plot_sbc_histograms(self: Self, post_samples, validation_sims, prior: Prior):
        f = diag.plot_sbc_histograms(
            post_samples, validation_sims["parameters"], param_names=prior.param_names)
        return f, "sbc_ecdf_histograms.png"

    @figure_wrapper
    def plot_recovery(self: Self,  post_samples, validation_sims, prior):
        f = diag.plot_recovery(
            post_samples, validation_sims["parameters"], param_names=prior.param_names)
        return f, "recovery.png"

    @figure_wrapper
    def plot_z_score_contraction(self: Self, post_samples, validation_sims, prior):
        f = diag.plot_z_score_contraction(
            post_samples, validation_sims["parameters"], param_names=prior.param_names)
        return f, "z_score_contraction.png"

    @figure_wrapper
    def plot_posterior_2d(self: Self, post_samples, prior: Prior, with_prior: bool = False):
        if with_prior:
            f = diag.plot_posterior_2d(post_samples, prior=prior)
            return f, "posterior_2d_with_prior.png"
        else:
            f = diag.plot_posterior_2d(
                post_samples, param_names=prior.param_names)
            return f, "posterior_2d.png"

    def plot_all(self: Self, history, config, prior, prior_scaler, simulator_function, generative_model, trainer, amortizer):
        self.plot_prior2d(prior)

        self.diagnose_latent2d(trainer)

        self.plot_losses(history)

        self.diagnose_sbc_histograms(trainer)

        # Generate some validation data
        validation_sims = trainer.configurator(
            generative_model(batch_size=3000))

        # Generate posterior draws for all simulations
        post_samples = amortizer.sample(validation_sims, n_samples=100)

        # Create ECDF plot
        self.plot_sbc_ecdf(
            post_samples, validation_sims, prior)

        self.plot_sbc_ecdf_stacked(
            post_samples, validation_sims,)

        self.plot_sbc_histograms(
            post_samples, validation_sims, prior)

        post_samples = amortizer.sample(validation_sims, n_samples=1000)
        self.plot_recovery(
            post_samples, validation_sims, prior)

        self.plot_z_score_contraction(
            post_samples, validation_sims, prior)

        # Format data into a 3D array of shape (1, n_time_steps, 1) and perform log transform
        obs_data = np.log1p(config["obs_data"])[
            np.newaxis, :, np.newaxis].astype(np.float32)

        # vec_num_obs = config["obs_data"].shape[0] * \
        #     np.ones((obs_data.shape[0], 1))
        # # Obtain 500 posterior draws given real data
        # post_samples = amortizer.sample({"summary_conditions": obs_data, "direct_conditions": np.sqrt(
        #     vec_num_obs).astype(np.float32)}, 5000)

        post_samples = amortizer.sample({"summary_conditions": obs_data}, 5000)
        # Undo standardization to get parameters on their original (unstandardized) scales
        post_samples = prior_scaler.inverse_transform(post_samples)

        self.plot_posterior_2d(post_samples, prior)

        self.plot_posterior_2d(post_samples, prior, with_prior=True)

        self.plot_ppc(config, post_samples, simulator_function)
