import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

from plotting_settings import *

model_colors = {"ODE": colors["Light green"], "LCT": colors["Rose"],
                "IDE": colors["Light blue"], "ABM": colors["Teal"]}

set_fontsize()

labelpad = 10


def plot_three_pdfs(exp_params, erlang_params, lognorm_params, xmax, xlabel, ylabel, save_dir=None):

    xmax = 12.

    x = np.linspace(0, xmax, 1000)
    plt.figure(figsize=(7, 4))

    # Exponential PDF
    plt.plot(
        x,
        scipy.stats.expon.pdf(x, scale=exp_params[0]),
        label='Exponential', color=model_colors["ODE"])

    # print("Mean Exponential: ", scipy.stats.expon.stats(
    #     scale=exp_params[0], moments='mv'))

    # Erlang PDF (special case of gamma: shape=k, scale=theta)
    plt.plot(
        x,
        scipy.stats.gamma.pdf(x, a=erlang_params[0], scale=erlang_params[1]),
        label='Erlang', color=model_colors["LCT"])

    # print("Mean Erlang: ", scipy.stats.gamma.stats(
    #     a=erlang_params[0], scale=erlang_params[1], moments='mv'))

    # Lognormal PDF
    plt.plot(
        x,
        scipy.stats.lognorm.pdf(
            x, s=lognorm_params[0], scale=lognorm_params[1]),
        label='Lognormal', color=model_colors["IDE"])

    # print("Mean Lognormal: ", scipy.stats.lognorm.stats(
    #     s=lognorm_params[0], scale=lognorm_params[1],  moments='mv'))

    plt.xlabel(xlabel, labelpad=labelpad)
    plt.ylabel(ylabel, labelpad=labelpad)
    plt.xlim(0, xmax)
    plt.ylim(0, None)
    plt.title("Scenario 2", pad=labelpad)
    plt.legend()
    plt.tight_layout()

    if save_dir:

        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)

        plt.savefig(os.path.join(save_dir, "all_pdfs" + ".png"), dpi=dpi)


def plot_exponential_pdf(exp_params, xmax, xlabel, ylabel, save_dir=None):

    x = np.linspace(0, xmax, 1000)
    plt.figure(figsize=(7, 4))

    # Exponential PDF
    plt.plot(
        x,
        scipy.stats.expon.pdf(x, scale=exp_params[0]),
        label='Exponential', color=model_colors["ODE"])

    # set_fontsize()

    plt.xlabel(xlabel, labelpad=labelpad)
    plt.ylabel(
        ylabel, labelpad=labelpad)
    plt.xlim(0, xmax)
    plt.ylim(0, None)
    plt.title("Scenario 1", pad=labelpad)
    plt.legend()
    plt.tight_layout()

    if save_dir:

        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)

        plt.savefig(os.path.join(
            save_dir, "exponential_pdf" + ".png"), dpi=dpi)


if __name__ == "__main__":

    # Path where plots will be stored.
    plot_dir = os.path.join(os.path.dirname(
        __file__), "../../..", "plots/pdfs/")

    # Plot distributions that are used for transition from Exposed to Carrier compartment.
    exp_params = [4.5]                          # scale
    erlang_params = [9, 4.5/9.]                 # shape, scale = mean/shape
    lognorm_params = [0.32459285, 4.26907484]   # shape, scale

    # Define some parameters for plot.
    xmax = 12.
    xlabel = "Time since transmission"
    ylabel = "PDF"

    # PDFs for Scenario 1
    plot_exponential_pdf(exp_params, xmax, xlabel, ylabel, save_dir=plot_dir)

    # PDFs for Scenario 2
    plot_three_pdfs(exp_params, erlang_params,
                    lognorm_params, xmax, xlabel, ylabel, save_dir=plot_dir)
