import matplotlib.pyplot as plt
import numpy as np


def plot_ppc(config, post_samples, ode_model, logscale=True, color="#8f2727", dummy=True, figsize=(12, 6), font_size=18):
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
    ax.plot(np.median(sims, axis=0), label="Median predicted cases", color=color)
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
    return f
