import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


def plot_column_means(csv_path: str, time_col: str = "Time", out_path: str | None = None):
    df = pd.read_csv(csv_path)

    # Drop the time column (and any other non-numeric columns) before averaging
    data_cols = [c for c in df.columns if c != time_col]
    means = df[data_cols].mean()

    fig, axes = plt.subplots(1, 1, figsize=(16, 8.5))
    axes.yaxis.grid(True, linestyle="--", alpha=0.4)
    axes.set_axisbelow(True)

    n_bars = len(means)
    colors = ["#332288", "#44AA99", "#882255", "#88CCEE",  "#2E8B57"]
    colors = colors[::-1]
    x_pos = np.arange(n_bars)
    bars = axes.bar(x_pos, means.values, color=colors, width=0.6,
                    edgecolor="white", linewidth=0.8, alpha=1.)

    # --- Value labels on top of each bar ---
    for rect, val in zip(bars, means.values):
        axes.annotate(
            f"{val:.2f}",
            xy=(rect.get_x() + rect.get_width() / 2, rect.get_height()),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=20,
        )

    # --- Remove x-tick labels; use a legend instead ---
    axes.set_xticks(x_pos)
    axes.set_xticklabels([])
    axes.tick_params(axis="x", length=0)  # hide tick marks too

    axes.tick_params(axis='y', labelsize=16)

    axes.set_title("Run time", fontsize=32, pad=12)
    axes.set_ylabel("Time [seconds]", fontsize=32)

    # Light horizontal gridlines behind the bars, no top/right spines

    # axes.spines["top"].set_visible(False)
    # axes.spines["right"].set_visible(False)

    # Legend on the right-hand side of the plot, one entry per bar/column

    labels = ["Smoother cosine", r"$\mathcal{C^1}$ smoothstep", r"$\mathcal{C^2}$ smoothstep",
              r"$\mathcal{C^3}$ smoothstep", r"$\mathcal{C^4}$ smoothstep"]
    axes.legend(
        bars,
        labels,
        title=None,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        # frameon=False,
        fontsize=28,
        ncols=1
    )

    fig.tight_layout()

    if out_path:
        if not os.path.isdir(out_path):
            os.makedirs(out_path)
        fig.savefig(os.path.join(out_path, "benchmark_results.png"),
                    dpi=150, bbox_inches="tight")
        print(f"Saved plot to {out_path}")
    else:
        plt.show()

    # Also print the numeric means for reference
    print("\nColumn means:")
    print(means.to_string())

    return means


if __name__ == "__main__":

    csv_path = os.path.join(os.path.dirname(__file__),
                            f"../simulation_results/2026-07-30/benchmark_results/smoother_times_20warmupruns_100runs_t0=0_tmax=100_numsteps=100000001.csv")
    plot_dir = os.path.join(os.path.dirname(__file__),
                            f"../plots/2026-07-30/benchmark_results/")

    time_col = "Time"

    plot_column_means(csv_path, time_col, plot_dir)
