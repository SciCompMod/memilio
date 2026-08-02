import os
import pandas as pd
import matplotlib.pyplot as plt


def plot_column_means(csv_path: str, time_col: str = "Time", out_path: str | None = None):
    df = pd.read_csv(csv_path)

    # Drop the time column (and any other non-numeric columns) before averaging
    data_cols = [c for c in df.columns if c != time_col]
    means = df[data_cols].mean()

    fig, axes = plt.subplots(1, 1, figsize=(12, 5))

    # --- Bar chart of means ---
    axes.bar(means.index, means.values, color="steelblue")
    axes.set_title("Mean per column")
    axes.set_ylabel("Mean value")
    axes.tick_params(axis="x", rotation=45)
    for label in axes.get_xticklabels():
        label.set_ha("right")

    # # --- Box plot: mean + spread per column (alternative view) ---
    # axes[1].boxplot([df[c].dropna().values for c in data_cols], tick_labels=data_cols, showmeans=True)
    # axes[1].set_title("Distribution per column (mean shown as green triangle)")
    # axes[1].tick_params(axis="x", rotation=45)
    # for label in axes[1].get_xticklabels():
    #     label.set_ha("right")

    fig.tight_layout()

    if out_path:
        if not os.path.isdir(out_path):
            os.makedirs(out_path)
        fig.savefig(os.path.join(out_path, "benchmark_results.png"), dpi=150)
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
                                    f"../plots/2026-07-30/benachmark_results/")
   
    time_col = "Time"


    plot_column_means(csv_path, time_col, plot_dir)