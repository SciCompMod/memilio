import pandas as pd
import matplotlib.patches as mpatches
import numpy as np

from plotting_settings import *

def plot_ts_output(axis, filename, percentiles, color, compartments):
    while len(percentiles) > 0:
        if len(percentiles) == 1:
            df = pd.read_csv(filename + f"_p{percentiles[0]}.csv")
            axis.plot(df["Time"], df[compartments].sum(axis=1), color=color)
            del percentiles[0]
        else:
            df_low = pd.read_csv(filename + f"_p{percentiles[0]}.csv")
            df_high = pd.read_csv(filename + f"_p{percentiles[-1]}.csv")
            axis.fill_between(df_low["Time"], df_low[compartments].sum(axis=1), df_high[compartments].sum(axis=1), alpha=0.3, color=color)
            del percentiles[0]
            del percentiles[-1]
    return axis

result_dir = "/home/bick_ju//Documents/MEmilioPaper/HybridApplication/"
colors = {"dabm": "#20A398", "hybrid": "#BA421E", "osecir": "#8200AE", "init": "#E5945A", "sim": "#155489"}
set_fontsize()

# Plot time series
panel = [0.18, 0.18, 0.8, 0.8]
figsize = (5, 3.5)
fig = plt.figure(figsize=figsize)
ax = fig.add_axes(panel)
ax = plot_ts_output(ax, result_dir + "dabm", ["05", "50", "95"], colors["dabm"], ["E", "Ins", "InsC", "Isy", "IsyC", "Isev", "Icri"])
ax = plot_ts_output(ax, result_dir + "osecir", ["05", "50", "95"], colors["osecir"], ["E", "Ins", "InsC", "Isy", "IsyC", "Isev", "Icri"])
ax = plot_ts_output(ax, result_dir + "hybrid", ["05", "50", "95"], colors["hybrid"], ["E", "Ins", "InsC", "Isy", "IsyC", "Isev", "Icri"])
ax.set_ylabel("Infected [#]")
ax.set_xlabel("Time [days]")
fig.savefig(result_dir + "hybrid_application_sim_results.png", dpi=300)

fig = plt.figure(figsize=figsize)
fig.legend(handles=[mpatches.Patch(color=colors["dabm"], label="ABM"), mpatches.Patch(color=colors["osecir"], label="PBM"), mpatches.Patch(color=colors["hybrid"], label="Hybrid")], loc='center', ncol=1)
fig.savefig(result_dir + "hybrid_application_legend.png", dpi=300)

# Runtimes
sim_time = []
init_time = []

for key in colors.keys():
    if key in ["init", "sim"]:
        continue
    df_init = pd.read_csv(result_dir + key + "_init_time.csv")
    df_sim = pd.read_csv(result_dir + key + "_sim_time.csv")
    sim_time.append(df_sim["C1"].mean())
    init_time.append(df_init["C1"].mean())
    
fig = plt.figure(figsize=figsize)
ax = fig.add_axes(panel)
labels = ["ABM", "Hybrid", "PBM"]
x = np.arange(len(labels))
width = 0.35

ax.bar(x - width/2, init_time, width, label='Initialization', color = colors["init"])
ax.bar(x + width/2, sim_time, width, label='Simulation', color = colors["sim"])

ax.set_ylabel('Runtime [s]')
ax.set_yscale('log')
ax.set_xticks(x)
ax.set_xticklabels(labels)
fig.savefig(result_dir + "hybrid_application_runtimes.png", dpi=300)

fig = plt.figure(figsize=figsize)
fig.legend(handles=[mpatches.Patch(color=colors["init"], label="Initialization"), mpatches.Patch(color=colors["sim"], label="Simulation")], loc='center', ncol=1)
fig.savefig(result_dir + "hybrid_application_runtimes_legend.png", dpi=300)
