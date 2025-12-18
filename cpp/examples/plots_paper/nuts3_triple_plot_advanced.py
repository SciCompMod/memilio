from plotting_settings import plotting_dir, set_fontsize, colors, dpi
from memilio.simulation.osecir import InfectionState

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd
import os

# params

color_infected_line = colors["Teal"]

damping_day = 7
change_day = 60
state = "InfectedSymptoms"
scale = 100_000

results_dir = os.path.join(plotting_dir, "../simulation_paper/results")
save_dir = os.path.join(plotting_dir, "plots")

def preprocess_data(filepath):
    h5file = h5py.File(filepath, 'r')
    # time_array = np.array(h5file['1001']["Time"])

    region_names = list(h5file.keys())

    # Get shape info from first region
    first_region = h5file[region_names[0]]
    time = first_region["Time"][:]                   # 1D (n_time,)
    data_shape = first_region["Total"].shape         # (n_time, n_state)

    # Initialize accumulator
    total_sum = np.zeros(data_shape, dtype=np.float64)

    # --- Vectorized accumulation across all regions ---
    for region in region_names:
        region_data = h5file[region]["Total"][:]  # or ["group_1"]
        total_sum += region_data
    
    # --- Convert to DataFrame --- 
    infection_state_names = []
    for name, value in InfectionState.__members__.items():
        infection_state_names.append(f"{name}")

    df = pd.DataFrame(total_sum, index=time, columns=infection_state_names)
    df.index.name = "time"

    return df


df_inference = preprocess_data(os.path.join(plotting_dir, "../../../data/results_run0.h5"))

# ================================================= PLOT =================================================

df_open = preprocess_data(os.path.join(results_dir, "result_open.h5"))
df_same = preprocess_data(os.path.join(results_dir, "result_same.h5"))
df_lockdown = preprocess_data(os.path.join(results_dir, "result_lockdown.h5"))
df_dynamic = preprocess_data(os.path.join(results_dir, "result_dynamic.h5"))

set_fontsize()        
figsize = (5, 3.5)
fig = plt.figure(figsize=figsize)
left = 0.17
width = 0.79
# bottom axis
ax_low = fig.add_axes((left, 0.18, width, 0.44))
# top axis
ax_high = fig.add_axes((left, 0.63, width, 0.30), sharex=ax_low)


def plot_all(ax):
    df_inference[state].divide(scale).plot(
        ax=ax, lw=2, alpha=0.6,
        color=color_infected_line, linestyle="--", label='_nolegend_'
    )

    ax.plot(df_open[state].index + change_day, df_open[state].values / scale,
            color=colors["Green"], label='No NPIs')
    ax.plot(df_same[state].index + change_day, df_same[state].values / scale,
            color=color_infected_line, label='Static NPIs')
    ax.plot(df_lockdown[state].index + change_day, df_lockdown[state].values / scale,
            color=colors["Red"], label='Strict NPIs')
    ax.plot(df_dynamic[state].index + change_day, df_dynamic[state].values / scale,
            color=colors["Orange"], label='Dynamic NPIs')

plot_all(ax_low)
plot_all(ax_high)

ax_low.set_ylim(0, 2e6 / scale)       # low-level dynamics
ax_high.set_ylim(2e6 / scale, 2e7 / scale)    # exploding peak

# hide spines between axes
ax_high.spines['bottom'].set_visible(False)
ax_low.spines['top'].set_visible(False)

# tick_step = 1_000_000
# ax_low.yaxis.set_major_locator(MultipleLocator(tick_step))
# ax_high.yaxis.set_major_locator(MultipleLocator(tick_step))

ax_high.tick_params(labelbottom=False)
ax_high.tick_params(
    axis='x',
    which='both',
    bottom=False,
    top=False,
    labelbottom=False
)
ax_low.set_yscale("linear")
ax_high.set_yscale('symlog', linthresh=1)
yticks = [50, 100, 200]
ax_high.set_yticks(yticks)
ax_high.set_yticklabels([str(t) for t in yticks])

ax_high.set_xlabel("")   
ax_low.set_xlabel("Time [days]")
ax_low.set_ylabel("Infected [#]")
ax_low.yaxis.set_label_coords(-0.1, 1)  

ax_high.grid(False)
ax_low.grid(False)

# Add diagonal break marks (important for readability)
d = 0.015
kwargs = dict(transform=ax_high.transAxes, color='k', clip_on=False)
ax_high.plot((-d, +d), (-d, +d), **kwargs)
ax_high.plot((1 - d, 1 + d), (-d, +d), **kwargs)

scale_d = 0.7
kwargs.update(transform=ax_low.transAxes, color='k', clip_on=False)
ax_low.plot((-d, +d), (1 - scale_d*d, 1 + scale_d*d), **kwargs)
ax_low.plot((1 - d, 1 + d), (1 - scale_d*d, 1 + scale_d*d), **kwargs)

ax_high.legend()
ax_low.set_xlim(0, 120)

plt.legend()
plt.savefig(os.path.join(save_dir, 'nuts3_triple_plot_advanced.png'), dpi=dpi)