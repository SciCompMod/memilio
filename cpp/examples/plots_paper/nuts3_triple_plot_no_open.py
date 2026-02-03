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

# df_open = preprocess_data(os.path.join(results_dir, "result_open.h5"))
df_same = preprocess_data(os.path.join(results_dir, "result_same.h5"))
df_lockdown = preprocess_data(os.path.join(results_dir, "result_lockdown.h5"))
df_dynamic = preprocess_data(os.path.join(results_dir, "result_dynamic.h5"))

set_fontsize()        
figsize = (5.5, 3.5)
panel = (0.14, 0.26, 0.82, 0.68)
fig = plt.figure(figsize=figsize)
ax = fig.add_axes(panel)

# Plot first 60 days with dotted line
df_inference[state].plot(ax=ax, lw=2, alpha=0.6, color=color_infected_line, linestyle="--", label='_nolegend_')

# Plot the rest normally
# ax.plot(df_open[state].index + change_day, df_open[state].values, color=colors["Green"], alpha=1.0, label='No NPIs')
ax.plot(df_same[state].index + change_day, df_same[state].values, color=color_infected_line, linewidth=2.5, alpha=1.0, label='Static NPIs')
ax.plot(df_lockdown[state].index + change_day, df_lockdown[state].values, color=colors["Red"], linewidth=2.5, alpha=1.0, label='Strict NPIs')
ax.plot(df_dynamic[state].index + change_day, df_dynamic[state].values, color=colors["Orange"], linewidth=2.5, alpha=1.0, label='Dynamic NPIs')

# Plot X-Axis
start_date = pd.Timestamp("2020-10-01")
date_range = pd.date_range(start=start_date, periods=121, freq="D")
weekly_ticks = np.arange(0, 121, 14)  # Tick every 14 days
ax.set_xticks(weekly_ticks)
ax.set_xticklabels(date_range[weekly_ticks].strftime("%Y-%m-%d"), rotation=45)

ax.yaxis.set_major_locator(MultipleLocator(0.5e6))

# Optional: style
# plt.xlabel("Time [days]")
plt.ylabel("Infected [#]")
plt.grid(False)

ax.set_xlim(0, 120)
ax.set_ylim(0, 2100000)
# plt.legend() # no legend, as we plot it seperatly
plt.savefig(os.path.join(save_dir, 'nuts3_triple_plot_no_open.png'), dpi=dpi)