from plotting_settings import plotting_dir, set_fontsize, colors, dpi
from memilio.simulation.osecir import InfectionState

import h5py
import numpy as np
import matplotlib.pyplot as plt
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

df_open = preprocess_data(os.path.join(results_dir, "result_open.h5"))
df_same = preprocess_data(os.path.join(results_dir, "result_same.h5"))
df_lockdown = preprocess_data(os.path.join(results_dir, "result_lockdown.h5"))
df_dynamic = preprocess_data(os.path.join(results_dir, "result_dynamic.h5"))

set_fontsize()        
figsize = (5, 3.5)
panel = (0.17, 0.18, 0.79, 0.76)
fig = plt.figure(figsize=figsize)
ax = fig.add_axes(panel)

# Plot first 60 days with dotted line
df_inference[state].plot(ax=ax, lw=2, alpha=0.6, color=color_infected_line, linestyle="--", label='_nolegend_')

# Plot the rest normally
ax.plot(df_open[state].index + change_day, df_open[state].values, color=colors["Green"], alpha=1.0, label='No NPIs')
ax.plot(df_same[state].index + change_day, df_same[state].values, color=color_infected_line, alpha=1.0, label='Static NPIs')
ax.plot(df_lockdown[state].index + change_day, df_lockdown[state].values, color=colors["Red"], alpha=1.0, label='Strict NPIs')
ax.plot(df_dynamic[state].index + change_day, df_dynamic[state].values, color=colors["Orange"], alpha=1.0, label='Dynamic NPIs')

# Optional: style
plt.xlabel("Time [days]")
plt.ylabel("Infected [#]")
plt.grid(False)

ax.set_xlim(0, 120)
# ax.set_ylim(bottom=None, top=df[state].max()*1.05)
plt.legend()
plt.savefig(os.path.join(save_dir, 'nuts3_triple_plot_simple.png'), dpi=dpi)