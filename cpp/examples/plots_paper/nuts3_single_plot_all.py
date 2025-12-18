from plotting_settings import plotting_dir, set_fontsize, colors, dpi
from memilio.simulation.osecir import InfectionState

import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# params

color_infected_line = colors["Teal"]
color_damping_background = colors["Red"]

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

def plot_curve(df_inference, df, color=color_infected_line):

    set_fontsize()        
    figsize = (5, 3.5)
    panel = (0.21, 0.18, 0.75, 0.78)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(panel)

    # Plot first 60 days with dotted line
    df_inference[state].plot(ax=ax, color=color, linestyle="--", alpha=1.0, label='_nolegend_')

    # Plot the rest normally
    ax.plot(df[state].index + change_day, df[state].values, color=color, alpha=1.0, label='_nolegend_')

    # Damping background
    plt.axvspan(damping_day, df.index.max(), color=color_damping_background, alpha=0.1, label="Damping active", linewidth=0)

    # Optional: style
    plt.xlabel("Time [days]")
    plt.ylabel("Infected [#]")
    plt.grid(False)

    ax.set_xlim(0, 120)
    # ax.set_ylim(bottom=None, top=df[state].max()*1.05)
    plt.legend()

df_inference = preprocess_data(os.path.join(plotting_dir, "../../../data/results_run0.h5"))

# ================================================= PLOT =================================================

df = preprocess_data(os.path.join(results_dir, "result_open.h5"))
plot_curve(df_inference, df)

plt.savefig(os.path.join(save_dir, 'germany_nuts3_open.png'), dpi=dpi)

# ================================================= PLOT =================================================

df = preprocess_data(os.path.join(results_dir, "result_same.h5"))
plot_curve(df_inference, df)

plt.axvspan(change_day, df.index.max() + change_day, color=color_damping_background, alpha=0.1, label="Damping active", linewidth=0)

plt.savefig(os.path.join(save_dir, 'germany_nuts3_same.png'), dpi=dpi)

# ================================================= PLOT =================================================

df = preprocess_data(os.path.join(results_dir, "result_lockdown.h5"))
plot_curve(df_inference, df)

plt.axvspan(change_day, df.index.max() + change_day, color=color_damping_background, alpha=0.2, label="_nolegend_", linewidth=0)

plt.savefig(os.path.join(save_dir, 'germany_nuts3_lockdown.png'), dpi=dpi)