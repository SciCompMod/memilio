from plotting_settings import plotting_dir, set_fontsize, colors, dpi
from memilio.simulation.osecir import InfectionState

import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

infection_state_names = []
for name, value in InfectionState.__members__.items():
    infection_state_names.append(f"{name}")


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

    df = pd.DataFrame(total_sum, index=time, columns=infection_state_names)
    df.index.name = "time"

    return df

df_inference = preprocess_data(os.path.join(plotting_dir, "../../../data/results_run0.h5"))

damping_day = 7
change_day = 60
state = "InfectedSymptoms"

results_dir = os.path.join(plotting_dir, "../simulation_paper/results")
save_dir = os.path.join(plotting_dir, "plots")

# First Plot
df = preprocess_data(os.path.join(results_dir, "result_open.h5"))

fig, ax = plt.subplots(figsize=(5, 5))

# Plot first 60 days with dotted line
df_inference[state].plot(ax=ax, lw=2, alpha=0.6, color='tab:blue', linestyle="--", label='_nolegend_')

# Plot the rest normally
ax.plot(df[state].index + change_day, df[state].values, lw=2, color='tab:blue', alpha=1.0, label='_nolegend_')

#df["InfectedNoSymptoms"].plot(ax=ax, lw=2, label="InfectedNoSymptoms")

plt.xlabel("Time [days]")
plt.ylabel("Infected [#]")
plt.grid(False)

# plt.axvline(damping_day, color="red", linestyle="--", lw=2)
# plt.annotate("Damping introduced",
#              xy=(damping_day, df[state].max()*0.8),
#              xytext=(damping_day+1, df[state].max()*0.85),
#              arrowprops=dict(arrowstyle="->", color="red"),
#              color="red")

# plt.axvline(change_day, color="gray", linestyle="--", lw=2)

plt.axvspan(damping_day, df_inference.index.max(), color="red", alpha=0.1, label="Damping active")
# plt.axvspan(change_day, df.index.max(), color="red", alpha=0.2)
ax.set_xlim(0, 120)
ax.ticklabel_format(style='plain', axis='y')

plt.legend(loc='upper left')
plt.savefig(os.path.join(save_dir, 'germany_nuts3_open.png'), dpi=dpi)

# Second Plot

df = preprocess_data(os.path.join(results_dir, "result_same.h5"))

fig, ax = plt.subplots(figsize=(5, 5))

# Plot first 60 days with dotted line
df_inference[state].plot(ax=ax, lw=2, alpha=0.6, color='tab:blue', linestyle="--", label='_nolegend_')

# Plot the rest normally
ax.plot(df[state].index + change_day, df[state].values, lw=2, color='tab:blue', alpha=1.0, label='_nolegend_')

# Optional: style
plt.xlabel("Time [days]")
plt.ylabel("Infected [#]")
plt.grid(False)
plt.yscale('log')

plt.axvspan(damping_day, df.index.max() + change_day, color="red", alpha=0.1, label="Damping active")

ax.set_xlim(0, 120)
# ax.set_ylim(bottom=None, top=df[state].max()*1.05)
plt.legend(loc='upper left')
plt.savefig(os.path.join(save_dir, 'germany_nuts3_same.png'), dpi=dpi)

# Third Plot

df = preprocess_data(os.path.join(results_dir, "result_lockdown.h5"))

fig, ax = plt.subplots(figsize=(5, 5))

# Plot first 60 days with dotted line
df_inference[state].plot(ax=ax, lw=2, alpha=0.6, color='tab:blue', linestyle="--", label='_nolegend_')

# Plot the rest normally
ax.plot(df[state].index + change_day, df[state].values, color='tab:blue', alpha=1.0, label='_nolegend_')

# Optional: style
plt.xlabel("Time [days]")
plt.ylabel("Infected [#]")
plt.grid(False)

plt.axvspan(damping_day, df.index.max(), color="red", alpha=0.1, label="Damping active")
plt.axvspan(change_day, df.index.max() + change_day, color="red", alpha=0.2, label="_nolegend_")

ax.set_xlim(0, 120)
# ax.set_ylim(bottom=None, top=df[state].max()*1.05)
plt.legend(loc='upper left')
plt.savefig(os.path.join(save_dir, 'germany_nuts3_lockdown.png'), dpi=dpi)

# Fourth Plot

df = preprocess_data(os.path.join(results_dir, "result_dynamic.h5"))

threshold1 = 50 * 830   # your threshold value
threshold2 = 250 * 830   # your threshold value

fig, ax = plt.subplots(figsize=(5, 5))

# Plot first 60 days with dotted line
df_inference[state].plot(ax=ax, lw=2, alpha=0.6, color='tab:blue', linestyle="--", label='_nolegend_')

# Plot the rest normally
ax.plot(df[state].index + change_day, df[state].values, lw=2, color='tab:blue', alpha=1.0, label='_nolegend_')

# Add horizontal threshold line
ax.axhline(y=threshold1, color='red', linestyle='--', lw=2, label=f"Threshold Low = {threshold1}")
ax.axhline(y=threshold2, color='red', linestyle='--', lw=2, label=f"Threshold High = {threshold2}")

# ax.fill_between(
#     df.index,
#     threshold,                 # start shading at threshold
#     df[state].max()*1.05,      # end slightly above max for clarity
#     color='red',
#     alpha=0.1
# )

plt.axvspan(damping_day, change_day, color="red", alpha=0.1, label="Damping active")

# Optional: style
ax.set_xlabel("Time [days]")
plt.ylabel("Infected [#]")
plt.grid(False)

ax.legend(loc='upper left')
ax.set_xlim(0, 120)
ax.set_ylim(bottom=None, top=df[state].max()*1.05)
plt.savefig(os.path.join(save_dir, 'germany_nuts3_dynamic.png'), dpi=dpi)