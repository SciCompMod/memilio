from plotting_settings import plotting_dir, set_fontsize, colors, dpi
from memilio.simulation.osecir import InfectionState

import h5py
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_rgba
import matplotlib.pyplot as plt
import pandas as pd
import os

# params

color_infected_line = colors["Teal"]

damping_day = 7
change_day = 60
state = "InfectedSymptoms"
regions_to_plot = [
# "11000",
# "6437",
# "5766",
# "4011",
# "14511",
# "11000", # Berlin
"5334", # Aachen
"8125", # Heilbronn
# "1001", # Flensburg
# "12054", # Potsdam
# "3457" # Leer
]
highlighted_region_to_plot = "1001"

results_dir = os.path.join(plotting_dir, "../simulation_paper/results")
save_dir = os.path.join(plotting_dir, "plots")

# index: (region, time)
# columns: infection states

def preprocess_data_multiregion(filepath):
    h5file = h5py.File(filepath, 'r')

    region_names = list(h5file.keys())
    all_dfs = []

    infection_state_names = list(InfectionState.__members__.keys())

    for region in region_names:
        grp = h5file[region]
        time = grp["Time"][:]
        data = grp["Total"][:]

        df = pd.DataFrame(
            data,
            index=time,
            columns=infection_state_names
        )
        df["region"] = region
        df["time"] = time

        all_dfs.append(df)

    df_all = pd.concat(all_dfs)
    df_all.set_index(["region", "time"], inplace=True)

    return df_all

df_inference = preprocess_data_multiregion(os.path.join(plotting_dir, "../../../data/results_run0.h5"))
df_same = preprocess_data_multiregion(os.path.join(results_dir, "result_same.h5"))
df_dynamic = preprocess_data_multiregion(os.path.join(results_dir, "result_dynamic.h5"))


population = df_inference.sum(axis=1).round()

df_inference["Incidence"] = (
    df_inference["InfectedSymptoms"] / population * 100_000
)
df_same["Incidence"] = (
    df_same["InfectedSymptoms"] / population * 100_000
)
df_dynamic["Incidence"] = (
    df_dynamic["InfectedSymptoms"] / population * 100_000
)

set_fontsize()
figsize = (9.5, 4.5)
fig, ax = plt.subplots(figsize=figsize)

# Add horizontal threshold line
threshold1 = 100
threshold2 = 500
plt.axhline(y=threshold2, color=colors["Blue"], linestyle='-', lw=1, alpha = 0.8, label=f"Threshold High = {threshold2}")
plt.axhline(y=threshold1, color=colors["Blue"], linestyle='-', lw=1, alpha = 0.5, label=f"Threshold Low = {threshold1}")

# Make fade
fade = np.linspace(0, 1, 256).reshape(256, 1)

base_fade_color = to_rgba(colors["Blue"])

custom_red_cmap = LinearSegmentedColormap.from_list(
    "custom_fade",
    [
        (*base_fade_color[:3], 0.0),   # fully transparent
        (*base_fade_color[:3], 1.0),   # fully opaque
    ],
)

ax.imshow(
    fade,
    extent=[0, 120, threshold2 + 3, threshold2 + 500],
    origin="upper",
    aspect="auto",
    cmap=custom_red_cmap,
    alpha=0.5,
    zorder=0
)
ax.imshow(
    fade,
    extent=[0, 120, threshold1 + 3, threshold1 + 300],
    origin="upper",
    aspect="auto",
    cmap=custom_red_cmap,
    alpha=0.3,
    zorder=0
)

for region in regions_to_plot:
    # --- First phase: df_inference (days 0–60)
    inf_series = df_inference.loc[
        (region, slice(0, 60)), "Incidence"
    ]

    plt.plot(
        inf_series.index.get_level_values("time"),
        inf_series.values,
        linestyle="--",
        color=colors["Grey"],
        alpha = 1,
        linewidth=1,
        label='_nolegend_',
    )

    # --- Third phase: df_dynamic (days 60–120)
    dyn_series = df_dynamic.loc[(region, slice(None)), "Incidence"]

    # shift time axis
    shifted_time = (
        dyn_series.index.get_level_values("time") + 60
    )

    plt.plot(
        shifted_time,
        dyn_series.values,
        linestyle="-",
        color=colors["Dark grey"],
        alpha = 1,
        linewidth=1,
        label='_nolegend_'   
    )

    # --- Second phase: df_same (days 60–120)
    same_series = df_same.loc[
        (region, slice(0, 60)), "Incidence"
    ]

    # shift time axis
    shifted_time = (
        same_series.index.get_level_values("time") + 60
    )

    plt.plot(
        shifted_time,
        same_series.values,
        linestyle="-",      # different style
        color=colors["Grey"],
        alpha = 1,
        linewidth=1,
        label='_nolegend_',
    )

# Plot highlighted region
inf_series = df_inference.loc[
    (highlighted_region_to_plot, slice(0, 60)), "Incidence"
]

plt.plot(
    inf_series.index.get_level_values("time"),
    inf_series.values,
    linestyle="--",
    color=colors["Teal"],
    linewidth=2.5,
    label='_nolegend_',
)

# --- Third phase: df_dynamic (days 60–120)
dyn_series = df_dynamic.loc[(highlighted_region_to_plot, slice(None)), "Incidence"]

# shift time axis
shifted_time = (
    dyn_series.index.get_level_values("time") + 60
)

plt.plot(
    shifted_time,
    dyn_series.values,
    linestyle="-",
    color=colors["Orange"],
    linewidth=2.5,
    label=f"{highlighted_region_to_plot} with Dynamic NPIs"
)

# --- Second phase: df_same (days 60–120)
same_series = df_same.loc[
    (highlighted_region_to_plot, slice(0, 60)), "Incidence"
]

# shift time axis
shifted_time = (
    same_series.index.get_level_values("time") + 60
)

plt.plot(
    shifted_time,
    same_series.values,
    linestyle="-",      # different style
    color=colors["Teal"],
    linewidth=2.5,
    label=f"{highlighted_region_to_plot}",
)

# Plot X-Axis
start_date = pd.Timestamp("2020-10-01")
date_range = pd.date_range(start=start_date, periods=121, freq="D")
weekly_ticks = np.arange(0, 121, 14)  # Tick every 14 days
ax.set_xticks(weekly_ticks)
ax.set_xticklabels(date_range[weekly_ticks].strftime("%Y-%m-%d"), rotation=45)

ymin, ymax = ax.get_ylim()
ax.set_ylim(0, ymax)
ax.set_xlim(0, 120)
# plt.xlabel("Time [days]")
plt.ylabel("Infected [# per 100k]")
# plt.legend() # no legend, as we plot it seperatly

plt.tight_layout()
plt.savefig(os.path.join(save_dir, 'nuts3_dynamic_regional.png'), dpi=dpi)