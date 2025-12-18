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
regions_to_plot = [
"11000",
"6437",
"5766",
# "5770",
# "5774",
# "3251",
# "3352",
# "3356",
# "3357",
# "3361",
# "3401",
# "3402",
# "3403",
# "3405",
# "3451",
# "3452",
# "3453",
# "3455",
"3457",
"3458",
"3461",
# "3462",
# "4011",
# "4012",
# "9464",
# "9475",
# "14511",
"7313"
]
highlighted_region_to_plot = "3457"

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
fig, ax = plt.subplots(figsize=(10, 6))

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
        linewidth=0.5,
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
        linewidth=0.5,
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
        linewidth=0.5,
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

# Add horizontal threshold line
threshold1 = 250
threshold2 = 500
plt.axhline(y=threshold2, color=colors["Red"], linestyle='-', lw=1, alpha = 0.8, label=f"Threshold High = {threshold2}")
plt.axhline(y=threshold1, color=colors["Red"], linestyle='-', lw=1, alpha = 0.5, label=f"Threshold Low = {threshold1}")

ax.set_xlim(0, 120)
plt.xlabel("Time [days]")
plt.ylabel("Infected per 100,000 [#]")
plt.legend()

plt.savefig(os.path.join(save_dir, 'nuts3_dynamic_regional.png'), dpi=dpi)