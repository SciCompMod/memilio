import os
import h5py
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# configuration
# -----------------------------
DATASET = "lha_data_2026-03-24"
SUBFOLDER = "2026-03-24"
BASE_DIR = f"./data/results/{SUBFOLDER}/{DATASET}"
SCENARIOS = ["original", "altered_vaccinations"]
PERCENTILES = ["p25", "p50", "p75"]

OUTPUT_DIR = f"plots/{SUBFOLDER}/{DATASET}/infcri"
os.makedirs(OUTPUT_DIR, exist_ok=True)

AGE_GROUPS = ["Group1", "Group2", "Group3", "Group4", "Group5", "Group6"]

COMPARTMENTS = {
    0: 'SusceptibleNaive',
    1: 'SusceptiblePartialImmunity',
    2: 'ExposedNaive',
    3: 'ExposedPartialImmunity',
    4: 'ExposedImprovedImmunity',
    5: 'InfectedNoSymptomsNaive',
    6: 'InfectedNoSymptomsPartialImmunity',
    7: 'InfectedNoSymptomsImprovedImmunity',
    8: 'InfectedNoSymptomsNaiveConfirmed',
    9: 'InfectedNoSymptomsPartialImmunityConfirmed',
    10: 'InfectedNoSymptomsImprovedImmunityConfirmed',
    11: 'InfectedSymptomsNaive',
    12: 'InfectedSymptomsPartialImmunity',
    13: 'InfectedSymptomsImprovedImmunity',
    14: 'InfectedSymptomsNaiveConfirmed',
    15: 'InfectedSymptomsPartialImmunityConfirmed',
    16: 'InfectedSymptomsImprovedImmunityConfirmed',
    17: 'InfectedSevereNaive',
    18: 'InfectedSeverePartialImmunity',
    19: 'InfectedSevereImprovedImmunity',
    20: 'InfectedCriticalNaive',
    21: 'InfectedCriticalPartialImmunity',
    22: 'InfectedCriticalImprovedImmunity',
    23: 'SusceptibleImprovedImmunity',
    24: 'DeadNaive',
    25: 'DeadPartialImmunity',
    26: 'DeadImprovedImmunity'
}


# -----------------------------
# load and aggregate
# -----------------------------

def load_aggregated_results(filepath):

    with h5py.File(filepath, "r") as f:

        sim_key = list(f.keys())[0]
        g = f[sim_key]

        time = g["Time"][:]

        aggregated = None

        for group in AGE_GROUPS:

            data = g[group][:]

            if aggregated is None:
                aggregated = data.copy()
            else:
                aggregated += data

    return time, aggregated


# -----------------------------
# load all datasets
# -----------------------------

results = {}

for scenario in SCENARIOS:

    results[scenario] = {}

    for p in PERCENTILES:

        path = os.path.join(
            BASE_DIR,
            scenario,
            "non_aggregated_results",
            p,
            "Results.h5"
        )

        if os.path.exists(path):

            time, data = load_aggregated_results(path)

            results[scenario][p] = {
                "time": time,
                "data": data
            }


# -----------------------------
# plotting
# -----------------------------

colors = {
    "original": "tab:blue",
    "altered_vaccinations": "tab:red"
}

infcri_comps = ["InfectedCriticalNaive",
                "InfectedCriticalPartialImmunity", "InfectedCriticalImprovedImmunity"]

infcri_labels = ["ICU (Naive)", "ICU (Partial immunity)",
                 "ICU (Improved immunity)"]

# find indices for the three InfectedCritical compartments
comp_indices = []
for comp in infcri_comps:
    found = None
    for idx, name in COMPARTMENTS.items():
        if name == comp:
            found = idx
            break
    if found is not None:
        comp_indices.append((found, comp))
# if any missing, proceed with those found

labels = ["Fully local data (median)", "Fully local data (p25–p75)", "Reduced information data (median)",
          "Reduced information data (p25–p75)"]

if comp_indices:
    n = len(comp_indices)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 6), sharex=True)

    # ensure axes is iterable
    if n == 1:
        axes = [axes]

    for i, (comp_idx, comp_name) in enumerate(comp_indices):
        ax = axes[i]

        for scenario in SCENARIOS:

            if not all(p in results[scenario] for p in PERCENTILES):
                continue

            time = results[scenario]["p50"]["time"]

            p25 = results[scenario]["p25"]["data"][:, comp_idx]
            p50 = results[scenario]["p50"]["data"][:, comp_idx]
            p75 = results[scenario]["p75"]["data"][:, comp_idx]

            color = colors.get(scenario, None)

            if scenario == "original":
                label = "Fully local data"
            else:
                label = "Reduced information data"

            ax.plot(
                time,
                p50,
                color=color,
                linewidth=2,
                # label=f"{label} median"
            )

            ax.fill_between(
                time,
                p25,
                p75,
                alpha=0.25,
                color=color,
                # label=f"{label} p25–p75"
            )

        ax.set_title(infcri_labels[i] if i < len(
            infcri_labels) else comp_name, fontsize=20, fontweight="bold", pad=10)
        # ax.set_xlabel("Time")
        if i == 0:
            ax.set_ylabel("Population", labelpad=10, fontsize=16)
        if i == 1:
            ax.set_xlabel("Time", fontsize=16, labelpad=10)
        ax.grid(True)

    # fig.supxlabel("        Time", fontsize=14)
    fig.legend(labels, ncol=2, fontsize=18, loc="lower center",
               bbox_to_anchor=(0.52, -0.15), bbox_transform=fig.transFigure)  #

    # plt.tight_layout()
    plt.tight_layout()  # rect=(0., 0., 0.9, 1.)

    outfile = os.path.join(
        OUTPUT_DIR,
        f"InfectedCritical_all.png"
    )

    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close()

print("All plots saved to:", OUTPUT_DIR)
