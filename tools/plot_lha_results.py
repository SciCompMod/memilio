import os
import h5py
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# configuration
# -----------------------------
DATASET = "lha_data_2026-03-12_full_pop"
BASE_DIR = f"./data/results/{DATASET}"
SCENARIOS = ["original", "altered_vaccinations"]
PERCENTILES = ["p25", "p50", "p75"]

OUTPUT_DIR = f"plots/{DATASET}"
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

for comp_idx, comp_name in COMPARTMENTS.items():

    plt.figure(figsize=(10, 6))

    for scenario in SCENARIOS:

        if not all(p in results[scenario] for p in PERCENTILES):
            continue

        time = results[scenario]["p50"]["time"]

        p25 = results[scenario]["p25"]["data"][:, comp_idx]
        p50 = results[scenario]["p50"]["data"][:, comp_idx]
        p75 = results[scenario]["p75"]["data"][:, comp_idx]

        color = colors[scenario]

        # percentile band
        plt.fill_between(
            time,
            p25,
            p75,
            alpha=0.25,
            color=color,
            label=f"{scenario} p25–p75"
        )

        # median
        plt.plot(
            time,
            p50,
            color=color,
            linewidth=2,
            label=f"{scenario} median"
        )

    plt.title(comp_name)
    plt.xlabel("Time")
    plt.ylabel("Population")
    plt.grid(True)
    plt.legend()

    plt.tight_layout()

    outfile = os.path.join(
        OUTPUT_DIR,
        f"{comp_idx:02d}_{comp_name}.png"
    )

    plt.savefig(outfile, dpi=300)
    plt.close()

print("All plots saved to:", OUTPUT_DIR)
