# import pandas as pd
# import numpy as np
# import glob
# import matplotlib.pyplot as plt

# files = sorted(glob.glob("population_time_series_run_*.csv"))
# dfs = [pd.read_csv(f) for f in files]

# time = dfs[0]["Time"]

# def compute_aggregates(df):
#     out = pd.DataFrame()
#     out["Time"] = df["Time"]
#     out["Total_Infected"] = df[
#         [
#             "InfectedNoSymptomsNaive",
#             "InfectedNoSymptomsPartialImmunity",
#             "InfectedNoSymptomsImprovedImmunity",
#             "InfectedSymptomsNaive",
#             "InfectedSymptomsPartialImmunity",
#             "InfectedSymptomsImprovedImmunity",
#         ]
#     ].sum(axis=1)
#     out["Total_Severe"] = df[
#         [
#             "InfectedSevereNaive",
#             "InfectedSeverePartialImmunity",
#             "InfectedSevereImprovedImmunity",
#         ]
#     ].sum(axis=1)
#     out["Total_Critical"] = df[
#         [
#             "InfectedCriticalNaive",
#             "InfectedCriticalPartialImmunity",
#             "InfectedCriticalImprovedImmunity",
#         ]
#     ].sum(axis=1)
#     out["Total_Dead"] = (
#         df[
#             ["DeadNaive", "DeadPartialImmunity", "DeadImprovedImmunity"]
#         ].sum(axis=1)
#         - df[["DeadNaive", "DeadPartialImmunity", "DeadImprovedImmunity"]]
#         .sum(axis=1)
#         .iloc[0]
#     )
#     return out

# agg = [compute_aggregates(df) for df in dfs]

# stacked = {var: np.vstack([a[var].values for a in agg]) 
#            for var in ["Total_Infected", "Total_Severe", "Total_Critical", "Total_Dead"]}

# percentiles = {
#     name: {
#         "p5": np.percentile(stacked[name], 5, axis=0),
#         "p50": np.percentile(stacked[name], 50, axis=0),
#         "p95": np.percentile(stacked[name], 95, axis=0),
#     }
#     for name in stacked
# }

# plt.figure(figsize=(12, 7))
# plt.axhline(1000000, linestyle="--", label="Path Constraint")

# for name, color in zip(
#     ["Total_Infected", "Total_Severe", "Total_Critical", "Total_Dead"],
#     ["C0", "C1", "C2", "C3"]
# ):
#     p = percentiles[name]
#     plt.plot(time, p["p50"], label=f"{name} (median)", color=color)
#     plt.fill_between(time, p["p5"], p["p95"], alpha=0.25, color=color)

# plt.yscale("log")
# plt.xlabel("Time (days)")
# plt.ylabel("Population")
# plt.title("Percentile Ranges Over All Simulation Runs")
# plt.legend()
# plt.tight_layout()
# plt.savefig("percentile_plot.png", dpi=300)
# plt.show()
import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

files = sorted(glob.glob("population_time_series_run_*.csv"))
dfs = [pd.read_csv(f) for f in files]

print("Number of runs:", len(files))

time = dfs[0]["Time"]

def compute_aggregates(df):
    out = pd.DataFrame()
    out["Time"] = df["Time"]
    out["Infected"] = df[
        [
            "InfectedNoSymptomsNaive",
            "InfectedNoSymptomsPartialImmunity",
            "InfectedNoSymptomsImprovedImmunity",
            "InfectedSymptomsNaive",
            "InfectedSymptomsPartialImmunity",
            "InfectedSymptomsImprovedImmunity",
        ]
    ].sum(axis=1)
    out["Severe"] = df[
        [
            "InfectedSevereNaive",
            "InfectedSeverePartialImmunity",
            "InfectedSevereImprovedImmunity",
        ]
    ].sum(axis=1)
    out["Critical"] = df[
        [
            "InfectedCriticalNaive",
            "InfectedCriticalPartialImmunity",
            "InfectedCriticalImprovedImmunity",
        ]
    ].sum(axis=1)
    out["Dead"] = (
        df[["DeadNaive", "DeadPartialImmunity", "DeadImprovedImmunity"]].sum(axis=1)
        - df[["DeadNaive", "DeadPartialImmunity", "DeadImprovedImmunity"]].sum(axis=1).iloc[0]
    )
    return out

agg = [compute_aggregates(df) for df in dfs]

plt.figure(figsize=(12, 7))
plt.axhline(100_000, linestyle="--", label="Path Constraint")

# plt.axhline(70, linestyle="--", label="Path Constraint", color="red")

variables = ["Infected", "Severe", "Critical", "Dead"]
colors = ["C0", "C1", "C2", "C3"]

# Plot all runs, same linewidth
for name, color in zip(variables, colors):
    for a in agg:
        plt.plot(
            time,
            a[name].values,
            color=color,
            alpha=0.25,
            linewidth=0.8
        )

# Legend handles without plotting extra lines
legend_elements = [
    Line2D([0], [0], color=color, lw=2, label=name)
    for name, color in zip(variables, colors)
]

plt.yscale("log")
plt.xlabel("Time (days)")
plt.ylabel("Population")
plt.title("Optimization Under Uncertainty")
plt.legend(handles=legend_elements)
plt.tight_layout()
plt.savefig("all_runs_plot.png", dpi=300)
plt.show()
