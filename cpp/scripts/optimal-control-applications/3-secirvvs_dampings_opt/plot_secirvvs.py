# import pandas as pd
# import matplotlib.pyplot as plt

# # --- Configuration ---
# csv_file = 'population_time_series.csv'

# # --- Load data ---
# df = pd.read_csv(csv_file)

# # --- Aggregate columns into broader categories ---
# df['Total_Susceptible'] = df[['SusceptibleNaive',
#                               'SusceptiblePartialImmunity',
#                               'SusceptibleImprovedImmunity']].sum(axis=1)

# df['Total_Exposed'] = df[['ExposedNaive',
#                           'ExposedPartialImmunity',
#                           'ExposedImprovedImmunity']].sum(axis=1)

# df['Total_Infected_Asymptomatic'] = df[['InfectedNoSymptomsNaive',
#                                         'InfectedNoSymptomsPartialImmunity',
#                                         'InfectedNoSymptomsImprovedImmunity']].sum(axis=1)

# df['Total_Infected_Symptomatic'] = df[['InfectedSymptomsNaive',
#                                        'InfectedSymptomsPartialImmunity',
#                                        'InfectedSymptomsImprovedImmunity']].sum(axis=1)

# df['Total_Infected_Asymptomatic_Confirmed'] = df[['InfectedNoSymptomsNaiveConfirmed',
#                                         'InfectedNoSymptomsPartialImmunityConfirmed',
#                                         'InfectedNoSymptomsImprovedImmunityConfirmed']].sum(axis=1)

# df['Total_Infected_Symptomatic_Confirmed'] = df[['InfectedSymptomsNaiveConfirmed',
#                                        'InfectedSymptomsPartialImmunityConfirmed',
#                                        'InfectedSymptomsImprovedImmunityConfirmed']].sum(axis=1)

# df['Total_Infected_Severe'] = df[['InfectedSevereNaive',
#                                    'InfectedSeverePartialImmunity',
#                                    'InfectedSevereImprovedImmunity']].sum(axis=1)

# df['Total_Infected_Critical'] = df[['InfectedCriticalNaive',
#                                      'InfectedCriticalPartialImmunity',
#                                      'InfectedCriticalImprovedImmunity']].sum(axis=1)

# df['Total_Dead'] = df[['DeadNaive',
#                        'DeadPartialImmunity',
#                        'DeadImprovedImmunity']].sum(axis=1)

# # --- Plotting ---
# plt.figure(figsize=(12, 8))

# for col in ['Total_Susceptible',
#             'Total_Exposed',
#             'Total_Infected_Asymptomatic',
#             'Total_Infected_Symptomatic',
#             'Total_Infected_Asymptomatic_Confirmed',
#             'Total_Infected_Symptomatic_Confirmed',
#             'Total_Infected_Severe',
#             'Total_Infected_Critical',
#             'Total_Dead']:
#     plt.plot(df['Time'], df[col], label=col.replace('_', ' '))

# plt.xlabel('Time')
# plt.ylabel('Population')
# plt.title('Disease Progression Over Time')
# plt.legend(loc='upper right', fontsize='small')
# plt.tight_layout()

# plt.yscale('log')

# plt.savefig('disease_progression.png', dpi=300)
# plt.show()


import pandas as pd
import matplotlib.pyplot as plt

# --- Configuration ---
csv_file = 'population_time_series.csv'   # <-- replace with your CSV filename/path

# --- Load data ---
df = pd.read_csv(csv_file)

# --- Compute aggregate infection states ---
df['Total_Infected'] = (
    df[['InfectedNoSymptomsNaive',
        'InfectedNoSymptomsPartialImmunity',
        'InfectedNoSymptomsImprovedImmunity',
        'InfectedSymptomsNaive',
        'InfectedSymptomsPartialImmunity',
        'InfectedSymptomsImprovedImmunity']].sum(axis=1)
)

df['Total_Severe'] = (
    df[['InfectedSevereNaive',
        'InfectedSeverePartialImmunity',
        'InfectedSevereImprovedImmunity']].sum(axis=1)
)

df['Total_Critical'] = (
    df[['InfectedCriticalNaive',
        'InfectedCriticalPartialImmunity',
        'InfectedCriticalImprovedImmunity']].sum(axis=1)
)

df['Total_Dead'] = (
    df[['DeadNaive',
        'DeadPartialImmunity',
        'DeadImprovedImmunity']].sum(axis=1)
)

# --- Adjust deaths so they start at 0 ---
df['Total_Dead'] = df['Total_Dead'] - df['Total_Dead'].iloc[0]

# --- Plotting ---
plt.figure(figsize=(10, 6))

plt.axhline(y=125000, color='orange', linestyle='--', linewidth=1.5, label='Path Constraint (Severe)')

plt.plot(df['Time'], df['Total_Infected'], label='Infected')
plt.plot(df['Time'], df['Total_Severe'], label='Severe')
plt.plot(df['Time'], df['Total_Critical'], label='Critical')
plt.plot(df['Time'], df['Total_Dead'], label='Dead')

# plt.axhline(y=100, color='green', linestyle='--', linewidth=1.5, label='ICU-Capacity')

plt.xlabel('Time (days)')
plt.ylabel('Population')
plt.title('Infected, Severe, Critical, and Dead Over Time')
plt.legend()
plt.tight_layout()

# Uncomment for logarithmic scale if desired:
plt.yscale('log')

# Optional: save plot
plt.savefig('disease_progression.png', dpi=300)
plt.show()


# --- Load data ---
df = pd.read_csv("control_parameters.csv")  # replace with your filename if needed

# --- Create plot ---
plt.figure(figsize=(10, 6))

# Columns to plot, HomeOffice first
cols = [
    "SchoolClosure",
    "HomeOffice",
    "PhysicalDistancingSchool",
    "PhysicalDistancingWork",
    "PhysicalDistancingOther"
]

lw = 3

plt.step(df["Time"], df["HomeOffice"], where="post", label="HomeOffice", color="tab:red", linewidth=lw)
plt.step(df["Time"], df["SchoolClosure"], where="post", label="SchoolClosure", color="tab:brown", linewidth=lw)
plt.step(df["Time"], df["PhysicalDistancingSchool"], where="post", label="PhysicalDistancingSchool", color="tab:blue", linewidth=lw)
plt.step(df["Time"], df["PhysicalDistancingWork"], where="post", label="PhysicalDistancingWork", color="tab:orange", linewidth=lw, linestyle='--',
    dashes=(2, 2))
plt.step(df["Time"], df["PhysicalDistancingOther"], where="post", label="PhysicalDistancingOther", color="tab:green", linewidth=lw, linestyle='--',
    dashes=(2, 2))

# --- Styling ---
plt.xlabel("Time (days)")
plt.ylabel("Control intensity")
plt.title("Weekly Control Measures")
# plt.grid(True, linestyle="--", alpha=0.6)
plt.legend(loc="best")
plt.tight_layout()

# --- Show plot ---
plt.show()
plt.savefig("control_measures.png", dpi=300)
