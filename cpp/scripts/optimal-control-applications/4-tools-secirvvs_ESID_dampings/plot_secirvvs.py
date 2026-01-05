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

plt.axhline(y=1000000, color='blue', linestyle='--', linewidth=1.5, label='Path Constraint')

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
