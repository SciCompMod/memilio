import pandas as pd
import matplotlib.pyplot as plt

# --- Configuration ---
csv_file = 'population_time_series.csv'   # <-- replace with your CSV filename/path

# --- Load data ---
df = pd.read_csv(csv_file)

# --- Aggregate columns into broader categories ---
df['Total_Susceptible'] = df[['SusceptibleNaive',
                              'SusceptiblePartialImmunity',
                              'SusceptibleImprovedImmunity']].sum(axis=1)

df['Total_Exposed'] = df[['ExposedNaive',
                          'ExposedPartialImmunity',
                          'ExposedImprovedImmunity']].sum(axis=1)

df['Total_Infected_Asymptomatic'] = df[['InfectedNoSymptomsNaive',
                                        'InfectedNoSymptomsPartialImmunity',
                                        'InfectedNoSymptomsImprovedImmunity']].sum(axis=1)

df['Total_Infected_Symptomatic'] = df[['InfectedSymptomsNaive',
                                       'InfectedSymptomsPartialImmunity',
                                       'InfectedSymptomsImprovedImmunity']].sum(axis=1)

df['Total_Infected_Asymptomatic_Confirmed'] = df[['InfectedNoSymptomsNaiveConfirmed',
                                        'InfectedNoSymptomsPartialImmunityConfirmed',
                                        'InfectedNoSymptomsImprovedImmunityConfirmed']].sum(axis=1)

df['Total_Infected_Symptomatic_Confirmed'] = df[['InfectedSymptomsNaiveConfirmed',
                                       'InfectedSymptomsPartialImmunityConfirmed',
                                       'InfectedSymptomsImprovedImmunityConfirmed']].sum(axis=1)

df['Total_Infected_Severe'] = df[['InfectedSevereNaive',
                                   'InfectedSeverePartialImmunity',
                                   'InfectedSevereImprovedImmunity']].sum(axis=1)

df['Total_Infected_Critical'] = df[['InfectedCriticalNaive',
                                     'InfectedCriticalPartialImmunity',
                                     'InfectedCriticalImprovedImmunity']].sum(axis=1)

df['Total_Dead'] = df[['DeadNaive',
                       'DeadPartialImmunity',
                       'DeadImprovedImmunity']].sum(axis=1)

# --- Plotting ---
plt.figure(figsize=(12, 8))

for col in ['Total_Susceptible',
            'Total_Exposed',
            'Total_Infected_Asymptomatic',
            'Total_Infected_Symptomatic',
            'Total_Infected_Asymptomatic_Confirmed',
            'Total_Infected_Symptomatic_Confirmed',
            'Total_Infected_Severe',
            'Total_Infected_Critical',
            'Total_Dead']:
    plt.plot(df['Time'], df[col], label=col.replace('_', ' '))

plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Disease Progression Over Time')
plt.legend(loc='upper right', fontsize='small')
plt.tight_layout()

# plt.yscale('log')

# Optional: save to file
plt.savefig('disease_progression.png', dpi=300)

plt.show()
