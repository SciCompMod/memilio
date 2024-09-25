import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.dates as mdates
import os

opacity = 0.15
lineWidth = 3
fontsize = 28
tickssize = 20  # Size of tick labels
week = 1

cwd = os.getcwd()
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# Load the CSV data
path_csv = "SARS-CoV-2-PCR-Testungen_in_Deutschland.csv"
df = pd.read_csv(path_csv, sep=",", decimal=",")

path_csv_bayern = 'PATH_TO//export_2024-09-12.csv'
df_bayern = pd.read_csv(path_csv_bayern, sep=";", decimal=",")

# eintrag für die erste woche von 2023
df_bayern_first_week_23 = df_bayern[df_bayern['Jahr'] == '2023'].iloc[0]


# lösche alle spalten mit 'Jahr' != 2022
df_bayern = df_bayern[df_bayern['Jahr'] == '2022']
# spalte Meldewoche als int
df_bayern['Meldewoche'] = df_bayern['Meldewoche'].astype(int)
df_bayern = df_bayern[df_bayern['Meldewoche'] >= week].reset_index(drop=True)
# Spalte 'Positivrate' ist in Prozent angegeben. Als float umwandeln
df_bayern['Positivrate'] = df_bayern['Positivrate'].str.replace(
    '%', '').str.replace(',', '.').astype(float) / 100

df_bayern_positivrate = df_bayern['Positivrate'].astype(float)
df_bayern_tests = df_bayern['Gesamtanzahl'].astype(float)


# Filter entries that start with 2022 or the first week of 2023
df['date'] = df['date'].str.replace('2023-W1$', '2023-W01', regex=True)
df = df[df['date'].str.startswith(('2022', '2023-W01'))].reset_index(drop=True)

# Convert 'weeks' to datetime format representing the end of each week
df['date'] = pd.to_datetime(df['date'] + '-0', format='%Y-W%U-%w')
weeks = df['date']
df_positive_ratio = df['tests_positive_ratio'].astype(float)
df_tests_total = df['tests_total'].astype(float)

# ab 39 woche
# df_positive_ratio = df_positive_ratio[39:]
# df_tests_total = df_tests_total[39:]
weeks = weeks[week:]

# Plot into one figure
fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.set_xlabel('Date', fontsize=fontsize)
ax1.set_ylabel('Positive Rate', fontsize=fontsize, color=colors[0])
ax1.plot(weeks, df_bayern_positivrate, color=colors[0], linewidth=lineWidth)
ax1.tick_params(axis='y', labelcolor=colors[0], labelsize=tickssize)
ax1.grid(True)

# Format x-axis to show "Month Year"
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%B %Y'))
ax1.xaxis.set_major_locator(mdates.MonthLocator())
plt.setp(ax1.get_xticklabels(), rotation=45, ha="right", fontsize=tickssize)

# Plot the secondary axis
ax2 = ax1.twinx()
ax2.set_ylabel('Total Tests', fontsize=fontsize, color=colors[1])
ax2.plot(weeks, df_bayern_tests, color=colors[1], linewidth=lineWidth)
ax2.tick_params(axis='y', labelcolor=colors[1], labelsize=tickssize)

fig.tight_layout()
plt.savefig(
    "PATH_TO//tests_germany.png")
plt.show()
