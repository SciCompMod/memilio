import pandas as pd
import matplotlib.pyplot as plt
import os

# download data from
# https://github.com/robert-koch-institut/Intensivkapazitaeten_und_COVID-19-Intensivbettenbelegung_in_Deutschland/blob/main/Intensivregister_Deutschland_Kapazitaeten.csv


# get path of current file
dir = os.path.dirname(os.path.abspath(__file__))
path_csv = os.path.join(dir, 'Intensivregister_Deutschland_Kapazitaeten.csv')

df = pd.read_csv(path_csv, sep=',')

icu_occ = df['intensivbetten_belegt']
icu_avail = df['intensivbetten_frei']

# ersetze in beiden daten NAN werte durch 0
icu_occ = icu_occ.fillna(0)
icu_avail = icu_avail.fillna(0)

icu_cap = icu_occ + icu_avail

dates = df['datum']

# plot data
plt.figure(figsize=(10, 5))
plt.plot(dates, icu_cap, label='ICU capacity')
# plt.plot(dates, icu_avail, label='ICU available')
plt.plot(dates, icu_occ, label='ICU Occupancy')
plt.xticks(dates[::360], rotation=45)
plt.legend()
plt.title('ICU beds in Germany')
plt.ylabel('Number of ICU beds')
plt.xlabel('Date')
plt.tight_layout()
plt.savefig(dir + "/icu_cap.png")
# plt.show()
