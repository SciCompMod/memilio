import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os


# from the first file 'datensatzbeschreibung_massnahmen' we want to extract the following information:
# on the site DSB BL we have the following information:
#first column is the short name of the measure
# third column is the explanation of the measure
#we want a map from the short name to the explanation

# get the data from the first file
csv_file = pd.read_csv(
    r'C:\Users\korf_sa\Documents\rep\cpp\simulations\datensatzbeschreibung_massnahmen.csv', header=None, skiprows=1)




csv_file = pd.read_csv(
    r'C:\Users\korf_sa\Documents\rep\cpp\simulations\germany_counties_npi_maincat.csv', header=None, skiprows=1)
# csv_file = pd.read_csv(
#     r'C:\Users\korf_sa\Documents\rep\cpp\simulations\germany_counties_npi_subcat.csv', header=None, skiprows=1)

zip_code_brunswick = 3101
start_date = '2021-03-01'
end_date = '2021-05-31'

#mapping from column entry to column name:
dict_leisure = {1: 'ContactPriv', 2: 'Schools', 3: 'Kita', 4: 'IndoorEvents',
                5: 'OutdoorEvents', 6: 'Culture', 7: 'Retail', 8: 'Gastronomy'}

#delete all files with different zip code
csv_file = csv_file[csv_file[2] == zip_code_brunswick]
# also delete all files with different date
csv_file = csv_file[(csv_file[1] >= start_date) & (csv_file[1] <= end_date)]
#delete all the columns that are zero
# csv_file = csv_file.loc[:, (csv_file != 0).any(axis=0)]
#delete all the columns that are not zero but keep header number 0

#plot a heatmap of the resulting 0/1 matrix within the dataframe

sns.heatmap(csv_file.iloc[:, 3:], cmap='coolwarm', cbar=False)
plt.show()

# we want 

x=1