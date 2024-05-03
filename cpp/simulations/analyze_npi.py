import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import openpyxl

# # from the first file 'datensatzbeschreibung_massnahmen' we want to extract the following information:
# # on the site DSB BL we have the following information:
# #first column is the short name of the measure
# # third column is the explanation of the measure
# #we want a map from the short name to the explanation

# # get the data from the first file
# csv_file = pd.read_csv(
#     r'C:\Users\korf_sa\Documents\rep\cpp\simulations\datensatzbeschreibung_massnahmen.csv', header=None, skiprows=1)




# csv_file = pd.read_csv(
#     r'C:\Users\korf_sa\Documents\rep\cpp\simulations\germany_counties_npi_maincat.csv', header=None, skiprows=1)
# # csv_file = pd.read_csv(
# #     r'C:\Users\korf_sa\Documents\rep\cpp\simulations\germany_counties_npi_subcat.csv', header=None, skiprows=1)



# #mapping from column entry to column name:
# dict_leisure = {1: 'ContactPriv', 2: 'Schools', 3: 'Kita', 4: 'IndoorEvents',
#                 5: 'OutdoorEvents', 6: 'Culture', 7: 'Retail', 8: 'Gastronomy'}

# #delete all files with different zip code
# csv_file = csv_file[csv_file[2] == zip_code_brunswick]
# # also delete all files with different date
# csv_file = csv_file[(csv_file[1] >= start_date) & (csv_file[1] <= end_date)]
# #delete all the columns that are zero
# # csv_file = csv_file.loc[:, (csv_file != 0).any(axis=0)]
# #delete all the columns that are not zero but keep header number 0

# #plot a heatmap of the resulting 0/1 matrix within the dataframe

# sns.heatmap(csv_file.iloc[:, 3:], cmap='coolwarm', cbar=False)
# plt.show()

# # we want 

# x=1



#import the first excel file with the measures and ther abbreviations 
df_abb = pd.read_excel('/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/cpp/simulations/datensatzbeschreibung_massnahmen.xlsx', sheet_name='DSB BL')

#read in the matrix which tells us which measure is active on which day
df_measure_matrix = pd.read_csv(r'/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/cpp/simulations/germany_counties_npi_subcat.csv', header=None, skiprows=1)

zip_code_brunswick = 3101
start_date = '2021-03-01'
end_date = '2021-05-31'

#now we want to shorten the matrix to the zip code of Brunswick and the time period of interest
#delete all files with different zip code
df_measure_matrix = df_measure_matrix[df_measure_matrix[2] == zip_code_brunswick]
# also delete all files with different date
df_measure_matrix = df_measure_matrix[(df_measure_matrix[1] >= start_date) & (df_measure_matrix[1] <= end_date)]

# now we want to do the following:
# 1. for every measure we want to count how many days it was active
# 2. We want to sort the measures by the number of days they were active
# 3. We want to plot the measures in a bar plot with the number of days they were active
# 4. We want a txt file which has for each measure the full explanation and the days it was active

# begin with 1.
# get the number of days a measure was active
measure_by_days = df_measure_matrix.iloc[:, 3:]
# count the number of days a measure was active
measure_by_days = measure_by_days.sum()
# sort the measures by the number of days they were active
measure_by_days = measure_by_days.sort_values(ascending=False)
#delete all the columns that are zero
measure_by_days = measure_by_days[measure_by_days != 0]

# now we need create a txt file with the measures and the number of days they were active and the full explanation
# create a txt file
f = open('measures_active_days.txt', 'w')
#get first row of df_abb


# for every measure
for measure in measure_by_days.axes[0]:
    # get the full explanation
    full_explanation = df_abb['Beschreibung'][measure-3]
    # write the full explanation and the number of days the measure was active to the txt file
    f.write(full_explanation + ' ' + str(measure_by_days[measure]) + '\n')

x=1
