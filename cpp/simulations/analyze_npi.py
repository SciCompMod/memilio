import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import openpyxl


#define a function which takes in a list of dates in the form yyyy-mm-dd and summarizes them in a range of dates where the dates are consecutive, e.g. [2021-01-01, 2021-01-02, 2021-01-03, 2021-01-05] -> [2021-01-01 to 2021-01-03, 2021-01-05]
def summarize_dates(dates):
    #sort the list of dates
    dates.sort()
    #initialize the list of ranges
    ranges = []
    #initialize the start of the range
    start = dates[0]
    #initialize the end of the range
    end = dates[0]
    #iterate through the list of dates
    for i in range(1, len(dates)):
        #if the date is the next day of the previous date
        if pd.to_datetime(dates[i]) == pd.to_datetime(dates[i-1]) + pd.DateOffset(days=1):
            #set the end of the range to the current date
            end = dates[i]
        #if the date is not the next day of the previous date
        else:
            #append the range to the list of ranges
            ranges.append(str(start) + ' to ' + str(end))
            #set the start of the range to the current date
            start = dates[i]
            #set the end of the range to the current date
            end = dates[i]
    #append the last range to the list of ranges
    ranges.append(str(start) + ' to ' + str(end))
    return ranges



#import the first excel file with the measures and ther abbreviations 
df_abb = pd.read_excel(
    '/Users/david/Documents/HZI/memilio/data/npi_data/datensatzbeschreibung_massnahmen.xlsx', sheet_name='DSB BL')

#read in the matrix which tells us which measure is active on which day
df_measure_matrix = pd.read_csv(
    r'/Users/david/Documents/HZI/memilio/data/npi_data/germany_counties_npi_subcat.csv', header=None, skiprows=1)

zip_code_brunswick = 3101
start_date = '2021-03-01'
end_date = '2021-05-31'


one_day_after = pd.to_datetime(start_date) + pd.DateOffset(days=1)
#now we want to shorten the matrix to the zip code of Brunswick and the time period of interest
#delete all files with different zip code
df_measure_matrix = df_measure_matrix[df_measure_matrix['ID_County'] == zip_code_brunswick]
#also delete all files with different date
df_measure_matrix = df_measure_matrix[(df_measure_matrix['Date'] >= start_date) & (df_measure_matrix['Date'] <= end_date)]
#delete all the columns that just contain zeros
df_measure_matrix = df_measure_matrix.loc[:, (df_measure_matrix != 0).any(axis=0)]


#rename the first column from the abbreviation to the full explanation
# create a map from the short name to the explanation
dict_measure = {}
for i in range(len(df_abb)):
    dict_measure[df_abb['Variablenname'][i]] = df_abb['Beschreibung'][i]
df_measure_matrix = df_measure_matrix.rename(columns=dict_measure)


# now we need create a txt file with the measures and the number of days they were active and the full explanation
f = open('measures_active_days.txt', 'w')
# for every measure
for measure in df_measure_matrix.axes[1][3:]:
    # count the days it was active
    days = df_measure_matrix[measure].sum()
    # get the dates it was active on in one list
    dates = df_measure_matrix[df_measure_matrix[measure] == 1]['Date'].to_list()
    # summarize the dates in ranges
    date_ranges = summarize_dates(dates)
    # write the measure and the number of days it was active in the txt file also write down each date the measure was active
    f.write(measure + ': ' + str(days) + ' days ' + 'Dates: ' + str(date_ranges) + '\n\n')
    

f.close()

df_measure_matrix = df_measure_matrix.iloc[:, measure_by_days.axes[0]]
sns.heatmap(df_measure_matrix.iloc[:, 3:], cmap='coolwarm', cbar=False)
plt.show()

x = 1
