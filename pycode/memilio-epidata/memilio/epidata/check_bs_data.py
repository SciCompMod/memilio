import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


bd = pd.read_csv(r'C:\Users\korf_sa\Documents\rep\pycode\memilio-epidata\memilio\epidata\bs.csv', header=None)
print(bd.dtypes)

bd.rename(
    columns={0: 'idTrafficZone', 1: 'tripID', 2: 'personID', 3: 'tripChain', 4: 'startZone', 5: 'destZone', 6: 'loc_id_start', 7: 'loc_id_end',
             8: 'countyStart', 9: 'countyEnd', 10: 'hhID', 11: 'TripID', 12: 'tripDistance', 13: 'startTime', 14: 'travelTime', 19: 'vehicleChoice', 20:
             'ActivityBefore', 21: 'ActivityAfter', 15: 'loCs', 16: 'laCs', 17: 'loCe', 18: 'laCe', 22: 'age'},
    inplace=True)

dict_leisure = {1: 'work', 2: 'education', 3: 'Shopping', 4: 'free time', 5: 'private matters', 6: 'others', 7: 'home', 0: 'not specified'}

#### Counting people ####
bd_persons = bd.groupby(['personID']).size().reset_index(name='counts')
# Get the frequency of each number in the 'numbers' column
bd_personss = bd_persons['counts'].value_counts().sort_index()

# Plot the frequency of each number
bd_personss.plot(kind='bar')

# # Add labels and title to the plot
# plt.xlabel('Trip count')
# plt.ylabel('Number of persons')
# plt.title('Trip count per person (Braunschweig)' + ' (n = ' + str(len(bd_persons)) + ' for ' + str(len(bd)) + ' trips)')

# # Show the plot
plt.show()

# Count traffic zones where person is residing and longitude and latuitude of this in "EPSG:3035" format
bd_traffic_zones_persons = bd.groupby(['idTrafficZone']).size().reset_index(name='counts').sort_values(by=['counts'], ascending=False, ignore_index=True)
for row in bd_traffic_zones_persons.iterrows():
    bd_traffic_zones_persons.loc[row[0], ['longitude', 'latitude']
                                 ] = bd.loc[bd_traffic_zones_persons.loc[row[0], 'idTrafficZone'] == bd['startZone'], ['loCs', 'laCs']].values[0]
print(bd_traffic_zones_persons)

# Count traffic zones where people are starting their trips
bd_traffic_zones_start = bd.groupby(['startZone']).size().reset_index(name='counts').sort_values(by=['counts'], ascending=False, ignore_index=True)
for row in bd_traffic_zones_start.iterrows():
    bd_traffic_zones_start.loc[row[0], ['longitude', 'latitude']
                               ] = bd.loc[bd_traffic_zones_start.loc[row[0], 'startZone'] == bd['startZone'], ['loCs', 'laCs']].values[0]
print(bd_traffic_zones_start)

# Count traffic zones where people are going and leaving
bd_traffic_zones_end = bd.groupby(['destZone']).size().reset_index(name='counts').sort_values(by=['counts'], ascending=False, ignore_index=True)
for row in bd_traffic_zones_end.iterrows():
    bd_traffic_zones_end.loc[row[0], ['longitude', 'latitude']
                             ] = bd.loc[bd_traffic_zones_end.loc[row[0], 'destZone'] == bd['destZone'], ['loCs', 'laCs']].values[0]
print(bd_traffic_zones_end)


# Time of day
bd_time = bd.groupby(['startTime']).size().reset_index(name='counts')
print(bd_time)
# bd_time.plot(kind='line', x='startTime', y='counts')
# plt.xlabel('Time of day')
# plt.ylabel('Number of trips')
# plt.title('Time of day (Braunschweig)' + ' (n = ' + str(len(bd)) + ' trips)')
# plt.show()


# time of day rolling average
bd_time['rolling_mean'] = bd_time['counts'].rolling(window=7).mean()
# bd_time.plot(kind='line', x='startTime', y='rolling_mean')
# plt.xlabel('Time of day')
# plt.ylabel('Number of trips')
# plt.title('Time of day (Braunschweig)' + ' (n = ' + str(len(bd)) + ' trips), rolling mean witg window size 7')
# plt.show()

# Frequency matrix of trips between traffic zones
matrix_freq = pd.crosstab(bd['startZone'], bd['destZone'])

nonzero_matrix_freq = matrix_freq[matrix_freq > 0]
nonzero_matrix_freq_greater = matrix_freq[matrix_freq > 100]
nnz_matrix_freq = np.count_nonzero(matrix_freq)
subset_matrix_freq_gretaer = nonzero_matrix_freq_greater.iloc[0:30, 0:30]

sns.heatmap(nonzero_matrix_freq_greater,  vmin=0)
sns.heatmap(subset_matrix_freq_gretaer,  vmin=0)


# frequency matrix of trips from which leisure activity to which leisure activity

# check if all activities of bd activity before are in the dictionary
if bd['ActivityBefore'].isin(dict_leisure.keys()).all():
    print('All activities in dictionary')
# check if there is a activity missing that is in the dictionary
if not set(dict_leisure.keys()).issubset(set(bd['ActivityBefore'].unique())):
    print('There is an activity missing that is in the dictionary:')
    for key in (set(dict_leisure.keys()).difference(set(bd['ActivityBefore'].unique()))):
        print(key, dict_leisure[key])


matrix_frqu_leisure = pd.crosstab(
    bd['ActivityBefore'],
    bd['ActivityAfter'],
    normalize='index')
sns.heatmap(matrix_frqu_leisure,  vmin=0, xticklabels=dict_leisure.values(),
            yticklabels=dict_leisure.values())


# quick check if a data point in the matrix is correct
check = bd.loc[(bd['startZone'] == 31010011) & (bd['destZone'] == 31010011), 'startZone'].count()


# # plotting the duration of trips
# fig, axs = plt.subplots(2)
# axs[0].hist(bd['travelTime'], bins=100, range=(0, max(bd['travelTime'])*1.2))
# axs[0].set_title('Trip duration (Braunschweig)' + ' (n = ' + str(len(bd)) + ' trips)')
# axs[0].set_xlabel('Trip duration in seconds')
# axs[0].set_ylabel('Number of trips')
# axs[1].hist(bd['travelTime'], bins=100, range=(0, max(bd['travelTime'])*1.2))
# axs[1].set_yscale('log')
# axs[1].set_title('Trip duration (Braunschweig)' + ' (n = ' + str(len(bd)) + ' trips) log scale')
# axs[1].set_xlabel('Trip duration in seconds')
# axs[1].set_ylabel('Number of trips')

# # same thing with the distance of the trips
# fig2, axs2 = plt.subplots(2)
# axs2[0].hist(bd['tripDistance'], bins=100, range=(0, max(bd['tripDistance'])*1.2))
# axs2[0].set_title('Trip distance (Braunschweig)' + ' (n = ' + str(len(bd)) + ' trips)')
# axs2[0].set_xlabel('Trip distance in meters')
# axs2[0].set_ylabel('Number of trips')
# axs2[1].hist(bd['tripDistance'], bins=100, range=(0, max(bd['tripDistance'])*1.2))
# axs2[1].set_yscale('log')
# axs2[1].set_title('Trip distance (Braunschweig)' + ' (n = ' + str(len(bd)) + ' trips) log scale')
# axs2[1].set_xlabel('Trip distance in meters')
# axs2[1].set_ylabel('Number of trips')

# analyze age distribution
# bd_age = bd.groupby(['age']).size().reset_index(name='counts')
bd_age = bd.groupby(['age']).size().reset_index(name='counts').rolling(window=7, min_periods=1, step=5).mean()
bd_age.plot(kind='bar', x='age', y='counts')


# relation between age and trip duration
bd_age_duration = bd.groupby(['age']).mean().reset_index()
bd_age_duration.plot(kind='bar', x='age', y='travelTime')

# also do this for trip distance
bd_age_distance = bd.groupby(['age']).mean().reset_index()
bd_age_distance.plot(kind='bar', x='age', y='tripDistance')

# also do a scatter plot of trip duration and trip distance with a regression line
bd_age_distance.plot(kind='scatter', x='age', y='tripDistance')
plt.plot(np.unique(bd_age_distance['age']), np.poly1d(np.polyfit(bd_age_distance['age'],
         bd_age_distance['tripDistance'], 1))(np.unique(bd_age_distance['age'])))
bd_age_duration.plot(kind='scatter', x='age', y='travelTime')
plt.plot(np.unique(bd_age_duration['age']), np.poly1d(np.polyfit(bd_age_duration['age'], bd_age_duration['travelTime'], 1))(np.unique(bd_age_duration['age'])))
plt.show()


x = 42
