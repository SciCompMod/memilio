import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os


def get_trip_chain_activity_after(person_id):
    bd_persons_trip_chain_activity_after = bd.loc[bd['personID'] == person_id, ['tripChain', 'ActivityAfter']]
    bd_persons_trip_chain_activity_after = bd_persons_trip_chain_activity_after.sort_values(by=['tripChain'], ascending=True, ignore_index=True)
    bd_persons_trip_chain_activity_after['ActivityAfter'] = bd_persons_trip_chain_activity_after['ActivityAfter'].map(dict_leisure)
    # and vehicle choice
    print("\n")
    bd_persons_trip_chain_vehicle_choice = bd.loc[bd['personID'] == person_id, ['tripChain', 'vehicleChoice']]
    bd_persons_trip_chain_vehicle_choice = bd_persons_trip_chain_vehicle_choice.sort_values(by=['tripChain'], ascending=True, ignore_index=True)
    bd_persons_trip_chain_vehicle_choice['vehicleChoice'] = bd_persons_trip_chain_vehicle_choice['vehicleChoice'].map(dict_vehicle)
    return bd_persons_trip_chain_activity_after, bd_persons_trip_chain_vehicle_choice


bd = pd.read_csv(r'/Users/saschakorf/Documents/Arbeit/memilio.nosynch/memilio/pycode/memilio-epidata/memilio/epidata/bs.csv', header=None)
figs_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figs_bs_data')

bd.rename(
    columns={0: 'idTrafficZone', 1: 'tripID', 2: 'personID', 3: 'tripChain', 4: 'startZone', 5: 'destZone', 6: 'loc_id_start', 7: 'loc_id_end',
             8: 'countyStart', 9: 'countyEnd', 10: 'hhID', 11: 'TripID', 12: 'tripDistance', 13: 'startTime', 14: 'travelTime', 19: 'vehicleChoice', 20:
             'ActivityBefore', 21: 'ActivityAfter', 15: 'loCs', 16: 'laCs', 17: 'loCe', 18: 'laCe', 22: 'age'},
    inplace=True)

dict_leisure = {1: 'work', 2: 'education', 3: 'Shopping', 4: 'free time', 5: 'private matters', 6: 'others', 7: 'home', 0: 'not specified'}
dict_vehicle = {1: 'bicyle', 2: 'car_driver', 3: 'car_codriver', 4: 'public transport', 5: 'walk'}

#### Counting people ####
bd_persons = bd.groupby(['personID']).size().reset_index(name='counts')
# Get the frequency of each number in the 'numbers' column
bd_personss = bd_persons['counts'].value_counts().sort_index()

# probably the traffic zones that are in braunschweig
bd_tz_bs = bd.groupby(['idTrafficZone']).size().reset_index(name='counts').sort_values(by=['counts'], ascending=False, ignore_index=True)

# count who does how many trips and sort them
bd_persons = bd_persons.sort_values(by=['counts'], ascending=False, ignore_index=True)
# Plot the frequency of each number
bd_personss.plot(kind='bar')
plt.xlabel('Number of trips a person does')
plt.ylabel('Number of persons who do this number of trips')
plt.title('Number of trips per person tripcount (Braunschweig)')

# plot the amount of trips of the highest 200 persons
bd_persons[0:4000].plot(kind='bar', x='personID', y='counts')
plt.xlabel('Person ID')
plt.ylabel('Number of trips')
plt.title('Number of trips per person (Braunschweig)' + ' (n = ' + str(len(bd_persons)) + ' persons)')
plt.savefig(os.path.join(figs_path, 'tripcount.png'), dpi=300)

# get the id with the persons with the highest number of trips
[id_person_max_trips, id_person_max_1_trips] = [bd_persons['personID'][0], bd_persons['personID'][1]]
print(get_trip_chain_activity_after(id_person_max_trips))

# select rows where the start county is in id traffic zone
bd_persons_inside_bs = bd.loc[bd['startZone'].isin(bd_tz_bs['idTrafficZone']) & bd['destZone'].isin(bd_tz_bs['idTrafficZone'])]


# see which age takes which mode of transport
bd_persons_age_vehicle_choice = bd.loc[:, ['personID', 'age', 'vehicleChoice']]
bd_persons_age_vehicle_choice = bd_persons_age_vehicle_choice.sort_values(by=['age'], ascending=True, ignore_index=True)
bd_persons_age_vehicle_choice['vehicleChoice'] = bd_persons_age_vehicle_choice['vehicleChoice'].map(dict_vehicle)
# accumulate the number of trips per age and vehicle choice
bd_persons_age_vehicle_choice = bd_persons_age_vehicle_choice.groupby(
    ['age', 'vehicleChoice']).size().reset_index(
    name='counts').sort_values(
    by=['age'],
    ascending=True, ignore_index=True)
# assign each age to a age cohort
bd_persons_age_vehicle_choice['ageCohort'] = pd.cut(bd_persons_age_vehicle_choice['age'], bins=[0, 18, 25, 35, 45, 55, 65, 75, 85, 95, 105], labels=[
                                                    '0-18', '19-25', '26-35', '36-45', '46-55', '56-65', '66-75', '76-85', '86-95', '96-105'])
# plot a cake chart for each age cohort
bd_persons_age_vehicle_choice_cake = bd_persons_age_vehicle_choice.groupby(['ageCohort', 'vehicleChoice']).sum().reset_index()
bd_persons_age_vehicle_choice_cake = bd_persons_age_vehicle_choice_cake.pivot(index='ageCohort', columns='vehicleChoice', values='counts')
bd_persons_age_vehicle_choice_cake.plot.pie(subplots=True, figsize=(20, 10), title='Vehicle choice per age cohort (Braunschweig)')
plt.savefig(os.path.join(figs_path, 'vehicle_choice.png'), dpi=300)

# switch age cohort and vehicle choice
# bd_persons_age_vehicle_choice_cake = bd_persons_age_vehicle_choice_cake.transpose()
# bd_persons_age_vehicle_choice_cake.plot.pie(subplots=True, figsize=(20, 10), title='Age cohort per vehicle choice (Braunschweig)')
# plt.show()


# plot which traffic zone is the most popular for each age cohort
bd_persons_age_traffic_zone = bd.loc[:, ['personID', 'age', 'idTrafficZone']]
bd_persons_age_traffic_zone = bd_persons_age_traffic_zone.sort_values(by=['age'], ascending=True, ignore_index=True)
# accumulate the number of trips per age and traffic zone
bd_persons_age_traffic_zone = bd_persons_age_traffic_zone.groupby(
    ['age', 'idTrafficZone']).size().reset_index(
    name='counts').sort_values(
    by=['age'],
    ascending=True, ignore_index=True)
# assign each age to a age cohort
bd_persons_age_traffic_zone['ageCohort'] = pd.cut(bd_persons_age_traffic_zone['age'], bins=[0, 18, 25, 35, 45, 55, 65, 75, 85, 95, 105], labels=[
    '0-18', '19-25', '26-35', '36-45', '46-55', '56-65', '66-75', '76-85', '86-95', '96-105'])
# plot the number of trips per age and traffic zone
bd_persons_age_traffic_zone.plot(x='ageCohort', y='counts', kind='bar', figsize=(20, 10), title='Number of trips per age and traffic zone (Braunschweig)')


# Add labels and title to the plot
plt.xlabel('Age cohort')
plt.ylabel('Number of persons')
plt.title('Trip count per person (Braunschweig)' + ' (n = ' + str(len(bd_persons)) + ' for ' + str(len(bd)) + ' trips)')

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
# plot longitude and latitude of traffic zones
bd_traffic_zones_start.plot(
    kind='scatter', x='longitude', y='latitude', s=bd_traffic_zones_start['counts'] / 100, figsize=(20, 10),
    title='Traffic zones where people are starting their trips (Braunschweig)')

# Count traffic zones where people are going and leaving
bd_traffic_zones_end = bd.groupby(['destZone']).size().reset_index(name='counts').sort_values(by=['counts'], ascending=False, ignore_index=True)
for row in bd_traffic_zones_end.iterrows():
    bd_traffic_zones_end.loc[row[0], ['longitude', 'latitude']
                             ] = bd.loc[bd_traffic_zones_end.loc[row[0], 'destZone'] == bd['destZone'], ['loCs', 'laCs']].values[0]
print(bd_traffic_zones_end)


# Time of day
bd_time = bd.groupby(['startTime']).size().reset_index(name='counts')
print(bd_time)
bd_time.plot(kind='line', x='startTime', y='counts')
plt.xlabel('Time of day')
plt.ylabel('Number of trips')
plt.title('Time of day (Braunschweig)' + ' (n = ' + str(len(bd)) + ' trips)')
plt.show()

# analyze diff between each time points

# time of day rolling average
bd_time['rolling_mean'] = bd_time['counts'].rolling(window=7).mean()
bd_time.plot(kind='line', x='startTime', y='rolling_mean')
plt.xlabel('Time of day')
plt.ylabel('Number of trips')
plt.title('Time of day (Braunschweig)' + ' (n = ' + str(len(bd)) + ' trips), rolling mean witg window size 7')
plt.show()

# Frequency matrix of trips between traffic zones
matrix_freq = pd.crosstab(bd['startZone'], bd['destZone'])
nnz_matrix_freq = np.count_nonzero(matrix_freq)
subset_matrix_freq = matrix_freq.iloc[0:30, 0:30]
nnz_sub_matrix_freq = np.count_nonzero(subset_matrix_freq)
sns.heatmap(subset_matrix_freq, cmap="Blues",  vmin=0)
plt.xlabel('Destination traffic zone')
plt.ylabel('Start traffic zone')
plt.title('Number of non-zero entries in frequency matrix: ' + str(nnz_sub_matrix_freq) + ' out of ' + str(subset_matrix_freq.size) + ' entries (' +
          str(round(nnz_sub_matrix_freq / subset_matrix_freq.size * 100, 2)) + '%)')
plt.savefig(os.path.join(figs_path, 'heatmap.png'), dpi=300)

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


# plotting the duration of trips
fig, axs = plt.subplots(2)
axs[0].hist(bd['travelTime'], bins=100, range=(0, max(bd['travelTime'])*1.2))
axs[0].set_title('Trip duration (Braunschweig)' + ' (n = ' + str(len(bd)) + ' trips)')
axs[0].set_xlabel('Trip duration in seconds')
axs[0].set_ylabel('Number of trips')
axs[1].hist(bd['travelTime'], bins=100, range=(0, max(bd['travelTime'])*1.2))
axs[1].set_yscale('log')
axs[1].set_title('Trip duration (Braunschweig)' + ' (n = ' + str(len(bd)) + ' trips) log scale')
axs[1].set_xlabel('Trip duration in seconds')
axs[1].set_ylabel('Number of trips')

# same thing with the distance of the trips
fig2, axs2 = plt.subplots(2)
axs2[0].hist(bd['tripDistance'], bins=100, range=(0, max(bd['tripDistance'])*1.2))
axs2[0].set_title('Trip distance (Braunschweig)' + ' (n = ' + str(len(bd)) + ' trips)')
axs2[0].set_xlabel('Trip distance in meters')
axs2[0].set_ylabel('Number of trips')
axs2[1].hist(bd['tripDistance'], bins=100, range=(0, max(bd['tripDistance'])*1.2))
axs2[1].set_yscale('log')
axs2[1].set_title('Trip distance (Braunschweig)' + ' (n = ' + str(len(bd)) + ' trips) log scale')
axs2[1].set_xlabel('Trip distance in meters')
axs2[1].set_ylabel('Number of trips')

# analyze age distribution
bd_age = bd.groupby(['age']).size().reset_index(name='counts')
bd_age = bd.groupby(['age']).size().reset_index(name='counts').rolling(window=7, min_periods=1, step=5).mean()
bd_age.plot(kind='bar', x='age', y='counts')
# TODO: redo this in age cohorts
# analysze age distribution with age cohorts
pd.cut(bd_persons_age_traffic_zone['age'], bins=[0, 18, 25, 35, 45, 55, 65, 75, 85, 95, 105], labels=[
    '0-18', '19-25', '26-35', '36-45', '46-55', '56-65', '66-75', '76-85', '86-95', '96-105'])
bd_age_cohorts = bd.groupby(['age']).size().reset_index(name='counts')
bd_age_cohorts.plot(kind='bar', x='age_cohort', y='counts')


# relation between age and trip duration in minutes and hours
bd_age_duration = bd.groupby(['age']).mean().reset_index()
bd_age_duration['travelTime'] = bd_age_duration['travelTime'] / 60
bd_age_duration.plot(kind='bar', x='age', y='travelTime')

# also do this for trip distance #TODO: further investigate what tripdistance means
bd_age_distance = bd.groupby(['age']).mean().reset_index()
bd_age_distance.plot(kind='bar', x='age', y='tripDistance')

# also do a scatter plot of trip duration and trip distance with a regression line
bd_age_distance.plot(kind='scatter', x='age', y='tripDistance')
plt.plot(np.unique(bd_age_distance['age']), np.poly1d(np.polyfit(bd_age_distance['age'],
         bd_age_distance['tripDistance'], 1))(np.unique(bd_age_distance['age'])))
bd_age_duration.plot(kind='scatter', x='age', y='travelTime')
plt.plot(np.unique(bd_age_duration['age']), np.poly1d(np.polyfit(bd_age_duration['age'], bd_age_duration['travelTime'], 1))(np.unique(bd_age_duration['age'])))
plt.show()

# TODO: Vergleich der tripldistance mit vehicle choice
#

x = 42
