import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os


####### minimal sanity check on data #######
bd = pd.read_csv(r'data/mobility/bs.csv', header=None, skiprows=1)

# setup dictionary for the leisure activities, and vehicle choice and column names
bd.rename(
    columns={0: 'idTrafficZone', 1: 'tripID', 2: 'personID', 3: 'tripChain', 4: 'startZone', 5: 'destZone', 6: 'loc_id_start', 7: 'loc_id_end',
             8: 'countyStart', 9: 'countyEnd', 10: 'hhID', 11: 'tripChainID', 12: 'tripDistance', 13: 'startTime', 14: 'travelTime', 19: 'vehicleChoice', 20:
             'ActivityBefore', 21: 'ActivityAfter', 15: 'loCs', 16: 'laCs', 17: 'loCe', 18: 'laCe', 22: 'age'},
    inplace=True)

dict_leisure = {1: 'work', 2: 'education', 3: 'Shopping', 4: 'free time',
                5: 'private matters', 6: 'others', 7: 'home', 0: 'not specified'}
dict_vehicle = {1: 'bicyle', 2: 'car_driver',
                3: 'car_codriver', 4: 'public transport', 5: 'walk'}

# # check if people do the same trip more than once
# trips = bd[['personID', 'loc_id_start', 'loc_id_end', 'startTime']]
# duplicate_trips = trips.duplicated()
# if (duplicate_trips.any()):
#     print('Error: The Person does the same trip more than once. Number of duplicate trips: ' + str(duplicate_trips[duplicate_trips==True].size))
#     if (~bd[['tripID']].duplicated().any()):
#         print('There are no duplicate TripIDs. There are multiple TripIDs for the same Trip. \n')
# activities_after_duplicate_trips  = bd[['personID', 'loc_id_start', 'loc_id_end', 'startTime', 'ActivityAfter']].loc[trips.duplicated(keep=False)]
# if (~activities_after_duplicate_trips.duplicated(keep=False).all()):
#     print('Error: Multiple activities after the same trip. \n')

# # check if persons have more than one home
# person_homes = bd[['personID', 'loc_id_end', 'ActivityAfter']].loc[bd['ActivityAfter']==7]
# person_homes = person_homes.drop_duplicates().groupby(['personID']).size().reset_index(name='counts')
# person_homes.drop(person_homes.loc[person_homes['counts']<2].index, inplace=True)
# if(person_homes.size > 0):
#     print('Error: There are people that have more than one home. \n')

# # check if the proportion of single-person-households is too high
# households = bd[['personID', 'hhID']].drop_duplicates().groupby(['hhID']).size().reset_index(name='counts').sort_values(by=['counts'], ascending=False, ignore_index=True)

# if households.drop(households.loc[households['counts']>1].index).size / bd[['personID']].drop_duplicates().size >= 0.4:
#     print("Error: The proportion of single-person-households is too high: " + 
#           str(households.drop(households.loc[households['counts']>1].index).size / bd[['personID']].drop_duplicates().size))

# # check if there are invalid entries
# if not bd['idTrafficZone'].ge(30000000).all():
#     print('Error: There is an entry in "tripID" that is not assignable. \n')
#     # number of entries that are not assignable
#     print('Number of entries that are not assignable: ' + str(bd.loc[~bd['tripID'].ge(30000000)].size) + '. \n')
# if not bd['personID'].ge(100000000).all():
#     print('Error: There is an entry in "personID" that is not assignable. \n')
#     # number of entries that are not assignable
#     print('Number of entries that are not assignable: ' + str(bd.loc[~bd['personID'].ge(100000000)].size) + '. \n')
# if not bd['tripChain'].between(1, 100).all():
#     print('Error: There is an entry in "tripChain" that is not assignable. \n')
#     # number of entries that are not assignable
#     print('Number of entries that are not assignable: ' + str(bd.loc[~bd['tripChain'].between(1, 100)].size) + '. \n')
#     # max assigned value
#     print('Max assigned value: ' + str(bd['tripChain'].max()) + '. \n')
# if not bd['countyStart'].ge(30000000).all():
#     print('Error: There is an entry in "countyStart" that is not assignable. \n')
#     # which ones are not assignable
#     print('Number of entries that are not assignable: ' + str(bd.loc[~bd['countyStart'].ge(30000000)].size) + '. \n')
# if not bd['countyEnd'].ge(30000000).all():
#     print('Error: There is an entry in "countyEnd" that is not assignable. \n')
#     # which ones are not assignable
#     print('Number of entries that are not assignable: ' + str(bd.loc[~bd['countyEnd'].ge(30000000)].size) + '. \n')
# if not bd['hhID'].ge(100000000).all():
#     print('Error: There is an entry in "hhID" that is not assignable. \n')
# if not bd['tripChainID'].between(-1, -1).all():
#     print('Error: There is an entry in "tripChainID" that is not assignable. \n')
# if not bd['vehicleChoice'].between(1, 5).all():
#     print('Error: There is an entry in "vehicleChoice" that is not assignable. \n')
#     #print number of entries that are not assignable
#     print('Number of entries that are not assignable: ' + str(bd.loc[~bd['vehicleChoice'].between(1, 5)].size) + '. \n')
# if not bd['ActivityAfter'].between(0, 7).all():
#     print('Error: There is an entry in "ActivityAfter" that is not assignable. \n')
#     #print number of entries that are not assignable
#     print('Number of entries that are not assignable: ' + str(bd.loc[~bd['ActivityAfter'].between(0, 7)].size) + '. \n')

# # check if there are empty cells
# for header in bd.columns:
#     if (bd[header].isna().any()):
#         print('Error: ' + str(len(bd[bd[header].isna()])) + ' empty entries in column' + str(header) + '. \n')

# # check age groups with schools
# number_of_people = bd[['personID']].drop_duplicates().size
# print(str(number_of_people) + ' are people. \n')

# number_of_trips = bd[['tripID']].drop_duplicates().size
# print(str(number_of_trips) + ' trips. \n')

# students = bd[['personID', 'loc_id_end', 'age']].loc[bd['ActivityAfter']==2]
# print('Minimal age of people going to school: ' + str(students['age'].min()) + '. Maximal age of people going to school: ' + str(students['age'].max()) + '.\n')
# print(str(students.loc[students['age'] > 20].size) + ' persons of ' + str(students.size) + ' people going to school are in a higher age group than 10. \n')

# children = bd[['personID','ActivityAfter']].loc[bd['age']<=20]
# number_of_children = children['personID'].drop_duplicates().size
# number_of_children_school = children[children['ActivityAfter'] == 2]['personID'].drop_duplicates().size
# print(str(number_of_children - number_of_children_school) + ' of ' + str(number_of_children) + ' children do not go to school. \n')





###############################################





####### visual check on values in data frame #######

first_trip = bd[['loc_id_start']][bd['tripChain']==1]
last_trip = bd[['personID', 'tripChain', 'loc_id_end']].sort_values(by=['personID', 'tripChain'])
last_trip.drop_duplicates(subset=['personID'], inplace=True, keep='last')
first_trip.compare(last_trip)



# persons = bd[['personID']].drop_duplicates()
# for person in persons.iterrows():
#     trip = bd[['personID', 'tripChain', 'loc_id_end']][bd[['personID']]==person[1]].max()
# last_trip = bd[['personID', 'tripChain']].groupby('personID', group_keys=False).apply(lambda x: max(x))
# last_trip = bd[['loc_id_end']][~bd['personID'].duplicated(keep='last')]



def get_trip_chain_activity_after(person_id):
    bd_persons_trip_chain_activity_after = bd.loc[bd['personID'] == person_id, [
        'tripChain', 'ActivityAfter']]
    bd_persons_trip_chain_activity_after = bd_persons_trip_chain_activity_after.sort_values(
        by=['tripChain'], ascending=True, ignore_index=True)
    bd_persons_trip_chain_activity_after['ActivityAfter'] = bd_persons_trip_chain_activity_after['ActivityAfter'].map(
        dict_leisure)
    # and vehicle choice
    bd_persons_trip_chain_vehicle_choice = bd.loc[bd['personID'] == person_id, [
        'tripChain', 'vehicleChoice']]
    bd_persons_trip_chain_vehicle_choice = bd_persons_trip_chain_vehicle_choice.sort_values(
        by=['tripChain'], ascending=True, ignore_index=True)
    bd_persons_trip_chain_vehicle_choice['vehicleChoice'] = bd_persons_trip_chain_vehicle_choice['vehicleChoice'].map(
        dict_vehicle)
    return bd_persons_trip_chain_activity_after, bd_persons_trip_chain_vehicle_choice


# read in the data
if not os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figs_bs_data')):
    os.makedirs(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), 'figs_bs_data'))
figs_path = os.path.join(os.path.dirname(
    os.path.abspath(__file__)), 'figs_bs_data')

# probably the traffic zones that are in braunschweig because in idTrafficzones the people of bs are located
bd_tz_bs = bd.groupby(['idTrafficZone']).size().reset_index(
    name='counts').sort_values(by=['counts'], ascending=False, ignore_index=True)
# Count traffic zones where person is residing and longitude and latuitude of this in "EPSG:3035" format
bd_traffic_zones_persons = bd.groupby(['idTrafficZone']).size().reset_index(
    name='counts').sort_values(by=['counts'], ascending=False, ignore_index=True)
for row in bd_traffic_zones_persons.iterrows():
    bd_traffic_zones_persons.loc[row[0], ['longitude', 'latitude']
                                 ] = bd.loc[bd_traffic_zones_persons.loc[row[0], 'idTrafficZone'] == bd['startZone'], ['loCs', 'laCs']].values[0]
print(bd_traffic_zones_persons)

#### Counting people ####
bd_persons = bd.groupby(['personID']).size().reset_index(
    name='counts').sort_values(by=['counts'], ascending=False, ignore_index=True)
# Get the frequency of each number in the 'numbers' column
bd_personss = bd_persons['counts'].value_counts().sort_index()
# Plot the frequency of each number
bd_personss.plot(kind='bar')
plt.xlabel('Number of trips per day')
plt.ylabel('Number of persons')
plt.title('Number of trips per person tripcount (Braunschweig)')
plt.savefig(os.path.join(figs_path, 'number_of_trips_per_person.png'), dpi=300)

bd_persons_ages = bd[['personID', 'age']].drop_duplicates()
bd_persons_ages['AgeCohort'] = pd.cut(bd_persons_ages['age'], bins=[-1, 18, 25, 35, 45, 55, 65, 75, 85, 105], labels=[
    '0-18', '19-25', '26-35', '36-45', '46-55', '56-65', '66-75', '76-85', '86+',])
bd_persons_ages_cohorts = bd_persons_ages.groupby(
    ['AgeCohort']).size().reset_index(name='counts')
bd_persons_ages_cohorts.index = bd_persons_ages_cohorts['AgeCohort'].to_list()
bd_persons_ages_cohorts.plot.bar(figsize=(
    10, 7), title='Age distribution', xlabel='Age cohorts', ylabel='Number of persons', legend=None)
plt.savefig(os.path.join(figs_path, 'age_distribution.png'), dpi=300)

# get the id with the persons with the highest number and print them
[id_person_max_trips, id_person_max_1_trips] = [
    bd_persons['personID'][0], bd_persons['personID'][1]]
print(get_trip_chain_activity_after(id_person_max_trips))

# select internal and external trips (from/to Braunschweig)
bd_persons_inside_bs = bd.loc[bd['startZone'].isin(
    bd_tz_bs['idTrafficZone']) & bd['destZone'].isin(bd_tz_bs['idTrafficZone'])]
bd_persons_outside_bs = bd.loc[~bd['startZone'].isin(
    bd_tz_bs['idTrafficZone']) & ~bd['destZone'].isin(bd_tz_bs['idTrafficZone'])]

# check how many trips are in the same traffic zone
bd_taz_start = bd.groupby('startZone').size().reset_index(name='counts')
#bd_same_taz = bd.loc[bd['startZone'] == bd['destZone']]
bd_different_taz = bd.loc[bd['startZone'] != bd['destZone']]
#bd_same_taz = bd_same_taz.groupby(['startZone']).size().reset_index(name='counts')
bd_different_taz = bd_different_taz.groupby(
    ['startZone']).size().reset_index(name='counts')
#bd_same_taz['percentage'] = bd_same_taz['counts']/bd_taz_start['counts']
bd_different_taz['complement_counts'] = bd_taz_start['counts'] - \
    bd_different_taz['counts']
bd_different_taz['total_counts'] = bd_taz_start['counts']
bd_different_taz = (bd_different_taz.sort_values(
    by='total_counts', ascending=False, ignore_index=True))
plt.figure()
plt.plot(bd_different_taz['counts'])
plt.plot(bd_different_taz['complement_counts'])
plt.legend(['different TAZ', 'same TAZ'])
plt.xlabel('Traffic analysis zones')
plt.ylabel('Amount of trips inside/outside a TAZ')
plt.savefig(os.path.join(figs_path, 'TAZ_distribution.png'), dpi=300)

# which age takes which mode of transport
bd_persons_age_vehicle_choice = bd.loc[:, ['personID', 'age', 'vehicleChoice']].sort_values(
    by=['age'], ascending=True, ignore_index=True)
bd_persons_age_vehicle_choice['vehicleChoice'] = bd_persons_age_vehicle_choice['vehicleChoice'].map(
    dict_vehicle)
# accumulate the number of trips per age and vehicle choice
bd_persons_age_vehicle_choice = bd_persons_age_vehicle_choice.groupby(
    ['age', 'vehicleChoice']).size().reset_index(
    name='counts').sort_values(
    by=['age'],
    ascending=True, ignore_index=True)
# assign each age to an age cohort
bd_persons_age_vehicle_choice['ageCohort'] = pd.cut(bd_persons_age_vehicle_choice['age'], bins=[-1, 18, 25, 35, 45, 55, 65, 75, 85, 105], labels=[
                                                    '0-18', '19-25', '26-35', '36-45', '46-55', '56-65', '66-75', '76-85', '86+',])
# plot a cake chart for each age cohort
bd_persons_age_vehicle_choice_cake = bd_persons_age_vehicle_choice.groupby(
    ['ageCohort', 'vehicleChoice']).sum().reset_index()
bd_persons_age_vehicle_choice_cake_age_veh = bd_persons_age_vehicle_choice_cake.pivot(
    index='ageCohort', columns='vehicleChoice', values='counts')
bd_persons_age_vehicle_choice_cake_age_veh.plot.pie(subplots=True, layout=(3, 2), figsize=(
    10, 10), title=bd_persons_age_vehicle_choice_cake_age_veh.columns.to_list(), legend=False, ylabel='')
plt.savefig(os.path.join(figs_path, 'vehicle_choice_per_age.png'), dpi=300)

# switch age cohort and vehicle choice
bd_persons_age_vehicle_choice_cake_veh_age = bd_persons_age_vehicle_choice_cake.pivot(
    index='vehicleChoice', columns='ageCohort', values='counts')
bd_persons_age_vehicle_choice_cake_veh_age.plot.pie(subplots=True, layout=(3, 3), figsize=(
    20, 10), title=bd_persons_age_vehicle_choice_cake_veh_age.columns.to_list(), legend=False, ylabel='')
plt.savefig(os.path.join(figs_path, 'vehicle_choice_reversed.png'), dpi=300)

# Analyse trip purposes per age cohort
bd_persons_age_purpose = bd.loc[:, ['personID', 'age', 'ActivityAfter']]
bd_persons_age_purpose['ActivityAfter'] = bd_persons_age_purpose['ActivityAfter'].map(
    dict_leisure)
bd_persons_age_purpose['ageCohort'] = pd.cut(bd_persons_age_purpose['age'], bins=[-1, 18, 25, 35, 45, 55, 65, 75, 85, 105], labels=[
    '0-18', '19-25', '26-35', '36-45', '46-55', '56-65', '66-75', '76-85', '86+',])
bd_persons_age_purpose_acc = bd_persons_age_purpose.groupby(
    ['ActivityAfter', 'ageCohort']).size().reset_index(name='counts')
bd_persons_age_purpose_cake = bd_persons_age_purpose_acc.pivot(
    index='ageCohort', columns='ActivityAfter', values='counts')
bd_persons_age_purpose_cake.plot.pie(subplots=True, layout=(2, 4), figsize=(
    15, 5), title=bd_persons_age_purpose_cake.columns.to_list(), legend=False, ylabel='')
plt.savefig(os.path.join(figs_path, 'trip_purpose.png'), dpi=300)

# switch age cohorts and purpose
bd_persons_age_purpose_cake_reversed = bd_persons_age_purpose_acc.pivot(
    index='ActivityAfter', columns='ageCohort', values='counts')
bd_persons_age_purpose_cake_reversed.plot.pie(subplots=True, layout=(2, 5), figsize=(
    15, 5), title=bd_persons_age_purpose_cake_reversed.columns.to_list(), legend=False, ylabel='')
plt.savefig(os.path.join(figs_path, 'trip_purpose_reversed.png'), dpi=300)

# Count traffic zones where people are starting their trips
bd_traffic_zones_start = bd.groupby(['startZone']).size().reset_index(
    name='counts').sort_values(by=['counts'], ascending=False, ignore_index=True)
for row in bd_traffic_zones_start.iterrows():
    bd_traffic_zones_start.loc[row[0], ['longitude', 'latitude']
                               ] = bd.loc[bd_traffic_zones_start.loc[row[0], 'startZone'] == bd['startZone'], ['loCs', 'laCs']].values[0]
print(bd_traffic_zones_start)

# plot longitude and latitude of traffic zones
bd_traffic_zones_start.plot(
    kind='scatter', x='longitude', y='latitude', s=bd_traffic_zones_start['counts'] / 100, figsize=(10, 15),
    title='Traffic zones where people are starting their trips (Braunschweig)')
plt.savefig(os.path.join(figs_path, 'traffic_zones_start.png'), dpi=300)

# Count traffic zones where people are going and leaving
bd_traffic_zones_end = bd.groupby(['destZone']).size().reset_index(
    name='counts').sort_values(by=['counts'], ascending=False, ignore_index=True)
for row in bd_traffic_zones_end.iterrows():
    bd_traffic_zones_end.loc[row[0], ['longitude', 'latitude']
                             ] = bd.loc[bd_traffic_zones_end.loc[row[0], 'destZone'] == bd['destZone'], ['loCs', 'laCs']].values[0]
print(bd_traffic_zones_end)

# plot longitude and latitude of traffic zones
bd_traffic_zones_end.plot(
    kind='scatter', x='longitude', y='latitude', s=bd_traffic_zones_end['counts'] / 100, figsize=(10, 10),
    title='Traffic zones where people are ending their trips (Braunschweig)')
plt.savefig(os.path.join(figs_path, 'traffic_zones_end.png'), dpi=300)

# Time analyzing
# Time of day
bd_time = bd.groupby(['startTime']).size().reset_index(name='counts')
print(bd_time)
bd_time.plot(kind='line', x='startTime', y='counts')
plt.xlabel('Time of day')
plt.ylabel('Number of trips')
plt.title('Time of day (Braunschweig)' + ' (n = ' + str(len(bd)) + ' trips)')
plt.savefig(os.path.join(figs_path, 'time_of_day.png'), dpi=300)

# time of day rolling average
bd_time['rolling_mean'] = bd_time['counts'].rolling(window=7).mean()
bd_time.plot(kind='line', x='startTime', y='rolling_mean')
plt.xlabel('Time of day')
plt.ylabel('Number of trips')
plt.title('Time of day (Braunschweig)' + ' (n = ' +
          str(len(bd)) + ' trips), rolling mean with window size 7')
plt.savefig(os.path.join(figs_path, 'time_of_day_rolling.png'), dpi=300)

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
plt.xlabel('Leisure activity after trip')
plt.ylabel('Leisure activity before trip')
plt.title('Frequency of trips from which leisure activity to which leisure activity')
plt.savefig(os.path.join(figs_path, 'heatmap_leisure.png'), dpi=300)


# quick check if a data point in the matrix is correct
check = bd.loc[(bd['startZone'] == 31010011) & (
    bd['destZone'] == 31010011), 'startZone'].count()


# plotting the duration of trips
fig, axs = plt.subplots(2)
axs[0].hist(bd['travelTime']/60/60, bins=100,
            range=(0, max(bd['travelTime'])*1.2/60/60))
axs[0].set_title('Trip duration (Braunschweig)' +
                 ' (n = ' + str(len(bd)) + ' trips)')
axs[0].set_xlabel('Trip duration in hours')
axs[0].set_ylabel('Number of trips')
axs[1].hist(bd['travelTime']/60/60, bins=100,
            range=(0, max(bd['travelTime'])*1.2/60/60))
axs[1].set_yscale('log')
axs[1].set_title('Trip duration (Braunschweig)' +
                 ' (n = ' + str(len(bd)) + ' trips) log scale')
axs[1].set_xlabel('Trip duration in hours')
axs[1].set_ylabel('Number of trips')

# same thing with the distance of the trips
fig2, axs2 = plt.subplots(2)
axs2[0].hist(bd['tripDistance'], bins=100,
             range=(0, max(bd['tripDistance'])*1.2,))
axs2[0].set_title('Trip distance (Braunschweig)' +
                  ' (n = ' + str(len(bd)) + ' trips)')
axs2[0].set_xlabel('Trip distance in Kilometers')
axs2[0].set_ylabel('Number of trips')
axs2[1].hist(bd['tripDistance'], bins=100,
             range=(0, max(bd['tripDistance'])*1.2,))
axs2[1].set_yscale('log')
axs2[1].set_title('Trip distance (Braunschweig)' +
                  ' (n = ' + str(len(bd)) + ' trips) log scale')
axs2[1].set_xlabel('Trip distance in Kilometers')
axs2[1].set_ylabel('Number of trips')


# analyze age distribution with age cohorts
bd_persons_id_and_age = bd[['personID', 'age']].drop_duplicates()
bd_age_cohorts = pd.cut(
    bd_persons_id_and_age['age'],
    bins=[-1, 18, 25, 35, 45, 55, 65, 75, 85, 95, 106],
    labels=['0-18', '19-25', '26-35', '36-45', '46-55', '56-65', '66-75', '76-85', '86-95', '96-105']).reset_index(
    name='age_cohort').groupby(
    ['age_cohort']).size().reset_index(
    name='counts')
bd_age_cohorts.plot(kind='bar', x='age_cohort', y='counts')
plt.xlabel('Age cohort')
plt.ylabel('Number of persons')
plt.title('Age distribution of persons in Braunschweig in age cohorts, number ')
plt.savefig(os.path.join(figs_path, 'age_cohorts.png'), dpi=300)

# relation between age and trip duration in minutes and hours
bd_age_duration = bd.groupby(['age']).mean().reset_index()
bd_age_duration['travelTime'] = bd_age_duration['travelTime'] / 60
bd_age_duration.plot(kind='bar', x='age', y='travelTime')

# also do this for trip distance
bd_age_distance = bd.groupby(['age']).mean().reset_index()
bd_age_distance.plot(kind='bar', x='age', y='tripDistance')

# analyze trip distance in distance cohort
bd_trip_distance_cohorts = pd.cut(
    bd['tripDistance'],
    bins=[-1, 1, 2, 5, 10, 20, 50, 423789798324],
    labels=['0-1', '1-2', '2-5', '5-10', '10-20', '20-50', '50+']).reset_index(
    name='trip_distance_cohort').groupby(
    ['trip_distance_cohort']).size().reset_index(
    name='counts')
bd_trip_distance_cohorts.plot(kind='bar', x='trip_distance_cohort', y='counts')
plt.xlabel('Trip distance cohort in km')
plt.ylabel('Number of trips')
plt.title('Trip distance distribution in Braunschweig in trip distance cohorts, in km. Number of trips: ' + str(len(bd)))
plt.savefig(os.path.join(figs_path, 'trip_distance_cohorts.png'), dpi=300)

# pie diagram of trip distance cohorts for each age cohort
bd_trip_distance_cohorts_vehicle = pd.cut(
    bd['tripDistance'],
    bins=[-1, 1, 2, 5, 10, 20, 50, 423789798324],
    labels=['0-1', '1-2', '2-5', '5-10', '10-20', '20-50', '50+']).reset_index(
    name='trip_distance_cohort')
bd_age_cohorts = pd.cut(bd['age'],
                        bins=[-1, 18, 25, 35, 45, 55, 65, 75, 85, 95],
                        labels=['0-18', '19-25', '26-35', '36-45', '46-55', '56-65', '66-75', '76-85', '86-95']).reset_index(
    name='age_cohort')
# attach the age cohort to the dataframe
bd_trip_distance_cohorts_vehicle = pd.concat(
    [bd_trip_distance_cohorts_vehicle['trip_distance_cohort'], bd_age_cohorts['age_cohort']], axis=1)
# group by age cohort and distance cohort
bd_trip_distance_cohorts_vehicle = bd_trip_distance_cohorts_vehicle.groupby(
    ['trip_distance_cohort', 'age_cohort']).size().reset_index(name='counts')
# plot a pie chart for each age cohort
for age_cohort in bd_trip_distance_cohorts_vehicle['age_cohort'].unique():
    bd_trip_distance_cohorts_vehicle[bd_trip_distance_cohorts_vehicle['age_cohort'] == age_cohort].plot(
        kind='pie', y='counts', labels=bd_trip_distance_cohorts_vehicle[bd_trip_distance_cohorts_vehicle['age_cohort'] == age_cohort]['trip_distance_cohort'], autopct='%1.1f%%', startangle=90)
    plt.title('Trip distance distribution in Braunschweig in trip distance cohorts, in km. Number of trips: ' +
              str(len(bd)) + ' for age cohort ' + age_cohort)
    plt.savefig(os.path.join(
        figs_path, 'trip_distance_cohorts_' + age_cohort + '.png'), dpi=300)

# now the same for trip duration vs vehicle type
bd_trip_duration_cohorts_vehicle = pd.cut(
    bd['travelTime'],
    bins=[-1, 500, 1000, 5000, 10000, 423789798324],
    labels=['0-500', '500-1000', '1000-5000', '5000-10000', '10000+']).reset_index(
    name='trip_duration_cohort')
bd_vehicle_type = bd['vehicleChoice'].reset_index(name='vehicleChoice')
# attach the age cohort to the dataframe
bd_trip_duration_cohorts_vehicle = pd.concat(
    [bd_trip_duration_cohorts_vehicle['trip_duration_cohort'], bd_vehicle_type['vehicleChoice']], axis=1)
# group by age cohort and distance cohort
bd_trip_duration_cohorts_vehicle = bd_trip_duration_cohorts_vehicle.groupby(
    ['trip_duration_cohort', 'vehicleChoice']).size().reset_index(name='counts')
# plot a pie chart for each age cohort
for vehicle_type in bd_trip_duration_cohorts_vehicle['vehicleChoice'].unique():
    bd_trip_duration_cohorts_vehicle[bd_trip_duration_cohorts_vehicle['vehicleChoice'] == vehicle_type].plot(
        kind='pie', y='counts', labels=bd_trip_duration_cohorts_vehicle[bd_trip_duration_cohorts_vehicle['vehicleChoice'] == vehicle_type]['trip_duration_cohort'], autopct='%1.1f%%', startangle=90)
    plt.title('Trip duration distribution in Braunschweig in trip duration cohorts, in minutes. Number of trips: ' + str(len(bd)
                                                                                                                         ) + ' for vehicle type ' + str(dict_vehicle[vehicle_type]))
    plt.savefig(os.path.join(figs_path, 'trip_duration_cohorts_' +
                str(dict_vehicle[vehicle_type]) + '.png'), dpi=300)


# also do a scatter plot of trip duration and trip distance with a regression line
bd_age_distance.plot(kind='scatter', x='age', y='tripDistance')
plt.plot(np.unique(bd_age_distance['age']), np.poly1d(np.polyfit(bd_age_distance['age'],
         bd_age_distance['tripDistance'], 1))(np.unique(bd_age_distance['age'])))
bd_age_duration.plot(kind='scatter', x='age', y='travelTime')
plt.plot(np.unique(bd_age_duration['age']), np.poly1d(np.polyfit(
    bd_age_duration['age'], bd_age_duration['travelTime'], 1))(np.unique(bd_age_duration['age'])))
plt.show()




x = 42


### draft to compute number of locations of type 6?
x = bd[bd.ActivityAfter==6].loc[:, 'loc_id_end'].nunique()
# downscaling of locations of type 6 to 2000 agents
x / (bd['personID'].nunique()/2000)