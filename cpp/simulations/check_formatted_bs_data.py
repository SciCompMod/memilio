import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from datetime import datetime, date


####### minimal sanity check on data #######
bd = pd.read_csv(
    '/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/mobility/braunschweig_result_ffa8_modified4.csv', header=None, skiprows=1)


bd.rename(
    columns={0: 'personID', 1: 'locID', 2: 'home_id', 3: 'startTime', 4: 'lon_end',5:'lat_end', 10:'home_in_bs', 11:'has_home_trip', 12:'location_type',
             8:'map_feature_key', 9:'map_feat_value',7:'age', 6:'activity'},
    inplace=True)

def calc_minutes(time):
    #calculate minutes from string with this format hh:mm
    return int(time.split(':')[0])*60 + int(time.split(':')[1])

location_type_dict = {'Home': 0, 'School': 2,
                        'Work': 1, 'BasicShop': 3, 'SocialEvent': 4}

minutes_long_trip = 420

bd=bd[['personID', 'startTime', 'map_feature_key','map_feat_value','has_home_trip','location_type','activity','locID','age','home_id']]
#for debug purposes we just use the first 20000 entries
# bd = bd.head(100000)

# find out if we have a locID with different location types
locID_location_type = bd[['locID', 'location_type']].drop_duplicates()
locID_location_type_grouped = locID_location_type.groupby(['locID']).agg({'location_type': 'nunique'}).reset_index()
print('Number of locID with different location types: ', len(locID_location_type_grouped[locID_location_type_grouped['location_type'] > 1]))
# print which locID have different location types
print('locID with different location types: ', locID_location_type_grouped[locID_location_type_grouped['location_type'] > 1], 'percentually: ', len(locID_location_type_grouped[locID_location_type_grouped['location_type'] > 1])/len(locID_location_type_grouped))

#take out trips which go to locIDs with different location types
# bd = bd[~bd['locID'].isin(locID_location_type_grouped[locID_location_type_grouped['location_type'] > 1]['locID'])]
# locID_location_type = bd[['locID', 'location_type']].drop_duplicates()
# locID_location_type_grouped = locID_location_type.groupby(['locID']).agg({'location_type': 'nunique'}).reset_index()
# print('Number of locID with different location types: ', len(locID_location_type_grouped[locID_location_type_grouped['location_type'] > 1]))
# # print which locID have different location types
# print('locID with different location types: ', locID_location_type_grouped[locID_location_type_grouped['location_type'] > 1], 'percentually: ', len(locID_location_type_grouped[locID_location_type_grouped['location_type'] > 1])/len(locID_location_type_grouped))

# now we want to see who is where at each timestep with all the trips and see how many persons are at each location type
# for this we fill a matrix, where on the y-axis represent unique locations and on the x-axis the time for each timestep (one hour) and we fill the matrix with a vector the amount of persons at each location at each timestep
# we need to iterate over the data and fill the matrix
bd_who_where = bd.copy()[['personID', 'startTime', 'locID', 'location_type']]
bd_who_where['startTime'] = pd.to_datetime(bd_who_where['startTime'], format='%H:%M')
bd_who_where['startTime'] = bd_who_where['startTime'].dt.hour
# we need to get the unique locations
unique_locations = bd_who_where['locID'].unique()
# we need to get the unique persons
unique_persons = bd_who_where['personID'].unique()

# we need to fill the matrix with the amount of persons at each location at each timestep
matrix = np.zeros((len(unique_locations), 24))
for i, row in bd_who_where.iterrows():
    next_trip_time = bd_who_where.iloc[i + 1]['startTime'] if i < len(bd_who_where) - 1 and bd_who_where.iloc[i + 1]['personID'] == row['personID'] else bd_who_where.iloc[bd_who_where[bd_who_where['personID'] == row['personID']].index[0]]['startTime']
    for hour in range(row['startTime'], next_trip_time):
        matrix[np.where(unique_locations == row['locID'])[0][0], hour % 24] += 1
    if next_trip_time < row['startTime']:
        for hour in range(row['startTime'] , 24):
            matrix[np.where(unique_locations == row['locID'])[0][0], hour % 24] += 1
        for hour in range(0, next_trip_time):
            matrix[np.where(unique_locations == row['locID'])[0][0], hour] += 1
    # if  next_trip_time == row['startTime'] and its the last trip of the person we need to fill the matrix for the last trip for the rest of the day
    if  next_trip_time == row['startTime'] &(i == len(bd_who_where) - 1 or bd_who_where.iloc[i + 1]['personID'] != row['personID']):
        for hour in range(row['startTime'], 24):
            matrix[np.where(unique_locations == row['locID'])[0][0], hour % 24] += 1
        for hour in range(0, row['startTime']):
            matrix[np.where(unique_locations == row['locID'])[0][0], hour] += 1


# quick check: for every hour we we need the sum of the matrix column to be the amount of persons
print('Sum of the matrix columns: ', matrix.sum(axis=0))
print('Should be the amount of persons: ', len(unique_persons))

# we assign each unique location the most amount of persons at this location at any timestep
unique_locations_with_location_type = bd_who_where[['locID']].drop_duplicates()
#vector with max for each location
max_persons = matrix.max(axis=1)
#assign eacj location the max amount of persons
unique_locations_with_location_type['max_persons'] = max_persons
#add column with location type
unique_locations_with_location_type=unique_locations_with_location_type.sort_values(by=['max_persons'], ascending=False, ignore_index=True).head(1000)



# find the first trip of the 10 locations with the most persons


#print top ten locations with the most persons
print('Top ten locations with the most persons: ',unique_locations_with_location_type)
#find the location type of the top ten locations
unique_locations_with_location_type['location_type'] = 0
for i, row in unique_locations_with_location_type.iterrows():
    unique_locations_with_location_type.at[i, 'location_type'] = bd_who_where[bd_who_where['locID'] == row['locID']].iloc[0]['location_type']
print('Top ten locations with the most persons and their location type: ',unique_locations_with_location_type)
##count the location type of these locaitons
unique_locations_with_location_type = unique_locations_with_location_type['location_type'].value_counts().reset_index(name='counts').sort_values(by=['counts'], ascending=False, ignore_index=True)
print('Top ten locations with the most persons and their location type: ',unique_locations_with_location_type)



# bd_hour_same=bd.copy()
# # how many trips haappen in the same hour as the previous trip of the same person the time before the forst trip is the last trip of the person
# bd_hour_same['startTime'] = pd.to_datetime(bd_hour_same['startTime'], format='%H:%M')
# bd_hour_same['startTime_prev'] = bd_hour_same.groupby(['personID'])['startTime'].shift(1)
# bd_last_trip = bd_hour_same.groupby(['personID']).last()
# # to fill NaN we take those and fill them with the last trip of the person with fillna and bd_last_trip['startTime']
# bd_hour_same['startTime_prev'] = bd_hour_same['startTime_prev'].fillna(bd_last_trip['startTime'].values[0])
# bd_hour_same['same_hour_as_prev'] = bd_hour_same['startTime'].dt.hour == bd_hour_same['startTime_prev'].dt.hour
# same_hour_as_prev = bd_hour_same['same_hour_as_prev'].value_counts(normalize=False).reset_index(name='counts')
# print('Percentage of trips that happen in the same hour as the previous trip of the same person: ', same_hour_as_prev)
# #also get the amount of trips that are hometrips and happen before another trip in the same hour

# # we iterate over the data and save the trip before a trip in the same hour
# bd_hour_same['trip_before_same_hour'] = 0
# for i, row in bd_hour_same.iterrows():
#     if row['same_hour_as_prev']:
#         bd_hour_same.at[i-1, 'trip_before_same_hour'] = 1

# # we want to see how many trips are home trips and happen before another trip in the same hour
# home_trip_before_same_hour = bd_hour_same[(bd_hour_same['trip_before_same_hour'] == 1) & (bd_hour_same['location_type'] == location_type_dict['Home'])]
# print('Number of home trips that happen before another trip in the same hour: ', len(home_trip_before_same_hour), '(which is percentually from home trips: ', len(home_trip_before_same_hour)/(bd_hour_same['location_type'] == location_type_dict['Home']).sum(), ')')




# # we meed to add a home trip for each person if they have no home trip at 17:00 into the data
# # ths easiest done when we go though it iteratively and add a trip when we see that the person has no home trip and it is either the last trip and before 17 or not the last trip but the first after 17
# # the trip has the same personID, location_type = 0, startTime = 17:00, map_feature_key = 0, map_feat_value = 0, has_home_trip = 1
# # the bd is sorted so we can go through it iteratively
# bd_addendum = pd.DataFrame(columns=['personID', 'startTime', 'map_feature_key','map_feat_value','has_home_trip','location_type'])
# for i, row in bd.iterrows():
#     is_new_person = i == 0 or bd.iloc[i-1]['personID'] != row['personID']
#     if is_new_person:
#         home_trip_added = False
#     is_last_trip_of_person = i == len(bd) - 1 or bd.iloc[i+1]['personID'] != row['personID']
#     is_first_trip_after_17 = row['startTime'] > '17:00' and bd.iloc[i-1]['startTime'] < '17:00'
#     if row['has_home_trip'] == 0 and (is_last_trip_of_person or is_first_trip_after_17):
#         add_trip = {'personID': row['personID'], 'startTime': '17:00', 'has_home_trip': 0, 'location_type': 0, 'map_feature_key': 'home', 'map_feat_value': 'Nan', 'activity': 'home'}
#         #we save this trip in a new dataframe and append it to the bd at the end and sort it again
#         bd_addendum = pd.concat([bd_addendum, pd.DataFrame([add_trip])], ignore_index=True)
#         home_trip_added = True
    
# bd = pd.concat([bd, bd_addendum], ignore_index=True)
# bd = bd.sort_values(by=['personID', 'startTime'], ignore_index=True)

# # for a quick check: see how many trips where added and how many unique persons are there with no home trip
# print('Number of trips added: ', len(bd_addendum))
# print('Number of unique persons with no home trip: ', len(bd[bd['has_home_trip'] == 0]['personID'].unique()))


# # now we want to get a dataframe with every trip of a person which is above 6 hours. for this we go through it again and save the length of each trip in a new column
# # add column for trip length
# bd['trip_length'] = 0
# # if the trip is longer than 6 hours we save it in a new dataframe
# long_trips = pd.DataFrame(columns=['personID', 'startTime', 'map_feature_key','map_feat_value','has_home_trip','location_type', 'trip_length','activity'])

# for i, row in bd.iterrows():
#     last_trip_of_person = i == len(bd) - 1 or bd.iloc[i+1]['personID'] != row['personID']
#     first_trip_of_person = i == 0 or bd.iloc[i-1]['personID'] != row['personID']
#     if first_trip_of_person:
#         time_of_first_trip = calc_minutes(row['startTime'])
#     if not last_trip_of_person:
#         duration = calc_minutes(bd.iloc[i+1]['startTime']) - calc_minutes(row['startTime'])
#     else:
#         duration = 1440 - calc_minutes(row['startTime']) + time_of_first_trip
#     bd.at[i, 'trip_length'] = duration
#     if duration > minutes_long_trip:
#         long_trips = pd.concat([long_trips, pd.DataFrame(bd.iloc[i]).T], ignore_index=True)

# # add the value in trip_length for every person in a seperate datafram
# bd_trip_length_cum = bd.groupby(['personID']).agg({'trip_length': 'sum'}).reset_index()
# # quick check: every person has 24 hours in minutes in the trip_length column
# print('Number of persons without 24 hours in trip_length: ', len(bd_trip_length_cum[
#     bd_trip_length_cum['trip_length'] != 1440]))


# # now we want to discard the trips which are longer than 6 hours and are home trips
# long_trips_wo_homes = long_trips[long_trips['location_type'] != location_type_dict['Home']]
# # quick check: how many long trips are home trips
# print('Number of long trips which are home trips: ', len(long_trips[long_trips['location_type'] == location_type_dict['Home']]))
# # quick check: how many long trips are not home trips
# print('Number of long trips which are not home trips: ', len(long_trips_wo_homes), 'which is percentually: ', len(long_trips_wo_homes)/len(long_trips))
# # quick check: how many long trips are there in total
# print('Number of long trips: ', len(long_trips))

# # now we want to get the location types of the long trips and sort it by location type
# long_trips_location_type = long_trips['location_type'].value_counts().reset_index(name='counts').sort_values(by=['location_type'], ascending=True, ignore_index=True)
# # quick check: how many long trips are there for each location type
# print('Number of long trips for each location type: ', long_trips_location_type, 'where the location types are: ', location_type_dict)

# # now we want to see the avarage amount of the long trips for each location type
# long_trips_location_type_avg = long_trips['location_type'].value_counts(normalize=True).reset_index(name='counts').sort_values(by=['location_type'], ascending=True, ignore_index=True)
# # quick check: how many long trips are there for each location type
# print('Avarage amount of long trips for each location type: ', long_trips_location_type_avg)

# # now we want to substract 6 hours from the long trips length
# long_trips_substracted = long_trips.copy()
# long_trips_substracted['trip_length'] = long_trips_substracted['trip_length'] - minutes_long_trip
# # sum over the trip lengths for each location type
# long_trips_substracted_sum = long_trips_substracted.groupby(['location_type']).agg({'trip_length': 'sum'}).reset_index()
# long_trips_substracted_avg = long_trips_substracted_sum.copy()
# long_trips_substracted_avg['trip_length'] = long_trips_substracted_avg['trip_length'] / long_trips_location_type['counts']

# print('Avarage amount of long trips for each location type after substracting long_trip minutes: ', long_trips_substracted_avg)

# # show me the long trips of basicShop
# long_trips_SE = long_trips[long_trips['location_type'] == location_type_dict['SocialEvent']]
# #group by the map_feature_key and map_feat_value
# long_trips_SE_grouped = long_trips_SE.groupby(['map_feature_key', 'map_feat_value']).agg({'trip_length': 'sum'}).reset_index().sort_values(by=['trip_length'], ascending=False, ignore_index=True)
# #note: most time is spend in commercial, supermarket and reatil which is classified as BS but maybe is work (short test around 40% is work i guess)

# # show me the long trips of basicShop
# long_trips_BS = long_trips[long_trips['location_type'] == location_type_dict['BasicShop']]
# #group by the map_feature_key and map_feat_value
# long_trips_BS_grouped = long_trips_BS.groupby(['map_feature_key', 'map_feat_value','activity']).agg({'trip_length': 'sum'}).reset_index().sort_values(by=['trip_length'], ascending=False, ignore_index=True)


# print('Long trips of SocialEvent: ', long_trips_SE_grouped.head(10))
# print('Long trips of BasicShop: ', long_trips_BS_grouped.head(10))








#############################



# #### Counting people ####
# bd_persons = bd.groupby(['personID']).size().reset_index(
#     name='counts').sort_values(by=['counts'], ascending=False, ignore_index=True)
# # Get the frequency of each number in the 'numbers' column
# # bd_personss = bd_persons['counts'].value_counts().sort_index()
# # # Plot the frequency of each number
# # bd_personss.plot(kind='bar')
# # plt.xlabel('Number of trips per day')
# # plt.ylabel('Number of persons')
# # plt.title('Number of trips per person tripcount (Braunschweig)')
# # plt.savefig(os.path.join(figs_path, 'number_of_trips_per_person.png'), dpi=300)

# # we want to calculate a matrix: the rows are the persons and the columns are the location types they visited in each cell we want to save the longest amount of time they spent in this location type
# # the time they spent in the location is the the startTime of the next trip minus the startTime of the current trip
# # we want to calculate the time they spent in each location type#

# def calc_minutes(time):
#     #calculate minutes from string with this format hh:mm
#     return int(time.split(':')[0])*60 + int(time.split(':')[1])

# #chekc how many have their last trip to home and how many not
# bd_last_trip = bd.groupby(['personID']).last()
# not_home = bd_last_trip[bd_last_trip['location_type'] != 0]
# home = bd_last_trip[bd_last_trip['location_type'] == 0]
# not_home_and_no_home_trip = bd_last_trip[(bd_last_trip['location_type'] != 0) & (bd_last_trip['has_home_trip'] == 0)]
# print('Number of last trips that are home trips: ', len(home))
# print('Number of last trips that are not home trips: ', len(not_home))

# #check how many persons have no home trip
# no_home_trip = bd[bd['has_home_trip'] == 0]
# print('Number of persons with no home trip: ', len(no_home_trip))

# # check how many persons habe no home trip and the last trip is not home
# no_home_trip_and_not_home = no_home_trip[no_home_trip['location_type'] != 0]
# print('Number of persons with no home trip and the last trip is not home: ', len(no_home_trip_and_not_home))

# #check how many of these are after 17:00
# no_home_trip_and_not_home_after_17 = no_home_trip_and_not_home[no_home_trip_and_not_home['startTime'] > '17:00']
# print('Number of persons with no home trip and the last trip is not home and after 17:00: ', len(no_home_trip_and_not_home_after_17))

# #show a distribution of the location types of the last trips of persons with no home trip and the last trip is not home and after 17:00 
# no_home_trip_and_not_home_after_17['location_type'].value_counts().plot(kind='bar')
# plt.xlabel('Location type')
# plt.ylabel('Number of persons')
# plt.title('Location type of the last trips of persons with no home trip and the last trip is not home and after 17:00')
# plt.show()





# # number_of_persons times number_of_location_types matrix
# matrix = np.zeros((len(bd_persons), 8))
# long_trips_matrix = []  # Matrix to store persons with long trips
# bd =  bd[bd['personID'] == 900971]

# # iterate over the persons
# for i, person in enumerate(bd_persons['personID']):
#     # get the trips of the person
#     trips_person = bd.loc[bd['personID'] == person]
    
#     #if person has no home trip add one with starttime 17:00
#     if len(trips_person[trips_person['has_home_trip'] == 1]) == 0:
#         add_trip = {'personID': 
#             person, 'loc_d_end': 0, 'home_id': 0, 'startTime': '17:00', 'lon_end': 0, 'lat_end': 0, 'lon_end_2': 0, 'home_in_bs': 0, 'has_home_trip': 1, 'location_type': 0, 'map_feature_key': 0, 'map_feat_value': 0}
#         trips_person.loc[len(trips_person)] = add_trip
        
#     # iterate over the trips
#     for k, trip in enumerate(trips_person.iterrows()):
#         # if it is the last trip
#         if k == len(trips_person) - 1:
#             duration = 1440 - calc_minutes(trips_person.iloc[k]['startTime']) + calc_minutes(trips_person.iloc[0]['startTime'])
#         else:
#             duration = calc_minutes(trips_person.iloc[k+1]['startTime']) - calc_minutes(trips_person.iloc[k]['startTime'])
   
#         matrix[i, trips_person.iloc[k]['location_type']] = max(duration, matrix[i, trips_person.iloc[k]['location_type']])
#         if duration > 480 and k != len(trips_person) - 1:  # 8 hours and not the last trip
#             long_trips_matrix.append(trips_person.iloc[k])

# # Convert list of long trips to DataFrame
# long_trips_df = pd.DataFrame(long_trips_matrix)

# long_trips_df = long_trips_df[['personID', 'startTime', 'location_type']]

# # amount of trips for each location type
# location_type_counts = long_trips_df.groupby(['location_type']).size().reset_index(name='counts').sort_values(by=['counts'], ascending=False, ignore_index=True)

# print(f"Number of persons with trips over 8 hours: {len(long_trips_matrix)}")

# # how many long trips have their last trip to home
# long_trips_df_last_trip = long_trips_df.groupby(['personID']).last()
# long_trips_df_last_trip_home = long_trips_df_last_trip[long_trips_df_last_trip['location_type'] == 0]
# print('Number of last trips that are home trips: ', len(long_trips_df_last_trip_home))
# long_trips_df_last_trip_not_home = long_trips_df_last_trip[long_trips_df_last_trip['location_type'] != 0]
# print('Number of last trips that are not home trips: ', len(long_trips_df_last_trip_not_home))

# # how many persons have no home trip
# long_trips_df_no_home_trip = long_trips_df[long_trips_df['has_home_trip'] == 0]
# print('Number of persons with no home trip: ', len(long_trips_df_no_home_trip))

# # group the trips which are longer than 8 hours by the tuple (map_feature_key, map_feat_value)
# grouped = long_trips_df.groupby(['personID', 'location_type']).agg({'startTime': 'count'}).reset_index()




# # save matrix as heatmap
# sns.heatmap(matrix)
# plt.show()
# # save matrix as png
# plt.savefig(os.path.join("/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/cpp/models/abm", 'matrix.png'), dpi=300)


# matrix