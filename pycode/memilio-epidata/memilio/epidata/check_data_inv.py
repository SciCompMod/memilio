import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import os
import math


# df_data_demo_1 = pd.read_csv(r'C:\Users\korf_sa\Documents\rep\pycode\memilio-epidata\memilio\epidata\Datenlieferung_20231214\Datenlieferung_20231214\home-location\demo\home-location_demo_2020-04-06.csv')
# df_data_tapas_1 = pd.read_csv(r'C:\Users\korf_sa\Documents\rep\pycode\memilio-epidata\memilio\epidata\Datenlieferung_20231214\Datenlieferung_20231214\home-location\tapas\home-location_tapas_2020-04-01.csv')
#columns:
# 0: area_id  
# 1: gender
# 2: age_group
# 3: home_locations_count
path = r"C:\Users\korf_sa\Documents\rep\pycode\memilio-epidata\memilio\epidata\Datenlieferung_20231214\Datenlieferung_20231214\home-location\demo"
dir = os.listdir( path )



#check if unique IDs are the same in all files and sorted
count_unique_id_demo = ['31010011', '31010021', '31010031', '31010041', '31010051',
       '31010061', '31010071', '31010081', '31010091', '31010101',
       '31010111', '31010121', '31010131', '31010141', '31010151',
       '31010161', '31010171', '31010181', '31010191', '31010201',
       '31010211']

count_unique_id_tapas = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '2',
       '21', '22', '23', '24', '25', '26', '27', '28', '3', '34', '35',
       '36', '37', '38', '4', '42', '43', '44', '45', '46', '47', '48',
       '49', '5', '53', '54', '55', '56', '59', '6', '60', '61', '64',
       '65', '66', '67', '68', '69', '7', '71', '72', '73', '74', '8',
       '9', '9115', '9116', '9180', '9238', '9250', '9251', '9930',
       '9932', '9940', '9970']



len_unique_id_1 = len(count_unique_id_demo)
#count unique IDs
for file in dir:
    file = pd.read_csv(path + '\\' + file)
    file = file[~file['area_id'].str.startswith('DE')]
    count_unique_id = file['area_id'].unique()
    #check if unique IDs are the same in all files
    if not np.array_equal(count_unique_id, count_unique_id_demo):
        print('not equal')

## check number of rows with NaN values per unique ID
max_number_of_rows_with_no_nan = 14 # because male + female and 7 age groups
abc = []
#percentage of rows with any NaN values
abc_total = []
for file in dir:
    #unique IDs in file
    file = pd.read_csv(path + '\\' + file)
    #delete every row where ID begins with DE
    file = file[~file['area_id'].str.startswith('DE')]
    #count unique IDs
    count_unique_id = file['area_id'].unique()  
    # for unique in count_unique_id:
    abc = []
    for unique in count_unique_id:
        # if unique == '9250':
        # file with just rows that have the unique ID
        file_id = file[file['area_id']==unique]
        #count how many rows have two NaN values
        count_two_nan = file_id[['gender', 'age_group']].isnull().all(axis=1).sum()
        #count how many rows have any NaN values
        count_any_nan = file_id[['gender', 'age_group']].isnull().any(axis=1).sum()
        #count how many rows have no NaN values
        count_no_nan = len(file_id)-count_any_nan
        #count how many rows have just one NaN values
        count_one_nan = count_any_nan-count_two_nan
        #count which percent of 14 rows have no NaN values
        perc_no_nan_from_14 = count_no_nan/max_number_of_rows_with_no_nan
        print('count_no_nan: ', count_no_nan, 'count_one_nan: ', count_one_nan, 'count_two_nan: ', count_two_nan, 'count_any_nan: ', count_any_nan, 'perc_no_nan_from_14: ', perc_no_nan_from_14)
        #append to vector
        abc.append([count_no_nan, count_one_nan, count_two_nan, count_any_nan, perc_no_nan_from_14])
    
    # for every unique ID, add the abc vector to the old vector
    abc_total.append(abc)

#convert abc_total to numpy array
abc_total = np.array(abc_total)
#mean percent of unqiue ids with no NaN values from 14 rows
mean_perc_no_nan_from_14 = abc_total[:,:,4].mean(axis=0)
#then max and min:
print('min_perc_no_nan_from_14: ', min(mean_perc_no_nan_from_14))
print('max_perc_no_nan_from_14: ', max(mean_perc_no_nan_from_14))
print('mean_perc_no_nan_from_14: ', mean_perc_no_nan_from_14)
#sum up all the rows with the same unique ID
abc_total = abc_total.sum(axis=0)
#calculate the percentage of rows with any NaN values
perc_any_nan = abc_total[:,3]/(abc_total[:,0]+abc_total[:,1]+abc_total[:,2])
#calculate the percentage of rows with two NaN values
perc_two_nan = abc_total[:,2]/(abc_total[:,0]+abc_total[:,1]+abc_total[:,2])
#calculate the percentage of rows with one NaN values
perc_one_nan = abc_total[:,1]/(abc_total[:,0]+abc_total[:,1]+abc_total[:,2])
#calculate the percentage of rows with no NaN values
perc_no_nan = abc_total[:,0]/(abc_total[:,0]+abc_total[:,1]+abc_total[:,2])
#print the percentages for every unique ID
for i in range(len(count_unique_id_tapas)):
    print('ID: ', count_unique_id_tapas[i], 'perc_any_nan: ', perc_any_nan[i], 'perc_two_nan: ', perc_two_nan[i], 'perc_one_nan: ', perc_one_nan[i], 'perc_no_nan: ', perc_no_nan[i])
# also print min and max
print('min perc_any_nan: ', min(perc_any_nan), 'max perc_any_nan: ', max(perc_any_nan))
        
#count how many home locations are there with two NaN values
# abc=[]
# rel_na = []
# for file in dir:
#     #unique IDs in file
#     file = pd.read_csv(path + '\\' + file)
#     #delete every row where ID begins with DE
#     file = file[~file['area_id'].str.startswith('DE')]
#     #count unique IDs
#     count_unique_id = file['area_id'].unique()  
#     # for unique in count_unique_id:
#     abc_unique = []
#     for unique in count_unique_id:
#         ff = file[file['area_id']==unique]
#         #number of home locations with two NaN values
#         avail_sum = ff[~ff[['gender', 'age_group']].isna().all(axis=1)]['home_locations_count'].sum()
#         na_sum = ff[ff[['gender', 'age_group']].isna().any(axis=1)]['home_locations_count'].sum()
#         abc_unique.append(avail_sum/(avail_sum+na_sum))
#     abc.append(abc_unique)

# #convert abc to numpy array
# abc = np.array(abc)
# #calculate mean, min and max of every unique ID
# mean = abc.mean(axis=0)
# min = abc.min(axis=0)
# max = abc.max(axis=0)
# #print the mean, min and max for every unique ID
# for i in range(len(mean)):
#     print('ID: ', count_unique_id_1[i], 'mean: ', mean[i], 'min: ', min[i], 'max: ', max[i])
