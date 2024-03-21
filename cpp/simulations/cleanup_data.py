import pandas as pd
import numpy as np
import os
import sys


def add_school_ids(pd):
    # Schools are anonymized with loc_id_start and loc_id_end = -2
    # for every row where loc_id_start or loc_id_end is -2 we assume a school in the county and change the id to something unique (somewhere  betwenn-100 to -10000)
    pd_data = pd['person_id', 'loc_id_start', 'loc_id_end', 'countyStart', 'countyEnd', 'activity_start' , 'activity_end']
    #list of schools
    schools = np.array([])
    # we go through every row and check if the loc_id_start or loc_id_end is -2 and the activity_start or activity_end is 2
    # if so we change the loc_id_start or loc_id_end to a unique id which represents a school in that county
    # if in one unique school are more than 600 persons we change the id to a new unique id which represents a new school in that county
    # one thing to consider is that a usually goes to and leaves the same school so we have to check if the school is already in the list and the person already is assigned to that school
    
    # array to save how many schools in each country and how many persons in each school
    schools_in_county = np.array([])
    persons_in_schools = np.array([])

  
    
        
    





def add_home_ids(pd):
    # when there is a -1 in loc_id_start or loc_id_end, the person is at home (due to anonymization we dont know exactly if the person goes to their home or another we assume their own home)
    pd_data = pd['loc_id_start', 'loc_id_end', 'huid', 'activity_start', 'activity_end']
    #every row where loc_id_start or loc_id_end is -1 and activity_start or activity_end is 7 is a person who goes /comes from home
    # there we change the loc_id_start and loc_id_end to huid
    for index, row in pd.iterrows():
        if(row['loc_id_start'] == -1 and row['activity_start'] == 7):
            pd.at[index, 'loc_id_start'] = row['huid']
        if(row['loc_id_end'] == -1 and row['activity_end'] == 7):
            pd.at[index, 'loc_id_end'] = row['huid']
    return pd




def add_home_is_in_bs_column(pd):
    # check if csv file has column home_in_bs
    if 'home_in_bs' in pd.columns:
        print("Column 'home_in_bs' is present in DataFrame.")
        # nothing tbd
    else:
        print("Column 'home_in_bs' is not present in DataFrame.")
        # add column home_in_bs
        pd['home_in_bs'] = 0
        #list of persons who are in braunschweig
        list_person_in_bs = np.array([])
        for index, row in pd.iterrows():
            if((row['countyStart']==3101 and row['ActivityBefore'] == 7) or (row['countyEnd']==3101 and row['ActivityAfter'] == 7)):
                list_person_in_bs = np.append(list_person_in_bs, row['personID'])
        #drop duplicates
        list_person_in_bs = np.unique(list_person_in_bs)
        for index, row in pd.iterrows():
            if row['personID'] in list_person_in_bs:
                pd.at[index, 'home_in_bs'] = 1
    return pd


bd = pd.read_csv(r'C:\\Users\\korf_sa\\Documents\\rep\\data\\mobility\\braunschweig_januar.csv', header=None, skiprows=1)

# read in column headers from first line of csv file
with open(r'C:\\Users\\korf_sa\\Documents\\rep\\data\\mobility\\braunschweig_januar.csv', 'r') as file:
    first_line = file.readline()
    column_headers = first_line.split(',')
    column_headers = [header.strip() for header in column_headers]

#check if anythiong changed since writing this script
if(column_headers[0:19] == ['puid', 'start_zone', 'end_zone', 'loc_id_start', 'loc_id_end', 'start_county', 'end_county', 'huid', 'trip_distance', 'start_time', 'travel_time_sec', 'lon_start', 'lat_start', 'lon_end', 'lat_end', 'travel_mode', 'activity_start', 'activity_end', 'age']):
    print("Column headers are correct")
else:
    print("Column headers are not correct")


# add column headers to DataFrame
pd_new = add_home_is_in_bs_column(bd)
pd_new = add_home_ids(pd_new)


    





