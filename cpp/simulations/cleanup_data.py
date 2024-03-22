import os
import sys
import numpy as np
import pandas as pd


def add_external_ids(pd):
    # Schools are anonymized with loc_id_start and loc_id_end = -2
    # for every row where loc_id_start or loc_id_end is -2 we assume a school in the zone and change the id to something unique (somewhere  betwenn-100 to -10000)

    # list of external locations
    work_places_dict = {}  # map zone to work place id and size
    person_to_work_place_dict = {}  # maps persons to work places to track persons
    # start of work place ids, decrement by 1 for each new work place
    work_place_ids = -10100
    max_persons_per_work_place = 40

    shops_dict = {}
    person_to_shop_dict = {}
    shop_ids = -20100
    max_persons_per_shop = 200

    social_event_dict = {}
    person_to_social_event_dict = {}
    social_event_ids = -30100
    max_persons_per_social_event = 100

    # we go through every row and check if the loc_id_start or loc_id_end is -3
    # if so we change the loc_id_start or loc_id_end to a unique id which represents an activity in that zone
    # if in one unique location are more than X (see above) persons we change the id to a new unique id which represents a new location in that zone

    start_or_end_loc = ''
    start_or_end_zone = ''
    for index, row in pd.iterrows():
        if ((row['loc_id_start'] == -3) | (row['loc_id_end'] == -3)):
            if row['loc_id_start'] == -3:
                start_or_end_loc = 'loc_id_start'
                start_or_end_zone = 'start_zone'
                start_or_end_activity = 'activity_start'
            else:
                start_or_end_loc = 'loc_id_end'
                start_or_end_zone = 'end_zone'
                start_or_end_activity = 'activity_end'

            found = False
            if (row[start_or_end_activity] == 1):  # work places
                if (not row['puid'] in person_to_work_place_dict.keys()):
                    # check if a work place is already available and not full
                    if (row[start_or_end_zone] in work_places_dict.keys()):
                        if (work_places_dict[row[start_or_end_zone]][1] < max_persons_per_work_place):
                            # add persons to existing work place
                            person_to_work_place_dict[row['puid']
                                                      ] = work_places_dict[row[start_or_end_zone]][0]
                            work_places_dict[row[start_or_end_zone]][1] += 1
                            found = True
                    if (not found):
                        # make new work place
                        work_places_dict[row[start_or_end_zone]] = [
                            work_place_ids, 1]
                        work_place_ids -= 1
                        person_to_work_place_dict[row['puid']
                                                  ] = work_places_dict[row[start_or_end_zone]][0]
                # overwrite "-3" in location_id
                pd.at[index, start_or_end_loc] = person_to_work_place_dict[row['puid']]

            # shops and private matters
            elif (row[start_or_end_activity] == 3 | row[start_or_end_activity] == 5):
                if (not row['puid'] in person_to_shop_dict.keys()):
                    # check if a shop is already available and not full
                    if (row[start_or_end_zone] in shops_dict.keys()):
                        if (shops_dict[row[start_or_end_zone]][1] < max_persons_per_shop):
                            # add persons to existing shop
                            person_to_shop_dict[row['puid']
                                                ] = shops_dict[row[start_or_end_zone]][0]
                            shops_dict[row[start_or_end_zone]][1] += 1
                            found = True
                    if (not found):
                        # make new shop
                        shops_dict[row[start_or_end_zone]] = [
                            shop_ids, 1]
                        shop_ids -= 1
                        person_to_shop_dict[row['puid']
                                            ] = shops_dict[row[start_or_end_zone]][0]
                # overwrite "-3" in location_id
                pd.at[index, start_or_end_loc] = person_to_shop_dict[row['puid']]

            # social events and 'other'
            elif (row[start_or_end_activity] == 4 | row[start_or_end_activity] == 6):
                if (not row['puid'] in person_to_social_event_dict.keys()):
                    # check if a social event is already available and not full
                    if (row[start_or_end_zone] in social_event_dict.keys()):
                        if (social_event_dict[row[start_or_end_zone]][1] < max_persons_per_social_event):
                            # add persons to existing social event
                            person_to_social_event_dict[row['puid']
                                                        ] = social_event_dict[row[start_or_end_zone]][0]
                            social_event_dict[row[start_or_end_zone]][1] += 1
                            found = True
                    if (not found):
                        # make new social event
                        social_event_dict[row[start_or_end_zone]] = [
                            social_event_ids, 1]
                        social_event_ids -= 1
                        person_to_social_event_dict[row['puid']
                                                    ] = social_event_dict[row[start_or_end_zone]][0]
                # overwrite "-3" in location_id
                pd.at[index, start_or_end_loc] = person_to_social_event_dict[row['puid']]

    return pd


def add_school_ids(pd):
    # Schools are anonymized with loc_id_start and loc_id_end = -2
    # for every row where loc_id_start or loc_id_end is -2 we assume a school in the zone and change the id to something unique (somewhere  betwenn-100 to -10000)

    # list of schools
    schools_dict = {}  # map zone to school id and size
    person_to_school_dict = {}  # maps persons to schools to track persons
    school_ids = -100  # start of school ids, decrement by 1 for each new school
    max_persons_per_school = 600

    # we go through every row and check if the loc_id_start or loc_id_end is -2
    # if so we change the loc_id_start or loc_id_end to a unique id which represents a school in that zone
    # if in one unique school are more than 600 persons we change the id to a new unique id which represents a new school in that zone

    start_or_end_loc = ''
    start_or_end_zone = ''
    for index, row in pd.iterrows():
        if (((row['loc_id_start'] == -2) | (row['loc_id_end'] == -2)) & (not row['puid'] in person_to_school_dict.keys())):
            if row['loc_id_start'] == -2:
                start_or_end_loc = 'loc_id_start'
                start_or_end_zone = 'start_zone'
            else:
                start_or_end_loc = 'loc_id_end'
                start_or_end_zone = 'end_zone'

            found = False
            # check if a school is already available and not full
            if (row[start_or_end_zone] in schools_dict.keys()):
                if (schools_dict[row[start_or_end_zone]][1] < max_persons_per_school):
                    # add persons to existing school
                    person_to_school_dict[row['puid']
                                          ] = schools_dict[row[start_or_end_zone]][0]
                    schools_dict[row[start_or_end_zone]][1] += 1
                    found = True
            if (not found):
                # make new school
                schools_dict[row[start_or_end_zone]] = [school_ids, 1]
                school_ids -= 1
                person_to_school_dict[row['puid']
                                      ] = schools_dict[row[start_or_end_zone]][0]
        # overwrite "-2" in location_id
        pd.at[index, start_or_end_loc] = person_to_school_dict[row['puid']]

    return pd


def add_home_ids(pd):
    # when there is a -1 in loc_id_start or loc_id_end, the person is at home (due to anonymization we dont know exactly if the person goes to their home or another we assume their own home)
    # every row where loc_id_start or loc_id_end is -1 and activity_start or activity_end is 7 is a person who goes /comes from home
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
            if ((row['start_county'] == 3101 and row['activity_start'] == 7) or (row['end_county'] == 3101 and row['activity_end'] == 7)):
                list_person_in_bs = np.append(list_person_in_bs, row['puid'])
        #drop duplicates
        list_person_in_bs = np.unique(list_person_in_bs)
        for index, row in pd.iterrows():
            if row['puid'] in list_person_in_bs:
                pd.at[index, 'home_in_bs'] = 1
    return pd


PATH = "/Users/david/Documents/HZI/memilio/data/mobility/braunschweig_result.csv"

bd = pd.read_csv(PATH)

# read in column headers from first line of csv file
with open(PATH, 'r') as file:
    first_line = file.readline()
    column_headers = first_line.split(',')
    column_headers = [header.strip() for header in column_headers]

# check if anythiong changed since writing this script
if(column_headers[0:19] == ['puid', 'start_zone', 'end_zone', 'loc_id_start', 'loc_id_end', 'start_county', 'end_county', 'huid', 'trip_distance', 'start_time', 'travel_time_sec', 'lon_start', 'lat_start', 'lon_end', 'lat_end', 'travel_mode', 'activity_start', 'activity_end', 'age']):
    print("Column headers are correct")
else:
    print("Column headers are not correct")

bd_new = add_home_is_in_bs_column(bd)
print('Added flag for persons who live in Braunschweig.')
bd_new = add_home_ids(bd_new)
print('Home IDs set.')
bd_new = add_school_ids(bd_new)
print('School IDs set.')
bd_new = add_external_ids(bd_new)
print('External location IDs set.')

# Write data back to disk
bd_new.to_csv('modified_braunschweig_result.csv', index=False)
