import pandas as pd
import numpy as np
import os
import memilio.epidata.geoModificationGermany as gMG
import matplotlib.pyplot as plt
import memilio.epidata.getPopulationData as gPD

print(os.listdir())
county22 = pd.read_csv('pycode/memilio-epidata/memilio/epidata/trip_count_matrix_kreis_by_age_activity_2022.csv', header=None)
# ask for data types
print(county22.dtypes)
# rename columns
county22.rename(columns={0: 'idCountyFrom', 1: 'idCountyTo', 2: 'ageGroup', 3: 'activity', 4: 'commuters'}, inplace=True)

for col in ['idCountyFrom', 'idCountyTo']:
    print(county22[col].nunique() == 400)
    print(1001 in county22[col])
    print((16056 in county22[col].values) == False) # Eisenach shouldnt be there...
    print(len(county22[county22[col]==16056]) == 0)
    print(16063 in county22[col].values) # Wartburgkreis

colFrom = 'idCountyFrom'
countyIDs = gMG.get_county_ids()
print((np.sort(county22[colFrom].unique()) == countyIDs).all())
print(county22['commuters'].min() >= 0)

for id in countyIDs:
    if not (county22[county22[colFrom]==id]['ageGroup'].nunique() == 5):
        print('Error in county ' + str(id))
    if not (county22[county22[colFrom]==id]['activity'].nunique()==7):
        print('Error in county ' + str(id))

    if not (np.sort(county22[county22[colFrom]==id]['ageGroup'].unique()) == [1,2,3,4,5]).all():
        print('Error in county ' + str(id))
    if not (np.sort(county22[county22[colFrom]==id]['activity'].unique()) == [1,2,3,4,5,6,7]).all():
        print('Error in county ' + str(id))



# county = pd.read_csv('pycode/memilio-epidata/memilio/epidata/trip_count_matrix_kreis_by_age_activity.csv', header=None)
# # ask for data types
# print(county.dtypes)
# # rename columns
# county.rename(columns={0: 'idCountyFrom', 1: 'idCountyTo', 2: 'ageGroup', 3: 'activity', 4: 'commuters'}, inplace=True)

# for col in ['idCountyFrom', 'idCountyTo']:
#     print(county[col].nunique() == 401)
#     print(1001 in county[col])
#     print((16056 in county[col].values) == True) # Eisenach should be there...
#     print(len(county[county[col]==16056]) > 0)
#     print(16063 in county[col].values) # Wartburgkreis

# colTo = 'idCountyFrom'
# countyIDs = gMG.get_county_ids(merge_eisenach=False)
# print((np.sort(county[colTo].unique()) == countyIDs).all())
# print(county['commuters'].min() >= 0)

# for id in countyIDs:
#     if not (county[county[colTo]==id]['ageGroup'].nunique() == 5):
#         print('Error in county ' + str(id))
#     if not (county[county[colTo]==id]['activity'].nunique()==7):
#         print('Error in county ' + str(id))

#     if not (np.sort(county[county[colTo]==id]['ageGroup'].unique()) == [1,2,3,4,5]).all():
#         print('Error in county ' + str(id))
#     if not (np.sort(county[county[colTo]==id]['activity'].unique()) == [1,2,3,4,5,6,7]).all():
#         print('Error in county ' + str(id))


states = gMG.get_state_ids()
states_to_counties = gMG.get_stateid_to_countyids_map()

countyIDtoIdx = {i : countyIDs.index(i) for i in countyIDs}

mat_commuters_total = np.zeros((400, 400))

population = gPD.get_population_data()[['ID_County', 'Population']]
population['id_idx'] = population['ID_County'].map(countyIDtoIdx)

county22['idxFrom'] = county22['idCountyFrom'].map(countyIDtoIdx)
county22['idxTo'] = county22['idCountyTo'].map(countyIDtoIdx)

mobility_agg = county22.groupby(['idxFrom', 'idxTo']).agg({'commuters': sum}).reset_index()

# total
perc_trips_inside = np.zeros(400)
perc_trips_outside = np.zeros(400)

for i in range(population['id_idx'].min(), population['id_idx'].max()+1):

    perc_trips_inside[i] = mobility_agg[(mobility_agg['idxFrom']==i) & (mobility_agg['idxTo']==i)]['commuters'].sum()/population[population['id_idx']==i]['Population'].sum()

    a=[k for k in range(0,400)]
    a.remove(i)

    perc_trips_outside[i] = mobility_agg[(mobility_agg['idxFrom']==i) & (mobility_agg['idxTo']).isin(a)]['commuters'].sum()/population[population['id_idx']==i]['Population'].sum()

plt.plot([k for k in range(0,400)], perc_trips_inside, label='inside')

plt.plot([k for k in range(0,400)], perc_trips_outside, label='outside')
plt.title('Trips in-/outside county relative to county population')
plt.legend()

# work
county22_work = county22.copy()
county22_work = county22_work[county22_work['activity']==1]

mobility_agg = county22_work.groupby(['idxFrom', 'idxTo']).agg({'commuters': sum}).reset_index()

perc_trips_inside = np.zeros(400)
perc_trips_outside = np.zeros(400)

for i in range(population['id_idx'].min(), population['id_idx'].max()+1):

    perc_trips_inside[i] = mobility_agg[(mobility_agg['idxFrom']==i) & (mobility_agg['idxTo']==i)]['commuters'].sum()/population[population['id_idx']==i]['Population'].sum()

    a=[k for k in range(0,400)]
    a.remove(i)

    perc_trips_outside[i] = mobility_agg[(mobility_agg['idxFrom']==i) & (mobility_agg['idxTo']).isin(a)]['commuters'].sum()/population[population['id_idx']==i]['Population'].sum()

plt.plot([k for k in range(0,400)], perc_trips_inside, label='inside')

plt.plot([k for k in range(0,400)], perc_trips_outside, label='outside')
plt.title('Trips in-/outside county relative to county population')
plt.legend()

for i in states:
    counties = states_to_counties[i]
    if i == 1:
        counties.append(states_to_counties[2])

    mat_commuters_fedstate = county22[(county22['idCountyFrom'].isin(counties)) & (county22['idCountyTo'].isin(counties))][['idxFrom', 'idxTo', 'commuters']]

    im = plt.imshow(county22)
    




x=15