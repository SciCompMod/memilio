from pyproj import Proj, Transformer
import pandas as pd
import numpy as np
import os


def transform_coordinates(x, y, w, z, a):
    transformer = Transformer.from_crs(
        "EPSG:3035", "EPSG:4326", always_xy=True)
    # returns a tuple * expands inner tuples
    return (*transformer.transform(x, y), *transformer.transform(w, z), a)


# read in the data
bd = pd.read_csv(r'~/Documents/HZI/memilio/data/mobility/bs.csv', header=None)
#number_of_rows_to_transform = 100

bd.rename(
    columns={0: 'idTrafficZone', 1: 'tripID', 2: 'personID', 3: 'tripChain', 4: 'startZone', 5: 'destZone', 6: 'loc_id_start', 7: 'loc_id_end',
             8: 'countyStart', 9: 'countyEnd', 10: 'hhID', 11: 'TripID', 12: 'tripTime', 13: 'startTime', 14: 'travelTime', 19: 'vehicleChoice', 20:
             'ActivityBefore', 21: 'ActivityAfter', 15: 'loCs', 16: 'laCs', 17: 'loCe', 18: 'laCe', 22: 'age'},
    inplace=True)

# bd_long_lat = bd[['loCs', 'laCs', 'loCe', 'laCe', 'vehicleChoice']].head(number_of_rows_to_transform)
bd_long_lat = bd[['loCs', 'laCs', 'loCe', 'laCe']]
# aggregate the same trips
bd_long_lat = bd_long_lat.groupby(['loCs', 'laCs', 'loCe', 'laCe']).size(
).reset_index(name='counts').sort_values('counts', ascending=False)
# take out the trips with less than 11 occurences
bd_long_lat = bd_long_lat[bd_long_lat['counts'] > 10]

bd_long_lat = bd_long_lat.apply(lambda x: transform_coordinates(
    x['loCs'], x['laCs'], x['loCe'], x['laCe'], x['counts']), axis=1, result_type='expand')

bd_long_lat.rename(columns={0: 'longitude_start', 1: 'latitude_start',
                   2: 'longitude_end', 3: 'latitude_end', 4: 'counts'}, inplace=True)

print(bd_long_lat.head(100))
bd_long_lat.to_csv(
    r'~/Documents/HZI/memilio/data/mobility/bs_lola.csv', index=False)
