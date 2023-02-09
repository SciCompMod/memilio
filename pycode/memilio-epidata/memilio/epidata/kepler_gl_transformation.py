from pyproj import Proj, Transformer
import pandas as pd
import numpy as np
import os

def transform_coordinates(x, y,w,z,a):
    transformer = Transformer.from_crs("EPSG:3035", "EPSG:4326", always_xy=True)
    return (*transformer.transform(x, y), *transformer.transform(w, z),a) # returns a tuple * expands inner tuples

# read in the data
bd = pd.read_csv('C:/Users/korf_sa/Documents/rep/pycode/memilio-epidata/memilio/epidata/bs.csv', header=None);


bd.rename(columns={0: 'idTrafficZone', 1: 'tripID', 2: 'personID', 3: 'tripChain', 4: 'startZone', 5:'destZone', 6:'loc_id_start', 7:'loc_id_end', 8:'countyStart', 9:'countyEnd', 10:'hhID', 11: 'TripID', 12: 'tripTime', 13: 'startTime', 14: 'travelTime', 19: 'vehicleChoice', 20: 'ActivityBefore', 21: 'ActivityAfter', 15: 'loCs', 16: 'laCs', 17: 'loCe', 18: 'laCe', 22:'age'}, inplace=True);

bd_long_lat = bd[['loCs', 'laCs', 'loCe', 'laCe','vehicleChoice']].head(10000)

bd_long_lat=bd_long_lat.apply(lambda x: transform_coordinates(x['loCs'], x['laCs'],x['loCe'], x['laCe'],x['vehicleChoice'] ), axis=1, result_type='expand')

bd_long_lat.rename(columns={0: 'longitude_start', 1: 'latitude_start', 2: 'longitude_end', 3: 'latitude_end',4:'vehicle_choice'}, inplace=True)

print(bd_long_lat.head(100))
bd_long_lat.to_csv('C:/Users/korf_sa/Documents/rep/pycode/memilio-epidata/memilio/epidata/bs_long_lat.csv', index=False)



