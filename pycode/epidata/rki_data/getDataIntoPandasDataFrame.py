import sys
from urllib.request import urlopen
import json
import pandas
import matplotlib.pyplot as plt

def loadGeojson( itemId, apiUrl = 'https://opendata.arcgis.com/datasets/', 
                 extension = 'geojson' ):
    """ Loads ArcGIS data sets in GeoJSON format. (pandas DataFrame)

    This routine loads ArcGIS data sets in GeoJSON format of the given public 
    data item ID into a pandas DataFrame and returns the DataFrame. Trivial 
    information gets removed by JSON normalization and dropping of always 
    constant data fields.

    Keyword arguments:
    itemId -- ArcGIS public data item ID (string)
    apiUrl -- ArcGIS data sets API URL (string, default
              'https://opendata.arcgis.com/datasets/')
    extension -- Data format extension (string, default 'geojson')

    """
    url = apiUrl + itemId + '_0.' + extension

    with urlopen( url ) as res:
        data = json.loads( res.read().decode() )

    # Shape data:
    df = pandas.json_normalize( data, 'features' )

    # Make-up (removing trivial information):
    df.drop( columns = ['type', 'geometry'], inplace = True )
    df.rename( columns = lambda s: s[11:], inplace = True )

    return df

def loadCsv( itemId, apiUrl = 'https://opendata.arcgis.com/datasets/', 
                 extension = 'csv' ):
    """ Loads ArcGIS data sets in CSV format. (pandas DataFrame)

    This routine loads ArcGIS data sets in CSV format of the given public data 
    item ID into a pandas DataFrame and returns the DataFrame. 

    Keyword arguments:
    itemId -- ArcGIS public data item ID (string)
    apiUrl -- ArcGIS data sets API URL (string, default
              'https://opendata.arcgis.com/datasets/')
    extension -- Data format extension (string, default 'csv')

    """
    url = apiUrl + itemId + '.' + extension

    df = pandas.read_csv( url )

    return df

