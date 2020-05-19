import sys
from urllib.request import urlopen
import json
import pandas
import argparse


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

def cli(description):
  
   out_path_default = ""
 
   parser = argparse.ArgumentParser(description=description)
   
   parser.add_argument('-r',  '--read-from-disk',
                       help='Reads the data from file "json" instead of downloading it.',
                       action='store_true')
   parser.add_argument('-p', '--plot', help='Plots the data.',
                       action='store_true')
   parser.add_argument('-h5', '--hdf5', help='Changes output format from json to hdf5.',
                       action='store_true')
   parser.add_argument('-o', '--out_path', type=str, default=out_path_default, help='Defines folder for output.')
                       #action='store_true')

   args = parser.parse_args()

   READ_DATA = args.read_from_disk
   MAKE_PLOT = args.plot
   OUT_FORM = "hdf5" if args.hdf5 else "json"

   
   return [READ_DATA, MAKE_PLOT, OUT_FORM, args.out_path]

