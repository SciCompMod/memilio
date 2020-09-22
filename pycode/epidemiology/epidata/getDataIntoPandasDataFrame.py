import os
import sys
from urllib.request import urlopen
import json
import pandas
import argparse

from epidemiology.epidata import defaultDict as dd

def loadGeojson( targetFileName, apiUrl = 'https://opendata.arcgis.com/datasets/', 
                 extension = 'geojson' ):

    """ Loads data default: ArcGIS data sets in GeoJSON format. (pandas DataFrame)

    This routine loads datasets default: ArcGIS data sets in GeoJSON format of the given public
    data item ID into a pandas DataFrame and returns the DataFrame. Trivial 
    information gets removed by JSON normalization and dropping of always 
    constant data fields.

    Keyword arguments:
    targetFileName -- File name without ending, for ArcGIS: public data item ID (string)
    apiUrl -- URL to file, default: ArcGIS data sets API URL (string, default
              'https://opendata.arcgis.com/datasets/')
    extension -- Data format extension (string, default 'geojson')

    """
    url = apiUrl + targetFileName + '.' + extension

    try:
       with urlopen( url ) as res:
          data = json.loads( res.read().decode() )
    except OSError as e:
        exit_string = "ERROR: URL " + url + " could not be opened."
        sys.exit(exit_string)

    # Shape data:
    df = pandas.json_normalize( data, 'features' )

    # Make-up (removing trivial information):
    df.drop( columns = ['type', 'geometry'], inplace = True )
    df.rename( columns = lambda s: s[11:], inplace = True )

    return df


def loadCsv( targetFileName, apiUrl = 'https://opendata.arcgis.com/datasets/', 
                 extension = 'csv' ):
    """ Loads ArcGIS data sets in CSV format. (pandas DataFrame)

    This routine loads data sets (default from ArcGIS) in CSV format of the given public data 
    item ID into a pandas DataFrame and returns the DataFrame. 

    Keyword arguments:
    targerFileName -- file name which should be downloaded, for ArcGIS it should be public data item ID (string)
    apiUrl -- API URL (string, default
              'https://opendata.arcgis.com/datasets/')
    extension -- Data format extension (string, default 'csv')

    """
    url = apiUrl + targetFileName + '.' + extension

    try:
       df = pandas.read_csv( url )
    except OSError as e:
        exit_string = "ERROR: URL " + url + " could not be opened."
        sys.exit(exit_string)

    return df


def cli(description):
  
   out_path_default = dd.defaultDict['out_folder']
   out_path_default = os.path.join(out_path_default, 'pydata')
   check_dir(out_path_default)

   parser = argparse.ArgumentParser(description=description)
   
   parser.add_argument('-r',  '--read-from-disk',
                       help='Reads the data from file "json" instead of downloading it.',
                       action='store_true')
   parser.add_argument('-p', '--plot', help='Plots the data.',
                       action='store_true')
   parser.add_argument('-h5', '--hdf5', help='Changes output format from json to hdf5.',
                       action='store_true')
   parser.add_argument('-o', '--out-path', type=str, default=out_path_default, help='Defines folder for output.')

   args = parser.parse_args()

   READ_DATA = args.read_from_disk
   MAKE_PLOT = args.plot
   OUT_FORM = "hdf5" if args.hdf5 else dd.defaultDict['out_form']

   return [READ_DATA, MAKE_PLOT, OUT_FORM, args.out_path]

def check_dir(directory):
    # check directory exists or create it
    if not os.path.exists(directory):
        os.makedirs(directory)

def write_dataframe(df, directory, file_prefix, file_type):

   outForm = {
       'json': [".json", {"orient": "records"}],
       'hdf5': [".h5", {"key": "data"}]
   }

   try:
      outFormEnd = outForm[file_type][0]
      outFormSpec = outForm[file_type][1]
   except KeyError:
       exit_string = "Error: The file format: " + file_type + " does not exist. Use another one."
       sys.exit(exit_string)

   if file_type == "json":
       df.to_json(os.path.join(directory, file_prefix + outFormEnd), **outFormSpec)
   elif file_type == "hdf5":
       df.to_hdf(os.path.join(directory, file_prefix + outFormEnd), **outFormSpec)
