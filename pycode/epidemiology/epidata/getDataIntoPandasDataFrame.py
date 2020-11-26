import os
import sys
from urllib.request import urlopen
import json
import pandas
import argparse
import datetime

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
                 extension = '.csv' ):
    """ Loads ArcGIS data sets in CSV format. (pandas DataFrame)

    This routine loads data sets (default from ArcGIS) in CSV format of the given public data 
    item ID into a pandas DataFrame and returns the DataFrame. 

    Keyword arguments:
    targerFileName -- file name which should be downloaded, for ArcGIS it should be public data item ID (string)
    apiUrl -- API URL (string, default
              'https://opendata.arcgis.com/datasets/')
    extension -- Data format extension (string, default 'csv')

    """
    url = apiUrl + targetFileName  + extension

    try:
       df = pandas.read_csv( url )
    except OSError as e:
        exit_string = "ERROR: URL " + url + " could not be opened."
        sys.exit(exit_string)

    return df


# function to return list of keys for any value
# def get_key(val, my_dict):
#    key_list = []
#    for key, value_list in my_dict.items():
#        if val in value_list:
#            key_list.append(key)

#    return key_list

def cli(what):

   # TODO: may it would be easier to make a dict like the following one together with a function to get key:
   # TODO: all should automatically do everything
   # cli_dict2 = {"end_date": ['divi'],
   #                "plot": ['rki'],
   #                "start_date": ['divi'],
   #                "update": ['divi']                 }

   cli_dict = {"divi": ['Downloads data from DIVI', 'start_date', 'end_date', 'update'],
               "rki": ['Download data from RKI', 'plot', 'split_berlin'],
               "rkiest": ['Download data from RKI and JH and estimate recovered and deaths', 'plot'],
               "spain": ['Download of spain data'],
               "population": ['Download population data'],
               "jh" : ['Downloads data from JH'],
               "all": ['Download all possible data', 'plot','start_date', 'end_date', 'update']}

   try:
      what_list = cli_dict[what]
   except KeyError:
      exit_string = "Wrong key or cli_dict."
      sys.exit(exit_string)

   out_path_default = dd.defaultDict['out_folder']
   out_path_default = os.path.join(out_path_default, 'pydata')

   check_dir(out_path_default)

   parser = argparse.ArgumentParser(description=what_list[0])
   
   parser.add_argument('-r',  '--read-from-disk',
                       help='Reads the data from file "json" instead of downloading it.',
                       action='store_true')
   parser.add_argument('-ff', '--file-format', type=str, default=dd.defaultDict['out_form'],
                       choices=['json', 'hdf5', 'json_timeasstring'],
                       help='Defines output format for data files. Default is \"' + str(dd.defaultDict['out_form']+ "\"."))
   parser.add_argument('-o', '--out-path', type=str, default=out_path_default, help='Defines folder for output.')

   if 'split_berlin' in what_list:
       parser.add_argument('-cb', '--split_berlin',
                           help='Berlin data is split into different counties,'
                                ' instead of having only one county for Berlin.',
                           action='store_true')
   if 'end_date' in what_list:
       parser.add_argument('-ed', '--end-date',
                           help='Defines date after which data download is stopped.'
                                'Should have form: YYYY-mm-dd. Default is today',
                           type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d').date(),
                           default=dd.defaultDict['end_date'])
   if 'plot' in what_list:
      parser.add_argument('-p', '--plot', help='Plots the data.',
                          action='store_true')
   if 'start_date' in what_list:
      parser.add_argument('-sd',  '--start-date',
                          help='Defines start date for data download. Should have form: YYYY-mm-dd.'
                               'Default is 2020-04-24',
                          type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d').date(),
                          default = dd.defaultDict['start_date'])

   if 'update' in what_list:
      parser.add_argument('-u',  '--update',
                          help='Reads the data from file "json", downloads and adds data from today.',
                          action='store_true')

   args = parser.parse_args()

   arg_list = []

   READ_DATA = args.read_from_disk
   arg_list.append(READ_DATA)
   arg_list.append(args.file_format)
   arg_list.append(args.out_path)

   # add additional arguments in alphabetical order
   # TODO: check if it is possible to automatically generate this
   if 'split_berlin' in what_list:
       arg_list.append(args.split_berlin)
   if 'end_date' in what_list:
       arg_list.append(args.end_date)
   if 'plot' in what_list:
      arg_list.append(args.plot)
   if 'start_date' in what_list:
      arg_list.append(args.start_date)
   if 'update' in what_list:
      UPDATE_DATA = args.update
      arg_list.append(UPDATE_DATA)

      # TODO: Change arguments such that one argument + parameter can be either read_data or update
      if READ_DATA:
         exit_string = "You called the program with '--read-from-disk' and '--update'." \
                       "Please choose just one. Both together is not possible."
         sys.exit(exit_string)

   return arg_list

def check_dir(directory):
    # check if directory exists or create it
    if not os.path.exists(directory):
        os.makedirs(directory)

def write_dataframe(df, directory, file_prefix, file_type):

   outForm = {
       'json': [".json", {"orient": "records"}],
       'json_timeasstring': [".json", {"orient": "records"}],
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
   elif file_type == "json_timeasstring":
       if dd.EngEng['date'] in df.columns:
            if type(df.Date.values[0]) != type("string"):
                 df.Date = df.Date.dt.strftime('%Y-%m-%d')
       df.to_json(os.path.join(directory, file_prefix + outFormEnd), **outFormSpec)
   elif file_type == "hdf5":
       df.to_hdf(os.path.join(directory, file_prefix + outFormEnd), **outFormSpec)
