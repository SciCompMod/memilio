"""
@file getDataIntoPandasDataFrame.py
@brief Tools to convert data to pandas dataframes

This tool contains
- load geojson format
- load csv format
- organizes the command line interface
- check if directory exists and if not creates it
- writes pandas dataframe to file of three different formats
"""

import os
import sys
from urllib.request import urlopen
import json
import argparse
import datetime
import pandas

from epidemiology.epidata import defaultDict as dd


def loadGeojson(targetFileName, apiUrl='https://opendata.arcgis.com/datasets/', extension='geojson'):
    """! Loads data default: ArcGIS data sets in GeoJSON format. (pandas DataFrame)
    This routine loads datasets default: ArcGIS data sets in GeoJSON format of the given public
    data item ID into a pandas DataFrame and returns the DataFrame. Trivial
    information gets removed by JSON normalization and dropping of always
    constant data fields.

    @param targetFileName -- File name without ending, for ArcGIS use public data item ID (string)
    @param apiUrl -- URL to file (default: ArcGIS data sets API URL (string, default
              'https://opendata.arcgis.com/datasets/'))
    @param extension -- Data format extension (string, default 'geojson')
    return dataframe
    """
    url = apiUrl + targetFileName + '.' + extension

    try:
        with urlopen(url) as res:
            data = json.loads(res.read().decode())
    except OSError:
        exit_string = "ERROR: URL " + url + " could not be opened."
        sys.exit(exit_string)

    # Shape data:
    df = pandas.json_normalize(data, 'features')

    # Make-up (removing trivial information):
    df.drop(columns=['type', 'geometry'], inplace=True)
    df.rename(columns=lambda s: s[11:], inplace=True)

    return df


def loadCsv(targetFileName, apiUrl='https://opendata.arcgis.com/datasets/', extension='.csv'):
    """! Loads data sets in CSV format. (pandas DataFrame)
    This routine loads data sets (default from ArcGIS) in CSV format of the given public data
    item ID into a pandas DataFrame and returns the DataFrame.

    @param targetFileName -- file name which should be downloaded, for ArcGIS it should be public data item ID (string)
    @param apiUrl -- API URL (string, default
              'https://opendata.arcgis.com/datasets/')
    @param extension -- Data format extension (string, default 'csv')
    return dataframe
    """

    url = apiUrl + targetFileName + extension

    try:
        df = pandas.read_csv(url)
    except OSError:
        exit_string = "ERROR: URL " + url + " could not be opened."
        sys.exit(exit_string)

    return df


def loadExcel(targetFileName, apiUrl='https://opendata.arcgis.com/datasets/',
              extension='.xls', sheet_name=0, header=0, engine='openpyxl'):
    """ Loads ArcGIS data sets in Excel formats (xls,xlsx,xlsm,xlsb,odf,ods,odt). (pandas DataFrame)

    This routine loads data sets (default from ArcGIS) in Excel format of the given public data
    item ID into a pandas DataFrame and returns the DataFrame.

    Keyword arguments:
    targetFileName -- file name which should be downloaded, for ArcGIS it should be public data item ID (string)
    apiUrl -- API URL (string, default
              'https://opendata.arcgis.com/datasets/')
    extension -- Data format extension (string, default '.xls')
    sheet -- sheet from Excel file which should be in DataFrame
            (string (sheetname) or integer (zero-indexed sheet position), default 0)
    header -- row to use for column labels (Use None if there is no header) (int, default 0)
    """
    url = apiUrl + targetFileName + extension

    try:
        df = pandas.read_excel(url, sheet_name=sheet_name, header=header, engine=engine)
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
    """! Defines command line interface

    The function parameter "what" is used as a dictionary key.
    The return of the dictionary is either a list of a string and a list of keywords.
    The string is the message that should be printed when working on the specific package.
    The further list, contains all special command line arguments which are needed for this package.

    If the key is nor part of the dictionary the program is stopped.

    Three default arguments are added to the parser:
    - read-from-disk, Default = False
    - file-format, Default = json_timeasstring, choices = ['json', 'hdf5', 'json_timeasstring']
    - out_path Default = data/pydata/

    Depending on what following parser can be added:
    - end_date
    - plot
    - split_berlin
    - moving-average
    - start_date
    - update

    @param what Defines what packages calls and thus what kind of command line arguments should be defined.
    """

    # TODO: may it would be easier to make a dict like the following one together with a function to get key:
    # TODO: all should automatically do everything
    # cli_dict2 = {"end_date": ['divi'],
    #                "plot": ['rki'],
    #                "start_date": ['divi'],
    #                "update": ['divi']                 }

    cli_dict = {"divi": ['Downloads data from DIVI', 'start_date', 'end_date', 'update'],
                "rki": ['Download data from RKI', 'fill_dates', 'make_plot', 'moving_average', 'split_berlin'],
                "rkiest": ['Download data from RKI and JH and estimate recovered and deaths', 'make_plot'],
                "spain": ['Download of spain data'],
                "population": ['Download population data'],
                "vaccine": ['Download vaccine data'],
                "jh" : ['Downloads data from JH'],
                "sim": ['Download all data needed for simulations', 'start_date', 'end_date', 'update',
                        'fill_dates', 'make_plot', 'moving_average', 'split_berlin']}

    try:
        what_list = cli_dict[what]
    except KeyError:
        exit_string = "Wrong key or cli_dict."
        sys.exit(exit_string)

    out_path_default = dd.defaultDict['out_folder']
    out_path_default = os.path.join(out_path_default, 'pydata')

    check_dir(out_path_default)

    parser = argparse.ArgumentParser(description=what_list[0])

    parser.add_argument('-r', '--read-from-disk',
                        help='Reads the data from file "json" instead of downloading it.',
                        action='store_true')
    parser.add_argument('-f', '--file-format', type=str, default=dd.defaultDict['out_form'],
                        choices=['json', 'hdf5', 'json_timeasstring'],
                        help='Defines output format for data files. Default is \"' + str(
                            dd.defaultDict['out_form'] + "\"."))
    parser.add_argument('-o', '--out-path', type=str, default=out_path_default, help='Defines folder for output.')
    parser.add_argument('-n', '--no-raw', default=dd.defaultDict['no_raw'],
                        help='Defines if raw data fill be stored for further use.',
                        action='store_true')

    if 'end_date' in what_list:
        parser.add_argument('-e', '--end-date',
                            help='Defines date after which data download is stopped.'
                                 'Should have form: YYYY-mm-dd. Default is today',
                            type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d').date(),
                            default=dd.defaultDict['end_date'])
    if 'fill_dates' in what_list:
        parser.add_argument('-d', '--fill-dates',
                            help='the resulting dfs contain all dates instead of'
                                 ' omitting dates where no new cases were reported',
                            action='store_true')
    if 'make_plot' in what_list:
        parser.add_argument('-p', '--plot', help='Plots the data.',
                            action='store_true')
    if 'moving_average' in what_list:
        parser.add_argument('-m', '--moving-average',
                            help='The moving average is computed instead of the real values',
                            action='store_true')
    if 'split_berlin' in what_list:
        parser.add_argument('-b', '--split-berlin',
                            help='Berlin data is split into different counties,'
                                 ' instead of having only one county for Berlin.',
                            action='store_true')
    if 'start_date' in what_list:
        parser.add_argument('-s', '--start-date',
                            help='Defines start date for data download. Should have form: YYYY-mm-dd.'
                                 'Default is 2020-04-24',
                            type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d').date(),
                            default=dd.defaultDict['start_date'])
    if 'update' in what_list:
        parser.add_argument('-u', '--update',
                            help='Reads the data from file "json", downloads and adds data from today.',
                            action='store_true')

    args = parser.parse_args()

    arg_list = []

    read_data = args.read_from_disk
    arg_list.append(read_data)
    arg_list.append(args.file_format)
    arg_list.append(args.out_path)
    arg_list.append(args.no_raw)

    # add additional arguments in alphabetical order
    # TODO: check if it is possible to automatically generate this
    if 'end_date' in what_list:
        arg_list.append(args.end_date)
    if 'fill_dates' in what_list:
        arg_list.append(args.fill_dates)
    if 'make_plot' in what_list:
        arg_list.append(args.plot)
    if 'moving_average' in what_list:
        arg_list.append(args.moving_average)
    if 'split_berlin' in what_list:
        arg_list.append(args.split_berlin)
    if 'start_date' in what_list:
        arg_list.append(args.start_date)
    if 'update' in what_list:
        update_data = args.update
        arg_list.append(update_data)

        # TODO: Change arguments such that one argument + parameter can be either read_data or update
        if read_data and update_data:
            exit_string = "You called the program with '--read-from-disk' and '--update'." \
                          "Please choose just one. Both together is not possible."
            sys.exit(exit_string)

    return arg_list


def check_dir(directory):
    """! Checks existence and creates folder

    It is checked if the folder given in the parameter "directory" exists.
    If it does not exist it is created.

    @param directory directory which should exist
    """

    # check if directory exists or create it
    if not os.path.exists(directory):
        os.makedirs(directory)


def write_dataframe(df, directory, file_prefix, file_type):
    """! Writes pandas dataframe to file

    This routine writes a pandas dataframe to a file in a given format.
    The filename is given without ending.
    A file_type can be
    - json
    - json_timeasstring [Default]
    - hdf5
    The file_type defines the file format and thus also the file ending.
    The file format can be json or hdf5.
    For this option the column Date is converted from datetime to string.

    @param df pandas dataframe (pandas DataFrame)
    @param directory directory where to safe (string)
    @param file_prefix filename without ending (string)
    @param file_type defines ending (string)

    """

    outForm = {'json': [".json", {"orient": "records"}],
               'json_timeasstring': [".json", {"orient": "records"}],
               'hdf5': [".h5", {"key": "data"}]}

    try:
        outFormEnd = outForm[file_type][0]
        outFormSpec = outForm[file_type][1]
    except KeyError:
        exit_string = "Error: The file format: " + file_type + " does not exist. Use another one."
        sys.exit(exit_string)

    out_path = os.path.join(directory, file_prefix + outFormEnd)

    if file_type == "json":
        df.to_json(out_path, **outFormSpec)
    elif file_type == "json_timeasstring":
        if dd.EngEng['date'] in df.columns:
            if not isinstance(df.Date.values[0], type("string")):
                df.Date = df.Date.dt.strftime('%Y-%m-%d')
        df.to_json(out_path, **outFormSpec)
    elif file_type == "hdf5":
        df.to_hdf(out_path, **outFormSpec)

    print("Information: Data has been written to", out_path)
