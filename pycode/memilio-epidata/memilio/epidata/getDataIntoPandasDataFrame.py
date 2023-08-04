#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Kathrin Rack
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################
"""
@file getDataIntoPandasDataFrame.py
@brief Tools to convert data to pandas dataframes

This tool contains
- load excel format
- load csv format
- organizes the command line interface
- check if directory exists and if not creates it
- writes pandas dataframe to file of three different formats
"""

import os
import argparse
import datetime
import requests
import magic
import urllib3
from io import BytesIO
from zipfile import ZipFile
from warnings import warn

import pandas as pd

from memilio.epidata import defaultDict as dd
from memilio.epidata import progress_indicator


def user_choice(message, default=False):
    while True:
        user_input = input(message + " (y/n): ")[0].lower()
        if user_input == "y":
            return True
        elif user_input == "n":
            return False
        else:
            print("Please answer with y (yes) or n (no)")


def download_file(
        url, chunk_size=1024, timeout=None, progress_function=None,
        verify=True):
    """! Download a file using GET over HTTP.

    @param url Full url of the file to download.
    @param chunk_size Number of Bytes downloaded at once. Only used when a
        progress_function is specified. For a good display of progress, this
        size should be about the speed of your internet connection in Bytes/s.
        Can be set to None to let the server decide the chunk size (may be
        equal to the file size).
    @param timeout Timeout in seconds for the GET request.
    @param progress_function Function called regularly, with the current
        download progress in [0,1] as a float argument.
    @param verify bool or "interactive". If False, ignores the connection's
        security. If True, only starts downloads from secure connections, and
        insecure connections raise a FileNotFoundError. If "interactive",
        prompts the user whether or not to allow insecure connections.
    @return File as BytesIO
    """
    if verify not in [True, False, "interactive"]:
        warn('Invalid input for argument verify. Expected True, False, or'
             ' "interactive", got ' + str(verify) + '.'
             ' Proceeding with "verify=True".', category=RuntimeWarning)
        verify = True
    # send GET request as stream so the content is not downloaded at once
    try:
        req = requests.get(
            url, stream=True, timeout=timeout,
            verify=verify in [True, "interactive"])
    except OSError:
        if verify == "interactive" and user_choice(
            url +
            " could not be opened due to an insecure connection. "
                "Do you want to open it anyways?\n"):
            urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
            req = requests.get(url, stream=True, timeout=timeout, verify=False)
        else:
            raise FileNotFoundError(
                "Error: URL " + url + " could not be opened.")
    if req.status_code != 200:  # e.g. 404
        raise requests.exceptions.HTTPError("HTTPError: "+str(req.status_code))
    # get file size from http header
    # this is only the number of bytes downloaded, the size of the actual file
    # may be larger (e.g. when 'content-encoding' is gzip; decoding is handled
    # by iter_content)
    file_size = int(req.headers.get('content-length'))
    file = bytearray()  # file to be downloaded
    if progress_function:
        progress = 0
        # download file as bytes via iter_content
        for chunk in req.iter_content(chunk_size=chunk_size):
            file += chunk  # append chunk to file
            # note: len(chunk) may be larger (e.g. encoding)
            # or smaller (for the last chunk) than chunk_size
            progress = min(progress+chunk_size, file_size)
            progress_function(progress/file_size)
    else:  # download without tracking progress
        for chunk in req.iter_content(chunk_size=None):
            file += chunk  # append chunk to file
    # return the downloaded content as file like object
    return BytesIO(file)


def extract_zip(file, **param_dict):
    """! reads a zip file and returns a list of dataframes for every file in the zip folder.
    If only one file is readable for func_to_use a single dataframe is returned instead of a list with one entry.

    @param file String. Path to Zipfile to read.
    @param param_dict Dict. Additional information for download functions (e.g. engine, sheet_name, header...)

    @return list od all dataframes (one for each file).
    """
    with ZipFile(file, 'r') as zipObj:
        names = zipObj.namelist()
        all_dfs = [[] for i in range(len(names))]
        for i in range(len(names)):
            with zipObj.open(names[i]) as file2:
                try:
                    all_dfs[i] = pd.read_excel(file2.read(), **param_dict)
                except:
                    pass
    if len(all_dfs) == 1:
        all_dfs = all_dfs[0]
    return all_dfs


def get_file(
        filepath='', url='', read_data=dd.defaultDict['read_data'],
        param_dict={},
        interactive=False):
    """! Loads data from filepath and stores it in a pandas dataframe.
    If data can't be read from given filepath the user is asked whether the file should be downloaded from the given url or not.
    Uses the progress indicator to give feedback.

    @param filepath String. Filepath from where the data is read.
    @param url String. URL to download the dataset.
    @param read_data True or False. Defines if item is opened from directory (True) or downloaded (False).
    @param param_dct Dict. Additional information for download functions (e.g. engine, sheet_name, header...)
    @param interactive bool. Whether to ask for user input. If False, raises Errors instead.

    @return pandas dataframe
    """
    param_dict_excel = {"sheet_name": 0, "header": 0, "engine": 'openpyxl'}
    param_dict_csv = {"sep": ',', "header": 0, "encoding": None, 'dtype': None}
    param_dict_zip = {}

    filetype_dict = {
        'text': pd.read_csv,
        'Composite Document File V2 Document': pd.read_excel,
        'Excel': pd.read_excel, 'Zip': extract_zip}
    param_dict_dict = {
        pd.read_csv: param_dict_csv, pd.read_excel: param_dict_excel,
        extract_zip: param_dict_zip}

    if read_data:
        try:
            df = pd.read_json(filepath)
        except FileNotFoundError:
            if interactive and user_choice(
                "Warning: The file: " + filepath +
                " does not exist in the directory. Do you want to download "
                    "the file from " + url + " instead?\n"):
                df = get_file(filepath=filepath, url=url,
                              read_data=False, param_dict={})
            else:
                error_message = "Error: The file from " + filepath + \
                    " does not exist. Call program without -r flag to get it."
                raise FileNotFoundError(error_message)
    else:
        try:  # to download file from url and show download progress
            with progress_indicator.Percentage(message="Downloading " + url) as p:
                file = download_file(
                    url, 1024, None, p.set_progress,
                    verify="interactive" if interactive else True)
                # read first 2048 bytes to find file type
                ftype = magic.from_buffer(file.read(2048))
                # set pointer back to starting position
                file.seek(0)
                # find file type in dict and use function to read
                func_to_use = [
                    val for key, val in filetype_dict.items()
                    if key in ftype]
                # use different default dict for different functions
                dict_to_use = param_dict_dict[func_to_use[0]]
                # adjust dict
                for k in dict_to_use:
                    if k not in param_dict:
                        param_dict[k] = dict_to_use[k]
                # create dataframe
                df = func_to_use[0](file, **param_dict)
        except OSError:
            raise FileNotFoundError(
                "Error: URL " + url + " could not be opened.")
    try:
        if df.empty:
            raise DataError("Error: Dataframe is empty.")
    except AttributeError:
        if isinstance(df, list) or isinstance(df, dict):
            for i in df:
                if df[i].empty:
                    raise DataError("Error: Dataframe is empty.")
        else:
            raise DataError("Could not catch type of df: " + str(type(df)))
    return df


def cli(what):
    """! Defines command line interface

    The function parameter "what" is used as a dictionary key.
    The return of the dictionary is either a list of a string and a list of keywords.
    The string is the message that should be printed when working on the specific package.
    The further list, contains all special command line arguments which are needed for this package.

    If the key is not part of the dictionary the program is stopped.

    The following default arguments are added to the parser:
    - read-from-disk
    - file-format, choices = ['json', 'hdf5', 'json_timeasstring']
    - out_path
    - no_raw
    - no_progress_indicators (excluded from dict)
    The default values are defined in default dict.

    Depending on what following parser can be added:
    - start_date
    - end_date
    - impute_dates
    - moving_average
    - make_plot
    - split_berlin
    - rep_date
    - sanitize_data

    @param what Defines what packages calls and thus what kind of command line arguments should be defined.
    """

    # TODO: may it would be easier to make a dict like the following one together with a function to get key:
    # TODO: all should automatically do everything
    # cli_dict2 = {"end_date": ['divi'],
    #                "plot": ['cases'],
    #                "start_date": ['divi']                 }

    cli_dict = {"divi": ['Downloads data from DIVI', 'start_date', 'end_date', 'impute_dates', 'moving_average', 'make_plot'],
                "cases": ['Download case data from RKI', 'start_date', 'end_date', 'impute_dates', 'moving_average', 'make_plot', 'split_berlin', 'rep_date'],
                "cases_est": ['Download case data from RKI and JHU and estimate recovered and deaths', 'start_date', 'end_date', 'impute_dates', 'moving_average', 'make_plot', 'split_berlin', 'rep_date'],
                "population": ['Download population data from official sources', 'username'],
                "commuter_official": ['Download commuter data from official sources', 'make_plot'],
                "vaccination": ['Download vaccination data', 'start_date', 'end_date', 'impute_dates', 'moving_average', 'make_plot', 'sanitize_data'],
                "testing": ['Download testing data', 'start_date', 'end_date', 'impute_dates', 'moving_average', 'make_plot'],
                "jh": ['Downloads data from Johns Hopkins University', 'start_date', 'end_date', 'impute_dates', 'moving_average', 'make_plot'],
                "hospitalization": ['Download hospitalization data', 'start_date', 'end_date', 'impute_dates', 'moving_average', 'make_plot'],
                "sim": ['Download all data needed for simulations', 'start_date', 'end_date', 'impute_dates', 'moving_average', 'make_plot', 'split_berlin', 'rep_date', 'sanitize_data']}

    try:
        what_list = cli_dict[what]
    except KeyError:
        raise ValueError("Wrong key or cli_dict.")

    out_path_default = dd.defaultDict['out_folder']

    check_dir(out_path_default)

    parser = argparse.ArgumentParser(description=what_list[0])
    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        '-r', '--read-data',
        help='Reads the data from file "json" instead of downloading it.',
        action='store_true')

    parser.add_argument(
        '-f', '--file-format', type=str, default=dd.defaultDict
        ['file_format'],
        choices=['json', 'hdf5', 'json_timeasstring'],
        help='Defines output format for data files. Default is \"' +
        str(dd.defaultDict['file_format'] + "\"."))
    parser.add_argument('-o', '--out-folder', type=str,
                        default=out_path_default,
                        help='Defines folder for output.')
    parser.add_argument(
        '-n', '--no-raw', default=dd.defaultDict['no_raw'],
        help='Defines if raw data will be stored for further use.',
        action='store_true')

    if 'start_date' in what_list:
        parser.add_argument(
            '-s', '--start-date',
            help='Defines start date for data download. Should have form: YYYY-mm-dd.'
            'Default is 2020-04-24',
            type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d').date(),
            default=dd.defaultDict['start_date'])
    if 'end_date' in what_list:
        parser.add_argument(
            '-e', '--end-date',
            help='Defines date after which data download is stopped.'
            'Should have form: YYYY-mm-dd. Default is today',
            type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d').date(),
            default=dd.defaultDict['end_date'])
    if 'impute_dates' in what_list:
        parser.add_argument(
            '-i', '--impute-dates',
            help='the resulting dfs contain all dates instead of'
            ' omitting dates where no data was reported', action='store_true')
    if 'moving_average' in what_list:
        parser.add_argument(
            '-m', '--moving-average', type=int, default=0,
            help='Compute a moving average of N days over the time series')
    if 'make_plot' in what_list:
        parser.add_argument('-p', '--make-plot', help='Plots the data.',
                            action='store_true')
    if 'split_berlin' in what_list:
        parser.add_argument(
            '-b', '--split-berlin',
            help='Berlin data is split into different counties,'
            ' instead of having only one county for Berlin.',
            action='store_true')
    if 'rep_date' in what_list:
        parser.add_argument(
            '--rep-date', default=False,
            help='If reporting date is activated, the reporting date'
            'will be prefered over possibly given dates of disease onset.',
            action='store_true')
    if 'sanitize_data' in what_list:
        parser.add_argument(
            '-sd', '--sanitize_data', type=int, default=1,
            help='Redistributes cases of every county either based on regions ratios or on thresholds and population'
        )

    parser.add_argument(
        '--no-progress-indicators',
        help='Disables all progress indicators (used for downloads etc.).',
        action='store_true')

    if 'username' in what_list:
        parser.add_argument(
            '--username', type=str
        )

        parser.add_argument(
            '--password', type=str
        )
    args = vars(parser.parse_args())
    # disable progress indicators globally, if the argument --no-progress-indicators was specified
    progress_indicator.ProgressIndicator.disable_indicators(
        args["no_progress_indicators"])
    # remove the no_progress_indicators entry from the dict
    # (after disabling indicators, its value is no longer usefull)
    args.pop("no_progress_indicators")

    return args


def append_filename(
        filename='', impute_dates=False, moving_average=0, split_berlin=False,
        rep_date=False):
    """! Creates consistent file names for all output.
    """
    # split_berlin and repdate especially for case data
    if split_berlin:
        filename = filename + '_split_berlin'
    if rep_date:
        filename = filename + '_repdate'

    if moving_average > 0:
        filename = filename + '_ma' + str(moving_average)
    elif impute_dates:
        filename = filename + '_all_dates'

    return filename


def check_dir(directory):
    """! Checks existence and creates folder

    It is checked if the folder given in the parameter "directory" exists.
    If it does not exist it is created.

    @param directory directory which should exist
    """

    # check if directory exists or create it
    if not os.path.exists(directory):
        os.makedirs(directory)


def write_dataframe(df, directory, file_prefix, file_type, param_dict={}):
    """! Writes pandas dataframe to file

    This routine writes a pandas dataframe to a file in a given format.
    The filename is given without ending.
    A file_type can be
    - json
    - json_timeasstring [Default]
    - hdf5
    - txt
    The file_type defines the file format and thus also the file ending.
    The file format can be json, hdf5 or txt.
    For this option the column Date is converted from datetime to string.

    @param df pandas dataframe (pandas DataFrame)
    @param directory directory where to safe (string)
    @param file_prefix filename without ending (string)
    @param file_type defines ending (string)
    @param param_dict defines parameters for to_csv/txt(dictionary)

    """

    outForm = {'json': [".json", {"orient": "records"}],
               'json_timeasstring': [".json", {"orient": "records"}],
               'hdf5': [".h5", {"key": "data"}],
               'txt': [".txt", param_dict]}

    try:
        outFormEnd = outForm[file_type][0]
        outFormSpec = outForm[file_type][1]
    except KeyError:
        raise ValueError(
            "Error: The file format: " + file_type +
            " does not exist. Use json, json_timeasstring, hdf5 or txt.")

    out_path = os.path.join(directory, file_prefix + outFormEnd)

    if file_type == "json":
        df.to_json(out_path, **outFormSpec)
    elif file_type == "json_timeasstring":
        if dd.EngEng['date'] in df.columns:
            if not isinstance(df.Date.values[0], str):
                df.Date = df.Date.dt.strftime('%Y-%m-%d')
        df.to_json(out_path, **outFormSpec)
    elif file_type == "hdf5":
        df.to_hdf(out_path, **outFormSpec)
    elif file_type == "txt":
        df.to_csv(out_path, **outFormSpec)

    print("Information: Data has been written to", out_path)


class DataError(Exception):
    """ Error for handling incomplete or unexpected Data """
    pass
