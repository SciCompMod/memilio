#############################################################################
# Copyright (C) 2020-2026 MEmilio
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
:strong:`getDataIntoPandasDataFrame.py`
Tools to convert data to pandas dataframes

This tool contains

- load excel format
- load csv format
- organizes the command line interface
- check if directory exists and if not creates it
- writes pandas dataframe to file of three different formats
"""
import sys
import os
import argparse
import configparser
import datetime
import requests
import magic
import urllib3
import warnings
import matplotlib
from io import BytesIO
from zipfile import ZipFile
from enum import Enum
from pkg_resources import parse_version

import pandas as pd

from memilio.epidata import defaultDict as dd
from memilio.epidata import progress_indicator


class VerbosityLevel(Enum):
    """ """
    Off = 0
    Critical = 1
    Error = 2
    Warning = 3
    Info = 4
    Debug = 5
    Trace = 6


class Conf:
    """Configures all relevant download outputs etc."""

    v_level = 'Info'
    show_progr = False
    if parse_version(pd.__version__) < parse_version('2.2'):
        excel_engine = 'openpyxl'
    else:
        # calamine is faster, but cannot be used for pandas < 2.2
        # also there are issues with pd >= 2.2 and openpyxl engine
        excel_engine = 'calamine'

    def __init__(self, out_folder, **kwargs):

        # change v_level from int to str
        if 'verbosity_level' in kwargs.keys():
            if isinstance(kwargs['verbosity_level'], int):
                kwargs['verbosity_level'] = VerbosityLevel(
                    kwargs['verbosity_level']).name

        path = os.path.join(os.path.dirname(
            os.path.abspath(__file__)), 'download_config.conf')

        # activate CoW for more predictable behaviour of pandas DataFrames
        pd.options.mode.copy_on_write = True
        # read in config file
        # if no config file is given, use default values
        if os.path.exists(path):
            parser = configparser.ConfigParser()
            parser.read(path)
            # all values will be read in as string

            if parser['SETTINGS']['path_to_use'] == 'default':
                self.path_to_use = out_folder
            else:
                self.path_to_use = parser['SETTINGS']['path_to_use']

            matplotlib.use(str(parser['SETTINGS']['mpl_backend']))

            # merge kwargs with config data
            # Do not overwrite kwargs, just add from parser
            for key in parser['SETTINGS']:
                if key not in kwargs:
                    kwargs.update({key: parser['SETTINGS'][key]})

            Conf.show_progr = True if str(
                kwargs['show_progress']) == 'True' else False
            Conf.v_level = str(kwargs['verbosity_level'])
            self.checks = True if str(
                kwargs['run_checks']) == 'True' else False
            self.interactive = True if str(
                kwargs['interactive']) == 'True' else False
            self.plot = True if str(kwargs['make_plot']) == 'True' else False
            self.no_raw = True if str(kwargs['no_raw']) == 'True' else False
            self.to_dataset = True if str(
                kwargs['to_dataset']) == 'True' else False
        else:
            # default values:
            Conf.show_progr = kwargs['show_progress'] if 'show_progress' in kwargs.keys(
            ) else Conf.show_progr
            Conf.v_level = kwargs['verbosity_level'] if 'verbosity_level' in kwargs.keys(
            ) else Conf.v_level
            self.checks = kwargs['run_checks'] if 'run_checks' in kwargs.keys(
            ) else True
            self.interactive = kwargs['interactive'] if 'interactive' in kwargs.keys(
            ) else False
            self.plot = kwargs['make_plot'] if 'make_plot' in kwargs.keys(
            ) else dd.defaultDict['make_plot']
            self.no_raw = kwargs['no_raw'] if 'no_raw' in kwargs.keys(
            ) else dd.defaultDict['no_raw']
            self.path_to_use = out_folder
            self.to_dataset = kwargs['to_dataset'] if 'to_dataset' in kwargs.keys(
            ) else False

        # suppress Future & DepricationWarnings
        if VerbosityLevel[Conf.v_level].value <= 2:
            warnings.simplefilter(action='ignore', category=FutureWarning)
            warnings.simplefilter(action='ignore', category=DeprecationWarning)
        # deactivate (or activate progress indicator)
        if Conf.show_progr == True:
            progress_indicator.ProgressIndicator.disable_indicators(False)
        else:
            progress_indicator.ProgressIndicator.disable_indicators(True)


def default_print(verbosity_level, message):
    """

    :param verbosity_level: 
    :param message: 

    """
    if VerbosityLevel[verbosity_level].value <= VerbosityLevel[Conf.v_level].value:
        print(verbosity_level + ": " + message)


def user_choice(message, default=False):
    """

    :param message: 
    :param default:  (Default value = False)

    """
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
    """ Download a file using GET over HTTP.

    :param url: Full url of the file to download.
    :param chunk_size: Number of Bytes downloaded at once. Only used when a
        progress_function is specified. For a good display of progress, this
        size should be about the speed of your internet connection in Bytes/s.
        Can be set to None to let the server decide the chunk size (may be
        equal to the file size). (Default value = 1024)
    :param timeout: Timeout in seconds for the GET request. (Default value = None)
    :param progress_function: Function called regularly, with the current
        download progress in [0,1] as a float argument. (Default value = None)
    :param verify: bool or "interactive". If False, ignores the connection's
        security. If True, only starts downloads from secure connections, and
        insecure connections raise a FileNotFoundError. If "interactive",
        prompts the user whether or not to allow insecure connections. (Default value = True)
    :returns: File as BytesIO

    """
    if verify not in [True, False, "interactive"]:
        warnings.warn('Invalid input for argument verify. Expected True, False, or'
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
    if ('content-length' in req.headers) and progress_function:
        # get file size from http header
        # this is only the number of bytes downloaded, the size of the actual file
        # may be larger (e.g. when 'content-encoding' is gzip; decoding is handled
        # by iter_content)
        # this is only needed for the progress indicator
        file_size = int(req.headers.get('content-length'))
        # if content length is not known, a progress cant be set.
        set_progr = True
    else:
        set_progr = False
    file = bytearray()  # file to be downloaded
    if set_progr:
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
    """ reads a zip file and returns a list of dataframes for every file in the zip folder.
    If only one file is readable for func_to_use a single dataframe is returned instead of a list with one entry.

    :param file: String. Path to Zipfile to read.
    :param param_dict: Dict. Additional information for download functions (e.g. engine, sheet_name, header...)
    :returns: list od all dataframes (one for each file).

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
    """ Loads data from filepath and stores it in a pandas dataframe.
    If data can't be read from given filepath the user is asked whether the file should be downloaded from the given url or not.
    Uses the progress indicator to give feedback.

    :param filepath: String. Filepath from where the data is read. (Default value = '')
    :param url: String. URL to download the dataset. (Default value = '')
    :param read_data: True or False. Defines if item is opened from directory (True) or downloaded (False). (Default value = dd.defaultDict['read_data'])
    :param param_dict: Dict. Additional information for download functions (e.g. engine, sheet_name, header...) (Default value = {})
    :param interactive: bool. Whether to ask for user input. If False, raises Errors instead. (Default value = False)
    :returns: pandas dataframe

    """
    param_dict_excel = {"sheet_name": 0,
                        "header": 0, "engine": Conf.excel_engine}
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
    """ Defines command line interface

    The function parameter "what" is used as a dictionary key.
    The return of the dictionary is either a list of a string and a list of keywords.
    The string is the message that should be printed when working on the specific package.
    The further list, contains all special command line arguments which are needed for this package.

    If the key is not part of the dictionary the program is stopped.

    The following default arguments are added to the parser:
    - read-file
    - file-format, choices = ['json', 'hdf5', 'json_timeasstring']
    - out_path
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
    - no_progress_indicator
    - interactive
    - verbose
    - skip_checks
    - no_raw
    - to_dataset

    :param what: Defines what packages calls and thus what kind of command line arguments should be defined.

    """

    # TODO: may it would be easier to make a dict like the following one together with a function to get key:
    # TODO: all should automatically do everything
    # cli_dict2 = {"end_date": ['divi'],
    #                "plot": ['cases'],
    #                "start_date": ['divi']                 }

    cli_dict = {"divi": ['Downloads data from DIVI', 'start_date', 'end_date', 'impute_dates', 'moving_average'],
                "cases": ['Download case data from RKI', 'start_date', 'end_date', 'impute_dates', 'moving_average', 'split_berlin', 'rep_date', 'files'],
                "population": ['Download population data from official sources'],
                "commuter_official": ['Download commuter data from official sources'],
                "vaccination": ['Download vaccination data', 'start_date', 'end_date', 'impute_dates', 'moving_average', 'sanitize_data'],
                "testing": ['Download testing data', 'start_date', 'end_date', 'impute_dates', 'moving_average'],
                "jh": ['Downloads data from Johns Hopkins University', 'start_date', 'end_date', 'impute_dates', 'moving_average'],
                "hospitalization": ['Download hospitalization data', 'start_date', 'end_date', 'impute_dates', 'moving_average'],
                "sim": ['Download all data needed for simulations', 'start_date', 'end_date', 'impute_dates', 'moving_average', 'split_berlin', 'rep_date', 'sanitize_data']}

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

    if 'start_date' in what_list:
        if what == 'divi':
            start_date_default = datetime.date(2020, 4, 24)
        elif what == 'jh':
            start_date_default = datetime.date(2020, 1, 22)
        else:
            start_date_default = dd.defaultDict['start_date']
        parser.add_argument(
            '-s', '--start-date', default=start_date_default,
            help='Defines start date for data download. Should have form: YYYY-mm-dd.'
            'Default is ' +
            str(dd.defaultDict['start_date']) +
            ' (2020-04-24 for divi and 2020-01-22 for jh)',
            type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d').date())
    if 'end_date' in what_list:
        parser.add_argument(
            '-e', '--end-date', default=dd.defaultDict['end_date'],
            help='Defines date after which data download is stopped.'
            'Should have form: YYYY-mm-dd. Default is today',
            type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d').date())
    if 'impute_dates' in what_list:
        parser.add_argument(
            '-i', '--impute-dates', default=dd.defaultDict['impute_dates'],
            help='the resulting dfs contain all dates instead of'
            ' omitting dates where no data was reported', action='store_true')
    if 'moving_average' in what_list:
        parser.add_argument(
            '-m', '--moving-average', type=int, default=dd.defaultDict['moving_average'],
            help='Compute a moving average of N days over the time series. Default is ' + str(dd.defaultDict['moving_average']))
    if 'split_berlin' in what_list:
        parser.add_argument(
            '-b', '--split-berlin', default=dd.defaultDict['split_berlin'],
            help='Berlin data is split into different counties,'
            ' instead of having only one county for Berlin.',
            action='store_true')
    if 'rep_date' in what_list:
        parser.add_argument(
            '--rep-date', default=dd.defaultDict['rep_date'],
            help='If reporting date is activated, the reporting date'
            'will be prefered over possibly given dates of disease onset.',
            action='store_true')
    if 'sanitize_data' in what_list:
        parser.add_argument(
            '-sd', '--sanitize-data', type=int, default=dd.defaultDict['sanitize_data'], dest='sanitize_data',
            help='Redistributes cases of every county either based on regions ratios or on thresholds and population'
        )
    if 'files' in what_list:
        parser.add_argument(
            '--files', nargs="*", default='All'
        )
    if 'ref_year' in what_list:
        parser.add_argument(
            '--ref-year', default='newest',
            help='Considered year.'
        )

    # add optional download options
    if '--no-progress-indicators' in sys.argv:
        parser.add_argument(
            '--no-progress-indicators', dest='show_progress',
            help='Disables all progress indicators (used for downloads etc.).',
            action='store_false')

    if not {'--no-raw', '-n'}.isdisjoint(sys.argv):
        parser.add_argument(
            '-n', '--no-raw',
            help='Defines if raw data will be stored for further use.',
            action='store_true')

    if not {'--make_plot', '-p'}.isdisjoint(sys.argv):
        parser.add_argument('-p', '--make-plot',
                            help='Plots the data.', action='store_true')

    if '--interactive' in sys.argv:
        parser.add_argument(
            '--interactive',
            help='Interactive download (Handle warnings, passwords etc.).', action='store_true')

    if not {'--verbose', '-v', '-vv', '-vvv', '-vvvv', '-vvvvv', '-vvvvvv'}.isdisjoint(sys.argv):
        parser.add_argument(
            '-v', '--verbose', dest='verbosity_level',
            help='Increases verbosity level (Trace, Debug, Info, Warning, Error, Critical, Off).',
            action='count', default=0)

    if '--skip-checks' in sys.argv:
        parser.add_argument(
            '--skip-checks', dest='run_checks', action='store_false',
            help='Skips sanity checks etc.')

    if '--to-dataset' in sys.argv:
        parser.add_argument(
            '--to-dataset', dest='to_dataset',
            help="To return saved dataframes as objects.",
            action='store_true'
        )

    args = vars(parser.parse_args())

    return args


def append_filename(
        filename='', impute_dates=False, moving_average=0, split_berlin=False,
        rep_date=False):
    """ Creates consistent file names for all output.

    :param filename:  (Default value = '')
    :param impute_dates:  (Default value = False)
    :param moving_average:  (Default value = 0)
    :param split_berlin:  (Default value = False)
    :param rep_date:  (Default value = False)

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
    """ Checks existence and creates folder

    It is checked if the folder given in the parameter "directory" exists.
    If it does not exist it is created.

    :param directory: directory which should exist

    """

    # check if directory exists or create it
    if not os.path.exists(directory):
        os.makedirs(directory)


def write_dataframe(df, directory, file_prefix, file_type, param_dict={}):
    """ Writes pandas dataframe to file

    This routine writes a pandas dataframe to a file in a given format.
    The filename is given without ending.
    A file_type can be
    - json
    - json_timeasstring [Default]
    - hdf5
    - csv
    - txt
    The file_type defines the file format and thus also the file ending.
    The file format can be json, hdf5, csv or txt.
    For this option the column Date is converted from datetime to string.

    :param df: pandas dataframe (pandas DataFrame)
    :param directory: directory where to safe (string)
    :param file_prefix: filename without ending (string)
    :param file_type: defines ending (string)
    :param param_dict: defines parameters for to_csv/txt(dictionary) (Default value = {})

    """

    outForm = {'json': [".json", {"orient": "records"}],
               'json_timeasstring': [".json", {"orient": "records"}],
               'hdf5': [".h5", {"key": "data"}],
               'csv': [".csv", {}],
               'txt': [".txt", param_dict]}

    try:
        outFormEnd = outForm[file_type][0]
        outFormSpec = outForm[file_type][1]
    except KeyError:
        raise ValueError(
            "Error: The file format: " + file_type +
            " does not exist. Use json, json_timeasstring, hdf5, csv or txt.")

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
    elif file_type == 'csv':
        df.to_csv(out_path)
    elif file_type == "txt":
        df.to_csv(out_path, **outFormSpec)

    default_print('Info', "Data has been written to " + out_path)


class DataError(Exception):
    """Error for handling incomplete or unexpected Data"""
    pass
