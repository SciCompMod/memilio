#############################################################################
# Copyright (C) 2020-2024 MEmilio
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
:strong:`getJHData.py`

Download data from John Hopkins University
"""

import os
from datetime import date

import pandas as pd

from memilio.epidata import defaultDict as dd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import modifyDataframeSeries as mdfs

# activate CoW for more predictable behaviour of pandas DataFrames
pd.options.mode.copy_on_write = True


def get_jh_data(read_data=dd.defaultDict['read_data'],
                file_format=dd.defaultDict['file_format'],
                out_folder=dd.defaultDict['out_folder'],
                start_date=date(2020, 1, 22),
                end_date=dd.defaultDict['end_date'],
                impute_dates=dd.defaultDict['impute_dates'],
                **kwargs):
    """ Download data from John Hopkins University

       Data is either downloaded and afterwards stored or loaded from a stored filed.
       The file is "FullData_JohnHopkins.json"

       Working with the data includes
       - rename columns such that "/" is deleted, e.g Country/Region becomes CountryRegion
       - data of all countries together are written to a file
       - download the data from following countries in a separate file
       and are stored in the according folders with the country name
       - Germany, SouthKorea, Spain, France, Italy, US, China
       - furthermore, all countries, for which provinces are added, are written to a file

    :param read_data: True or False. Defines if data is read from file or downloaded. Default defined in defaultDict. (Default value = dd.defaultDict['read_data'])
    :param file_format: File format which is used for writing the data. Default defined in defaultDict. (Default value = dd.defaultDict['file_format'])
    :param out_folder: Folder where data is written to. Default defined in defaultDict. (Default value = dd.defaultDict['out_folder'])
    :param start_date: Date of first date in dataframe. Default defined in defaultDict. (Default value = date(2020, 1, 22))
    :param end_date: Date of last date in dataframe. Default defined in defaultDict. (Default value = dd.defaultDict['end_date'])
    :param impute_dates: Currently not used] True or False. Defines if values for dates without new information are imputed. Default defined in defaultDict. (Default value = dd.defaultDict['impute_dates'])
    :param **kwargs: 

    """
    conf = gd.Conf(out_folder, **kwargs)
    out_folder = conf.path_to_use
    # Generate folders if needed
    directory_ger = os.path.join(out_folder, 'Germany', 'pydata')
    directory_es = os.path.join(out_folder, 'Spain', 'pydata')
    directory_fr = os.path.join(out_folder, 'France', 'pydata')
    directory_it = os.path.join(out_folder, 'Italy', 'pydata')
    directory_us = os.path.join(out_folder, 'US', 'pydata')
    directory_rok = os.path.join(out_folder, 'SouthKorea', 'pydata')
    directory_prc = os.path.join(out_folder, 'China', 'pydata')
    directory_glb = os.path.join(out_folder, 'Global', 'pydata')
    # dictionary of countries
    countries = {
        "Germany": directory_ger,
        "Spain": directory_es,
        "France": directory_fr,
        "Italy": directory_it,
        "US": directory_us,
        "SouthKorea": directory_rok,
        "China": directory_prc,
    }
    gd.check_dir(directory_glb)

    no_raw = conf.no_raw
    if start_date < date(2020, 1, 22):
        gd.default_print("Warning", "First data available on 2020-01-22. "
                         "You asked for " + start_date.strftime("%Y-%m-%d") +
                         ". Changed it to 2020-01-22.")
        start_date = date(2020, 1, 22)

    filename = "FullData_JohnHopkins"
    url = "https://raw.githubusercontent.com/datasets/covid-19/master/data/time-series-19-covid-combined.csv"
    path = os.path.join(out_folder, "Global", "pydata", filename + ".json")
    df = gd.get_file(path, url, read_data, param_dict={},
                     interactive=conf.interactive)

    if not no_raw:
        gd.write_dataframe(df, directory_glb, filename, "json")

    df.rename({'Country/Region': 'CountryRegion',
              'Province/State': 'ProvinceState'}, axis=1, inplace=True)
    gd.default_print("Debug", "Available columns: " + df.columns)

    # extract subframe of dates
    df = mdfs.extract_subframe_based_on_dates(df, start_date, end_date)

    # Change "Korea, South" to SouthKorea
    df.loc[df['CountryRegion'] == "Korea, South",
           ['CountryRegion']] = 'SouthKorea'

    ########### Countries ##########################

    gb = df.groupby(['CountryRegion', 'Date']).agg(
        {"Confirmed": "sum", "Recovered": "sum", "Deaths": "sum"})

    gd.write_dataframe(gb.reset_index(), directory_glb,
                       "all_countries_jh", file_format)

    for key in countries:
        # get data for specific countries
        gb_country = gb.reset_index()[gb.reset_index()["CountryRegion"] == key]
        dir_country = countries[key]
        gd.check_dir(dir_country)
        gd.write_dataframe(gb_country, dir_country,
                           "whole_country_" + key + "_jh", file_format)

    # Check what about external provinces. Should they be added?

    ################ Countries with given States ################################

    # Get all countries where States are given
    dfD = df[~df["ProvinceState"].isnull()]

    gb = dfD.groupby(['CountryRegion', 'ProvinceState', 'Date']).agg(
        {"Confirmed": "sum", "Recovered": "sum", "Deaths": "sum"})

    gd.write_dataframe(gb.reset_index(), directory_glb,
                       "all_provincestate_jh", file_format)

    # TODO: How to handle empty values which become NaN in the beginnin but after woking on the data its just 0.0
    # One solution is to preserve them with : df['b'] = df['b'].astype(str)
    # However, what to do with the cases where after some times values occur? Do those cases exist?

    # TODO: US more detailes


def main():
    """ Main program entry."""

    arg_dict = gd.cli("jh")
    get_jh_data(**arg_dict)


if __name__ == "__main__":
    main()
