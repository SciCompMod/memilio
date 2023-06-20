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
@file getJHData.py

@brief Download data from John Hopkins University
"""

import os
from datetime import date

import pandas

from memilio.epidata import defaultDict as dd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import modifyDataframeSeries as mdfs


def get_jh_data(read_data=dd.defaultDict['read_data'],
                file_format=dd.defaultDict['file_format'],
                out_folder=dd.defaultDict['out_folder'],
                no_raw=dd.defaultDict['no_raw'],
                start_date=date(2020, 1, 22),
                end_date=dd.defaultDict['end_date'],
                impute_dates=dd.defaultDict['impute_dates'],
                moving_average=dd.defaultDict['moving_average'],
                make_plot=dd.defaultDict['make_plot']):
    """! Download data from John Hopkins University

   Data is either downloaded and afterwards stored or loaded from a stored filed.
   The file is "FullData_JohnHopkins.json"

   Working with the data includes
   - rename columns such that "/" is deleted, e.g Country/Region becomes CountryRegion
   - data of all countries together are written to a file
   - download the data from following countries in a separate file
   and are stored in the according folders with the country name
       - Germany, SouthKorea, Spain, France, Italy, US, China
   - furthermore, all countries, for which provinces are added, are written to a file

    @param read_data True or False. Defines if data is read from file or downloaded. Default defined in defaultDict.
    @param file_format File format which is used for writing the data. Default defined in defaultDict.
    @param out_folder Folder where data is written to. Default defined in defaultDict.
    @param no_raw True or False. Defines if unchanged raw data is saved or not. Default defined in defaultDict.
    @param start_date Date of first date in dataframe. Default defined in defaultDict.
    @param end_date Date of last date in dataframe. Default defined in defaultDict.
    @param impute_dates [Currently not used] True or False. Defines if values for dates without new information are imputed. Default defined in defaultDict.
    @param moving_average [Currently not used] Integers >=0. Applies an 'moving_average'-days moving average on all time series
        to smooth out effects of irregular reporting. Default defined in defaultDict.
    @param make_plot [Currently not used] True or False. Defines if plots are generated with matplotlib. Default defined in defaultDict.
   """

    filename = "FullData_JohnHopkins"
    url = "https://raw.githubusercontent.com/datasets/covid-19/master/data/time-series-19-covid-combined.csv"
    path = os.path.join(out_folder, filename + ".json")
    df = gd.get_file(path, url, read_data, param_dict={}, interactive=True)

    if not no_raw:
        gd.write_dataframe(df, out_folder, filename, "json")

    df.rename({'Country/Region': 'CountryRegion',
              'Province/State': 'ProvinceState'}, axis=1, inplace=True)
    print("Available columns:", df.columns)

    # extract subframe of dates
    df = mdfs.extract_subframe_based_on_dates(df, start_date, end_date)

    # Change "Korea, South" to SouthKorea
    df.loc[df['CountryRegion'] == "Korea, South",
           ['CountryRegion']] = 'SouthKorea'

    # Generate folders if needed
    directory_ger = os.path.join(out_folder, 'Germany/')
    directory_es = os.path.join(out_folder, 'Spain/')
    directory_fr = os.path.join(out_folder, 'France/')
    directory_it = os.path.join(out_folder, 'Italy/')
    directory_us = os.path.join(out_folder, 'US/')
    directory_rok = os.path.join(out_folder, 'SouthKorea/')
    directory_prc = os.path.join(out_folder, 'China/')

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

    ########### Countries ##########################

    gb = df.groupby(['CountryRegion', 'Date']).agg(
        {"Confirmed": sum, "Recovered": sum, "Deaths": sum})

    gd.write_dataframe(gb.reset_index(), out_folder,
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
        {"Confirmed": sum, "Recovered": sum, "Deaths": sum})

    gd.write_dataframe(gb.reset_index(), out_folder,
                       "all_provincestate_jh", file_format)

    # print(dfD[dfD.ProvinceState=="Saskatchewan"])
    # print(gb.reset_index()[gb.reset_index().ProvinceState=="Saskatchewan"])

    # TODO: How to handle empty values which become NaN in the beginnin but after woking on the data its just 0.0
    # One solution is to preserve them with : df['b'] = df['b'].astype(str)
    # However, what to do with the cases where after some times values occur? Do those cases exist?

    # TODO: US more detailes


def main():
    """! Main program entry."""

    arg_dict = gd.cli("jh")
    get_jh_data(**arg_dict)


if __name__ == "__main__":
    main()
