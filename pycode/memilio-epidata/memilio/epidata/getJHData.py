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
import pandas

from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import defaultDict as dd


def get_jh_data(data_folder,
                read_data=dd.defaultDict['read_data'],
                file_format=dd.defaultDict['file_format'],
                no_raw=dd.defaultDict['no_raw']):
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

   @param data_folder Path to directory where data is written in. Here data_folder/Germany.
   @param read_data False [Default] or True. Defines if data is read from file or downloaded.
   @param file_format File format which is used for writing the data. Default defined in defaultDict.
   @param no_raw True or False [Default]. Defines if unchanged raw data is saved or not.
   """

    filename = "FullData_JohnHopkins"

    if read_data:
        file_in = os.path.join(data_folder, filename + ".json")
        # if once dowloaded just read json file
        try:
            df = pandas.read_json(file_in)
        except ValueError as err:
            raise FileNotFoundError("Error: The file: " + file_in + \
                                  " does not exist." + \
                                  " Call program without -r flag to get it.") \
                from err
    else:
        # Get data:
        # https://raw.githubusercontent.com/datasets/covid-19/master/data/time-series-19-covid-combined.csv
        df = gd.loadCsv('time-series-19-covid-combined',
                        apiUrl='https://raw.githubusercontent.com/datasets/covid-19/master/data/')

        # convert "Datenstand" to real date:
        # df.Datenstand = pandas.to_datetime( df.Datenstand, format='%d.%m.%Y, %H:%M Uhr').dt.tz_localize
        # ('Europe/Berlin')

        # output data to not always download it
        if not no_raw:
            gd.write_dataframe(df, data_folder, filename, "json")

    df.rename({'Country/Region': 'CountryRegion', 'Province/State': 'ProvinceState'}, axis=1, inplace=True)
    print("Available columns:", df.columns)

    # Change "Korea, South" to SouthKorea
    df.loc[df['CountryRegion'] == "Korea, South", ['CountryRegion']] = 'SouthKorea'

    # Generate folders if needed
    directory_ger = os.path.join(data_folder, 'Germany/')
    directory_es = os.path.join(data_folder, 'Spain/')
    directory_fr = os.path.join(data_folder, 'France/')
    directory_it = os.path.join(data_folder, 'Italy/')
    directory_us = os.path.join(data_folder, 'US/')
    directory_rok = os.path.join(data_folder, 'SouthKorea/')
    directory_prc = os.path.join(data_folder, 'China/')

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

    gb = df.groupby(['CountryRegion', 'Date']).agg({"Confirmed": sum, "Recovered": sum, "Deaths": sum})

    gd.write_dataframe(gb.reset_index(), data_folder, "all_countries_jh", file_format)

    for key in countries:
        # get data for specific countries
        gb_country = gb.reset_index()[gb.reset_index()["CountryRegion"] == key]
        dir_country = countries[key]
        gd.check_dir(dir_country)
        gd.write_dataframe(gb_country, dir_country, "whole_country_" + key + "_jh", file_format)

    # Check what about external provinces. Should they be added?

    ################ Countries with given States ################################

    # Get all countries where States are given
    dfD = df[~df["ProvinceState"].isnull()]

    gb = dfD.groupby(['CountryRegion', 'ProvinceState', 'Date']).agg(
        {"Confirmed": sum, "Recovered": sum, "Deaths": sum})

    gd.write_dataframe(gb.reset_index(), data_folder, "all_provincestate_jh", file_format)

    # print(dfD[dfD.ProvinceState=="Saskatchewan"])
    # print(gb.reset_index()[gb.reset_index().ProvinceState=="Saskatchewan"])

    # TODO: How to handle empty values which become NaN in the beginnin but after woking on the data its just 0.0
    # One solution is to preserve them with : df['b'] = df['b'].astype(str)
    # However, what to do with the cases where after some times values occur? Do those cases exist?

    # TODO: US more detailes


def main():
    """! Main program entry."""

    path = os.path.join(os.getcwd(), 'data', 'pydata')
    arg_dict = gd.cli("jh")
    get_jh_data(path, **arg_dict)


if __name__ == "__main__":
    main()
