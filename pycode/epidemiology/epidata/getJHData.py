"""
@file getJHData.py

@brief Download data from John Hopkins University
"""

import os
import sys
import pandas

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd


def get_jh_data(read_data=dd.defaultDict['read_data'],
                out_form=dd.defaultDict['out_form'],
                out_folder=dd.defaultDict['out_folder']):
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

   @param read_data False [Default] or True. Defines if data is read from file or downloaded.
   @param out_form File format which is used for writing the data. Default defined in defaultDict.
   @param out_folder Path to folder where data is written in folder out_folder/Germany.
   """

    filename = "FullData_JohnHopkins"

    if read_data:
        file_in = os.path.join(out_folder, filename + ".json")
        # if once dowloaded just read json file
        try:
            df = pandas.read_json(file_in)
        except ValueError:
            exit_string = "Error: The file: " + file_in + " does not exist. Call program without -r flag to get it."
            sys.exit(exit_string)
    else:
        # Get data:
        # https://raw.githubusercontent.com/datasets/covid-19/master/data/time-series-19-covid-combined.csv
        df = gd.loadCsv('time-series-19-covid-combined',
                        apiUrl='https://raw.githubusercontent.com/datasets/covid-19/master/data/')

        # convert "Datenstand" to real date:
        # df.Datenstand = pandas.to_datetime( df.Datenstand, format='%d.%m.%Y, %H:%M Uhr').dt.tz_localize
        # ('Europe/Berlin')

        # output data to not always download it
        gd.write_dataframe(df, "", filename, "json")

    df.rename({'Country/Region': 'CountryRegion', 'Province/State': 'ProvinceState'}, axis=1, inplace=True)
    print("Available columns:", df.columns)

    # Change "Korea, South" to SouthKorea
    df.loc[df['CountryRegion'] == "Korea, South", ['CountryRegion']] = 'SouthKorea'

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

    gb = df.groupby(['CountryRegion', 'Date']).agg({"Confirmed": sum, "Recovered": sum, "Deaths": sum})

    gd.write_dataframe(gb.reset_index(), out_folder, "all_countries_jh", out_form)

    for key in countries:
        # get data for specific countries
        gb_country = gb.reset_index()[gb.reset_index()["CountryRegion"] == key]
        dir_country = countries[key]
        gd.check_dir(dir_country)
        gd.write_dataframe(gb_country, dir_country, "whole_country_" + key + "_jh", out_form)

    # Check what about external provinces. Should they be added?

    ################ Countries with given States ################################

    # Get all countries where States are given
    dfD = df[~df["ProvinceState"].isnull()]

    gb = dfD.groupby(['CountryRegion', 'ProvinceState', 'Date']).agg(
        {"Confirmed": sum, "Recovered": sum, "Deaths": sum})

    gd.write_dataframe(gb.reset_index(), out_folder, "all_provincestate_jh", out_form)

    # print(dfD[dfD.ProvinceState=="Saskatchewan"])
    # print(gb.reset_index()[gb.reset_index().ProvinceState=="Saskatchewan"])

    # TODO: How to handle empty values which become NaN in the beginnin but after woking on the data its just 0.0
    # One solution is to preserve them with : df['b'] = df['b'].astype(str)
    # However, what to do with the cases where after some times values occur? Do those cases exist?

    # TODO: US more detailes


def main():
    """! Main program entry."""

    [read_data, out_form, out_folder] = gd.cli("jh")
    get_jh_data(read_data, out_form, out_folder)


if __name__ == "__main__":
    main()
