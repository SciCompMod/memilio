#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Patrick Lenz
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
@file getHospitalizationData.py
@brief Downloads the hospitalization data of the Robert Koch-Institute (RKI) and provides it in different ways.

The raw hospitalization data we download can be found at
https://github.com/robert-koch-institut/COVID-19-Hospitalisierungen_in_Deutschland
"""

# Imports
import os
import pandas as pd
import numpy as np

from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import defaultDict as dd

def download_hospitalization_data():
    # RKI content from github
    url = 'https://raw.githubusercontent.com/robert-koch-institut/COVID-19-Hospitalisierungen_in_Deutschland/master/Aktuell_Deutschland_COVID-19-Hospitalisierungen.csv'
    # empty data frame to return if not read correctly
    df = pd.DataFrame()
    # try to read csv
    try:
        df = pd.read_csv(url)
    except Exception:
        print("Error in reading csv while downloading vaccination data.")
        raise
    sanity_checks(df)

    return df
def sanity_checks(df):
    #test if dataframe is empty
    if df.empty:
        raise gd.DataError("Download of Vaccination Data failed. File is empty.")

    actual_strings_list = df.columns.tolist()
    # check number of data categories
    if len(actual_strings_list) != 6:
        raise gd.DataError("Error: Number of data categories changed.")

    # these strings need to be in the header
    test_strings = {
        'Datum', 'Bundesland', 'Bundesland_Id', 'Altersgruppe',
        '7T_Hospitalisierung_Faelle', '7T_Hospitalisierung_Inzidenz'}
    # check if headers are those we want
    for name in test_strings:
        if(name not in actual_strings_list):
            raise gd.DataError("Error: Data categories have changed.")
def get_hospitalization_data(read_data=dd.defaultDict['read_data'],
                             file_format=dd.defaultDict['file_format'],
                             out_folder=dd.defaultDict['out_folder'],
                             no_raw=dd.defaultDict['no_raw'],
                             start_date=dd.defaultDict['start_date'],
                             end_date=dd.defaultDict['end_date'],
                             make_plot=dd.defaultDict['make_plot'],
                             moving_average=dd.defaultDict['moving_average']):
    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)

    # get raw dataframe
    filename = "RKIHospitFull"
    if read_data:
        # read json file for already downloaded data
        file_in = os.path.join(directory, filename + ".json")

        try:
            df_raw = pd.read_json(file_in)
        except ValueError:
            raise FileNotFoundError("Error: The file: " + file_in +
                                    " does not exist. Call program without"
                                    " -r flag to get it.")
    else:

        df_raw = download_hospitalization_data()
        if not no_raw:
            gd.write_dataframe(df_raw, directory, filename, file_format)
    df = df_raw.copy()
    # drop unwanted columns, rename and sort dataframe
    df.rename(dd.GerEng, axis = 1, inplace = True)
    df.rename(columns={'Datum': dd.EngEng['date']}, inplace=True)
    df = df.drop(columns=['State', '7T_Hospitalisierung_Inzidenz'])
    df['Age_RKI'] = np.where((df['Age_RKI'] =='00+'), 'Total', df['Age_RKI'])
    df.sort_values(by=['Date', 'ID_State', 'Age_RKI']).reset_index(drop=True)
    # get data for germany
    df_germany = df[['Date', 'ID_State', 'Age_RKI']].copy()
    df_germany = df_germany[df_germany['ID_State']==1]
    df_germany = df_germany.drop(columns=['ID_State'])
    hospit_ger = np.zeros((7*len(df['Date'].unique())))
    for state in df['ID_State'].unique():
        hospit_ger += df[df.ID_State==state]['7T_Hospitalisierung_Faelle'].values
    df_germany['7T_Hospitalisierung_Faelle'] = hospit_ger
    x=15

def main():
    """! Main program entry."""

    get_hospitalization_data()

if __name__ == "__main__":

    main()
