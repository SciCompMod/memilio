#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Martin J. Kuehn
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
from datetime import datetime, timedelta
import time
import os
import pandas as pd
import numpy as np
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata  import getDataIntoPandasDataFrame as gd
from memilio.epidata import defaultDict as dd
from memilio.epidata import modifyDataframeSeries


def transformWeatherData(data_folder,
                         read_data=dd.defaultDict['read_data'],
                         file_format=dd.defaultDict['file_format'],
                         start_date=dd.defaultDict['start_date'],
                         end_date=dd.defaultDict['end_date'],
                         make_plot=dd.defaultDict['make_plot'],
                         moving_average=dd.defaultDict['moving_average'],
                         merge_berlin=True,
                         merge_eisenach=False
                         ):
    """! ...
    @param data_folder Path to folder where data is written in folder 
        data_folder/Germany.
    @param file_format File format which is used for writing the data. 
        Default defined in defaultDict.
    @param start_date [Default = '', taken from read data] Start date
        of stored data frames.
    @param end_date [Default = '', taken from read data] End date of
        stored data frames.
    @param make_plot False [Default] or True. Defines if plots are
        generated with matplotlib.
    @param moving_average 0 [Default] or Number > 0. Defines the number of
        days for which a centered moving average is computed.
    """

    directory = os.path.join(data_folder, 'Germany/')
    gd.check_dir(directory)

    if not read_data:

        df_weather_old = gd.loadCsv('','wetterdaten','.csv')
        df_weather_old.rename(dd.GerEng, axis=1, inplace=True)
        col_old = ['kr_tamm_',  # average temperature
                   'kr_tadn_',  # temperature minimum
                   'kr_tadx_',  # temperature maximum
                   'kr_rsms_',  # rainfall height
                   'kr_sdms_',  # sunshine duration
                   'kr_dimm_']  # dryness index after de Martonne
        col_new = ['temp_av', 'temp_min', 'temp_max', 'rain', 'sun', 'dry']
        cold_old_to_new = dict(zip(col_old, col_new))

    else:  # read formatted file

        filename = 'germany_counties_weather'
        df_weather = pd.read_json(directory + filename + ".json")

    # transform data from original format to desired format
    if not read_data:
        # get county ids
        unique_geo_entities = geoger.get_county_ids(merge_berlin, merge_eisenach, zfill=False)

        start_weather_cols = list(
            df_weather_old.columns).index(
            dd.EngEng['county']) + 1

        if len(unique_geo_entities) < len(df_weather_old):
            raise gd.DataError('Error: County-IDs do not match with file')

        # create new data frame for all NPIs given in the columns, resolved by
        # county and day
        df_weather = pd.DataFrame(
            columns=[dd.EngEng['date']] + [dd.EngEng['idCounty']] +
            col_new)
        # convert NPI data from object to int such that correlations can be
        # computed
        df_weather = df_weather.astype(dict(
            zip(
                [dd.EngEng['date']] + [dd.EngEng['idCounty']] +
                list(cold_old_to_new.values()), ['str', 'int'] +
                ['int' for i in list(cold_old_to_new.values())])))

        # store string dates 'dYYYYMMDD' in list before parsing
        str_dates = sorted(
            set(
                [str_col[-6:] + '15'
                 for str_col in df_weather_old.iloc
                 [:, start_weather_cols:].columns]))
        # convert string dates into other format
        dates_new = [datetime.strptime(old_date, "%Y%m%d")
                     for old_date in str_dates]
        min_date = dates_new[0]+timedelta(days=-14)
        max_date = dates_new[-1]+timedelta(days=14)

        # get corresponding row ranges for all features given by col_old
        col_old_vars = [[] for i in col_old]
        for i in range(len(df_weather_old.columns[start_weather_cols:])):
            for j in range(len(col_old)):
                if col_old[j] in df_weather_old.columns[start_weather_cols+i]:
                    col_old_vars[j].append(i)
        # transform to numpy array to allow for adding integers onto the indices
        for j in range(len(col_old_vars)):
            col_old_vars[j] = np.array(col_old_vars[j])

        # iterate over countyIDs
        start_time = time.perf_counter()
        for countyID in unique_geo_entities:

            # get county-local data frame
            df_local_old = df_weather_old[df_weather_old[dd.EngEng['idCounty']]
                                          == countyID].copy()

            # create columns for date, county ID and NPI code
            df_local_new = pd.DataFrame(columns=df_weather.columns)

            # access weather values matrix and store it as floats
            weather_vals = df_local_old.iloc[:,
                                             start_weather_cols+col_old_vars[0]].astype(float).values
            for i in range(1, len(col_old_vars)):
                # add new rows
                weather_vals = np.append(
                    weather_vals, df_local_old.iloc[:,
                                                    start_weather_cols + col_old_vars[i]].astype(float).values, axis=0)

            # fill in NPI values by transposing from columns to rows
            df_local_new[dd.EngEng['date']] = dates_new
            df_local_new[dd.EngEng['idCounty']] = countyID
            # possible resorting of rows such that they are sorted according to
            # a literal sorting of the code strings
            df_local_new[col_new] = np.transpose(weather_vals)

            df_local_new = modifyDataframeSeries.impute_and_reduce_df(
                                                                      df_local_new,
                                                                      {},
                                                                      col_new,
                                                                      impute='forward',
                                                                      moving_average=moving_average,
                                                                      min_date=min_date,
                                                                      max_date=max_date,
                                                                      start_w_firstval=True)

            df_weather = df_weather.append(df_local_new.copy())

        # reset index and drop old index column
        df_weather.reset_index(inplace=True)
        try:
            df_weather = df_weather.drop(columns='index')
        except KeyError:
            pass
        try:
            df_weather = df_weather.drop(columns='level_0')
        except KeyError:
            pass

        print(
            "Time needed: " + str(time.perf_counter()-start_time) + " sec")

        #### start validation ####
        
        #### end validation ####

        if moving_average > 0:
            filename = 'germany_counties_weather_ma_' + str(moving_average)
        else:
            filename = 'germany_counties_weather'
        gd.write_dataframe(df_weather, directory, filename, file_format)


def main():
    """! Main program entry."""

    # arg_dict = gd.cli("testing")
    path = os.path.join(os.getcwd(), 'data', 'pydata')
    transformWeatherData(path, read_data=False, make_plot=True, moving_average=30)


if __name__ == "__main__":

   main()
