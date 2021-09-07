#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Kathrin Rack, Wadim Koslow
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
@file getRKIData.py
@brief Downloads the data of the Robert-Koch-Institut (RKI) and provides it in different ways.

The RKI data we download can be found at
https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0

Be careful: Recovered and deaths are not correct set in this case
"""

# Imports
import os
import sys
import itertools
import pandas
import numpy as np
import matplotlib.pyplot as plt

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd


def check_for_completeness(df):
    """! Checks if all states are mentioned

   This check had to be added, due to incomplete data downloads because of changes from RKI-site
   It is checked if all 16 federal states and all 401, resp. 412 (400+Berlin), counties are part of the data.
   If data is incomplete the data is downloaded from another source.
   Note: There is no check if data for every day and every county is available (which can happen).

   @param df pandas dataframe to check
   @return Boolean to say if data is complete or not
   """

    if not df.empty:
        # check for 400 counties outside Berlin plus 12 districts of Berlin
        if len(df["IdLandkreis"].unique()) < 412:
            print("Downloaded RKI data is not complete. Another option will be tested.")
            return False

        # if it seems complete
        return True

    # if it is empty
    return False


def fill_df(rki_old, group_by, avg_by, moving_average):
    """! Fills missing dates of df and optionally calculates the the 7 day moving average of the rki data

    @param rki_old pandas dataframe
    @param group_by List of columns which hold the specifications of the dataframe (Region, age, gender,...)
    @param avg_by List of columns for which the rolling average should be computed (Confirmed, recovered, deaths)
    @param moving_average Flag which indicates whether to compute the moving average
    @return dataframe with filled dates/moving average
    """
    # drop time from date
    try:
        rki_old.Date = rki_old.Date.dt.date
    except:
        rki_old['Date'] = pandas.to_datetime(rki_old['Date'])
        rki_old.Date = rki_old.Date.dt.date

    #create empty copy of the df
    rki_new = pandas.DataFrame(columns=rki_old.columns)

    # range of dates which should be filled
    idx = pandas.date_range(min(rki_old.Date), max(rki_old.Date))

    # create list of all possible groupby columns combinations
    unique_ids = []
    for group in group_by:
        unique_ids.append(list(rki_old[group].unique()))
    unique_ids_comb = list(itertools.product(*unique_ids))

    #loop over all regions/ages/gender
    for ids in unique_ids_comb:
        df = rki_old.copy()
        counter = 0
        # filter df
        while counter < len(ids):
            df = df[df[group_by[counter]] == ids[counter]]
            counter += 1

        # create missing dates
        df.index = df.Date
        df_new = df.reindex(idx)
        df_new.Date = idx
        df_new.index = (range(len(idx)))

        if len(df) > 0:
            # create values for first date
            values = {}
            for column in df.columns:
                values[column] = df[column][0]
            for avg in avg_by:
                values[avg] = 0

            # fill values of missing dates based on last entry
            df_new.fillna(method='ffill', inplace=True)
            # fill value of the first date, if it doesn't exist yet
            df_new.fillna(values, limit=1, inplace=True)
            # fill remaining values (between first date and first reported date of the df)
            df_new.fillna(method='ffill', inplace=True)

            # compute 7 day moving average
            if moving_average:
                for avg in avg_by:
                    df_new['MA' + avg] = df_new[avg].rolling(window=7, min_periods=1).mean()
                    df_new['MA' + avg] = df_new['MA' + avg].fillna(df_new[avg])
                    df_new[avg] = df_new['MA' + avg]
                    df_new.drop('MA' + avg, axis=1, inplace=True)

            rki_new = rki_new.append(df_new)
            rki_new.index = (range(len(rki_new)))

    return rki_new


def fuse_berlin(df):
    """! Concatenates the different districts of Berlin into one district

    The RKI data for Berlin is devided into 7 different districts.
    This does not correspond to the other datasets, which usually has only
    one entry for Berlin.
    This function is used to replace the entries of the 7 different
    districts with only one county, which is called 'Berlin'.

    @param df pandas dataframe
    @return dataframe with fused Berlin
    """

    berlin = df[(df['ID_County'].values / 1000).astype(int) == 11]
    berlin = berlin.groupby(['Date', 'Gender', 'ID_State', 'State', 'County', 'Age_RKI']).agg('sum').reset_index()

    berlin[dd.EngEng['idCounty']] = 11000
    berlin[dd.EngEng['county']] = dd.County[11000]

    new_df = df[(df[dd.EngEng['idCounty']].values / 1000).astype(int) != 11]
    new_df = pandas.concat([new_df, berlin], axis=0)

    dateToUse = 'Date'
    new_df.sort_values([dateToUse], inplace=True)

    return new_df


def get_rki_data(read_data=dd.defaultDict['read_data'],
                 file_format=dd.defaultDict['file_format'],
                 out_folder=dd.defaultDict['out_folder'],
                 no_raw=dd.defaultDict['no_raw'],
                 fill_dates=dd.defaultDict['fill_dates'],
                 make_plot=dd.defaultDict['make_plot'],
                 moving_average=dd.defaultDict['moving_average'],
                 split_berlin=dd.defaultDict['split_berlin']
                 ):
    """! Downloads the RKI data and provides different kind of structured data

    The data is read either from the internet or from a json file (FullDataRKI.json), stored in an earlier run.
    If the data is read from the internet, before changing anything the data is stored in FullDataRKI.json.
    If data should be downloaded, it is checked if data contains all 16 states.
    If not a different source is tried, in this case a column has to be renamed.
    The file is read in or stored at the folder "out_folder"/Germany/.
    To store and change the data we use pandas

    While working with the data
    - the column names are changed to english depending on defaultDict
    - a new column "Date" is defined.
    - we are only interested in the values where the parameter NeuerFall, NeuerTodesfall, NeuGenesen are larger than 0.
    The values, when these parameters are negative are just useful,
    if one would want to get the difference to the previous day.
    For details we refer to the above mentioned webpage.
    - For all different parameters and different columns the values are added up for whole germany for every date
    and the cumulative sum is calculated. Unless something else is mentioned.
    - For Berlin all counties can be fused to one [Default]. If Berlin is split, different file names are used.
    - Following data is generated and written to the mentioned filename
        - All infected (current and past) for whole germany are stored in "infected_rki"
        - All deaths whole germany are stored in "deaths_rki"
        - Infected, deaths and recovered for whole germany are stored in "all_germany_rki"
        - Infected split for states are stored in "infected_state_rki"
        - Infected, deaths and recovered split for states are stored in "all_state_rki"
        - Infected split for counties are stored in "infected_county_rki(_split_berlin)"
        - Infected, deaths and recovered split for county are stored in "all_county_rki(_split_berlin)"
        - Infected, deaths and recovered split for gender are stored in "all_gender_rki"
        - Infected, deaths and recovered split for state and gender are stored in "all_state_gender_rki"
        - Infected, deaths and recovered split for county and gender are stored in "all_county_gender_rki(_split_berlin)"
        - Infected, deaths and recovered split for age are stored in "all_age_rki"
        - Infected, deaths and recovered split for state and age are stored in "all_state_age_rki"
        - Infected, deaths and recovered split for county and age are stored in "all_county_age_rki(_split_berlin)"

    @param read_data False [Default] or True. Defines if data is read from file or downloaded.
    @param file_format File format which is used for writing the data. Default defined in defaultDict.
    @param out_folder Path to folder where data is written in folder out_folder/Germany.
    @param no_raw True or False [Default]. Defines if unchanged raw data is saved or not.
    @param fill_dates False [Default] or True. Defines if values for dates without new information are imputed.
    @param make_plot False [Default] or True. Defines if plots are generated with matplotlib.
    @param moving_average True or False [Default]. Defines if files for 7 day moving average should be created.
    @param split_berlin True or False [Default]. Defines if Berlin counties is fused to just on county.
    """

    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)

    filename = "FullDataRKI"

    if read_data:
        # if once dowloaded just read json file

        file_in = os.path.join(directory, filename + ".json")
        try:
            df = pandas.read_json(file_in)
        except ValueError:
            exit_string = "Error: The file: " + file_in + " does not exist. Call program without -r flag to get it."
            sys.exit(exit_string)
    else:

        # Supported data formats:
        load = {
            'csv': gd.loadCsv,
            'geojson': gd.loadGeojson
        }

        # ArcGIS public data item ID:
        itemId = 'dd4580c810204019a7b8eb3e0b329dd6_0'

        # Get data:
        df = load['csv'](itemId)

        complete = check_for_completeness(df)

        # try another possibility if df was empty or incomplete
        if not complete:
            print("Note: RKI data is incomplete. Trying another source.")

            df = load['csv']("", "https://npgeo-de.maps.arcgis.com/sharing/rest/content/items/"
                                 "f10774f1c63e40168479a1feb6c7ca74/data", "")

            df.rename(columns={'FID': "ObjectId"}, inplace=True)
            complete = check_for_completeness(df)

        if complete:
            # output data to not always download it
            if not no_raw:
                gd.write_dataframe(df, directory, filename, "json")
        else:
            print("Information: dataframe was incomplete for csv. Trying geojson.")
            df = load['geojson'](itemId)

            complete = check_for_completeness(df)

            if not df.empty and complete:
                if not no_raw:
                    gd.write_dataframe(df, directory, filename, "json")
            else:
                exit_string = "Something went wrong, dataframe is empty for csv and geojson!"
                sys.exit(exit_string)

    # store dict values in parameter to not always call dict itself
    Altersgruppe2 = dd.GerEng['Altersgruppe2']
    Altersgruppe = dd.GerEng['Altersgruppe']
    Geschlecht = dd.GerEng['Geschlecht']
    AnzahlFall = dd.GerEng['AnzahlFall']
    AnzahlGenesen = dd.GerEng['AnzahlGenesen']
    AnzahlTodesfall = dd.GerEng['AnzahlTodesfall']
    IdBundesland = dd.GerEng['IdBundesland']
    Bundesland = dd.GerEng['Bundesland']
    IdLandkreis = dd.GerEng['IdLandkreis']
    Landkreis = dd.GerEng['Landkreis']

    # translate column gender from German to English and standardize
    df.loc[df.Geschlecht == 'unbekannt', ['Geschlecht']] = dd.GerEng['unbekannt']
    df.loc[df.Geschlecht == 'W', ['Geschlecht']] = dd.GerEng['W']
    df.loc[df.Geschlecht == 'M', ['Geschlecht']] = dd.GerEng['M']
    df.loc[df.Altersgruppe == 'unbekannt', ['Altersgruppe']] = dd.GerEng['unbekannt']
    df.loc[df.Altersgruppe2 == 'unbekannt', ['Altersgruppe2']] = dd.GerEng['unbekannt']

    # change names of columns
    df.rename(dd.GerEng, axis=1, inplace=True)

    # Add column 'Date' with Date= Refadtum if IstErkrankungsbeginn = 1 else take Meldedatum
    df['Date'] = np.where(df['IstErkrankungsbeginn'] == 1, df['Refdatum'], df['Meldedatum'])

    # TODO: uncomment if ALtersgruppe2 will again be provided
    # Add new column with Age with range 10 as spain data
    # conditions = [
    #   (df[Altersgruppe2] == '0-4') & (df[Altersgruppe2] == '5-9'),
    #   (df[Altersgruppe2] == '10-14') & (df[Altersgruppe2] == '15-19'),
    #   (df[Altersgruppe2] == '20-24') & (df[Altersgruppe2] == '25-29'),
    #   (df[Altersgruppe2] == '30-34') & (df[Altersgruppe2] == '35-39'),
    #   (df[Altersgruppe2] == '40-44') & (df[Altersgruppe2] == '45-49'),
    #   (df[Altersgruppe2] == '50-54') & (df[Altersgruppe2] == '55-59'),
    #   (df[Altersgruppe2] == '60-64') & (df[Altersgruppe2] == '65-69'),
    #   (df[Altersgruppe2] == '70-74') & (df[Altersgruppe2] == '75-79'),
    # ]

    # choices = ['0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70-79']
    # df['Age10'] = np.select(conditions, choices, default=dd.GerEng['unbekannt'])

    # convert "Datenstand" to real date:
    df.Datenstand = pandas.to_datetime(df.Datenstand, format='%d.%m.%Y, %H:%M Uhr')

    # Correct Timestampes:
    for col in ['Meldedatum', 'Refdatum', 'Date']:
        df[col] = df[col].astype('datetime64[ns]')

    # Be careful "Refdatum" may act different to official describtion
    # on https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0,
    # sometimes big difference identified between "Refdatum" and "Meldedatum"
    # New possibility Date is either Refdatum or Meldedatum after column 'IstErkrankungsbeginn' has been added.

    dateToUse = 'Date'
    df.sort_values(dateToUse, inplace=True)

    # Manipulate data to get rid of conditions: df.NeuerFall >= 0, df.NeuerTodesfall >= 0, df.NeuGenesen >=0
    # There might be a better way
    dfF = df

    dfF.loc[dfF.NeuerFall < 0, [AnzahlFall]] = 0
    dfF.loc[dfF.NeuerTodesfall < 0, [AnzahlTodesfall]] = 0
    dfF.loc[dfF.NeuGenesen < 0, [AnzahlGenesen]] = 0

    # get rid of unnecessary columns
    dfF = dfF.drop(['NeuerFall', 'NeuerTodesfall', 'NeuGenesen', "IstErkrankungsbeginn", "ObjectId",
                    "Meldedatum", "Datenstand", "Refdatum", Altersgruppe2], 1)

    print("Available columns:", df.columns)

    ######## Data for whole Germany all ages ##########

    # NeuerFall: Infected (incl. recovered) over "dateToUse":

    # make sum for one "dateToUse"
    # old way:
    # gbNF = df[df.NeuerFall >= 0].groupby( dateToUse ).sum()
    gbNF = df[df.NeuerFall >= 0].groupby(dateToUse).agg({AnzahlFall: sum})

    # make cumulative sum of "AnzahlFall" for "dateToUse"
    # old way:
    # gbNF_cs = gbNF.AnzahlFall.cumsum()
    gbNF_cs = gbNF.cumsum()

    # output to json file
    gd.write_dataframe(gbNF_cs.reset_index(), directory, "infected_rki", file_format)
    if fill_dates or moving_average:
        gbNF_cs = fill_df(gbNF_cs.reset_index(), [], ['Confirmed'], moving_average)
        if moving_average:
            gd.write_dataframe(gbNF_cs, directory, "infected_ma_rki", file_format)
        elif fill_dates:
            gd.write_dataframe(gbNF_cs, directory, "infected_all_dates_rki", file_format)

    if make_plot:
        # make plot
        gbNF_cs.plot(title='COVID-19 infections', grid=True,
                     style='-o')
        plt.tight_layout()
        plt.show()

    # Dead over Date:
    gbNT = df[df.NeuerTodesfall >= 0].groupby(dateToUse).agg({AnzahlTodesfall: sum})
    gbNT_cs = gbNT.cumsum()

    # output
    gd.write_dataframe(gbNT_cs.reset_index(), directory, "deaths_rki", file_format)
    if fill_dates or moving_average:
      gbNT_cs = fill_df(gbNT_cs.reset_index(), [], ['Deaths'], moving_average)
      if moving_average:
         gd.write_dataframe(gbNT_cs.reset_index(), directory, "deaths_ma_rki", file_format)
      elif fill_dates:
         gd.write_dataframe(gbNT_cs.reset_index(), directory, "deaths_all_dates_rki", file_format)

    if make_plot:
        gbNT_cs.plot(title='COVID-19 deaths', grid=True,
                     style='-o')
        plt.tight_layout()
        plt.show()

        dfF.agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum}) \
            .plot(title='COVID-19 infections, deaths, recovered', grid=True,
                  kind='bar')
        plt.tight_layout()
        plt.show()

    gbNF = df.groupby(dateToUse).agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    gbNF_cs = gbNF.cumsum()

    gd.write_dataframe(gbNF_cs.reset_index(), directory, "all_germany_rki", file_format)
    if fill_dates or moving_average:
        gbNF_cs = fill_df(gbNF_cs.reset_index(), [], ['Confirmed', 'Deaths', 'Recovered'], moving_average)
        if moving_average:
            gd.write_dataframe(gbNF_cs, directory, "all_germany_ma_rki", file_format)
        elif fill_dates:
            gd.write_dataframe(gbNF_cs, directory, "all_germany_all_dates_rki", file_format)

    ############## Data for states all ages ################

    # NeuerFall: Infected (incl. recovered) over "dateToUse" for every state ("Bundesland"):
    # gbNFst = df[df.NeuerFall >= 0].groupby( [IdBundesland','Bundesland', dateToUse]).AnzahlFall.sum()
    gbNFst = df[df.NeuerFall >= 0].groupby([IdBundesland, Bundesland, dateToUse])\
        .agg({AnzahlFall: sum})

    gbNFst_cs = gbNFst.groupby(level=1).cumsum().reset_index()

    # output
    gd.write_dataframe(gbNFst_cs, directory, "infected_state_rki", file_format)
    if fill_dates or moving_average:
        gbNFst_cs = fill_df(gbNFst_cs, ['ID_State'], ['Confirmed'], moving_average)
        if moving_average:
            gd.write_dataframe(gbNFst_cs, directory, "infected_state_ma_rki", file_format)
        elif fill_dates:
            gd.write_dataframe(gbNFst_cs, directory, "infected_state_all_dates_rki", file_format)

    # output nested json
    # gbNFst_cs.groupby(['IdBundesland', 'Bundesland'], as_index=False) \
    #            .apply(lambda x: x[[dateToUse,'AnzahlFall']].to_dict('r')) \
    #            .reset_index().rename(columns={0:'Dates'})\
    #            .to_json(directory + "gbNF_state_nested.json", orient='records')

    # infected (incl recovered), deaths and recovered together

    gbAllSt = dfF.groupby([IdBundesland, Bundesland, dateToUse]) \
        .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    gbAllSt_cs = gbAllSt.groupby(level=1).cumsum().reset_index()

    # output
    gd.write_dataframe(gbAllSt_cs, directory, "all_state_rki", file_format)
    if fill_dates or moving_average:
        gbAllSt_cs = fill_df(gbAllSt_cs, ['ID_State'], ['Confirmed', 'Deaths', 'Recovered'], moving_average)
        if moving_average:
            gd.write_dataframe(gbAllSt_cs, directory, "all_state_ma_rki", file_format)
        elif fill_dates:
            gd.write_dataframe(gbAllSt_cs, directory, "all_state_all_dates_rki", file_format)

    ############# Data for counties all ages ######################

    if not split_berlin:
        df = fuse_berlin(df)
    # NeuerFall: Infected (incl. recovered) over "dateToUse" for every county ("Landkreis"):
    gbNFc = df[df.NeuerFall >= 0].groupby([IdLandkreis, Landkreis, dateToUse]) \
        .agg({AnzahlFall: sum})
    gbNFc_cs = gbNFc.groupby(level=1).cumsum().reset_index()

    # output
    if split_berlin:
        gd.write_dataframe(gbNFc_cs, directory, "infected_county_split_berlin_rki", file_format)
        if fill_dates or moving_average:
            gbNFc_cs = fill_df(gbNFc_cs, ['ID_County'], ['Confirmed'], moving_average)
            if moving_average:
                gd.write_dataframe(gbNFc_cs, directory, "infected_county_split_berlin_ma_rki", file_format)
            elif fill_dates:
                gd.write_dataframe(gbNFc_cs, directory, "infected_county_split_berlin_all_dates_rki", file_format)
    else:
        gd.write_dataframe(gbNFc_cs, directory, "infected_county_rki", file_format)
        if fill_dates or moving_average:
            gbNFc_cs = fill_df(gbNFc_cs, ['ID_County'], ['Confirmed'], moving_average)
            if moving_average:
                gd.write_dataframe(gbNFc_cs, directory, "infected_county_ma_rki", file_format)
            elif fill_dates:
                gd.write_dataframe(gbNFc_cs, directory, "infected_county_all_dates_rki", file_format)

    # infected (incl recovered), deaths and recovered together

    if not split_berlin:
        dfF = fuse_berlin(dfF)
    gbAllC = dfF.groupby([IdLandkreis, Landkreis, dateToUse]).\
        agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    gbAllC_cs = gbAllC.groupby(level=1).cumsum().reset_index()

    # output
    if split_berlin:
        gd.write_dataframe(gbAllC_cs, directory, "all_county_split_berlin_rki", file_format)
        if fill_dates or moving_average:
            gbAllC_cs = fill_df(gbAllC_cs, ['ID_County'], ['Confirmed', 'Deaths', 'Recovered'], moving_average)
            if moving_average:
                gd.write_dataframe(gbAllC_cs, directory, "all_county_split_berlin_ma_rki", file_format)
            elif fill_dates:
                gd.write_dataframe(gbAllC_cs, directory, "all_county_split_berlin_all_dates_rki", file_format)
    else:
        gd.write_dataframe(gbAllC_cs, directory, "all_county_rki", file_format)
        if fill_dates or moving_average:
            gbAllC_cs = fill_df(gbAllC_cs, ['ID_County'], ['Confirmed', 'Deaths', 'Recovered'], moving_average)
            if moving_average:
                gd.write_dataframe(gbAllC_cs, directory, "all_county_ma_rki", file_format)
            elif fill_dates:
                gd.write_dataframe(gbAllC_cs, directory, "all_county_all_dates_rki", file_format)

    ######### Data whole Germany different gender ##################

    # infected (incl recovered), deaths and recovered together

    gbAllG = dfF.groupby([Geschlecht, dateToUse]) \
        .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})

    gbAllG_cs = gbAllG.groupby(level=0).cumsum().reset_index()

    # output
    gd.write_dataframe(gbAllG_cs, directory, "all_gender_rki", file_format)
    if fill_dates or moving_average:
        gbAllG_cs = fill_df(gbAllG_cs, ['Gender'], ['Confirmed', 'Deaths', 'Recovered'], moving_average)
        if moving_average:
            gd.write_dataframe(gbAllG_cs, directory, "all_gender_ma_rki", file_format)
        elif fill_dates:
            gd.write_dataframe(gbAllG_cs, directory, "all_gender_all_dates_rki", file_format)

    if make_plot:
        dfF.groupby(Geschlecht) \
            .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum}) \
            .plot(title='COVID-19 infections, deaths, recovered', grid=True,
                  kind='bar')
        plt.tight_layout()
        plt.show()

    ############################# Gender and State ######################################################

    # infected (incl recovered), deaths and recovered together

    gbAllGState = dfF.groupby([IdBundesland, Bundesland, Geschlecht, dateToUse]) \
        .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    gbAllGState_cs = gbAllGState.groupby(level=[1, 2]).cumsum().reset_index()

    # output
    gd.write_dataframe(gbAllGState_cs, directory, "all_state_gender_rki", file_format)
    if fill_dates or moving_average:
        gbAllGState_cs = fill_df(gbAllGState_cs, ['ID_State', 'Gender'], ['Confirmed', 'Deaths', 'Recovered'],
                                 moving_average)
        if moving_average:
            gd.write_dataframe(gbAllGState_cs, directory, "all_state_gender_ma_rki", file_format)
        elif fill_dates:
            gd.write_dataframe(gbAllGState_cs, directory, "all_state_gender_all_dates_rki", file_format)

    ############# Gender and County #####################

    gbAllGCounty = dfF.groupby([IdLandkreis, Landkreis, Geschlecht, dateToUse]) \
        .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    gbAllGCounty_cs = gbAllGCounty.groupby(level=[1, 2]).cumsum().reset_index()

    # output
    if split_berlin:
        gd.write_dataframe(gbAllGCounty_cs, directory, "all_county_gender_split_berlin_rki", file_format)
        if fill_dates or moving_average:
            gbAllGCounty_cs = fill_df(gbAllGCounty_cs, ['ID_County', 'Gender'], ['Confirmed', 'Deaths', 'Recovered'],
                                      moving_average)
            if moving_average:
                gd.write_dataframe(gbAllGCounty_cs, directory, "all_county_gender_split_berlin_ma_rki", file_format)
            elif fill_dates:
                gd.write_dataframe(gbAllGCounty_cs, directory, "all_county_gender_split_berlin_all_dates_rki", file_format)
    else:
        gd.write_dataframe(gbAllGCounty_cs, directory, "all_county_gender_rki", file_format)
        if fill_dates or moving_average:
            gbAllGCounty_cs = fill_df(gbAllGCounty_cs, ['ID_County', 'Gender'], ['Confirmed', 'Deaths', 'Recovered'],
                                      moving_average)
            if moving_average:
                gd.write_dataframe(gbAllGCounty_cs, directory, "all_county_gender_ma_rki", file_format)
            elif fill_dates:
                gd.write_dataframe(gbAllGCounty_cs, directory, "all_county_gender_all_dates_rki", file_format)

    ######### Data whole Germany different ages ####################

    # infected (incl recovered), deaths and recovered together

    gbAllA = dfF.groupby([Altersgruppe, dateToUse]) \
        .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    gbAllA_cs = gbAllA.groupby(level=0).cumsum().reset_index()

    # output
    gd.write_dataframe(gbAllA_cs, directory, "all_age_rki", file_format)
    if fill_dates or moving_average:
        gbAllA_cs = fill_df(gbAllA_cs, ['Age_RKI'],
                            ['Confirmed', 'Deaths', 'Recovered'], moving_average)
        if moving_average:
            gd.write_dataframe(gbAllA_cs, directory, "all_age_ma_rki", file_format)
        elif fill_dates:
            gd.write_dataframe(gbAllA_cs, directory, "all_age_all_dates_rki", file_format)

    if make_plot:
        dfF.groupby(Altersgruppe) \
            .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum}) \
            .plot(title='COVID-19 infections, deaths, recovered for diff ages', grid=True,
                  kind='bar')
        plt.tight_layout()
        plt.show()

        # Dead by "Altersgruppe":
        gbNTAG = df[df.NeuerTodesfall >= 0].groupby(Altersgruppe).agg({AnzahlTodesfall: sum})

        gbNTAG.plot(title='COVID-19 deaths', grid=True,
                    kind='bar')
        plt.tight_layout()
        plt.show()

    ############################# Age and State ######################################################

    ##### Age_RKI #####

    # infected (incl recovered), deaths and recovered together

    gbAllAgeState = dfF.groupby([IdBundesland, Bundesland, Altersgruppe, dateToUse]) \
        .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    gbAllAgeState_cs = gbAllAgeState.groupby(level=[1, 2]).cumsum().reset_index()

    # output
    gd.write_dataframe(gbAllAgeState_cs, directory, "all_state_age_rki", file_format)
    if fill_dates or moving_average:
        gbAllAgeState_cs = fill_df(gbAllAgeState_cs, ['ID_State', 'Age_RKI'],
                                   ['Confirmed', 'Deaths', 'Recovered'], moving_average)
        if moving_average:
            gd.write_dataframe(gbAllAgeState_cs, directory, "all_state_age_ma_rki", file_format)
        elif fill_dates:
            gd.write_dataframe(gbAllAgeState_cs, directory, "all_state_age_all_dates_rki", file_format)

    # TODO: uncomment if ALtersgruppe2 will again be provided
    ##### Age5 and Age10#####

    # infected (incl recovered), deaths and recovered together
    # gbAllAgeState = dfF.groupby([IdBundesland, Bundesland, dd.GerEng['Altersgruppe2'], dateToUse]) \
    #   .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})

    # gbAllAgeState_cs = gbAllAgeState.groupby(level=[1, 2]).cumsum().reset_index()

    # output
    # gd.write_dataframe(gbAllAgeState_cs, directory, "all_state_age5_rki", file_format)

    ##### Age10 #####

    # gbAllAgeState = dfF.groupby([IdBundesland, Bundesland, 'Age10', dateToUse]) \
    #   .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})

    # gbAllAgeState_cs = gbAllAgeState.groupby(level=[1, 2]).cumsum().reset_index()

    # output
    # gd.write_dataframe(gbAllAgeState_cs, directory, "all_state_age10_rki", file_format)

    ############# Age and County #####################

    gbAllAgeCounty = dfF.groupby([IdLandkreis, Landkreis, Altersgruppe, dateToUse]) \
        .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    gbAllAgeCounty_cs = gbAllAgeCounty.groupby(level=[1, 2]).cumsum().reset_index()

    # output
    if split_berlin:
        gd.write_dataframe(gbAllAgeCounty_cs, directory, "all_county_age_split_berlin_rki", file_format)
        if fill_dates or moving_average:
            gbAllAgeCounty_cs = fill_df(gbAllAgeCounty_cs, ['ID_County', 'Age_RKI'],
                                        ['Confirmed', 'Deaths', 'Recovered'], moving_average)
            if moving_average:
                gd.write_dataframe(gbAllAgeCounty_cs, directory, "all_county_age_split_berlin_ma_rki", file_format)
            elif fill_dates:
                gd.write_dataframe(gbAllAgeCounty_cs, directory, "all_county_age_split_berlin_all_dates_rki", file_format)
    else:
        gd.write_dataframe(gbAllAgeCounty_cs, directory, "all_county_age_rki", file_format)
        if fill_dates or moving_average:
            gbAllAgeCounty_cs = fill_df(gbAllAgeCounty_cs, ['ID_County', 'Age_RKI'],
                                        ['Confirmed', 'Deaths', 'Recovered'], moving_average)
            if moving_average:
                gd.write_dataframe(gbAllAgeCounty_cs, directory, "all_county_age_ma_rki", file_format)
            elif fill_dates:
                gd.write_dataframe(gbAllAgeCounty_cs, directory, "all_county_age_all_dates_rki", file_format)

    # TODO: uncomment if ALtersgruppe2 will again be provided
    #### age5 ####

    # gbAllAgeCounty = dfF.groupby([IdLandkreis, Landkreis, Altersgruppe2, dateToUse]) \
    #   .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    # gbAllAgeCounty_cs = gbAllAgeCounty.groupby(level=[1, 2]).cumsum().reset_index()

    # if split_berlin:
    #   gd.write_dataframe(gbAllAgeCounty_cs, directory, "all_county_age5_split_berlin_rki", file_format)
    # else:
    #   gd.write_dataframe(gbAllAgeCounty_cs, directory, "all_county_age5_rki", file_format)

    #### age10 ####

    # gbAllAgeCounty = dfF.groupby( [IdLandkreis, Landkreis, 'Age10', dateToUse])\
    #                  .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    # gbAllAgeCounty_cs = gbAllAgeCounty.groupby(level=[1,2]).cumsum().reset_index()

    # output

    # if split_berlin:
    #    gd.write_dataframe(gbAllAgeCounty_cs, directory, "all_county_age10_split_berlin_rki", file_format)
    # else:
    #    gd.write_dataframe(gbAllAgeCounty_cs, directory, "all_county_age10_rki", file_format)


def main():
    """! Main program entry."""

    arg_dict = gd.cli("rki")
    get_rki_data(**arg_dict)


if __name__ == "__main__":
    main()
