#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Kathrin Rack, Wadim Koslow, Martin J. Kuehn, Annette Lutz
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
@file getCaseData.py
@brief Downloads the case data of the Robert Koch-Institute (RKI) and provides it in different ways.

The raw case data we download can be found at
https://github.com/robert-koch-institut/SARS-CoV-2-Infektionen_in_Deutschland

Be careful: Date of deaths or recovery is not reported in original case data and will be
extrapolated in this script.
"""

# Imports
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import date

from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import defaultDict as dd
from memilio.epidata import modifyDataframeSeries as mdfs
from memilio.epidata import geoModificationGermany as geoger


def check_for_completeness(df, merge_berlin=False, merge_eisenach=True):
    """! Checks if all counties are mentioned in the case data set

   This check had to be added due to incomplete data downloads
   It is checked if all all counties are part of the data.
   If data is incomplete the data is downloaded from another source.
   Note: There is no check if data for every day and every county is available (which can happen).

   @param df pandas dataframe to check
   @return Boolean to say if data is complete or not
   """

    if not df.empty:
        return geoger.check_for_all_counties(
            df["IdLandkreis"].unique(),
            merge_berlin, merge_eisenach)
    # if it is empty
    return False


def get_case_data(read_data=dd.defaultDict['read_data'],
                  file_format=dd.defaultDict['file_format'],
                  out_folder=dd.defaultDict['out_folder'],
                  no_raw=dd.defaultDict['no_raw'],
                  start_date=date(2020, 1, 1),
                  end_date=dd.defaultDict['end_date'],
                  impute_dates=dd.defaultDict['impute_dates'],
                  moving_average=dd.defaultDict['moving_average'],
                  make_plot=dd.defaultDict['make_plot'],
                  split_berlin=dd.defaultDict['split_berlin'],
                  rep_date=dd.defaultDict['rep_date']
                  ):
    """! Downloads the case data and provides different kind of structured data

    The data is read either from the internet or from a json file (CaseDataFull.json), stored in an earlier run.
    If the data is read from the internet, before changing anything the data is stored in CaseDataFull.json.
    If data should be downloaded, it is checked if data contains all counties.
    If not a different source is tried.
    The file is read in or stored at the folder "out_folder"/Germany/.
    To store and change the data we use pandas.

    While working with the data
    - the column names are changed to english depending on defaultDict
    - a new column "Date" is defined.
    - we are only interested in the values where the parameter NeuerFall, NeuerTodesfall, NeuGenesen are larger than 0.
    The values, when these parameters are negative are just useful,
    if one would want to get the difference to the previous day.
    For details we refer to the above mentioned webpage.
    - For all different parameters and different columns the values are added up for whole germany for every date
    and the cumulative sum is calculated. Unless something else is mentioned.
    - For Berlin all districts can be merged into one [Default]. Otherwise, Berlin is divided into multiple districts and
        different file names are used.
    - Following data is generated and written to the mentioned filename
        - All infected (current and past) for whole germany are stored in "cases_infected"
        - All deaths whole germany are stored in "cases_deaths"
        - Infected, deaths and recovered for whole germany are stored in "cases_all_germany"
        - Infected split for states are stored in "cases_infected_state"
        - Infected, deaths and recovered split for states are stored in "cases_all_state"
        - Infected split for counties are stored in "cases_infected_county(_split_berlin)"
        - Infected, deaths and recovered split for county are stored in "cases_all_county(_split_berlin)"
        - Infected, deaths and recovered split for gender are stored in "cases_all_gender"
        - Infected, deaths and recovered split for state and gender are stored in "cases_all_state_gender"
        - Infected, deaths and recovered split for county and gender are stored in "cases_all_county_gender(_split_berlin)"
        - Infected, deaths and recovered split for age are stored in "cases_all_age"
        - Infected, deaths and recovered split for state and age are stored in "cases_all_state_age"
        - Infected, deaths and recovered split for county and age are stored in "cases_all_county_age(_split_berlin)"

    @param read_data True or False. Defines if data is read from file or downloaded. Default defined in defaultDict.
    @param file_format File format which is used for writing the data. Default defined in defaultDict.
    @param out_folder Folder where data is written to. Default defined in defaultDict.
    @param no_raw True or False. Defines if unchanged raw data is saved or not. Default defined in defaultDict.
    @param start_date Date of first date in dataframe. Default 2020-01-01.
    @param end_date Date of last date in dataframe. Default defined in defaultDict.
    @param impute_dates True or False. Defines if values for dates without new information are imputed. Default defined in defaultDict.
    @param moving_average Integers >=0. Applies an 'moving_average'-days moving average on all time series
        to smooth out effects of irregular reporting. Default defined in defaultDict.
    @param make_plot True or False. Defines if plots are generated with matplotlib. Default defined in defaultDict.
    @param split_berlin True or False. Defines if Berlin's disctricts are kept separated or get merged. Default defined in defaultDict.
    @param rep_date True or False. Defines if reporting date or reference date is taken into dataframe. Default defined in defaultDict.
    """

    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)
    filename = "CaseDataFull"

    complete = False
    filename = "CaseDataFull"
    path = os.path.join(directory + filename + ".csv")
    try:
        url = "https://media.githubusercontent.com/media/robert-koch-institut/" + \
            "SARS-CoV-2-Infektionen_in_Deutschland/main/Aktuell_Deutschland_SarsCov2_Infektionen.csv"
        df = gd.get_file(path, url, read_data, param_dict={})
        complete = check_for_completeness(df, merge_eisenach=True)
    except:
        pass
    if complete:
        if not read_data:
            # add column with state ids
            county_to_state_map = geoger.get_countyid_to_stateid_map(
                merge_berlin=False)
            df["IdBundesland"] = df["IdLandkreis"].map(county_to_state_map)
    else:
        # try another possibility if df was empty or incomplete
        print("Note: Case data is incomplete. Trying another source.")
        try:
            url="https://opendata.arcgis.com/datasets/66876b81065340a4a48710b062319336_0.csv"
            # if this file is encoded with utf-8 German umlauts are not displayed correctly because they take two bytes
            # utf_8_sig can identify those bytes as one sign and display it correctly
            df = gd.get_file(path, url, False, param_dict={"encoding": 'utf_8_sig'})
            complete = check_for_completeness(df, merge_eisenach=True)
        except:
            pass
        if not complete:
            print("Note: Case data is still incomplete. Trying a thrid source.")
            try:
                # If the data on github is not available we download the case data from rki from covid-19 datahub
                url = "https://npgeo-de.maps.arcgis.com/sharing/rest/content/" +\
                    "items/f10774f1c63e40168479a1feb6c7ca74/data"
                df = gd.get_file(path, url, False, param_dict={
                                 "encoding": 'utf_8_sig'})
                df.rename(columns={'FID': "OBJECTID"}, inplace=True)
                complete = check_for_completeness(df, merge_eisenach=True)
            except:
                pass
        if not complete:
            raise FileNotFoundError('Sorry')
        
        # drop columns that do not exist in data from github
        df = df.drop(["Altersgruppe2", "Datenstand", "OBJECTID",
                    "Bundesland", "Landkreis"], axis=1)

    df = df.convert_dtypes()

    # output data to not always download it
    if not no_raw:
        gd.write_dataframe(df, directory, filename, "json")

    # store dict values in parameter to not always call dict itself
    Altersgruppe = dd.GerEng['Altersgruppe']
    Geschlecht = dd.GerEng['Geschlecht']
    AnzahlFall = dd.GerEng['AnzahlFall']
    AnzahlGenesen = dd.GerEng['AnzahlGenesen']
    AnzahlTodesfall = dd.GerEng['AnzahlTodesfall']
    IdBundesland = dd.GerEng['IdBundesland']
    IdLandkreis = dd.GerEng['IdLandkreis']

    # translate column gender from German to English and standardize
    df.loc[df.Geschlecht == 'unbekannt', [
        'Geschlecht']] = dd.GerEng['unbekannt']
    df.loc[df.Geschlecht == 'W', ['Geschlecht']] = dd.GerEng['W']
    df.loc[df.Geschlecht == 'M', ['Geschlecht']] = dd.GerEng['M']
    df.loc[df.Altersgruppe == 'unbekannt', [
        'Altersgruppe']] = dd.GerEng['unbekannt']

    # change names of columns
    df.rename(dd.GerEng, axis=1, inplace=True)

    # Add column 'Date' with Date
    # = reporting date if rep_date is set
    # = reference date (date of disease onset) if IstErkrankungsbeginn = 1 else
    #       take Meldedatum (reporting date)
    if rep_date:
        df[dd.EngEng['date']] = df['Meldedatum'].astype('object')
    else:
        df[dd.EngEng['date']] = np.where(
            df['IstErkrankungsbeginn'] == 1, df['Refdatum'],
            df['Meldedatum'])

    df[dd.EngEng['date']] = pd.to_datetime(df[dd.EngEng['date']], format="%Y-%m-%d")

    # Date is either Refdatum or Meldedatum after column
    # 'IstErkrankungsbeginn' has been added. See also rep_date option.
    dateToUse = dd.EngEng['date']
    df.sort_values(dateToUse, inplace=True)

    # Manipulate data to get rid of conditions: df.NeuerFall >= 0, df.NeuerTodesfall >= 0, df.NeuGenesen >=0
    # There might be a better way
    df.loc[df.NeuerFall < 0, [AnzahlFall]] = 0
    df.loc[df.NeuerTodesfall < 0, [AnzahlTodesfall]] = 0
    df.loc[df.NeuGenesen < 0, [AnzahlGenesen]] = 0

    # get rid of unnecessary columns
    df = df.drop(['NeuerFall', 'NeuerTodesfall', 'NeuGenesen',
                  "IstErkrankungsbeginn", "Meldedatum", "Refdatum"], axis=1)

    print("Available columns:", df.columns)

    ######## Data for whole Germany all ages ##########

    # NeuerFall: Infected (incl. recovered) over "dateToUse":

    # make sum for one "dateToUse"
    gbNF = df.groupby(dateToUse).agg({AnzahlFall: sum})

    # make cumulative sum of "AnzahlFall" for "dateToUse"
    gbNF_cs = gbNF.cumsum()

    # output to json file
    # prefix for all files about cases data (apart from raw data)
    prefix = "cases_"
    filename = 'infected'
    if rep_date:
        filename_orig = filename + '_repdate'
    else:
        filename_orig = filename
    gd.write_dataframe(
        mdfs.extract_subframe_based_on_dates(
            gbNF_cs.reset_index(),
            start_date, end_date),
        directory, prefix + filename_orig, file_format)
    if impute_dates or moving_average > 0:
        gbNF_cs = mdfs.impute_and_reduce_df(
            gbNF_cs.reset_index(),
            {},
            ['Confirmed'],
            impute='forward', moving_average=moving_average,
            min_date=start_date, max_date=end_date)
        filename = gd.append_filename(filename, impute_dates, moving_average)
        if rep_date:
            filename = filename + '_repdate'
        gd.write_dataframe(
            mdfs.extract_subframe_based_on_dates(
                gbNF_cs, start_date, end_date),
            directory, prefix + filename, file_format)

    if make_plot:
        # make plot
        gbNF_cs.plot(title='COVID-19 infections', grid=True,
                     style='-o')
        plt.tight_layout()
        plt.show()

    # Dead over Date:
    gbNT = df.groupby(dateToUse).agg({AnzahlTodesfall: sum})
    gbNT = gbNT[gbNT[AnzahlTodesfall] != 0]
    gbNT_cs = gbNT.cumsum()

    # output
    filename = 'deaths'
    if rep_date:
        filename_orig = filename + '_repdate'
    else:
        filename_orig = filename
    gd.write_dataframe(
        mdfs.extract_subframe_based_on_dates(
            gbNT_cs.reset_index(),
            start_date, end_date),
        directory, prefix + filename_orig, file_format)
    if impute_dates or moving_average > 0:
        gbNT_cs = mdfs.impute_and_reduce_df(
            gbNT_cs.reset_index(),
            {},
            ['Deaths'],
            impute='forward', moving_average=moving_average,
            min_date=start_date, max_date=end_date)
        filename = gd.append_filename(filename, impute_dates, moving_average)
        if rep_date:
            filename = filename + '_repdate'
        gd.write_dataframe(
            mdfs.extract_subframe_based_on_dates(
                gbNT_cs.reset_index(),
                start_date, end_date),
            directory, prefix + filename, file_format)

    if make_plot:
        gbNT_cs.plot(title='COVID-19 deaths', grid=True,
                     style='-o')
        plt.tight_layout()
        plt.show()

        df.agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum}) \
            .plot(title='COVID-19 infections, deaths, recovered', grid=True,
                  kind='bar')
        plt.tight_layout()
        plt.show()

    gbNF = df.groupby(dateToUse).agg(
        {AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    gbNF_cs = gbNF.cumsum()

    filename = 'all_germany'
    if rep_date:
        filename_orig = filename + '_repdate'
    else:
        filename_orig = filename
    gd.write_dataframe(
        mdfs.extract_subframe_based_on_dates(
            gbNF_cs.reset_index(),
            start_date, end_date),
        directory, prefix + filename_orig, file_format)
    if impute_dates or moving_average > 0:
        gbNF_cs = mdfs.impute_and_reduce_df(
            gbNF_cs.reset_index(),
            {},
            ['Confirmed', 'Deaths', 'Recovered'],
            impute='forward', moving_average=moving_average,
            min_date=start_date, max_date=end_date)
        filename = gd.append_filename(filename, impute_dates, moving_average)
        if rep_date:
            filename = filename + '_repdate'
        gd.write_dataframe(
            mdfs.extract_subframe_based_on_dates(
                gbNF_cs, start_date, end_date),
            directory, prefix + filename, file_format)

    ############## Data for states all ages ################

    # NeuerFall: Infected (incl. recovered) over "dateToUse" for every state ("Bundesland"):
    groupby_list = [dateToUse, IdBundesland]
    gbNFst = df.groupby(groupby_list).agg({AnzahlFall: sum})
    gbNFst_cs = gbNFst.groupby(
        level=groupby_list.index(IdBundesland)).cumsum().reset_index()

    # output
    filename = 'infected_state'
    if rep_date:
        filename_orig = filename + '_repdate'
    else:
        filename_orig = filename
    gd.write_dataframe(
        mdfs.extract_subframe_based_on_dates(
            gbNFst_cs, start_date, end_date),
        directory, prefix + filename_orig, file_format)
    if impute_dates or moving_average > 0:
        gbNFst_cs = mdfs.impute_and_reduce_df(
            gbNFst_cs,
            {dd.EngEng["idState"]: [k for k, v in dd.State.items()]},
            ['Confirmed'],
            impute='forward', moving_average=moving_average,
            min_date=start_date, max_date=end_date)
        filename = gd.append_filename(filename, impute_dates, moving_average)
        if rep_date:
            filename = filename + '_repdate'
        gd.write_dataframe(
            mdfs.extract_subframe_based_on_dates(
                gbNFst_cs, start_date, end_date),
            directory, prefix + filename, file_format)

    # output nested json
    # gbNFst_cs.groupby(['IdBundesland'], as_index=False) \
    #            .apply(lambda x: x[[dateToUse,'AnzahlFall']].to_dict('r')) \
    #            .reset_index().rename(columns={0:'Dates'})\
    #            .to_json(directory + "gbNF_state_nested.json", orient='records')

    # infected (incl recovered), deaths and recovered together
    groupby_list = [dateToUse, IdBundesland]
    gbAllSt = df.groupby(groupby_list) \
        .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    gbAllSt_cs = gbAllSt.groupby(
        level=groupby_list.index(IdBundesland)).cumsum().reset_index()

    # output
    filename = 'all_state'
    if rep_date:
        filename_orig = filename + '_repdate'
    else:
        filename_orig = filename
    gd.write_dataframe(
        mdfs.extract_subframe_based_on_dates(
            gbAllSt_cs, start_date, end_date),
        directory, prefix + filename_orig, file_format)
    if impute_dates or moving_average > 0:
        gbAllSt_cs = mdfs.impute_and_reduce_df(
            gbAllSt_cs,
            {dd.EngEng["idState"]: [k for k, v in dd.State.items()]},
            ['Confirmed', 'Deaths', 'Recovered'],
            impute='forward', moving_average=moving_average,
            min_date=start_date, max_date=end_date)
        filename = gd.append_filename(filename, impute_dates, moving_average)
        if rep_date:
            filename = filename + '_repdate'
        gd.write_dataframe(
            mdfs.extract_subframe_based_on_dates(
                gbAllSt_cs, start_date, end_date),
            directory, prefix + filename, file_format)

    ############# Data for counties all ages ######################

    if not split_berlin:
        df = geoger.merge_df_counties(
            df, 11000, geoger.CountyMerging[11000],
            sorting=[dd.EngEng['date']],
            columns=[dd.EngEng['date'],
                     dd.EngEng['gender'],
                     dd.EngEng['idState'],
                     dd.EngEng['ageRKI']])

    # NeuerFall: Infected (incl. recovered) over "dateToUse" for every county ("Landkreis"):
    groupby_list = [dateToUse, IdLandkreis]
    gbNFc = df.groupby(groupby_list).agg({AnzahlFall: sum})
    gbNFc_cs = gbNFc.groupby(
        level=groupby_list.index(IdLandkreis)).cumsum().reset_index()

    # output
    if split_berlin:
        filename = 'infected_county_split_berlin'
    else:
        filename = 'infected_county'
    if rep_date:
        filename_orig = filename + '_repdate'
    else:
        filename_orig = filename
    gd.write_dataframe(mdfs.extract_subframe_based_on_dates(
        gbNFc_cs, start_date, end_date),
        directory, prefix + filename_orig, file_format)
    if impute_dates or moving_average > 0:
        gbNFc_cs = mdfs.impute_and_reduce_df(
            gbNFc_cs,
            {dd.EngEng["idCounty"]: sorted(set(df[dd.EngEng["idCounty"]].unique()))},
            ['Confirmed'],
            impute='forward', moving_average=moving_average,
            min_date=start_date, max_date=end_date)
        filename = gd.append_filename(filename, impute_dates, moving_average)
        if rep_date:
            filename = filename + '_repdate'
        gd.write_dataframe(
            mdfs.extract_subframe_based_on_dates(
                gbNFc_cs, start_date, end_date),
            directory, prefix + filename, file_format)

    # infected (incl recovered), deaths and recovered together
    groupby_list = [dateToUse, IdLandkreis]
    gbAllC = df.groupby(groupby_list).\
        agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    gbAllC_cs = gbAllC.groupby(
        level=groupby_list.index(IdLandkreis)).cumsum().reset_index()

    # output
    if split_berlin:
        filename = 'all_county_split_berlin'
    else:
        filename = 'all_county'
    if rep_date:
        filename_orig = filename + '_repdate'
    else:
        filename_orig = filename
    gd.write_dataframe(
        mdfs.extract_subframe_based_on_dates(
            gbAllC_cs, start_date, end_date),
        directory, prefix + filename_orig, file_format)
    if impute_dates or moving_average > 0:
        gbAllC_cs = mdfs.impute_and_reduce_df(
            gbAllC_cs,
            {dd.EngEng["idCounty"]: sorted(set(df[dd.EngEng["idCounty"]].unique()))},
            ['Confirmed', 'Deaths', 'Recovered'],
            impute='forward', moving_average=moving_average,
            min_date=start_date, max_date=end_date)
        filename = gd.append_filename(filename, impute_dates, moving_average)
        if rep_date:
            filename = filename + '_repdate'
        gd.write_dataframe(
            mdfs.extract_subframe_based_on_dates(
                gbAllC_cs, start_date, end_date),
            directory, prefix + filename, file_format)

    ######### Data whole Germany different gender ##################

    # infected (incl recovered), deaths and recovered together
    groupby_list = [dateToUse, Geschlecht]
    gbAllG = df.groupby(groupby_list) \
        .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})

    gbAllG_cs = gbAllG.groupby(
        level=groupby_list.index(Geschlecht)).cumsum().reset_index()

    # output
    filename = 'all_gender'
    if rep_date:
        filename_orig = filename + '_repdate'
    else:
        filename_orig = filename
    gd.write_dataframe(
        mdfs.extract_subframe_based_on_dates(
            gbAllG_cs, start_date, end_date),
        directory, prefix + filename_orig, file_format)
    if impute_dates or moving_average > 0:
        gbAllG_cs = mdfs.impute_and_reduce_df(
            gbAllG_cs,
            {dd.EngEng["gender"]: list(df[dd.EngEng["gender"]].unique())},
            ['Confirmed', 'Deaths', 'Recovered'],
            impute='forward', moving_average=moving_average,
            min_date=start_date, max_date=end_date)
        filename = gd.append_filename(filename, impute_dates, moving_average)
        if rep_date:
            filename = filename + '_repdate'
        gd.write_dataframe(
            mdfs.extract_subframe_based_on_dates(
                gbAllG_cs, start_date, end_date),
            directory, prefix + filename, file_format)

    if make_plot:
        df.groupby(Geschlecht) \
            .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum}) \
            .plot(title='COVID-19 infections, deaths, recovered', grid=True,
                  kind='bar')
        plt.tight_layout()
        plt.show()

    ############################# Gender and State ######################################################

    # infected (incl recovered), deaths and recovered together
    groupby_list = [dateToUse, IdBundesland, Geschlecht]
    gbAllGState = df.groupby(groupby_list) \
        .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    gbAllGState_cs = gbAllGState.groupby(
        level=[groupby_list.index(IdBundesland),
               groupby_list.index(Geschlecht)]).cumsum().reset_index()

    # output
    filename = 'all_state_gender'
    if rep_date:
        filename_orig = filename + '_repdate'
    else:
        filename_orig = filename
    gd.write_dataframe(
        mdfs.extract_subframe_based_on_dates(
            gbAllGState_cs, start_date, end_date),
        directory, prefix + filename_orig, file_format)
    if impute_dates or moving_average > 0:
        gbAllGState_cs = mdfs.impute_and_reduce_df(
            gbAllGState_cs,
            {dd.EngEng["idState"]: geoger.get_state_ids(),
             dd.EngEng["gender"]: list(df[dd.EngEng["gender"]].unique())},
            ['Confirmed', 'Deaths', 'Recovered'],
            impute='forward', moving_average=moving_average,
            min_date=start_date, max_date=end_date)
        filename = gd.append_filename(filename, impute_dates, moving_average)
        if rep_date:
            filename = filename + '_repdate'
        gd.write_dataframe(
            mdfs.extract_subframe_based_on_dates(
                gbAllGState_cs, start_date, end_date),
            directory, prefix + filename, file_format)

    ############# Gender and County #####################
    groupby_list = [dateToUse, IdLandkreis, Geschlecht]
    gbAllGCounty = df.groupby(groupby_list) \
        .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    gbAllGCounty_cs = gbAllGCounty.groupby(
        level=[groupby_list.index(IdLandkreis),
               groupby_list.index(Geschlecht)]).cumsum().reset_index()

    # output
    if split_berlin:
        filename = 'all_county_gender_split_berlin'
    else:
        filename = 'all_county_gender'
    if rep_date:
        filename_orig = filename + '_repdate'
    else:
        filename_orig = filename
    gd.write_dataframe(
        mdfs.extract_subframe_based_on_dates(
            gbAllGCounty_cs, start_date, end_date),
        directory, prefix + filename_orig, file_format)
    if impute_dates or moving_average > 0:
        gbAllGCounty_cs = mdfs.impute_and_reduce_df(
            gbAllGCounty_cs,
            {dd.EngEng["idCounty"]: sorted(set(df[dd.EngEng["idCounty"]].unique())),
             dd.EngEng["gender"]: list(df[dd.EngEng["gender"]].unique())},
            ['Confirmed', 'Deaths', 'Recovered'],
            impute='forward', moving_average=moving_average,
            min_date=start_date, max_date=end_date)
        filename = gd.append_filename(filename, impute_dates, moving_average)
        if rep_date:
            filename = filename + '_repdate'
        gd.write_dataframe(
            mdfs.extract_subframe_based_on_dates(
                gbAllGCounty_cs, start_date, end_date),
            directory, prefix + filename, file_format)

    ######### Data whole Germany different ages ####################

    # infected (incl recovered), deaths and recovered together
    groupby_list = [dateToUse, Altersgruppe]
    gbAllA = df.groupby(groupby_list) \
        .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    gbAllA_cs = gbAllA.groupby(
        level=groupby_list.index(Altersgruppe)).cumsum().reset_index()

    # output
    filename = 'all_age'
    if rep_date:
        filename_orig = filename + '_repdate'
    else:
        filename_orig = filename
    gd.write_dataframe(
        mdfs.extract_subframe_based_on_dates(
            gbAllA_cs, start_date, end_date),
        directory, prefix + filename_orig, file_format)
    if impute_dates or moving_average > 0:
        gbAllA_cs = mdfs.impute_and_reduce_df(
            gbAllA_cs,
            {dd.EngEng["ageRKI"]: sorted(set(df[dd.EngEng["ageRKI"]].unique()))},
            ['Confirmed', 'Deaths', 'Recovered'],
            impute='forward', moving_average=moving_average,
            min_date=start_date, max_date=end_date)
        filename = gd.append_filename(filename, impute_dates, moving_average)
        if rep_date:
            filename = filename + '_repdate'
        gd.write_dataframe(
            mdfs.extract_subframe_based_on_dates(
                gbAllA_cs, start_date, end_date),
            directory, prefix + filename, file_format)

    if make_plot:
        df.groupby(Altersgruppe).agg(
            {AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum}).plot(
            title='COVID-19 infections, deaths, recovered for diff ages',
            grid=True, kind='bar')
        plt.tight_layout()
        plt.show()

        # Dead by "Altersgruppe":
        gbNTAG = df.groupby(Altersgruppe).agg({AnzahlTodesfall: sum})

        gbNTAG.plot(title='COVID-19 deaths', grid=True,
                    kind='bar')
        plt.tight_layout()
        plt.show()

    ############################# Age and State ######################################################

    ##### Age_RKI #####

    # infected (incl recovered), deaths and recovered together
    groupby_list = [dateToUse, IdBundesland, Altersgruppe]
    gbAllAgeState = df.groupby(groupby_list) \
        .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    gbAllAgeState_cs = gbAllAgeState.groupby(
        level=[groupby_list.index(IdBundesland),
               groupby_list.index(Altersgruppe)]).cumsum().reset_index()

    # output
    filename = 'all_state_age'
    if rep_date:
        filename_orig = filename + '_repdate'
    else:
        filename_orig = filename
    gd.write_dataframe(
        mdfs.extract_subframe_based_on_dates(
            gbAllAgeState_cs, start_date, end_date),
        directory, prefix + filename_orig, file_format)
    if impute_dates or moving_average > 0:
        gbAllAgeState_cs = mdfs.impute_and_reduce_df(
            gbAllAgeState_cs,
            {dd.EngEng["idState"]: geoger.get_state_ids(),
             dd.EngEng["ageRKI"]: sorted(set(df[dd.EngEng["ageRKI"]].unique()))},
            ['Confirmed', 'Deaths', 'Recovered'],
            impute='forward', moving_average=moving_average,
            min_date=start_date, max_date=end_date)
        filename = gd.append_filename(filename, impute_dates, moving_average)
        if rep_date:
            filename = filename + '_repdate'
        gd.write_dataframe(
            mdfs.extract_subframe_based_on_dates(
                gbAllAgeState_cs, start_date, end_date),
            directory, prefix + filename, file_format)

    ############# Age and County #####################
    groupby_list = [dateToUse, IdLandkreis, Altersgruppe]
    gbAllAgeCounty = df.groupby(groupby_list) \
        .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
    gbAllAgeCounty_cs = gbAllAgeCounty.groupby(
        level=[groupby_list.index(IdLandkreis),
               groupby_list.index(Altersgruppe)]).cumsum().reset_index()

    # output
    if split_berlin:
        filename = 'all_county_age_split_berlin'
    else:
        filename = 'all_county_age'
    if rep_date:
        filename_orig = filename + '_repdate'
    else:
        filename_orig = filename
    gd.write_dataframe(
        mdfs.extract_subframe_based_on_dates(
            gbAllAgeCounty_cs, start_date, end_date),
        directory, prefix + filename_orig, file_format)
    if impute_dates or moving_average > 0:
        gbAllAgeCounty_cs = mdfs.impute_and_reduce_df(
            gbAllAgeCounty_cs,
            {dd.EngEng["idCounty"]: sorted(set(df[dd.EngEng["idCounty"]].unique())),
             dd.EngEng["ageRKI"]: sorted(set(df[dd.EngEng["ageRKI"]].unique()))},
            ['Confirmed', 'Deaths', 'Recovered'],
            impute='forward', moving_average=moving_average,
            min_date=start_date, max_date=end_date)
        filename = gd.append_filename(filename, impute_dates, moving_average)
        if rep_date:
            filename = filename + '_repdate'
        gd.write_dataframe(
            mdfs.extract_subframe_based_on_dates(
                gbAllAgeCounty_cs, start_date, end_date),
            directory, prefix + filename, file_format)


def main():
    """! Main program entry."""

    arg_dict = gd.cli("cases")
    get_case_data(**arg_dict)


if __name__ == "__main__":
    main()
