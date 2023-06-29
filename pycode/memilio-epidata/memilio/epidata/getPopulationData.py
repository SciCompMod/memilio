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
@file getPopulationData.py

@brief Downloads data about population statistic

"""

import os

import numpy as np
import pandas as pd

from memilio.epidata import defaultDict as dd
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getDataIntoPandasDataFrame as gd


def get_population_data(read_data=dd.defaultDict['read_data'],
                        file_format=dd.defaultDict['file_format'],
                        out_folder=dd.defaultDict['out_folder'],
                        no_raw=dd.defaultDict['no_raw'],
                        split_gender=False,
                        merge_eisenach=True):
    """! Download age-stratified population data for the German counties.

    There are two different data sources that can be transformed in a simple
    data file with age-resolved information per county and on which our other
    tools can continue to operate on.

    1.) Official 'Bevölkerungsfortschreibung' 12411-02-03-4:
    'Bevölkerung nach Geschlecht und Altersgruppen (17)' 
    of regionalstatistik.de. 
    ATTENTION: The raw file cannot be downloaded 
    automatically by our scripts without an Genesis Online account. In order to
    work on this dataset, please manually download it from:

    https://www.regionalstatistik.de/genesis/online -> "1: Gebiet, Bevölkerung,
    Arbeitsmarkt, Wahlen" -> "12: Bevölkerung" -> "12411 Fortschreibung des
    Bevölkerungsstandes" ->  "12411-02-03-4: Bevölkerung nach Geschlecht und 
    Altersgruppen (17) - Stichtag 31.12. - regionale Tiefe: Kreise und
    krfr. Städte". 

    Download the xlsx file and put it under dd.defaultDict['out_folder'], 
    this normally is Memilio/data/pydata/Germany. 
    The folders 'pydata/Germany' have to be created if they do not exist yet. 
    Then this script can be run.
    TODO: The following parameters have no effect for this source.

    2.) Combination of data from Federal Statistical Office of Germany and 
    Zensus2011 to compute very simple approximations of age-stratified number
    of inhabitants per county. Population from the Zensus data of all age groups is scaled to the
    total population of our more recent migration data by a factor which
    represents the relative increase/decrease in population size
    between 2011 and 2019 for each county"
    This data can either be downloaded automatically or read from
    "out_folder"/Germany/ if it was downloaded before.

    @param read_data False or True. Defines if data is read from file or
        downloaded. Default defined in defaultDict.
    @param file_format File format which is used for writing the data.
        Default defined in defaultDict.
    @param out_folder Path to folder where data is written in folder
        out_folder/Germany. Default defined in defaultDict.
    @param no_raw True or False. Defines if unchanged raw data is written or
        not. Default defined in defaultDict.
    @param split_gender [Default: False] or True. Defines whether data is
        splitted by gender
    @param merge_eisenach [Default: True] or False. Defines whether the
        counties 'Wartburgkreis' and 'Eisenach' are listed separately or
        combined as one entity 'Wartburgkreis'.
    @return DataFrame with adjusted population data for all ages to current level.
    """
    directory = os.path.join(dd.defaultDict['out_folder'], 'Germany')
    filename = '12411-02-03-4'  # '12411-09-01-4-B'
    new_data_file = os.path.join(directory, filename+'.xlsx')

    df_pop_raw = gd.get_file(
        new_data_file, url='https://www.regionalstatistik.de/genesis/online?operation=ergebnistabelleDownload&levelindex=3&levelid=1688024651637&option=xlsx', read_data=False,
        param_dict={"engine": "openpyxl",
                    "sheet_name": filename, "header": 4},
        interactive=False, auth=('',''))
    column_names = list(df_pop_raw.columns)
    # rename columns
    rename_columns = {
        column_names[0]: dd.EngEng['idCounty'],
        column_names[1]: dd.EngEng['county'],
        column_names[2]: dd.EngEng['ageRKI'],
        column_names[3]: dd.EngEng['number'],
        column_names[4]: dd.EngEng['male'],
        column_names[5]: dd.EngEng['female'],
    }
    df_pop_raw.rename(columns=rename_columns, inplace=True)
    # remove date and explanation rows at end of table
    # -33 for '12411-09-01-4-B', -38 for '12411-02-03-4'
    df_pop_raw = df_pop_raw[:-38].reset_index(drop=True)
    # get indices of counties first lines
    idCounty_idx = list(df_pop_raw.index[
        df_pop_raw[dd.EngEng['idCounty']].isna() == False])

    # read county list and create output data frame
    counties = np.array(geoger.get_county_names_and_ids(
        merge_berlin=True, merge_eisenach=False, zfill=True))
    new_cols = {dd.EngEng['idCounty']: counties[:, 1],
                dd.EngEng['county']: counties[:, 0]}
    age_cols = df_pop_raw.loc[
        idCounty_idx[0]: idCounty_idx[1] - 2,  # -1 for '12411-09-01-4-B'
        dd.EngEng['ageRKI']].copy().values
    nage = len(age_cols)
    for i in range(len(age_cols)):
        if i == 0:
            upper_bound = str(int(age_cols[i][
                age_cols[i].index('unter ')+6:].split(' ')[0])-1)
            age_cols[i] = '0-' + upper_bound
        elif i == len(age_cols)-1:
            lower_bound = age_cols[i].split(' ')[0]
            age_cols[i] = lower_bound + '-99'
        else:
            lower_bound = age_cols[i].split(' ')[0]
            upper_bound = str(int(age_cols[i][
                age_cols[i].index('unter ')+6:].split(' ')[0])-1)
            age_cols[i] = lower_bound + '-' + upper_bound

    # add age_cols with zero initilization to new_cols
    new_cols.update({age: 0 for age in age_cols})
    df_pop = pd.DataFrame(new_cols)

    empty_data = ['.', '-']
    for i in idCounty_idx:
        # county information needed
        if df_pop_raw.loc[i, dd.EngEng['idCounty']] in counties[:, 1]:
            # direct assignment of population data found
            df_pop.loc[df_pop[dd.EngEng['idCounty']] == df_pop_raw.loc
                        [i, dd.EngEng['idCounty']],
                        age_cols] = df_pop_raw.loc[i: i + nage - 1, dd.EngEng
                                                    ['number']].values
        # Berlin and Hamburg
        elif df_pop_raw.loc[i, dd.EngEng['idCounty']] + '000' in counties[:, 1]:
            # direct assignment of population data found
            df_pop.loc[df_pop[dd.EngEng['idCounty']] == df_pop_raw.loc
                        [i, dd.EngEng['idCounty']] + '000',
                        age_cols] = df_pop_raw.loc[i: i + nage - 1, dd.EngEng
                                                    ['number']].values
        # empty rows
        elif df_pop_raw.loc[i: i + nage - 1,
                            dd.EngEng['number']].values.any() in empty_data:
            if not df_pop_raw.loc[i: i + nage - 1,
                                    dd.EngEng['number']].values.all() in empty_data:
                raise gd.DataError(
                    'Error. Partially incomplete data for county ' +
                    df_pop_raw.loc[i, dd.EngEng['idCounty']])
        # additional information for local entities not needed
        elif df_pop_raw.loc[i, dd.EngEng['idCounty']] in ['03241001', '05334002', '10041100']:
            pass
        # Germany, federal states, and governing regions
        elif len(df_pop_raw.loc[i, dd.EngEng['idCounty']]) < 5:
            pass
        else:
            print('no data for ' + df_pop_raw.loc
                    [i, dd.EngEng['idCounty']])
            raise gd.DataError(
                'Error. County ID in input population data'
                'found which could not be assigned.')

    if df_pop[age_cols].sum().sum() != 83155031:
        raise gd.DataError('Wrong total population size')

    new_cols = [
        dd.EngEng['idCounty'],
        dd.EngEng['population'],
        '<3 years', '3-5 years', '6-14 years', '15-17 years',
        '18-24 years', '25-29 years', '30-39 years', '40-49 years',
        '50-64 years', '65-74 years', '>74 years']
    df_pop_export = pd.DataFrame(columns=new_cols)
    df_pop_export[df_pop.columns[0]] = df_pop[dd.EngEng['idCounty']]
    # <3 and 3-5
    df_pop_export[df_pop_export.columns[2:4]] = df_pop[df_pop.columns[2:4]]
    # 6-14
    df_pop_export[df_pop_export.columns[4]] = \
        df_pop[df_pop.columns[4:6]].sum(axis=1)
    # 15-17
    df_pop_export[df_pop_export.columns[5]] = df_pop[df_pop.columns[6]]
    # 18-24
    df_pop_export[df_pop_export.columns[6]] = \
        df_pop[df_pop.columns[7:9]].sum(axis=1)
    # 25-29
    df_pop_export[df_pop_export.columns[7]] = df_pop[df_pop.columns[9]]
    # 30-39
    df_pop_export[df_pop_export.columns[8]] = \
        df_pop[df_pop.columns[10:12]].sum(axis=1)
    # 40-49
    df_pop_export[df_pop_export.columns[9]] = \
        df_pop[df_pop.columns[12:14]].sum(axis=1)
    # 50-64
    df_pop_export[df_pop_export.columns[10]] = \
        df_pop[df_pop.columns[14:17]].sum(axis=1)
    # 65-74
    df_pop_export[df_pop_export.columns[11]] = df_pop[df_pop.columns[17]]
    # >74
    df_pop_export[df_pop_export.columns[12]] = df_pop[df_pop.columns[18]]

    df_pop_export[dd.EngEng['population']
                    ] = df_pop_export.iloc[:, 2:].sum(axis=1)

    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)

    if len(df_pop_export) == 401:
        filename = 'county_current_population_dim401'
        gd.write_dataframe(df_pop_export, directory, filename, file_format)

    if len(df_pop_export) == 400 or merge_eisenach == True:
        filename = 'county_current_population'

        # Merge Eisenach and Wartburgkreis
        df_pop_export = geoger.merge_df_counties_all(
            df_pop_export, sorting=[dd.EngEng["idCounty"]],
            columns=dd.EngEng["idCounty"])

        gd.write_dataframe(df_pop_export, directory, filename, file_format)

    return df_pop_export

def main():
    """! Main program entry."""

    arg_dict = gd.cli("population")
    get_population_data(**arg_dict)


if __name__ == "__main__":
    main()
