#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Kathrin Rack, Wadim Koslow, Patrick Lenz
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
:strong:`getPopulationData.py`

Downloads data about population statistic

"""
import warnings
import requests
import os
import io

import numpy as np
import pandas as pd

from memilio.epidata import defaultDict as dd
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getDataIntoPandasDataFrame as gd

# activate CoW for more predictable behaviour of pandas DataFrames
pd.options.mode.copy_on_write = True


def read_population_data(ref_year):
    """ Reads Population data from regionalstatistik.de

    A request is made to regionalstatistik.de and the StringIO is read in as a csv into the dataframe format.

    :param ref_year: Default: None or year (jjjj) convertible to str. Reference year.
    :returns: DataFrame

    """
    if ref_year is not None:
        try:
            download_url = 'https://www.regionalstatistik.de/genesis/online?operation=download&code=12411-02-03-4&option=csv&zeiten=' + \
                str(ref_year)
            req = requests.get(download_url)
            df_pop_raw = pd.read_csv(io.StringIO(req.text), sep=';', header=5)
        except pd.errors.ParserError:
            gd.default_print(
                'Warning', 'Data for year ' + str(ref_year) +
                ' is not available; downloading newest data instead.')
            ref_year = None
    if ref_year is None:
        download_url = 'https://www.regionalstatistik.de/genesis/online?operation=download&code=12411-02-03-4&option=csv'
        req = requests.get(download_url)
        df_pop_raw = pd.read_csv(io.StringIO(req.text), sep=';', header=5)

    return df_pop_raw, ref_year


def export_population_dataframe(
        df_pop: pd.DataFrame, directory: str, file_format: str,
        merge_eisenach: bool, ref_year):
    """ Writes population dataframe into directory with new column names and age groups

    :param df_pop: Population data DataFrame to be exported pd.DataFrame
    :param directory: Directory where data is written to. str
    :param file_format: File format which is used for writing the data. str
    :param merge_eisenach: Defines whether the counties 'Wartburgkreis'
        and 'Eisenach' are listed separately or
        combined as one entity 'Wartburgkreis'. bool
    :param ref_year: None or year (jjjj) convertible to str. Reference year.
    :returns: exported DataFrame

    """

    new_cols = [
        dd.EngEng['idCounty'],
        dd.EngEng['population'],
        '<3 years', '3-5 years', '6-14 years', '15-17 years',
        '18-24 years', '25-29 years', '30-39 years', '40-49 years',
        '50-64 years', '65-74 years', '>74 years']
    df_pop_export = pd.DataFrame(columns=new_cols)
    df_pop_export[new_cols[0]] = df_pop[dd.EngEng['idCounty']]
    # <3 and 3-5
    df_pop_export[new_cols[2:4]] = df_pop[df_pop.columns[2:4]]
    # 6-14
    df_pop_export[new_cols[4]] = \
        df_pop[df_pop.columns[4:6]].sum(axis=1)
    # 15-17
    df_pop_export[new_cols[5]] = df_pop[df_pop.columns[6]]
    # 18-24
    df_pop_export[new_cols[6]] = \
        df_pop[df_pop.columns[7:9]].sum(axis=1)
    # 25-29
    df_pop_export[new_cols[7]] = df_pop[df_pop.columns[9]]
    # 30-39
    df_pop_export[new_cols[8]] = \
        df_pop[df_pop.columns[10:12]].sum(axis=1)
    # 40-49
    df_pop_export[new_cols[9]] = \
        df_pop[df_pop.columns[12:14]].sum(axis=1)
    # 50-64
    df_pop_export[new_cols[10]] = \
        df_pop[df_pop.columns[14:17]].sum(axis=1)
    # 65-74
    df_pop_export[new_cols[11]] = df_pop[df_pop.columns[17]]
    # >74
    df_pop_export[new_cols[12]] = df_pop[df_pop.columns[18]]

    df_pop_export[dd.EngEng['population']
                  ] = df_pop_export.iloc[:, 2:].sum(axis=1)

    # merge eisenach if no data available
    if '16056' in df_pop_export[dd.EngEng['idCounty']].values:
        if df_pop_export[df_pop_export[dd.EngEng['idCounty']] == '16056'][dd.EngEng['population']].values[0] == 0:
            df_pop_export = geoger.merge_df_counties_all(
                df_pop_export, sorting=[dd.EngEng["idCounty"]],
                columns=dd.EngEng["idCounty"])

    gd.check_dir(directory)

    if ref_year is None:
        filename = 'county_current_population'
    else:
        filename = 'county_' + str(ref_year) + '_population'

    if len(df_pop_export) == 401:
        filename = filename + '_dim401'

    # Merge Eisenach and Wartburgkreis
    df_pop_export = geoger.merge_df_counties_all(
        df_pop_export, sorting=[dd.EngEng["idCounty"]],
        columns=dd.EngEng["idCounty"])

    gd.write_dataframe(df_pop_export, directory, filename, file_format)
    gd.write_dataframe(aggregate_to_state_level(df_pop_export),
                       directory, filename + '_states', file_format)
    gd.write_dataframe(aggregate_to_country_level(df_pop_export),
                       directory, filename + '_germany', file_format)

    return df_pop_export


def assign_population_data(df_pop_raw, counties, age_cols, idCounty_idx):
    """ Assigns population data of all counties of old dataframe in new created dataframe

    In df_pop_raw there might be additional information like federal states,
    governing regions etc. which is not necessary for the dataframe.
    Also checks for incomplete data.

    :param df_pop_raw: Raw Population DataFrame read from regionalstatistik.de
    :param counties: List of counties to be assigned in new DataFrame
    :param age_cols: Age groups in old DataFrame
    :param idCountyidx: indexes in old DataFrame where data of corresponding county starts
    :param idCounty_idx: 
    :returns: new DataFrame

    """

    new_cols = {dd.EngEng['idCounty']: counties[:, 1],
                dd.EngEng['county']: counties[:, 0]}

    # number of age_cols
    num_age_groups = len(age_cols)

    # add age_cols with zero initialization to new_cols
    new_cols.update({age: 0 for age in age_cols})

    df_pop = pd.DataFrame(new_cols)

    empty_data = ['.', '-']
    for start_idx in idCounty_idx:

        county_id = df_pop_raw.loc[start_idx, dd.EngEng['idCounty']]

        # check for empty rows
        if df_pop_raw.loc[start_idx: start_idx + num_age_groups - 1,
                          dd.EngEng['number']].values.any() in empty_data:
            if not df_pop_raw.loc[start_idx: start_idx + num_age_groups - 1,
                                  dd.EngEng['number']].values.all() in empty_data:
                raise gd.DataError(
                    'Error. Partially incomplete data for county ' +
                    county_id)

        # county information needed
        elif county_id in counties[:, 1]:
            # direct assignment of population data found
            df_pop.loc[df_pop[dd.EngEng['idCounty']] == df_pop_raw.loc
                       [start_idx, dd.EngEng['idCounty']],
                       age_cols] = df_pop_raw.loc[start_idx: start_idx +
                                                  num_age_groups - 1, dd.EngEng['number']].values.astype(int)
        # Berlin and Hamburg
        elif county_id + '000' in counties[:, 1]:
            # direct assignment of population data found
            df_pop.loc[df_pop[dd.EngEng['idCounty']] == df_pop_raw.loc
                       [start_idx, dd.EngEng['idCounty']] + '000',
                       age_cols] = df_pop_raw.loc[start_idx: start_idx + num_age_groups - 1, dd.EngEng
                                                  ['number']].values.astype(int)

        # additional information for local entities not needed
        elif county_id in ['03241001', '05334002', '10041100', "11001001", "11002002", "11003003",
                           "11004004", "11005005", "11006006", "11007007", "11008008", "11009009",
                           "11010010", "11011011", "11012012"]:
            pass
        # Germany, federal states, and governing regions
        elif len(county_id) < 5:
            pass
        else:
            raise gd.DataError(
                'No data for ' + df_pop_raw.loc
                [start_idx, dd.EngEng['idCounty']] +
                ' County ID in input population data '
                'found which could not be assigned.')

    return df_pop


def test_total_population(df_pop, age_cols):
    """ Tests if total population matches expectation

    :param df_pop: Population Dataframe with all counties
    :param age_cols: All age groups in DataFram

    """

    total_sum_expect = 84e6
    total_sum = df_pop[age_cols].sum().sum()

    if not isinstance(total_sum, (int, np.integer)):
        raise gd.DataError('Unexpected dtypes in Population Data.')
    # check if total population is +-5% accurate to 2024 population
    if (total_sum > 1.05*total_sum_expect) or (total_sum < 0.95*total_sum_expect):
        gd.default_print(
            'Warning', 'Total Population does not match expectation.')


def fetch_population_data(read_data: bool = dd.defaultDict['read_data'],
                          out_folder: str = dd.defaultDict['out_folder'],
                          ref_year=None,
                          **kwargs
                          ) -> pd.DataFrame:
    """ Downloads or reads the population data.
    If it does not already exist, the folder Germany is generated in the given out_folder.
    If read_data == True and the file "FullData_population.json" exists, the data is read form this file
    and stored in a pandas dataframe. If read_data = True and the file does not exist the program is stopped.
    The downloaded dataframe is written to the file "FullData_population".

    :param read_data: False or True. Defines if data is read from file or
        downloaded. Default defined in defaultDict. (Default value = dd.defaultDict['read_data'])
    :param out_folder: Path to folder where data is written in folder
        out_folder/Germany. Default defined in defaultDict. (Default value = dd.defaultDict['out_folder'])
    :param ref_year: (Default: None) or year (jjjj) convertible to str. Reference year.
    :param **kwargs: 
    :returns: DataFrame with adjusted population data for all ages to current level.

    """
    conf = gd.Conf(out_folder, **kwargs)
    out_folder = conf.path_to_use

    if read_data == True:
        gd.default_print(
            'Warning',
            'Read_data is not supportet for getPopulationData.py. Setting read_data = False')
        read_data = False

    directory = os.path.join(out_folder, 'Germany', 'pydata')
    gd.check_dir(directory)

    df_pop_raw, ref_year = read_population_data(ref_year)

    return df_pop_raw, ref_year


def preprocess_population_data(df_pop_raw: pd.DataFrame,
                               merge_eisenach: bool = True,
                               ) -> pd.DataFrame:
    """ Processing of the downloaded data
        * the columns are renamed to English and the state and county names are added.

    :param df_pop_raw: pd.DataFrame. A Dataframe containing input population data
    :param merge_eisenach: Default: True] or False. Defines whether the
     counties 'Wartburgkreis' and 'Eisenach' are listed separately or
     combined as one entity 'Wartburgkreis'. (Default value = True)
    :returns: df pd.DataFrame. Processed population data

    """
    column_names = list(df_pop_raw.columns)
    # rename columns
    if len(column_names) == 7:
        rename_columns = {
            column_names[0]: dd.EngEng['date'],
            column_names[1]: dd.EngEng['idCounty'],
            column_names[2]: dd.EngEng['county'],
            column_names[3]: dd.EngEng['ageRKI'],
            column_names[4]: dd.EngEng['number'],
            column_names[5]: dd.EngEng['male'],
            column_names[6]: dd.EngEng['female']
        }
    else:
        rename_columns = {
            column_names[0]: dd.EngEng['idCounty'],
            column_names[1]: dd.EngEng['county'],
            column_names[2]: dd.EngEng['ageRKI'],
            column_names[3]: dd.EngEng['number'],
            column_names[4]: dd.EngEng['male'],
            column_names[5]: dd.EngEng['female']
        }
    df_pop_raw.rename(columns=rename_columns, inplace=True)

    # remove date and explanation rows at end of table. If length of column_names is 7,
    # the explanation is in the date column, otherwise in the idCounty column
    if len(column_names) == 7:
        df_pop_raw = df_pop_raw[:np.where(df_pop_raw[dd.EngEng['date']].str.contains(
            '__') == True)[0][0]].reset_index(drop=True)
    else:
        df_pop_raw = df_pop_raw[:np.where(df_pop_raw[dd.EngEng['idCounty']].str.contains(
            '__') == True)[0][0]].reset_index(drop=True)
    # get indices of counties first lines
    idCounty_idx = df_pop_raw.groupby(
        dd.EngEng['idCounty']).head(1).index.tolist()

    # Delete column dd.EngEng[‘date’] as it was added to the data later and is not needed
    if len(column_names) == 7:
        df_pop_raw.drop(columns=[dd.EngEng['date']], inplace=True)

    # read county list and create output data frame
    counties = np.array(geoger.get_county_names_and_ids(
        merge_berlin=True, merge_eisenach=merge_eisenach, zfill=True))
    age_cols = df_pop_raw.loc[
        idCounty_idx[0]: idCounty_idx[1] - 2,
        dd.EngEng['ageRKI']].values.copy()
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

    df_pop = assign_population_data(
        df_pop_raw, counties, age_cols, idCounty_idx)
    test_total_population(df_pop, age_cols)
    return df_pop


def write_population_data(df_pop: pd.DataFrame,
                          out_folder: str = dd.defaultDict['out_folder'],
                          file_format: str = dd.defaultDict['file_format'],
                          merge_eisenach: bool = True,
                          ref_year=None
                          ) -> None or pd.DataFrame:
    """ Write the population data into json files
    Three kinds of structuring of the data are done.
    We obtain the chronological sequence of ICU and ICU_ventilated
    stored in the files "county_population".json", "state_population.json" and "germany_population.json"
    for counties, states and whole Germany, respectively.

    :param df_pop: pd.DataFrame. A Dataframe containing processed population data
    :param file_format: str. File format which is used for writing the data. Default defined in defaultDict. (Default value = dd.defaultDict['file_format'])
    :param out_folder: str. Folder where data is written to. Default defined in defaultDict. (Default value = dd.defaultDict['out_folder'])
    :param merge_eisenach: bool, Default: True. Defines whether the
        counties 'Wartburgkreis' and 'Eisenach' are listed separately or
        combined as one entity 'Wartburgkreis'. 
    :param ref_year: Default: None or year (jjjj) convertible to str. Reference year.
    :returns: None

    """
    directory = os.path.join(out_folder, 'Germany', 'pydata')
    df_pop_export = export_population_dataframe(
        df_pop, directory, file_format, merge_eisenach, ref_year)
    return df_pop_export


def get_population_data(read_data: bool = dd.defaultDict['read_data'],
                        file_format: str = dd.defaultDict['file_format'],
                        out_folder: str = dd.defaultDict['out_folder'],
                        merge_eisenach: bool = True,
                        ref_year=None,
                        **kwargs
                        ):
    """ Download age-stratified population data for the German counties.

    The data we use is:
    Official 'Bevölkerungsfortschreibung' 12411-02-03-4:
    'Bevölkerung nach Geschlecht und Altersgruppen (17)'
    of regionalstatistik.de.
    ATTENTION: The raw file cannot be downloaded
    automatically by our scripts without an Genesis Online account. In order to
    work on this dataset, please enter your username and password or manually download it from:

    https://www.regionalstatistik.de/genesis/online -> "1: Gebiet, Bevölkerung,
    Arbeitsmarkt, Wahlen" -> "12: Bevölkerung" -> "12411 Fortschreibung des
    Bevölkerungsstandes" ->  "12411-02-03-4: Bevölkerung nach Geschlecht und
    Altersgruppen (17) - Stichtag 31.12. - regionale Tiefe: Kreise und
    krfr. Städte".

    Download the xlsx or csv file and put it under dd.defaultDict['out_folder'],
    this normally is Memilio/data/pydata/Germany.
    The folders 'pydata/Germany' have to be created if they do not exist yet.
    Then this script can be run.

    :param read_data: False or True. Defines if data is read from file or
        downloaded. Default defined in defaultDict. (Default value = dd.defaultDict['read_data'])
    :param file_format: File format which is used for writing the data.
        Default defined in defaultDict. (Default value = dd.defaultDict['file_format'])
    :param out_folder: Path to folder where data is written in folder
        out_folder/Germany. Default defined in defaultDict. (Default value = dd.defaultDict['out_folder'])
    :param merge_eisenach: bool, Default: True]. Defines whether the
        counties 'Wartburgkreis' and 'Eisenach' are listed separately or
        combined as one entity 'Wartburgkreis'.
    :param ref_year: Default: None] or year (jjjj) convertible to str. Reference year.
    :param username: str. Username to sign in at regionalstatistik.de.
    :param password: str. Password to sign in at regionalstatistik.de.
    :param **kwargs: 
    :returns: DataFrame with adjusted population data for all ages to current level.

    """
    raw_df, ref_year = fetch_population_data(
        read_data=read_data,
        out_folder=out_folder,
        file_format=file_format,
        ref_year=ref_year,
        **kwargs
    )
    preprocess_df = preprocess_population_data(
        df_pop_raw=raw_df,
        merge_eisenach=merge_eisenach
    )
    df_pop_export = write_population_data(
        df_pop=preprocess_df,
        file_format=file_format,
        out_folder=out_folder,
        merge_eisenach=True,
        ref_year=ref_year
    )
    return df_pop_export


def aggregate_to_state_level(df_pop: pd.DataFrame):

    countyIDtostateID = geoger.get_countyid_to_stateid_map()

    df_pop['ID_State'] = df_pop[dd.EngEng['idCounty']].map(countyIDtostateID)
    df_pop = df_pop.drop(
        columns='ID_County').groupby(
        'ID_State', as_index=True).sum()
    df_pop['ID_State'] = df_pop.index
    return df_pop


def aggregate_to_country_level(df_pop: pd.DataFrame):

    df_pop['ID_Country'] = 0
    df_pop = df_pop.drop(
        columns=['ID_County', 'ID_State']).groupby(
        'ID_Country', as_index=True).sum()
    df_pop['ID_Country'] = df_pop.index
    return df_pop


def main():
    """ Main program entry."""

    arg_dict = gd.cli("population")
    get_population_data(**arg_dict)


if __name__ == "__main__":
    main()
