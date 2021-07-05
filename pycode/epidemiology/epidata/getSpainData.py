"""
@file getSpainData.py

@brief Downloads data for Spain

############################################################################################################
#                                                                                                          #
#        IMPORTANT NOTE: WHEN USING THIS DATA, WE HAVE TO CITE https://github.com/datadista/datasets       #
#                                                                                                          #
#                                                                                                          #
#        DO NOT USE DATA FROM THE FOLLOWING REGIONS SINCE THE COLUMNS HOSPITALIZED AND ICU                 #
#        ARE NOT CORRECTLY SUMMED TO TOTAL NUMBERS ! THE SAME APPLIES TO ALL AGE DATA AT THE MOMENT !      #
#                                                                                                          #
#               HOSPITALIZED                                   ICU                                         #
#               Castilla La Mancha (until 2020-04-11)          Castilla La Mancha (hasta 2020-04-12)       #
#               Comunidad Valenciana (hasta 2020-04-08)        Castilla y León (hasta 2020-04-17)          #
#               Madrid (hasta 2020-04-26)                      Comunidad Valenciana (hasta 2020-04-08)     #
#               Castilla y León (hasta 2020-04-06)             Galicia (hasta 2020-04-29)                  #
#               Madrid (hasta 2020-04-26)                                                                  #
#                                                                                                          #
############################################################################################################
"""

import os
import sys
import pandas
import numpy as np

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd


def get_spain_data(read_data=dd.defaultDict['read_data'],
                   file_format=dd.defaultDict['file_format'],
                   out_folder=dd.defaultDict['out_folder'],
                   no_raw=dd.defaultDict['no_raw']):
    """! Downloads data from Spain covid-19 data

   Source: https://github.com/datadista/datasets/tree/master/COVID%2019

   The data is read in, either from the internet or from two json files
   (raw_spain_all_age.json, raw_spain_all_state.json), stored in an earlier run.
   If the data is read from the internet, before changing anything the data is stored in above mentioned files.
   The file is read in or stored at the folder "out_folder"/Spain/.
   To store and change the data we use pandas

   Working with the data
   - Columns are renamed to english as defined in defaultDict
   - Spanish names of states, specific characters are removed.
   - For data splitted for states, it can be chosen if only positive PCR tests are taken into account
   or also positive antibody tests [Default].
   To change it, a parameter in the source code has t be changed.
   - Following files are written
      - whole spain together stored in spain
      - spain split for age stored in spain_all_age
      - spain split into states stored in spain_all_state

   @param read_data False [Default] or True. Defines if data is read from file or downloaded.
   @param file_format File format which is used for writing the data. Default defined in defaultDict.
   @param out_folder Path to folder where data is written in folder out_folder/Spain.
   @param no_raw True or False [Default]. Defines if unchanged raw data is saved or not.
   """

    directory = os.path.join(out_folder, 'Spain/')
    gd.check_dir(directory)

    _get_spain_all_age(read_data, file_format, no_raw, directory)
    _get_spain_all_state(read_data, file_format, no_raw, directory)


def _get_spain_all_age(read_data, file_format, no_raw, directory):

    ages_file = 'raw_spain_all_age'

    data_str = 'nacional_covid19_rango_edad'
    api_url = 'https://raw.githubusercontent.com/datadista/datasets/master/COVID%2019/old_series/'
    exit_string = "Files \'nacional_covid19_rango_edad\' are not available online. Check URL."

    df_age = _read_or_load_data(read_data, no_raw, directory, ages_file, data_str, api_url, exit_string)

    # Manipulate data
    if not df_age.empty:
        # standardization of column titles from Spanish to English
        # the stupid character in front of 'fecha' is correct here. There is a bug in the original file.
        df_age.rename(dd.EsEng, axis=1, inplace=True)

        print("Available age columns:", df_age.columns)

        # translate column gender from Spanish to English and standardize
        gender = dd.EngEng['gender']
        age10 = dd.EngEng['age10']
        df_age.loc[df_age[gender] == 'ambos', [gender]] = dd.EngEng['both']
        df_age.loc[df_age[gender] == 'mujeres', [gender]] = dd.EngEng['female']
        df_age.loc[df_age[gender] == 'hombres', [gender]] = dd.EngEng['male']
        df_age.loc[df_age[age10] == '80 y +', [age10]] = dd.EngEng['80+']
        df_age.loc[df_age[age10] == '90 y +', [age10]] = dd.EngEng['90+']
        df_age.loc[df_age[age10] == 'Total', [age10]] = dd.EngEng['all']

        # Correct Timestamps:
        date = dd.EngEng['date']
        df_age[date] = df_age[date].astype('datetime64[ns]')

    # only consider men AND women (through information on gender away)
    df_age = df_age.loc[df_age[dd.EngEng["gender"]] == dd.EngEng['both']]

    # write file for all age groups summed together
    df_agesum = df_age.loc[df_age[dd.EngEng["age10"]] == dd.EngEng['all']]

    gd.write_dataframe(df_agesum, directory, "spain", file_format)

    # write file with information on all age groups separately
    # age_groups = ['0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+']
    df_agesep = df_age.loc[df_age[dd.EngEng["age10"]] != dd.EngEng['all']]

    gd.write_dataframe(df_agesep, directory, "spain_all_age", file_format)


def _get_spain_all_state(read_data, file_format, no_raw, directory):

    stat_file = 'raw_spain_all_state'

    data_str = 'ccaa_covid19_datos_isciii'
    api_url = 'https://raw.githubusercontent.com/datadista/datasets/master/COVID%2019/old_series/'
    exit_string = "Files \'ccaa_covid19_datos_isciii\' are not available online. Check URL."

    df_state = _read_or_load_data(read_data, no_raw, directory,
                                  stat_file, data_str , api_url, exit_string)

    if not df_state.empty:
        # standardization of column titles from Spanish to English
        df_state.rename(dd.EsEng, axis=1, inplace=True)

        print("Available state columns:", df_state.columns)

        # fill empty cells (nan values) with zero
        df_state.replace(np.nan, 0, inplace=True)

        # remove special characters
        state = dd.EngEng['state']
        df_state.loc[df_state[state] == "Andalucía", [state]] = dd.EsEng["Andalucía"]
        df_state.loc[df_state[state] == "Castilla y León", [state]] = dd.EsEng["Castilla y León"]
        df_state.loc[df_state[state] == "Cataluña", [state]] = dd.EsEng["Cataluña"]
        df_state.loc[df_state[state] == "País Vasco", [state]] = dd.EsEng["País Vasco"]
        df_state.loc[df_state[state] == "Aragón", [state]] = dd.EsEng["Aragón"]

        # Correct Timestamps:
        df_state['Date'] = df_state['Date'].astype('datetime64[ns]')

    # Preparation for plotting/output:
    PCR_ONLY = False  # if pcr only is used
    # dff = df_state['state'].unique()

    # if PCR only is used, the number of confirmed cases is the number of people being tested positive with PCR test
    # otherwise, the number of positive antibody tests is also taken into account
    if PCR_ONLY:
        df_state.loc[df_state[dd.EngEng['confirmedTotal']] == 0, [dd.EngEng['confirmedTotal']]] \
            = df_state[dd.EngEng["confirmedPcr"]]
    else:
        df_state.loc[df_state[dd.EngEng['confirmedTotal']] == 0, [dd.EngEng['confirmedTotal']]] \
            = df_state[dd.EngEng["confirmedPcr"]] + df_state[dd.EngEng["confirmedAb"]]

    # output json
    gd.write_dataframe(df_state, directory, "spain_all_state", file_format)


def _read_or_load_data(read_data, no_raw, directory, file, data=None , api_url=None, exit_string=None):

    if read_data:
        JSONData = os.path.join(directory, file + '.json')

        # if once downloaded just read json file
        print("Read from local.")

        df = pandas.read_json(JSONData)

        if not df.empty:
            return df
        else:
            exit_string = "Something went wrong, dataframe is empty!"
            sys.exit(exit_string)

    else:

        print("Read Spanish data from online.")

        try:
            # Get data:
            df = gd.loadCsv(data, apiUrl=api_url)
        except:
            sys.exit(exit_string)

        if not df.empty:
            if not no_raw:
                # output data to not always download it
                gd.write_dataframe(df, directory, file, "json")
            return df
        else:
            exit_string = "Something went wrong, dataframe is empty!"
            sys.exit(exit_string)

def main():
    """! Main program entry."""

    arg_dict = gd.cli("spain")
    get_spain_data(**arg_dict)


if __name__ == "__main__":
    main()
