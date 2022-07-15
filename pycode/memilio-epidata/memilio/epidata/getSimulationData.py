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
@file getSimulationData.py

@brief Executes all data downloads which belong to the epidata package and downloads external data

The functions which are called are:
- getCaseData.get_case_data
- getPopulationData.get_population_data
- getVacccinationData.get_vaccination_data
- getDIVIData.get_divi_data
"""


import os
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import defaultDict as dd
from memilio.epidata import getVaccinationData
from memilio.epidata import getPopulationData
from memilio.epidata import getCaseData
from memilio.epidata import getDIVIData


def print_error(text):
    print('Error: Something went wrong while getting ' + text +
          ' data. This was likely caused by a changed file format'
          ' of the source material. Please report this as an issue. ' + text +
          ' data could not be stored correctly.')


def get_simulation_data(data_folder,
                        read_data=dd.defaultDict['read_data'],
                        file_format=dd.defaultDict['file_format'],
                        no_raw=dd.defaultDict['no_raw'],
                        end_date=dd.defaultDict['end_date'],
                        impute_dates=dd.defaultDict['impute_dates'],
                        make_plot=dd.defaultDict['make_plot'],
                        moving_average=dd.defaultDict['moving_average'],
                        split_berlin=dd.defaultDict['split_berlin'],
                        start_date=dd.defaultDict['start_date']
                        ):
    """! Downloads all data from external sources

    The functions which are called are:
    - getCaseData.get_case_data
    - getPopulationData.get_population_data
    - getVaccinationData.get_vaccination_data
    - getDIVIData.get_divi_data

    Keyword arguments:
    @param data_folder Path to folder where data is written in.
    @param read_data False [Default] or True. Defines if data is read from file or downloaded.
    @param file_format File format which is used for writing the data. Default defined in defaultDict.
    @param no_raw True or False [Default]. Defines if unchanged raw data is saved or not.
    @param end_date [Optional] Date to stop to download data [Default = today].
    @param impute_dates False [Default] or True. Defines if dates where nothing changed are added.
    @param make_plot False [Default] or True. Defines if plots are generated with matplotlib.
    @param moving_average True or False [Default]. Defines if files for 7 day moving average should be created
    @param split_berlin True [Default] or False. Defines if Berlin counties is fused to just on county.
    @param start_date [Optional] Date to start to download data [Default = 2020-04-24].
    """

    arg_dict_all = {
        "read_data": read_data, "file_format": file_format,
        "no_raw": no_raw}

    arg_dict_cases = {**arg_dict_all, "make_plot": make_plot,
                    "impute_dates": impute_dates,
                    "moving_average": moving_average,
                    "split_berlin": split_berlin}

    arg_dict_divi = {**arg_dict_all, "end_date": end_date,
                     "start_date": start_date, "moving_average": moving_average}

    arg_dict_vacc = {**arg_dict_all, "make_plot": make_plot,
                     "moving_average": moving_average}

    try:
        getCaseData.get_case_data(data_folder, **arg_dict_cases)
    except Exception as exp:
        print(str(type(exp).__name__) + ": " + str(exp))
        print_error('case')

    try:
        getPopulationData.get_population_data(data_folder, **arg_dict_all)
    except Exception as exp:
        print(str(type(exp).__name__) + ": " + str(exp))
        print_error('population')

    try:
        getDIVIData.get_divi_data(data_folder, **arg_dict_divi)
    except Exception as exp:
        print(str(type(exp).__name__) + ": " + str(exp))
        print_error('DIVI')

    try:
        getVaccinationData.get_vaccination_data(data_folder, **arg_dict_vacc)
    except Exception as exp:
        print(str(type(exp).__name__) + ": " + str(exp))
        print_error('vaccination')


def main():
    """! Main program entry."""

    path = os.path.join(os.getcwd(), 'data', 'pydata')
    arg_dict = gd.cli("sim")
    get_simulation_data(path, **arg_dict)


if __name__ == "__main__":
    main()
