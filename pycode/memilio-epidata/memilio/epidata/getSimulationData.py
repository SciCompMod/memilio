#############################################################################
# Copyright (C) 2020-2024 MEmilio
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
- getCommuterMobility.get_commuter_data
"""
import os

from memilio.epidata import defaultDict as dd
from memilio.epidata import getCaseData
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import getDIVIData, getPopulationData, getVaccinationData, getCommuterMobility, transformMobilityData


def print_error(text):
    gd.default_print('Error', 'Something went wrong while getting ' + text +
                     ' data. This was likely caused by a changed file format'
                     ' of the source material. Please report this as an issue. ' + text +
                     ' data could not be stored correctly.')


def get_simulation_data(read_data=dd.defaultDict['read_data'],
                        file_format=dd.defaultDict['file_format'],
                        out_folder=dd.defaultDict['out_folder'],
                        start_date=dd.defaultDict['start_date'],
                        end_date=dd.defaultDict['end_date'],
                        impute_dates=dd.defaultDict['impute_dates'],
                        moving_average=dd.defaultDict['moving_average'],
                        split_berlin=dd.defaultDict['split_berlin'],
                        rep_date=dd.defaultDict['rep_date'],
                        sanitize_data=dd.defaultDict['sanitize_data'],
                        ref_year = 2022,
                        **kwargs
                        ):
    """! Downloads all data from external sources

    The functions which are called are:
    - getCaseData.get_case_data
    - getPopulationData.get_population_data
    - getVaccinationData.get_vaccination_data
    - getDIVIData.get_divi_data
    - getCommuterMobility.get_commuter_data
    - transformMobilityData.updateMobility2022 (if ref_year < 2022)

    Keyword arguments:
    @param read_data True or False. Defines if data is read from file or downloaded. Default defined in defaultDict.
    @param file_format File format which is used for writing the data. Default defined in defaultDict.
    @param out_folder Folder where data is written to. Default defined in defaultDict.
    @param start_date Date of first date in dataframe. Default 2020-01-01.
    @param end_date Date of last date in dataframe. Default defined in defaultDict.
    @param impute_dates True or False. Defines if values for dates without new information are imputed. Default defined in defaultDict.
    @param moving_average Integers >=0. Applies an 'moving_average'-days moving average on all time series
        to smooth out effects of irregular reporting. Default defined in defaultDict.
    @param split_berlin True or False. Defines if Berlin's disctricts are kept separated or get merged. Default defined in defaultDict.
    @param rep_date True or False. Defines if reporting date or reference date is taken into dataframe. Default defined in defaultDict.
    @param sanitize_data Value in {0,1,2,3}. Redistributes cases of every county either based on regions' ratios or on thresholds and population.
    @param ref_year Year between 2013 and 2022 that specifies where the data should be taken from. Default value is 2022.
    """
    conf = gd.Conf(out_folder, **kwargs)
    out_folder = conf.path_to_use
    no_raw = conf.no_raw
    make_plot = conf.plot

    arg_dict_all = {
        "read_data": read_data, "file_format": file_format,
        "out_folder": out_folder, "no_raw": no_raw}

    arg_dict_data_download = {"start_date": start_date, "end_date": end_date,
                              "impute_dates": impute_dates, "moving_average": moving_average,
                              "make_plot": make_plot}

    arg_dict_cases = {**arg_dict_all, **arg_dict_data_download,
                      "split_berlin": split_berlin, "rep_date": rep_date}

    arg_dict_vacc = {**arg_dict_all, **arg_dict_data_download,
                     "sanitize_data": sanitize_data}

    arg_dict_divi = {**arg_dict_all, **arg_dict_data_download}

    arg_dict_mobility = {**arg_dict_all, **arg_dict_data_download,
                         "ref_year": ref_year}   
    
    try:
        getCaseData.get_case_data(**arg_dict_cases)
    except Exception as exp:
        gd.default_print('Error', str(type(exp).__name__) + ": " + str(exp))
        print_error('case')

    try:
        getPopulationData.get_population_data(**arg_dict_all)
    except Exception as exp:
        gd.default_print('Error', str(type(exp).__name__) + ": " + str(exp))
        print_error('population')

    try:
        getDIVIData.get_divi_data(**arg_dict_divi)
    except Exception as exp:
        gd.default_print('Error', str(type(exp).__name__) + ": " + str(exp))
        print_error('DIVI')

    try:
        getVaccinationData.get_vaccination_data(**arg_dict_vacc)
    except Exception as exp:
        gd.default_print('Error', str(type(exp).__name__) + ": " + str(exp))
        print_error('vaccination')

    try:
        getCommuterMobility.get_commuter_data(**arg_dict_mobility)
    except Exception as exp:
        gd.default_print('Error', str(type(exp).__name__) + ": " + str(exp))
        print_error('commuter mobility')
    
    try:
        transformMobilityData.updateMobility2022(out_folder, mobility_file='commuter_mobility')
        if(ref_year < 2022):
            transformMobilityData.updateMobility2022(out_folder, mobility_file='twitter_scaled_1252')
    except Exception as exp:
        gd.default_print('Error', str(type(exp).__name__) + ": " + str(exp))
        print_error('transform mobility')

    # rename commuter mobility file
    os.rename(f'{out_folder}/commuter_mobility_{ref_year}.txt', f'{out_folder}/commuter_mobility.txt')


def main():
    """! Main program entry."""

    arg_dict = gd.cli("sim")
    get_simulation_data(**arg_dict)


if __name__ == "__main__":
    main()
