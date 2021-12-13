#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors:
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
import unittest
from pyfakefs import fake_filesystem_unittest
from freezegun import freeze_time
from datetime import date, timedelta

import os
import pandas as pd

from epidemiology.epidata import getDIVIData as gdd
from epidemiology.epidata import getRKIData as grki
from epidemiology.epidata import getVaccinationData as gvd
from epidemiology.epidata import getPopulationData as gpd
from epidemiology.epidata import getSimulationData as gsd

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd
from unittest.mock import patch, call


class TestGetSimulationData(fake_filesystem_unittest.TestCase):
    # construct fake directory for testing
    maxDiff = None

    path = '/home/SumlationData'

    def setUp(self):
        self.setUpPyfakefs()

    @patch('epidemiology.epidata.getDIVIData.get_divi_data')
    @patch('epidemiology.epidata.getRKIData.get_rki_data')
    @patch('epidemiology.epidata.getPopulationData.get_population_data')
    @patch('epidemiology.epidata.getPopulationData.get_age_population_data')
    @patch('epidemiology.epidata.getVaccinationData.get_vaccination_data')
    def test_get_call_sub_functions(self, mock_vaccination, mock_agep,
                                    mock_popul, mock_rki, mock_divi):

        [read_data, file_format, out_folder, no_raw, end_date, impute_dates,
         make_plot, moving_average, split_berlin, start_date] = [False,
                                                                 "json_timeasstring", os.path.join(dd.defaultDict['out_folder']),
                                                                 False, dd.defaultDict['end_date'],
                                                                 dd.defaultDict['impute_dates'],
                                                                 dd.defaultDict['make_plot'],
                                                                 dd.defaultDict['moving_average'],
                                                                 dd.defaultDict['split_berlin'],
                                                                 dd.defaultDict['start_date']]

        gsd.get_simulation_data(
            read_data, file_format, out_folder, no_raw, end_date, impute_dates,
            make_plot, moving_average, split_berlin, start_date)

        arg_dict_all = {
            "read_data": dd.defaultDict['read_data'],
            "file_format": dd.defaultDict['file_format'],
            "out_folder": os.path.join(dd.defaultDict['out_folder']),
            'no_raw': dd.defaultDict["no_raw"]}

        arg_dict_rki = {
            **arg_dict_all, "make_plot": dd.defaultDict['make_plot'],
            "impute_dates": dd.defaultDict['impute_dates'],
            "moving_average": dd.defaultDict['moving_average'],
            "split_berlin": dd.defaultDict['split_berlin']}

        arg_dict_divi = {
            **arg_dict_all, "end_date": dd.defaultDict['end_date'],
            "start_date": dd.defaultDict['start_date'],
            "moving_average": dd.defaultDict['moving_average']}

        arg_dict_vaccination = {
            **arg_dict_all,
            "make_plot": dd.defaultDict['make_plot'],
            "moving_average": dd.defaultDict['moving_average']}

        mock_vaccination.assert_called()
        mock_vaccination.assert_called_with(**arg_dict_vaccination)

        mock_agep.assert_called()
        mock_agep.assert_called_with(**arg_dict_all)

        mock_popul.assert_called()
        mock_popul.assert_called_with(**arg_dict_all)

        mock_rki.assert_called()
        mock_rki.assert_called_with(**arg_dict_rki)

        mock_divi.assert_called()
        mock_divi.assert_called_with(**arg_dict_divi)

    @patch('builtins.print')
    @patch('epidemiology.epidata.getDIVIData.get_divi_data')
    @patch('epidemiology.epidata.getRKIData.get_rki_data')
    @patch('epidemiology.epidata.getPopulationData.get_population_data')
    @patch('epidemiology.epidata.getPopulationData.get_age_population_data')
    @patch('epidemiology.epidata.getVaccinationData.get_vaccination_data')
    def test_errors(
            self, mock_vaccination, mock_agep, mock_popul, mock_rki, mock_divi,
            mock_print):
        mock_vaccination.side_effect = Exception
        mock_agep.side_effect = Exception
        mock_popul.side_effect = Exception
        mock_rki.side_effect = Exception
        mock_divi.side_effect = Exception
        gsd.get_simulation_data()
        vaccprint = call(
            'Error: Something went wrong while getting ' + 'vaccination' +
            ' data. This was likely caused by a changed file format'
            ' of the source material. Please report this as an issue. ' +
            'vaccination' + ' data could not be stored correctly.')
        agepprint = call(
            'Error: Something went wrong while getting ' +
            'age-resolved population' +
            ' data. This was likely caused by a changed file format'
            ' of the source material. Please report this as an issue. ' +
            'age-resolved population' + ' data could not be stored correctly.')
        populprint = call(
            'Error: Something went wrong while getting ' + 'population' +
            ' data. This was likely caused by a changed file format'
            ' of the source material. Please report this as an issue. ' +
            'population' + ' data could not be stored correctly.')
        rkiprint = call('Error: Something went wrong while getting ' + 'RKI' +
                        ' data. This was likely caused by a changed file format'
                        ' of the source material. Please report this as an issue. ' + 'RKI' +
                        ' data could not be stored correctly.')
        diviprint = call(
            'Error: Something went wrong while getting ' + 'DIVI' +
            ' data. This was likely caused by a changed file format'
            ' of the source material. Please report this as an issue. ' +
            'DIVI' + ' data could not be stored correctly.')
        expected_calls = [rkiprint, populprint,
                          agepprint, diviprint, vaccprint]
        mock_print.assert_has_calls(expected_calls)


if __name__ == '__main__':
    unittest.main()
