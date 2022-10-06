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

from memilio.epidata import getSimulationData as gsd
from memilio.epidata import defaultDict as dd
from unittest.mock import patch, call


class TestGetSimulationData(fake_filesystem_unittest.TestCase):
    # construct fake directory for testing
    maxDiff = None

    path = '/home/SumlationData'

    def setUp(self):
        self.setUpPyfakefs()

    @patch('memilio.epidata.getVaccinationData.get_vaccination_data')
    @patch('memilio.epidata.getDIVIData.get_divi_data')
    @patch('memilio.epidata.getCaseData.get_case_data')
    @patch('memilio.epidata.getPopulationData.get_population_data')
    def test_get_call_sub_functions(self, mock_popul, mock_cases,
                                    mock_divi, mock_vaccination):

        [read_data, file_format, out_folder, no_raw, start_date, end_date,
         impute_dates, moving_average, make_plot, split_berlin, rep_date] = [
            False, "json_timeasstring", self.path, False,
            dd.defaultDict['start_date'],
            dd.defaultDict['end_date'],
            dd.defaultDict['impute_dates'],
            dd.defaultDict['moving_average'],
            dd.defaultDict['make_plot'],
            dd.defaultDict['split_berlin'],
            dd.defaultDict['rep_date']]

        gsd.get_simulation_data(
            read_data=read_data, file_format=file_format,
            out_folder=out_folder, no_raw=no_raw, start_date=start_date,
            end_date=end_date, impute_dates=impute_dates,
            moving_average=moving_average, make_plot=make_plot,
            split_berlin=split_berlin, rep_date=rep_date)

        arg_dict_all = {
            "read_data": dd.defaultDict['read_data'],
            "file_format": dd.defaultDict['file_format'],
            "out_folder": self.path,
            'no_raw': dd.defaultDict["no_raw"]}
        
        arg_dict_data_download = {
            "start_date": dd.defaultDict['start_date'],
            "end_date": dd.defaultDict['end_date'],
            "impute_dates": dd.defaultDict['impute_dates'],
            "moving_average": dd.defaultDict['moving_average'],
            "make_plot": dd.defaultDict['make_plot']}

        arg_dict_cases = {
            **arg_dict_all, **arg_dict_data_download,
            "rep_date": dd.defaultDict['rep_date'],
            "split_berlin": dd.defaultDict['split_berlin']}

        arg_dict_divi = {
            **arg_dict_all, **arg_dict_data_download}

        arg_dict_vaccination = {
            **arg_dict_all, **arg_dict_data_download,
            "sanitize_data": dd.defaultDict['sanitize_data']}

        mock_popul.assert_called()
        mock_popul.assert_called_with(**arg_dict_all)

        mock_cases.assert_called()
        mock_cases.assert_called_with(**arg_dict_cases)

        mock_divi.assert_called()
        mock_divi.assert_called_with(**arg_dict_divi)

        mock_vaccination.assert_called()
        mock_vaccination.assert_called_with(**arg_dict_vaccination)

    @patch('builtins.print')
    @patch('memilio.epidata.getVaccinationData.get_vaccination_data')
    @patch('memilio.epidata.getDIVIData.get_divi_data')
    @patch('memilio.epidata.getCaseData.get_case_data')
    @patch('memilio.epidata.getPopulationData.get_population_data')
    def test_errors(
            self, mock_popul, mock_cases, mock_divi, mock_vaccination,
            mock_print):
        mock_popul.side_effect = Exception
        mock_cases.side_effect = Exception
        mock_divi.side_effect = Exception
        mock_vaccination.side_effect = Exception
        gsd.get_simulation_data()
        populprint = call(
            'Error: Something went wrong while getting ' + 'population' +
            ' data. This was likely caused by a changed file format'
            ' of the source material. Please report this as an issue. ' +
            'population' + ' data could not be stored correctly.')
        casesprint = call(
            'Error: Something went wrong while getting ' + 'case' +
            ' data. This was likely caused by a changed file format'
            ' of the source material. Please report this as an issue. ' +
            'case' + ' data could not be stored correctly.')
        diviprint = call(
            'Error: Something went wrong while getting ' + 'DIVI' +
            ' data. This was likely caused by a changed file format'
            ' of the source material. Please report this as an issue. ' +
            'DIVI' + ' data could not be stored correctly.')
        vaccprint = call(
            'Error: Something went wrong while getting ' + 'vaccination' +
            ' data. This was likely caused by a changed file format'
            ' of the source material. Please report this as an issue. ' +
            'vaccination' + ' data could not be stored correctly.')

        exceptionprint = call('Exception: ')
        expected_calls = [
            exceptionprint, casesprint, exceptionprint, populprint,
            exceptionprint, diviprint, exceptionprint, vaccprint]
        mock_print.assert_has_calls(expected_calls)


if __name__ == '__main__':
    unittest.main()
