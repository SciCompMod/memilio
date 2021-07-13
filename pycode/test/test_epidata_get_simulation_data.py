import unittest
from pyfakefs import fake_filesystem_unittest
from freezegun import freeze_time
from datetime import date, timedelta

import os
import pandas as pd

from epidemiology.epidata import getDIVIData as gdd
from epidemiology.epidata import getRKIData as grki
from epidemiology.epidata import getVaccineData as gvd
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
    @patch('epidemiology.epidata.getVaccineData.get_vaccine_data')
    def test_get_call_sub_functions(self, mock_vaccine, mock_agep, mock_popul, mock_rki, mock_divi):

        [read_data, file_format, out_folder, no_raw, end_date, fill_dates, make_plot, moving_average, split_berlin,
         start_date, update_data] = [False, "json", self.path, False, dd.defaultDict['end_date'],
                                     dd.defaultDict['fill_dates'], dd.defaultDict['make_plot'],
                                     dd.defaultDict['moving_average'], dd.defaultDict['split_berlin'],
                                     dd.defaultDict['start_date'], dd.defaultDict['update_data']]

        gsd.get_simulation_data(read_data, file_format, out_folder, no_raw, end_date, fill_dates, make_plot,
                                moving_average, split_berlin, start_date, update_data)

        arg_dict_all = {"read_data": False, "file_format": "json", "out_folder": self.path, "no_raw": False}

        arg_dict_rki = {**arg_dict_all, "make_plot": dd.defaultDict['make_plot'],
                        "fill_dates": dd.defaultDict['fill_dates'], "moving_average": dd.defaultDict['moving_average'],
                        "split_berlin": dd.defaultDict['split_berlin']}

        arg_dict_divi = {**arg_dict_all, "end_date": dd.defaultDict['end_date'],
                         "start_date": dd.defaultDict['start_date'], "update_data": dd.defaultDict['update_data']}

        mock_vaccine.assert_called()
        mock_vaccine.assert_called_with(**arg_dict_all)

        mock_agep.assert_called()
        mock_agep.assert_called_with(**arg_dict_all)

        mock_popul.assert_called()
        mock_popul.assert_called_with(**arg_dict_all)

        mock_rki.assert_called()
        mock_rki.assert_called_with(**arg_dict_rki)

        mock_divi.assert_called()
        mock_divi.assert_called_with(**arg_dict_divi)


if __name__ == '__main__':
    unittest.main()
