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

        [read_data, out_form, out_folder, end_date, fill_dates, make_plot, moving_average, split_berlin, start_date,
         update_data] = [False, "json", self.path, dd.defaultDict['end_date'], dd.defaultDict['fill_dates'],
                         dd.defaultDict['make_plot'], dd.defaultDict['moving_average'], dd.defaultDict['split_berlin'],
                         dd.defaultDict['start_date'], dd.defaultDict['update_data']]

        gsd.get_simulation_data(read_data, out_form, out_folder, end_date, fill_dates, make_plot, moving_average,
                                split_berlin, start_date, update_data)

        mock_vaccine.assert_called()
        mock_vaccine.assert_called_with(False, "json", self.path)

        mock_agep.assert_called()
        mock_agep.assert_called_with(False, "json", self.path)

        mock_popul.assert_called()
        mock_popul.assert_called_with(False, "json", self.path)

        mock_rki.assert_called()
        mock_rki.assert_called_with(False, "json", self.path, dd.defaultDict['fill_dates'], dd.defaultDict['make_plot'],
                                    dd.defaultDict['moving_average'], dd.defaultDict['split_berlin'])

        mock_divi.assert_called()
        mock_divi.assert_called_with(False, "json", self.path, dd.defaultDict['end_date'], dd.defaultDict['start_date'],
                                     dd.defaultDict['update_data'])









