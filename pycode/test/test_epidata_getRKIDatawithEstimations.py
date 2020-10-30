import unittest
from pyfakefs import fake_filesystem_unittest
import numpy as np
from datetime import date, timedelta

import os
import pandas as pd

from epidemiology.epidata import getRKIDatawithEstimations as grdwd
from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd
from unittest.mock import patch, call

class TestGetRKIDatawithEstimations(fake_filesystem_unittest.TestCase):

    path = '/home/RKIEstimationData'

    str_whole_country_Germany_jh = ("""[{"CountryRegion":"Germany","Date":"2020-01-22","Confirmed":0,"Recovered":0.0,"Deaths":0.0},\
{"CountryRegion":"Germany","Date":"2020-01-23","Confirmed":0,"Recovered":0.0,"Deaths":0.0},\
{"CountryRegion":"Germany","Date":"2020-01-24","Confirmed":0,"Recovered":0.0,"Deaths":0.0},\
{"CountryRegion":"Germany","Date":"2020-01-25","Confirmed":0,"Recovered":0.0,"Deaths":0.0},\
{"CountryRegion":"Germany","Date":"2020-01-26","Confirmed":0,"Recovered":1.0,"Deaths":1.0},\
{"CountryRegion":"Germany","Date":"2020-01-27","Confirmed":1,"Recovered":2.0,"Deaths":1.0},\
{"CountryRegion":"Germany","Date":"2020-01-28","Confirmed":4,"Recovered":3.0,"Deaths":1.0},\
{"CountryRegion":"Germany","Date":"2020-01-29","Confirmed":4,"Recovered":3.0,"Deaths":1.0},\
{"CountryRegion":"Germany","Date":"2020-01-30","Confirmed":4,"Recovered":3.0,"Deaths":1.0},\
{"CountryRegion":"Germany","Date":"2020-01-31","Confirmed":5,"Recovered":3.0,"Deaths":2.0}]""")

    str_all_germany_rki = ("""[{"Date":1577836800000,"Confirmed":14,"Deaths":0,"Recovered":14},\
{"Date":1577923200000,"Confirmed":15,"Deaths":0,"Recovered":15},\
{"Date":1578009600000,"Confirmed":16,"Deaths":0,"Recovered":16},\
{"Date":1578096000000,"Confirmed":17,"Deaths":0,"Recovered":17},\
{"Date":1578441600000,"Confirmed":18,"Deaths":0,"Recovered":18},\
{"Date":1578528000000,"Confirmed":19,"Deaths":0,"Recovered":19},\
{"Date":1578614400000,"Confirmed":20,"Deaths":0,"Recovered":20},\
{"Date":1578700800000,"Confirmed":23,"Deaths":0,"Recovered":23},\
{"Date":1578873600000,"Confirmed":24,"Deaths":0,"Recovered":24},\
{"Date":1578960000000,"Confirmed":25,"Deaths":0,"Recovered":25},\
{"Date":1579046400000,"Confirmed":26,"Deaths":0,"Recovered":26},\
{"Date":1579132800000,"Confirmed":27,"Deaths":0,"Recovered":27},\
{"Date":1579219200000,"Confirmed":30,"Deaths":1,"Recovered":29},\
{"Date":1579305600000,"Confirmed":32,"Deaths":1,"Recovered":31},\
{"Date":1579392000000,"Confirmed":33,"Deaths":1,"Recovered":32},\
{"Date":1579478400000,"Confirmed":34,"Deaths":1,"Recovered":33},\
{"Date":1579564800000,"Confirmed":35,"Deaths":1,"Recovered":34},\
{"Date":1579651200000,"Confirmed":38,"Deaths":1,"Recovered":37},\
{"Date":1579737600000,"Confirmed":40,"Deaths":1,"Recovered":39},\
{"Date":1579824000000,"Confirmed":42,"Deaths":1,"Recovered":41},\
{"Date":1579910400000,"Confirmed":43,"Deaths":1,"Recovered":42},\
{"Date":1580083200000,"Confirmed":46,"Deaths":1,"Recovered":45},\
{"Date":1580169600000,"Confirmed":47,"Deaths":1,"Recovered":46},\
{"Date":1580256000000,"Confirmed":48,"Deaths":1,"Recovered":47},\
{"Date":1580342400000,"Confirmed":53,"Deaths":1,"Recovered":52},\
{"Date":1580428800000,"Confirmed":56,"Deaths":1,"Recovered":55}]""")

    rki_files_to_change = ["all_germany_rki", "all_gender_rki", "all_age_rki",
                           "all_state_rki", "all_state_gender_rki", "all_state_age_rki",
                           "all_county_rki", "all_county_gender_rki", "all_county_age_rki"]

    def setUp(self):
        self.setUpPyfakefs()


    def write_rki_data(self, out_folder):

        for file_to_change in self.rki_files_to_change:

            file_rki = file_to_change + ".json"
            file_rki_with_path = os.path.join(out_folder, file_rki)

            with open(file_rki_with_path, 'w') as f:
                f.write(self.str_all_germany_rki )

    def write_jh_data(self, out_folder):
        file_jh = "whole_country_Germany_jh.json"
        file_jh_with_path = os.path.join(out_folder, file_jh)

        with open(file_jh_with_path, 'w') as f:
            f.write(self.str_whole_country_Germany_jh)

    def test_get_rki_data_with_estimations(self):

        [read_data, make_plot, out_form, out_folder] \
            = [True, False, "json", self.path, ]

        # write files which should be read in by program

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        self.write_rki_data(directory)
        self.write_jh_data(directory)

        # check if expected files are written
        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(len(os.listdir(directory)), 1+len(self.rki_files_to_change))

        grdwd.get_rki_data_with_estimations(read_data, out_form, out_folder, make_plot)

        # check if expected files are written
        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(len(os.listdir(directory)), 1+ 2*len(self.rki_files_to_change))

        f_read = os.path.join(directory, "all_germany_rki_estimated.json")
        df = pd.read_json(f_read)

        confirmed = dd.EngEng['confirmed']
        recovered = dd.EngEng['recovered']
        deaths = dd.EngEng['deaths']
        date = dd.EngEng['date']
        recovered_estimated = recovered + "_estimated"
        deaths_estimated = deaths + "_estimated"

        data_list = df.columns.values.tolist()

        self.assertEqual(data_list, [date, confirmed, deaths, recovered, recovered_estimated, deaths_estimated])

        self.assertEqual(df[(df[date] == "2020-01-31")][recovered_estimated].item(), np.round(56*3./5.))

    @patch('epidemiology.epidata.getRKIDatawithEstimations.grd.get_rki_data')
    @patch('epidemiology.epidata.getRKIDatawithEstimations.gjd.get_jh_data')
    def test_get_rki_data_with_estimations_download(self,mock_get_jh_data, mock_get_rki_data):

        [read_data, make_plot, out_form, out_folder] \
            = [False, False, "json", self.path, ]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)


        mock_get_rki_data.side_effect =self.write_rki_data(directory)
        mock_get_jh_data.side_effect = self.write_jh_data(directory)

        # write files which should be read in by program

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        rki_files_to_change = ["all_germany_rki", "all_gender_rki", "all_age_rki",
                               "all_state_rki", "all_state_gender_rki", "all_state_age_rki",
                               "all_county_rki", "all_county_gender_rki", "all_county_age_rki"]

        for file_to_change in rki_files_to_change:

            file_rki = file_to_change + ".json"
            file_rki_with_path = os.path.join(directory, file_rki)

            with open(file_rki_with_path, 'w') as f:
                f.write(self.str_all_germany_rki )

        file_jh = "whole_country_Germany_jh.json"
        file_jh_with_path = os.path.join(directory, file_jh)

        with open(file_jh_with_path, 'w') as f:
            f.write(self.str_whole_country_Germany_jh)

        # check if expected files are written
        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(len(os.listdir(directory)), 1+len(rki_files_to_change))

        grdwd.get_rki_data_with_estimations(read_data, out_form, out_folder, make_plot)

        # check if expected files are written
        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(len(os.listdir(directory)), 1+ 2*len(rki_files_to_change))

        f_read = os.path.join(directory, "all_germany_rki_estimated.json")
        df = pd.read_json(f_read)

        confirmed = dd.EngEng['confirmed']
        recovered = dd.EngEng['recovered']
        deaths = dd.EngEng['deaths']
        date = dd.EngEng['date']
        recovered_estimated = recovered + "_estimated"
        deaths_estimated = deaths + "_estimated"

        data_list = df.columns.values.tolist()

        self.assertEqual(data_list, [date, confirmed, deaths, recovered, recovered_estimated, deaths_estimated])

        self.assertEqual(df[(df[date] == "2020-01-31")][recovered_estimated].item(), np.round(56*3./5.))



if __name__ == '__main__':
    unittest.main()
