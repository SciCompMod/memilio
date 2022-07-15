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
import os
import sys
from unittest.mock import patch

from memilio.epidata import cleanData as cd
from memilio.epidata import defaultDict as dd


class Test_cleanData(fake_filesystem_unittest.TestCase):

    path = '/home/x'
    maxDiff = None

    def setUp(self):
        self.setUpPyfakefs()

    def set_dirs_and_files(self, what):

        dir_dic_all = {
            'Germany':
            ["cases_a", "a_jh", "CaseDataFull", "PopulData",
             "county_population", "migration", "reg_key", "zensus", "FullVacc",
             "all_county_vacc", "all_state_vacc", "migration_bfa_2020_dim401",
             "states_testpos", "FullData_DIVI", "county_divi"],
            'Spain': ["b_jh"],
            'France': ["c_jh"],
            'Italy': ["d_jh"],
            'US': ["e_jh"],
            'SouthKorea': ["f_jh"],
            'China': ["g_jh"]}

        dir_dic_cases = {'Germany': ["cases_a", "CaseDataFull"]}

        dir_dic_popul = {'Germany': ["PopulData", "county_population", "migration", "reg_key", "zensus"]}

        dir_dic_jh = {'Germany': ["a_jh"],
                      'Spain': ["b_jh"],
                      'France': ["c_jh"],
                      'Italy': ["d_jh"],
                      'US': ["e_jh"],
                      'SouthKorea': ["f_jh"],
                      'China': ["g_jh"]}

        dir_dic_divi = {'Germany': ["FullData_DIVI", "county_divi"]}

        dir_dic_vacc = {'Germany': ["FullVacc", "all_county_vacc",
                        "all_state_vacc"]}

        dir_dic_commuter = {
            'Germany': ["migration_bfa_2020_dim401"]}

        dir_dic_testing = {'Germany': ["states_testpos"]}

        ending_all = [".json", ".h5"]
        ending_json = [".json"]

        dir_choose = {"all": dir_dic_all,
                      "cases": dir_dic_cases,
                      "jh": dir_dic_jh,
                      "popul": dir_dic_popul,
                      "divi": dir_dic_divi,
                      "vacc": dir_dic_vacc,
                      "commuter": dir_dic_commuter,
                      "testing": dir_dic_testing
                      }

        ending_choose = {"all": ending_all,
                         "cases": ending_json,
                         "jh": ending_json,
                         "popul": ending_json,
                         "divi": ending_json,
                         "vacc": ending_json,
                         "commuter": ending_json,
                         "testing": ending_json
                         }

        dir_dic = dir_choose[what]
        ending = ending_choose[what]

        file_list = ["all_jh", "FullJohnHopkins"]

        # make folders
        for key in dir_dic:
            dir_path = os.path.join(self.path, key)
            os.makedirs(dir_path)

            # make files
            for file in dir_dic[key]:
                for e in ending:
                    with open(os.path.join(dir_path, file + e), 'w') as f:
                        f.write('foo')

        if what == "all" or what == "jh":
            for file in file_list:
                for e in ending:
                    with open(os.path.join(self.path, file + e), 'w') as f:
                        f.write('foo')

    def test_set_dirs_and_files(self):

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        dir_list = ['Germany', 'Spain', 'France', 'Italy',
                    'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 11)

        self.assertEqual(
            os.listdir(self.path),
            dir_list +
            ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json',
             'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 30)

                fakefiles = [
                    "cases_a.h5", "a_jh.h5", "CaseDataFull.h5", "PopulData.h5",
                    "county_population.h5", "migration.h5", "reg_key.h5", "zensus.h5", "FullVacc.h5",
                    "all_county_vacc.h5", "all_state_vacc.h5",
                    "migration_bfa_2020_dim401.h5", "states_testpos.h5",
                    "FullData_DIVI.h5", "county_divi.h5", "cases_a.json",
                    "a_jh.json", "CaseDataFull.json", "PopulData.json",
                    "county_population.json", "migration.json", "reg_key.json", "zensus.json", "FullVacc.json",
                    "all_county_vacc.json", "all_state_vacc.json",
                    "migration_bfa_2020_dim401.json", "states_testpos.json",
                    "FullData_DIVI.json", "county_divi.json"]

                for file in fakefiles:
                    self.assertIn(file, os.listdir(dir_path))
                for file in os.listdir(dir_path):
                    self.assertIn(file, fakefiles)

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(
                    os.listdir(dir_path),
                    ["b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(
                        os.listdir(dir_path),
                        ["c_jh.json", "c_jh.h5"])

    # generate folder and files
    def test_clean_data_all_should_delete_all(self):

        self.set_dirs_and_files("all")

        cd.clean_data(True, False, False, False, False,
                      False, False, False, True, False, False, self.path)

        # Should delete everything
        self.assertEqual(len(os.listdir(self.path)), 0)

    def test_clean_data_all_should_not_delete_all(self):

        self.set_dirs_and_files("all")

        # add different files and folder
        os.makedirs(os.path.join(self.path, "ImportantDir"))

        with open(os.path.join(self.path, "wichtig.py"), 'w') as f:
            f.write('foo')

        dir_path = os.path.join(self.path, "China")
        with open(os.path.join(dir_path, "secret.txt"), 'w') as f:
            f.write('foo')

        cd.clean_data(True, False, False, False, False,
                      False, False, False, True, False, False, self.path)

        # Should delete everything
        self.assertEqual(len(os.listdir(self.path)), 3)
        self.assertEqual(
            os.listdir(self.path),
            ["China", "ImportantDir", "wichtig.py"])
        self.assertEqual(len(os.listdir(dir_path)), 1)
        self.assertEqual(os.listdir(dir_path), ["secret.txt"])

    def test_clean_data_cases(self):

        self.set_dirs_and_files("all")

        cd.clean_data(False, True, False, False, False,
                      False, False, False, True, False, False, self.path)

        dir_list = ['Germany', 'Spain', 'France', 'Italy',
                    'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 11)

        self.assertEqual(
            os.listdir(self.path),
            dir_list +
            ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json',
             'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 28)
                fakefiles = [
                    "cases_a.h5", "a_jh.h5", "CaseDataFull.h5", "PopulData.h5",
                    "county_population.h5", "migration.h5", "reg_key.h5", "zensus.h5", "FullVacc.h5",
                    "all_county_vacc.h5", "all_state_vacc.h5",
                    "migration_bfa_2020_dim401.h5", "states_testpos.h5",
                    "FullData_DIVI.h5", "county_divi.h5",
                    "a_jh.json", "PopulData.json",
                    "county_population.json", "migration.json", "reg_key.json", "zensus.json", "FullVacc.json",
                    "all_county_vacc.json", "all_state_vacc.json",
                    "migration_bfa_2020_dim401.json", "states_testpos.json",
                    "FullData_DIVI.json", "county_divi.json"]

                for file in fakefiles:
                    self.assertIn(file, os.listdir(dir_path))
                for file in os.listdir(dir_path):
                    self.assertIn(file, fakefiles)

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(
                    os.listdir(dir_path),
                    ["b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(
                        os.listdir(dir_path),
                        ["c_jh.json", "c_jh.h5"])

    def test_clean_data_cases_h5(self):

        self.set_dirs_and_files("all")

        cd.clean_data(False, True, False, False, False,
                      False, False, False, False, True, False, self.path)

        dir_list = ['Germany', 'Spain', 'France', 'Italy',
                    'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 11)

        self.assertEqual(
            os.listdir(self.path),
            dir_list +
            ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json',
             'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 28)
                fakefiles = [
                    "a_jh.h5", "PopulData.h5",
                    "county_population.h5", "migration.h5", "reg_key.h5", "zensus.h5", "FullVacc.h5",
                    "all_county_vacc.h5", "all_state_vacc.h5",
                    "migration_bfa_2020_dim401.h5", "states_testpos.h5",
                    "FullData_DIVI.h5", "county_divi.h5", "cases_a.json",
                    "a_jh.json", "CaseDataFull.json", "PopulData.json",
                    "county_population.json", "migration.json", "reg_key.json", "zensus.json", "FullVacc.json",
                    "all_county_vacc.json", "all_state_vacc.json",
                    "migration_bfa_2020_dim401.json", "states_testpos.json",
                    "FullData_DIVI.json", "county_divi.json"]

                for file in fakefiles:
                    self.assertIn(file, os.listdir(dir_path))
                for file in os.listdir(dir_path):
                    self.assertIn(file, fakefiles)

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(
                    os.listdir(dir_path),
                    ["b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(
                        os.listdir(dir_path),
                        ["c_jh.json", "c_jh.h5"])

    def test_clean_data_cases_del_dir(self):

        self.set_dirs_and_files("all")

        dir_path = os.path.join(self.path, "Germany")
        files = os.listdir(dir_path)

        # delete all files except which will be deleted
        for item in files:
            if item == "cases_a.json" or item == "CaseDataFull.json":
                continue
            else:
                os.remove(os.path.join(dir_path, item))

        self.assertEqual(len(os.listdir(dir_path)), 2)
        self.assertEqual(
            os.listdir(dir_path),
            ["cases_a.json", "CaseDataFull.json"])

        cd.clean_data(False, True, False, False, False,
                      False, False, False, True, False, False, self.path)

        dir_list = ['Spain', 'France', 'Italy',
                    'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure
        # Germany should be deleted, because no files where left after deletion
        self.assertEqual(len(os.listdir(self.path)), 10)

        self.assertEqual(
            os.listdir(self.path),
            dir_list +
            ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json',
             'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(
                    os.listdir(dir_path),
                    ["b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(
                        os.listdir(dir_path),
                        ["c_jh.json", "c_jh.h5"])

    def test_clean_data_population(self):

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        cd.clean_data(False, False, False, True, False,
                      False, False, False, True, False, False, self.path)

        dir_list = ['Germany', 'Spain', 'France', 'Italy',
                    'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 11)

        self.assertEqual(
            os.listdir(self.path),
            dir_list +
            ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json',
             'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 25)
                fakefiles = [
                    "cases_a.h5", "a_jh.h5", "CaseDataFull.h5", "PopulData.h5",
                    "county_population.h5", "migration.h5", "reg_key.h5", "zensus.h5", "FullVacc.h5",
                    "all_county_vacc.h5", "all_state_vacc.h5",
                    "migration_bfa_2020_dim401.h5", "states_testpos.h5",
                    "FullData_DIVI.h5", "county_divi.h5", "cases_a.json",
                    "a_jh.json", "CaseDataFull.json", "FullVacc.json",
                    "all_county_vacc.json", "all_state_vacc.json",
                    "migration_bfa_2020_dim401.json", "states_testpos.json",
                    "FullData_DIVI.json", "county_divi.json"]

                for file in fakefiles:
                    self.assertIn(file, os.listdir(dir_path))
                for file in os.listdir(dir_path):
                    self.assertIn(file, fakefiles)

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(
                    os.listdir(dir_path),
                    ["b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(
                        os.listdir(dir_path),
                        ["c_jh.json", "c_jh.h5"])

    def test_clean_data_population_hdf5(self):

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        cd.clean_data(False, False, False, True, False,
                      False, False, False, False, True, False, self.path)

        dir_list = ['Germany', 'Spain', 'France', 'Italy',
                    'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 11)

        self.assertEqual(
            os.listdir(self.path),
            dir_list +
            ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json',
             'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 25)
                fakefiles = [
                    "cases_a.h5", "a_jh.h5", "CaseDataFull.h5", "FullVacc.h5",
                    "all_county_vacc.h5", "all_state_vacc.h5",
                    "migration_bfa_2020_dim401.h5", "states_testpos.h5",
                    "FullData_DIVI.h5", "county_divi.h5", "cases_a.json",
                    "a_jh.json", "CaseDataFull.json", "PopulData.json",
                    "county_population.json", "migration.json", "reg_key.json", "zensus.json", "FullVacc.json",
                    "all_county_vacc.json", "all_state_vacc.json",
                    "migration_bfa_2020_dim401.json", "states_testpos.json",
                    "FullData_DIVI.json", "county_divi.json"]

                for file in fakefiles:
                    self.assertIn(file, os.listdir(dir_path))
                for file in os.listdir(dir_path):
                    self.assertIn(file, fakefiles)

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(
                    os.listdir(dir_path),
                    ["b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(
                        os.listdir(dir_path),
                        ["c_jh.json", "c_jh.h5"])

    def test_clean_data_population_del_dir(self):

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        dir_path = os.path.join(self.path, "Germany")
        files = os.listdir(dir_path)

        population_files = [
            "PopulData.json", "county_population.json", "migration.json",
            "reg_key.json", "zensus.json"]
        # delete all files except which will be deleted
        for item in files:
            if item in population_files:
                continue
            else:
                os.remove(os.path.join(dir_path, item))

        self.assertEqual(len(os.listdir(dir_path)), len(population_files))

        cd.clean_data(False, False, False, True, False,
                      False, False, False, True, False, False, self.path)

        dir_list = ['Spain', 'France', 'Italy',
                    'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 10)

        self.assertEqual(
            os.listdir(self.path),
            dir_list +
            ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json',
             'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(
                    os.listdir(dir_path),
                    ["b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(
                        os.listdir(dir_path),
                        ["c_jh.json", "c_jh.h5"])

    def test_all_false(self):

        cd.clean_data(False, False, False, False, False,
                      False, False, False, True, False, False, self.path)

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        dir_list = ['Germany', 'Spain', 'France', 'Italy',
                    'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 11)

        self.assertEqual(
            os.listdir(self.path),
            dir_list +
            ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json',
             'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 30)
                fakefiles = [
                    "cases_a.h5", "a_jh.h5", "CaseDataFull.h5", "PopulData.h5",
                    "county_population.h5", "migration.h5", "reg_key.h5", "zensus.h5", "FullVacc.h5",
                    "all_county_vacc.h5", "all_state_vacc.h5",
                    "migration_bfa_2020_dim401.h5", "states_testpos.h5",
                    "FullData_DIVI.h5", "county_divi.h5", "cases_a.json",
                    "a_jh.json", "CaseDataFull.json", "PopulData.json",
                    "county_population.json", "migration.json", "reg_key.json", "zensus.json", "FullVacc.json",
                    "all_county_vacc.json", "all_state_vacc.json",
                    "migration_bfa_2020_dim401.json", "states_testpos.json",
                    "FullData_DIVI.json", "county_divi.json"]

                for file in fakefiles:
                    self.assertIn(file, os.listdir(dir_path))
                for file in os.listdir(dir_path):
                    self.assertIn(file, fakefiles)

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(
                    os.listdir(dir_path),
                    ["b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(
                        os.listdir(dir_path),
                        ["c_jh.json", "c_jh.h5"])

    def test_wrong_path(self):

        self.set_dirs_and_files("all")

        dir1 = os.listdir(self.path)
        dir1a = os.listdir(os.path.join(self.path, 'Germany/'))

        cd.clean_data(True, False, False, False, False,
                      False, False, False, True, False, False, "/home/y")

        dir2 = os.listdir(self.path)
        dir2a = os.listdir(os.path.join(self.path, 'Germany/'))

        self.assertEqual(dir1, dir2)
        self.assertEqual(dir1a, dir2a)

    def test_clean_data_jh(self):

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        cd.clean_data(False, False, True, False, False,
                      False, False, False, True, False, False, self.path)

        dir_list = ['Germany', 'Spain', 'France', 'Italy',
                    'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 9)

        self.assertEqual(os.listdir(self.path),
                         dir_list + ['all_jh.h5', 'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 29)
                fakefiles = [
                    "cases_a.h5", "a_jh.h5", "CaseDataFull.h5", "PopulData.h5",
                    "county_population.h5", "migration.h5", "reg_key.h5", "zensus.h5", "FullVacc.h5",
                    "all_county_vacc.h5", "all_state_vacc.h5",
                    "migration_bfa_2020_dim401.h5", "states_testpos.h5",
                    "FullData_DIVI.h5", "county_divi.h5", "cases_a.json",
                    "CaseDataFull.json", "PopulData.json",
                    "county_population.json", "migration.json", "reg_key.json", "zensus.json", "FullVacc.json",
                    "all_county_vacc.json", "all_state_vacc.json",
                    "migration_bfa_2020_dim401.json", "states_testpos.json",
                    "FullData_DIVI.json", "county_divi.json"]

                for file in fakefiles:
                    self.assertIn(file, os.listdir(dir_path))
                for file in os.listdir(dir_path):
                    self.assertIn(file, fakefiles)

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 1)
                self.assertEqual(os.listdir(dir_path), ["b_jh.h5"])
            elif dir == "France":
                self.assertEqual(len(os.listdir(dir_path)), 1)
                self.assertEqual(os.listdir(dir_path), ["c_jh.h5"])

    def test_clean_data_jh_hdf5(self):

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        cd.clean_data(False, False, True, False, False,
                      False, False, False, False, True, False, self.path)

        dir_list = ['Germany', 'Spain', 'France', 'Italy',
                    'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 9)

        self.assertEqual(os.listdir(self.path),
                         dir_list + ['all_jh.json', 'FullJohnHopkins.json'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 29)
                fakefiles = [
                    "cases_a.h5", "CaseDataFull.h5", "PopulData.h5",
                    "county_population.h5", "migration.h5", "reg_key.h5", "zensus.h5", "FullVacc.h5",
                    "all_county_vacc.h5", "all_state_vacc.h5",
                    "migration_bfa_2020_dim401.h5", "states_testpos.h5",
                    "FullData_DIVI.h5", "county_divi.h5", "cases_a.json",
                    "a_jh.json", "CaseDataFull.json", "PopulData.json",
                    "county_population.json", "migration.json", "reg_key.json", "zensus.json", "FullVacc.json",
                    "all_county_vacc.json", "all_state_vacc.json",
                    "migration_bfa_2020_dim401.json", "states_testpos.json",
                    "FullData_DIVI.json", "county_divi.json"]

                for file in fakefiles:
                    self.assertIn(file, os.listdir(dir_path))
                for file in os.listdir(dir_path):
                    self.assertIn(file, fakefiles)

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 1)
                self.assertEqual(os.listdir(dir_path), ["b_jh.json"])
            elif dir == "France":
                self.assertEqual(len(os.listdir(dir_path)), 1)
                self.assertEqual(os.listdir(dir_path), ["c_jh.json"])

    def test_clean_data_jh_both_endings(self):

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        cd.clean_data(False, False, True, False, False,
                      False, False, False, True, False, False, self.path)
        cd.clean_data(False, False, True, False, False,
                      False, False, False, False, True, False, self.path)

        dir_list = ['Germany']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 1)

        self.assertEqual(os.listdir(self.path),
                         dir_list)

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 28)
                fakefiles = [
                    "cases_a.h5", "CaseDataFull.h5", "PopulData.h5",
                    "county_population.h5", "migration.h5", "reg_key.h5", "zensus.h5", "FullVacc.h5",
                    "all_county_vacc.h5", "all_state_vacc.h5",
                    "migration_bfa_2020_dim401.h5", "states_testpos.h5",
                    "FullData_DIVI.h5", "county_divi.h5", "cases_a.json",
                    "CaseDataFull.json", "PopulData.json",
                    "county_population.json", "migration.json", "reg_key.json", "zensus.json", "FullVacc.json",
                    "all_county_vacc.json", "all_state_vacc.json",
                    "migration_bfa_2020_dim401.json", "states_testpos.json",
                    "FullData_DIVI.json", "county_divi.json"]

                for file in fakefiles:
                    self.assertIn(file, os.listdir(dir_path))
                for file in os.listdir(dir_path):
                    self.assertIn(file, fakefiles)

    def test_clean_txt(self):
        # write file
        dir_path = os.path.join(self.path, 'Germany')
        os.makedirs(dir_path)
        with open(os.path.join(self.path, 'Germany', 'commuter_test.txt'), 'w') as f:
            f.write('foo')
        # check if file is written 
        self.assertEqual(os.listdir(os.path.join(self.path, 'Germany')), ['commuter_test.txt'])
        # delete file and folder
        cd.clean_data(False, False, False,False, False, False, True, False, False, False, True, self.path)
        #check if folder is deleted
        self.assertEqual(os.listdir(self.path), [])
        

    def test_file_not_found_cases(self):

        self.set_dirs_and_files("cases")

        # add different files and folder
        os.makedirs(os.path.join(self.path, "ImportantDir"))

        with open(os.path.join(self.path, "wichtig.py"), 'w') as f:
            f.write('foo')

        cd.clean_data(False, True, False, False, False,
                      False, False, False, True, False, False, self.path)

        self.assertEqual(len(os.listdir(self.path)), 2)
        self.assertEqual(os.listdir(self.path), ["ImportantDir", "wichtig.py"])

    def test_file_not_found_population(self):

        self.set_dirs_and_files("popul")

        # add different files and folder
        os.makedirs(os.path.join(self.path, "ImportantDir"))

        with open(os.path.join(self.path, "wichtig.py"), 'w') as f:
            f.write('foo')

        cd.clean_data(False, False, False, True, False,
                      False, False, False, True, False, False, self.path)

        self.assertEqual(len(os.listdir(self.path)), 2)
        self.assertEqual(os.listdir(self.path), ["ImportantDir", "wichtig.py"])

    def test_file_not_found_jh(self):

        self.set_dirs_and_files("jh")

        # add different files and folder
        os.makedirs(os.path.join(self.path, "ImportantDir"))

        with open(os.path.join(self.path, "wichtig.py"), 'w') as f:
            f.write('foo')

        cd.clean_data(False, False, True, False, False,
                      False, False, False, True, False, False, self.path)

        self.assertEqual(len(os.listdir(self.path)), 2)
        self.assertEqual(os.listdir(self.path), ["ImportantDir", "wichtig.py"])

    def test_file_not_found_divi(self):

        self.set_dirs_and_files("divi")

        # add different files and folder
        os.makedirs(os.path.join(self.path, "ImportantDir"))

        with open(os.path.join(self.path, "wichtig.py"), 'w') as f:
            f.write('foo')

        cd.clean_data(False, False, False, False, True,
                      False, False, False, True, False, False, self.path)

        self.assertEqual(len(os.listdir(self.path)), 2)
        self.assertEqual(os.listdir(self.path), ["ImportantDir", "wichtig.py"])

    def test_file_not_found_vacc(self):

        self.set_dirs_and_files("vacc")

        # add different files and folder
        os.makedirs(os.path.join(self.path, "ImportantDir"))

        with open(os.path.join(self.path, "wichtig.py"), 'w') as f:
            f.write('foo')

        cd.clean_data(False, False, False, False, False,
                      True, False, False, True, False, False, self.path)

        self.assertEqual(len(os.listdir(self.path)), 2)
        self.assertEqual(os.listdir(self.path), ["ImportantDir", "wichtig.py"])

    def test_file_not_found_commuter(self):

        self.set_dirs_and_files("commuter")

        # add different files and folder
        os.makedirs(os.path.join(self.path, "ImportantDir"))

        with open(os.path.join(self.path, "wichtig.py"), 'w') as f:
            f.write('foo')

        cd.clean_data(False, False, False, False, False,
                      False, True, False, True, False, False, self.path)

        self.assertEqual(len(os.listdir(self.path)), 2)
        self.assertEqual(os.listdir(self.path), ["ImportantDir", "wichtig.py"])

    def test_file_not_found_testing(self):

        self.set_dirs_and_files("testing")

        # add different files and folder
        os.makedirs(os.path.join(self.path, "ImportantDir"))

        with open(os.path.join(self.path, "wichtig.py"), 'w') as f:
            f.write('foo')

        cd.clean_data(False, False, False, False, False,
                      False, False, True, True, False, False, self.path)

        self.assertEqual(len(os.listdir(self.path)), 2)
        self.assertEqual(os.listdir(self.path), ["ImportantDir", "wichtig.py"])

    def test_no_files(self):

        # The following should run without any problem.
        # Every error should be cached and passed

        # no data in folder
        # jh
        cd.clean_data(False, False, True, False, False,
                      False, False, False, True, True, True, self.path)

        # population
        cd.clean_data(False, False, False, True, False,
                      False, False, False, True, True, True, self.path)

        # cases
        cd.clean_data(False, True, False, False, False,
                      False, False, False, True, True, True, self.path)

        # divi
        cd.clean_data(False, False, False, False, True,
                      False, False, False, True, True, True, self.path)

        # vaccination
        cd.clean_data(False, False, False, False, False,
                      True, False, False, True, True, True, self.path)

        # commuter
        cd.clean_data(False, False, False, False, False,
                      False, True, False, True, True, True, self.path)

        # testing
        cd.clean_data(False, False, False, False, False,
                      False, False, True, True, True, True, self.path)

    def test_cli_default(self):

        test_args = ["prog"]

        with patch.object(sys, 'argv', test_args):

            [all_data, cases, jh, popul, divi, vacc, commuter,
                testing, json, hdf5, txt] = cd.cli()
                
            print([all_data, cases, jh, popul, hdf5])

            self.assertEqual(all_data, False)
            self.assertEqual(cases, False)
            self.assertEqual(jh, False)
            self.assertEqual(popul, False)
            self.assertEqual(divi, False)
            self.assertEqual(vacc, False)
            self.assertEqual(commuter, False)
            self.assertEqual(testing, False)
            self.assertEqual(json, False)
            self.assertEqual(hdf5, False)
            self.assertEqual(txt, False)

    def test_cli_all(self):

        test_args = ["prog", '--all']

        with patch.object(sys, 'argv', test_args):

            [all_data, cases, jh, popul, divi, vacc, commuter,
                testing, json, hdf5, txt] = cd.cli()

            self.assertEqual(all_data, True)
            self.assertEqual(cases, False)
            self.assertEqual(jh, False)
            self.assertEqual(popul, False)
            self.assertEqual(divi, False)
            self.assertEqual(vacc, False)
            self.assertEqual(commuter, False)
            self.assertEqual(testing, False)
            self.assertEqual(json, False)
            self.assertEqual(hdf5, False)
            self.assertEqual(txt, False)

    def test_cli_cases(self):

        test_args = ["prog", '--cases', '--txt']

        with patch.object(sys, 'argv', test_args):

            [all_data, cases, jh, popul, divi, vacc, commuter,
                testing, json, hdf5, txt] = cd.cli()

            self.assertEqual(all_data, False)
            self.assertEqual(cases, True)
            self.assertEqual(jh, False)
            self.assertEqual(popul, False)
            self.assertEqual(divi, False)
            self.assertEqual(vacc, False)
            self.assertEqual(commuter, False)
            self.assertEqual(testing, False)
            self.assertEqual(json, False)
            self.assertEqual(hdf5, False)
            self.assertEqual(txt, True)

    def test_cli_jh(self):

        test_args = ["prog", '-j', '--json', '--hdf5']

        with patch.object(sys, 'argv', test_args):

            [all_data, cases, jh, popul, divi, vacc, commuter,
                testing, json, hdf5, txt] = cd.cli()

            self.assertEqual(all_data, False)
            self.assertEqual(cases, False)
            self.assertEqual(jh, True)
            self.assertEqual(popul, False)
            self.assertEqual(divi, False)
            self.assertEqual(vacc, False)
            self.assertEqual(commuter, False)
            self.assertEqual(testing, False)
            self.assertEqual(json, True)
            self.assertEqual(hdf5, True)
            self.assertEqual(txt, False)

    def test_cli_popul(self):

        test_args = ['prog', '--population', '-js', '-h5', '-tx']

        with patch.object(sys, 'argv', test_args):

            [all_data, cases, jh, popul, divi, vacc, commuter,
                testing, json, hdf5, txt] = cd.cli()

            self.assertEqual(all_data, False)
            self.assertEqual(cases, False)
            self.assertEqual(jh, False)
            self.assertEqual(popul, True)
            self.assertEqual(divi, False)
            self.assertEqual(vacc, False)
            self.assertEqual(commuter, False)
            self.assertEqual(testing, False)
            self.assertEqual(json, True)
            self.assertEqual(hdf5, True)
            self.assertEqual(txt, True)

    def test_cli_divi_vacc_commuter_testing(self):

        test_args = ['prog', '-d', '-v', '-co', '-t']

        with patch.object(sys, 'argv', test_args):

            [all_data, cases, jh, popul, divi, vacc, commuter,
                testing, json, hdf5, txt] = cd.cli()

            self.assertEqual(all_data, False)
            self.assertEqual(cases, False)
            self.assertEqual(jh, False)
            self.assertEqual(popul, False)
            self.assertEqual(divi, True)
            self.assertEqual(vacc, True)
            self.assertEqual(commuter, True)
            self.assertEqual(testing, True)
            self.assertEqual(json, False)
            self.assertEqual(hdf5, False)
            self.assertEqual(txt, False)

    def test_clean_divi_vacc_commuter_testing_json(self):

        # test cleaning of divi, vaccination, commuter & testing data
        self.set_dirs_and_files('all')

        cd.clean_data(False, False, False, False, True,
                      True, True, True, True, False, False, self.path)

        dir_list = ['Germany', 'Spain', 'France', 'Italy',
                    'US', 'SouthKorea', 'China']
        
        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 11)

        self.assertEqual(
            os.listdir(self.path),
            dir_list +
            ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json',
             'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 23)
                fakefiles = [
                    "cases_a.h5", "a_jh.h5", "CaseDataFull.h5", "PopulData.h5",
                    "county_population.h5", "migration.h5", "reg_key.h5",
                    "zensus.h5", "FullVacc.h5", "FullVacc.h5",
                    "all_county_vacc.h5", "all_state_vacc.h5",
                    "migration_bfa_2020_dim401.h5", "states_testpos.h5",
                    "FullData_DIVI.h5", "county_divi.h5", "a_jh.json",
                    "PopulData.json", "cases_a.json", "county_population.json",
                    "migration.json", "reg_key.json", "zensus.json",
                    "CaseDataFull.json"]

                for file in fakefiles:
                    self.assertIn(file, os.listdir(dir_path))
                for file in os.listdir(dir_path):
                    self.assertIn(file, fakefiles)

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(
                    os.listdir(dir_path),
                    ["b_jh.json", "b_jh.h5"])

            elif dir == "France":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(
                    os.listdir(dir_path),
                    ["c_jh.json", "c_jh.h5"])

if __name__ == '__main__':
    unittest.main()
