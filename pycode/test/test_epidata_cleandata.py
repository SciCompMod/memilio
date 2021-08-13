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

from epidemiology.epidata import cleanData as cd
from epidemiology.epidata import defaultDict as dd


class Test_cleanData(fake_filesystem_unittest.TestCase):

    path = '/home/x'

    def setUp(self):
        self.setUpPyfakefs()

    def set_dirs_and_files(self, what):

        dir_dic_all = { 'Germany' : ["a_rki", "a_jh", "FullRKI", "PopulData", "FullDataB", "FullDataL"],
                        'Spain': ["b_jh"],
                        'France': ["c_jh"],
                        'Italy': ["d_jh"],
                        'US' : ["e_jh"],
                        'SouthKorea' : ["f_jh"],
                        'China' : ["g_jh"]}

        dir_dic_rki = {'Germany': ["a_rki", "FullRKI"]}

        dir_dic_popul = {'Germany': ["PopulData", "FullDataB", "FullDataL"]}

        dir_dic_jh = {'Germany': ["a_jh"],
                      'Spain': ["b_jh"],
                      'France': ["c_jh"],
                      'Italy': ["d_jh"],
                      'US': ["e_jh"],
                      'SouthKorea': ["f_jh"],
                      'China': ["g_jh"]}

        ending_all = [".json", ".h5"]
        ending_json = [".json"]

        dir_choose = {"all": dir_dic_all,
                      "rki": dir_dic_rki,
                      "jh": dir_dic_jh,
                      "popul": dir_dic_popul
                      }

        ending_choose= {"all": ending_all,
                        "rki": ending_json,
                        "jh": ending_json,
                        "popul": ending_json
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

        dir_list = ['Germany', 'Spain', 'France', 'Italy', 'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 11)

        self.assertEqual(os.listdir(self.path),
                         dir_list + ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json', 'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 12)
                self.assertEqual(os.listdir(dir_path),
                                 ["a_rki.json", "a_rki.h5", "a_jh.json", "a_jh.h5", "FullRKI.json", "FullRKI.h5",
                                  "PopulData.json", "PopulData.h5", "FullDataB.json", "FullDataB.h5",
                                  "FullDataL.json", "FullDataL.h5"])

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(os.listdir(dir_path), ["b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    # generate folder and files
    def test_clean_data_all_should_delete_all(self):

        self.set_dirs_and_files("all")

        cd.clean_data(True, False, False, False, False, self.path)

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

        cd.clean_data(True, False, False, False, False, self.path)

        # Should delete everything
        self.assertEqual(len(os.listdir(self.path)), 3)
        self.assertEqual(os.listdir(self.path), ["China", "ImportantDir", "wichtig.py"])
        self.assertEqual(len(os.listdir(dir_path)), 1)
        self.assertEqual(os.listdir(dir_path), ["secret.txt"])

    def test_clean_data_rki(self):

        self.set_dirs_and_files("all")

        cd.clean_data(False, True, False, False, False, self.path)

        dir_list = ['Germany', 'Spain', 'France', 'Italy', 'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 11)

        self.assertEqual(os.listdir(self.path),
                         dir_list + ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json', 'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 10)
                self.assertEqual(os.listdir(dir_path),
                                 ["a_rki.h5", "a_jh.json", "a_jh.h5", "FullRKI.h5",
                                  "PopulData.json", "PopulData.h5", "FullDataB.json", "FullDataB.h5",
                                  "FullDataL.json", "FullDataL.h5"])

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(os.listdir(dir_path), ["b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    def test_clean_data_rki_h5(self):

        self.set_dirs_and_files("all")

        cd.clean_data(False, True, False, False, True, self.path)

        dir_list = ['Germany', 'Spain', 'France', 'Italy', 'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 11)

        self.assertEqual(os.listdir(self.path),
                         dir_list + ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json', 'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 10)
                self.assertEqual(os.listdir(dir_path),
                                 ["a_rki.json", "a_jh.json", "a_jh.h5", "FullRKI.json",
                                  "PopulData.json", "PopulData.h5", "FullDataB.json", "FullDataB.h5",
                                  "FullDataL.json", "FullDataL.h5"])

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(os.listdir(dir_path), ["b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    def test_clean_data_rki_del_dir(self):

        self.set_dirs_and_files("all")

        dir_path = os.path.join(self.path, "Germany")
        files = os.listdir(dir_path)

        # delete all files except which will be deleted
        for item in files:
            if item == "a_rki.json" or item == "FullRKI.json":
                continue
            else:
                os.remove(os.path.join(dir_path, item))

        self.assertEqual(len(os.listdir(dir_path)), 2)
        self.assertEqual(os.listdir(dir_path), ["a_rki.json", "FullRKI.json"])

        cd.clean_data(False, True, False, False, False, self.path)

        dir_list = ['Spain', 'France', 'Italy', 'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure
        # Germany should be deleted, because no files where left after deletion
        self.assertEqual(len(os.listdir(self.path)), 10)

        self.assertEqual(os.listdir(self.path),
                         dir_list + ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json', 'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(os.listdir(dir_path), ["b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    def test_clean_data_population(self):

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        cd.clean_data(False, False, False, True, False, self.path)

        dir_list = ['Germany', 'Spain', 'France', 'Italy', 'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 11)

        self.assertEqual(os.listdir(self.path),
                         dir_list + ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json', 'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 9)
                self.assertEqual(os.listdir(dir_path),
                                 ["a_rki.json", "a_rki.h5", "a_jh.json", "a_jh.h5", "FullRKI.json", "FullRKI.h5",
                                  "PopulData.h5", "FullDataB.h5", "FullDataL.h5"])

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(os.listdir(dir_path), ["b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    def test_clean_data_population_hdf5(self):

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        cd.clean_data(False, False, False, True, True, self.path)

        dir_list = ['Germany', 'Spain', 'France', 'Italy', 'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 11)

        self.assertEqual(os.listdir(self.path),
                         dir_list + ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json', 'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 9)
                self.assertEqual(os.listdir(dir_path),
                                 ["a_rki.json", "a_rki.h5", "a_jh.json", "a_jh.h5", "FullRKI.json", "FullRKI.h5",
                                  "PopulData.json", "FullDataB.json", "FullDataL.json"])

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(os.listdir(dir_path), ["b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    def test_clean_data_population_del_dir(self):

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        dir_path = os.path.join(self.path, "Germany")
        files = os.listdir(dir_path)

        # delete all files except which will be deleted
        for item in files:
            if item == "PopulData.json" or item == "FullDataB.json" or item == "FullDataL.json":
                continue
            else:
                os.remove(os.path.join(dir_path, item))

        self.assertEqual(len(os.listdir(dir_path)), 3)

        cd.clean_data(False, False, False, True, False, self.path)

        dir_list = ['Spain', 'France', 'Italy', 'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 10)

        self.assertEqual(os.listdir(self.path),
                         dir_list + ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json', 'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(os.listdir(dir_path), ["b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    def test_all_false(self):

        cd.clean_data(False, False, False, False, False, self.path)

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        dir_list = ['Germany', 'Spain', 'France', 'Italy', 'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 11)

        self.assertEqual(os.listdir(self.path),
                         dir_list + ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json', 'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 12)
                self.assertEqual(os.listdir(dir_path),
                                 ["a_rki.json", "a_rki.h5", "a_jh.json", "a_jh.h5", "FullRKI.json", "FullRKI.h5",
                                  "PopulData.json", "PopulData.h5", "FullDataB.json", "FullDataB.h5",
                                  "FullDataL.json", "FullDataL.h5"])

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(os.listdir(dir_path), ["b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    def test_wrong_path(self):

        cd.clean_data(True, False, False, False, False, "/home/y")

        # TODO add some test: but what? - nothing is happening in this case


    def test_clean_data_jh(self):

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        cd.clean_data(False, False, True, False, False, self.path)

        dir_list = ['Germany', 'Spain', 'France', 'Italy', 'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 9)

        self.assertEqual(os.listdir(self.path),
                         dir_list + ['all_jh.h5', 'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 11)
                self.assertEqual(os.listdir(dir_path),
                                 ["a_rki.json", "a_rki.h5", "a_jh.h5", "FullRKI.json", "FullRKI.h5",
                                  "PopulData.json", "PopulData.h5", "FullDataB.json", "FullDataB.h5",
                                  "FullDataL.json", "FullDataL.h5"])

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 1)
                self.assertEqual(os.listdir(dir_path), ["b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 1)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.h5"])

    def test_clean_data_jh_hdf5(self):

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        cd.clean_data(False, False, True, False, True, self.path)

        dir_list = ['Germany', 'Spain', 'France', 'Italy', 'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 9)

        self.assertEqual(os.listdir(self.path),
                         dir_list + ['all_jh.json', 'FullJohnHopkins.json'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 11)
                self.assertEqual(os.listdir(dir_path),
                                 ["a_rki.json", "a_rki.h5", "a_jh.json", "FullRKI.json", "FullRKI.h5",
                                  "PopulData.json", "PopulData.h5", "FullDataB.json", "FullDataB.h5",
                                  "FullDataL.json", "FullDataL.h5"])

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 1)
                self.assertEqual(os.listdir(dir_path), ["b_jh.json"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 1)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json"])

    def test_clean_data_jh_both_endings(self):

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        cd.clean_data(False, False, True, False, False, self.path)
        cd.clean_data(False, False, True, False, True, self.path)

        dir_list = ['Germany']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 1)

        self.assertEqual(os.listdir(self.path),
                         dir_list)

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 10)
                self.assertEqual(os.listdir(dir_path),
                                 ["a_rki.json", "a_rki.h5","FullRKI.json", "FullRKI.h5",
                                  "PopulData.json", "PopulData.h5", "FullDataB.json", "FullDataB.h5",
                                  "FullDataL.json", "FullDataL.h5"])

    def test_clean_data_rki_johns_hopkins(self):

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        cd.clean_data(False, True, True, False, False, self.path)

        dir_list = ['Germany', 'Spain', 'France', 'Italy', 'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 9)

        self.assertEqual(os.listdir(self.path),
                         dir_list + ['all_jh.h5', 'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 9)
                self.assertEqual(os.listdir(dir_path),
                                 ["a_rki.h5", "a_jh.h5", "FullRKI.h5",
                                  "PopulData.json", "PopulData.h5", "FullDataB.json", "FullDataB.h5",
                                  "FullDataL.json", "FullDataL.h5"])

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 1)
                self.assertEqual(os.listdir(dir_path), ["b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 1)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.h5"])

    def test_clean_data_rki_john_hokins_spain_population(self):

        # test if writte fct works as expected

        self.set_dirs_and_files("all")

        cd.clean_data(False, True, True, True, False, self.path)

        dir_list = ['Germany', 'Spain', 'France', 'Italy', 'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 9)

        self.assertEqual(os.listdir(self.path),
                         dir_list + ['all_jh.h5', 'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Germany":
                self.assertEqual(len(os.listdir(dir_path)), 6)
                self.assertEqual(os.listdir(dir_path),
                                 ["a_rki.h5", "a_jh.h5", "FullRKI.h5",
                                  "PopulData.h5", "FullDataB.h5", "FullDataL.h5"])

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 1)
                self.assertEqual(os.listdir(dir_path), ["b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 1)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.h5"])

    def test_file_not_found_rki(self):

        self.set_dirs_and_files("rki")

        # add different files and folder
        os.makedirs(os.path.join(self.path, "ImportantDir"))

        with open(os.path.join(self.path, "wichtig.py"), 'w') as f:
            f.write('foo')

        cd.clean_data(False, True, False, False, False, self.path)

        self.assertEqual(len(os.listdir(self.path)), 2)
        self.assertEqual(os.listdir(self.path), ["ImportantDir", "wichtig.py"])

    def test_file_not_found_population(self):

        self.set_dirs_and_files("popul")

        # add different files and folder
        os.makedirs(os.path.join(self.path, "ImportantDir"))

        with open(os.path.join(self.path, "wichtig.py"), 'w') as f:
            f.write('foo')

        cd.clean_data(False, False, False, True, False, self.path)

        self.assertEqual(len(os.listdir(self.path)), 2)
        self.assertEqual(os.listdir(self.path), ["ImportantDir", "wichtig.py"])

    def test_file_not_found_jh(self):

        self.set_dirs_and_files("jh")

        # add different files and folder
        os.makedirs(os.path.join(self.path, "ImportantDir"))

        with open(os.path.join(self.path, "wichtig.py"), 'w') as f:
            f.write('foo')

        cd.clean_data(False, False, True, False, False, self.path)

        self.assertEqual(len(os.listdir(self.path)), 2)
        self.assertEqual(os.listdir(self.path), ["ImportantDir", "wichtig.py"])

    def test_no_files(self):

        # The following should run without any problem.
        # Every error should be cached and passed

        # no data in folder
        cd.clean_data(False, False, True, False, False, self.path)

        # population
        cd.clean_data(False, False, False, True, False, self.path)

        # rki
        cd.clean_data(False, True, False, False, False, self.path)

    def test_cli_default(self):

        out_path_default = dd.defaultDict['out_folder']

        test_args =  ["prog"]

        with patch.object(sys, 'argv', test_args):

            [all_data, rki, jh, popul, hdf5, out_path] = cd.cli()

            print([all_data, rki, jh, popul, hdf5, out_path])

            self.assertEqual(all_data, False)
            self.assertEqual(rki, False)
            self.assertEqual(jh, False)
            self.assertEqual(popul, False)
            self.assertEqual(hdf5, False)
            self.assertEqual(out_path, out_path_default)

    def test_cli_folder(self):

        folder = "some_folder"
        test_args = ["prog", '--out_path', folder]

        with patch.object(sys, 'argv', test_args):

            [all_data, rki, jh, popul, hdf5, out_path] = cd.cli()

            self.assertEqual(all_data, False)
            self.assertEqual(rki, False)
            self.assertEqual(jh, False)
            self.assertEqual(popul, False)
            self.assertEqual(hdf5, False)
            self.assertEqual(out_path, folder)

    def test_cli_all(self):

        out_path_default = dd.defaultDict['out_folder']

        test_args =  ["prog", '--all']

        with patch.object(sys, 'argv', test_args):

            [all_data, rki, jh, popul, hdf5, out_path] = cd.cli()

            self.assertEqual(all_data, True)
            self.assertEqual(rki, False)
            self.assertEqual(jh, False)
            self.assertEqual(popul, False)
            self.assertEqual(hdf5, False)
            self.assertEqual(out_path, out_path_default)

    def test_cli_rki(self):

        out_path_default = dd.defaultDict['out_folder']

        test_args =  ["prog", '--rki']

        with patch.object(sys, 'argv', test_args):

            [all_data, rki, jh, popul, hdf5, out_path] = cd.cli()

            self.assertEqual(all_data, False)
            self.assertEqual(rki, True)
            self.assertEqual(jh, False)
            self.assertEqual(popul, False)
            self.assertEqual(hdf5, False)
            self.assertEqual(out_path, out_path_default)

    def test_cli_jh(self):

        out_path_default = dd.defaultDict['out_folder']

        test_args = ["prog", '-j', '--hdf5']

        with patch.object(sys, 'argv', test_args):

            [all_data, rki, jh, popul, hdf5, out_path] = cd.cli()

            self.assertEqual(all_data, False)
            self.assertEqual(rki, False)
            self.assertEqual(jh, True)
            self.assertEqual(popul, False)
            self.assertEqual(hdf5, True)
            self.assertEqual(out_path, out_path_default)


    def test_cli_popul(self):

        out_path_default = dd.defaultDict['out_folder']

        test_args = ['prog', '--population', '-h5']

        with patch.object(sys, 'argv', test_args):

            [all_data, rki, jh, popul, hdf5, out_path] = cd.cli()

            self.assertEqual(all_data, False)
            self.assertEqual(rki, False)
            self.assertEqual(jh, False)
            self.assertEqual(popul, True)
            self.assertEqual(hdf5, True)
            self.assertEqual(out_path, out_path_default)


if __name__ == '__main__':
    unittest.main()
