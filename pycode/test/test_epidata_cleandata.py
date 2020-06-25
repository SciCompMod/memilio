import unittest
from pyfakefs import fake_filesystem_unittest
import os

from epidemiology.epidata import cleanData as cd
#from unittest.mock import patch, mock_open

from epidemiology.epidata import getRKIData as grkid

class Test_cleanData(fake_filesystem_unittest.TestCase):

    path = '/home/x'

    def setUp(self):
        self.setUpPyfakefs()


    def set_dirs_and_files(self):

        dir_dic = { 'Germany' : ["a_rki", "a_jh", "FullRKI", "PopulData", "FullDataB", "FullDataL"],
                    'Spain': ["a_spain", "b_jh"],
                    'France': ["c_jh"],
                    'Italy': ["d_jh"],
                    'US' : ["e_jh"],
                    'SouthKorea' : ["f_jh"],
                    'China' : ["g_jh"],
                    }

        file_list = ["all_jh", "FullJohnHopkins"]
        ending = [".json", ".h5"]

        # make folders
        for key in dir_dic:
            dir_path = os.path.join(self.path, key)
            os.makedirs(dir_path)

            # make files
            for file in dir_dic[key]:
                for e in ending:
                    with open(os.path.join(dir_path, file + e), 'w') as f:
                        f.write('foo')

        for file in file_list:
            for e in ending:
                with open(os.path.join(self.path, file + e), 'w') as f:
                    f.write('foo')

    def test_set_dirs_and_files(self):

        # test if writte fct works as expected

        self.set_dirs_and_files()

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
                self.assertEqual(len(os.listdir(dir_path)), 4)
                self.assertEqual(os.listdir(dir_path), ["a_spain.json", "a_spain.h5", "b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    # generate folder and files
    def test_clean_data_all_should_delete_all(self):

        self.set_dirs_and_files()

        cd.clean_data(True, False, False, False, False, False, self.path)

        # Should delete everything
        self.assertEqual(len(os.listdir(self.path)), 0)

    def test_clean_data_all_should_not_delete_all(self):

        self.set_dirs_and_files()

        # add different files and folder
        os.makedirs(os.path.join(self.path, "ImportantDir"))

        with open(os.path.join(self.path, "wichtig.py"), 'w') as f:
            f.write('foo')

        dir_path = os.path.join(self.path, "China")
        with open(os.path.join(dir_path, "secret.txt"), 'w') as f:
            f.write('foo')

        cd.clean_data(True, False, False, False, False, False, self.path)

        # Should delete everything
        self.assertEqual(len(os.listdir(self.path)), 3)
        self.assertEqual(os.listdir(self.path), ["China", "ImportantDir", "wichtig.py"])
        self.assertEqual(len(os.listdir(dir_path)), 1)
        self.assertEqual(os.listdir(dir_path), ["secret.txt"])

    def test_clean_data_rki(self):

        self.set_dirs_and_files()

        cd.clean_data(False, True, False, False, False, False, self.path)

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
                self.assertEqual(len(os.listdir(dir_path)), 4)
                self.assertEqual(os.listdir(dir_path), ["a_spain.json", "a_spain.h5", "b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    def test_clean_data_rki_h5(self):

        self.set_dirs_and_files()

        cd.clean_data(False, True, False, False, False, True, self.path)

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
                self.assertEqual(len(os.listdir(dir_path)), 4)
                self.assertEqual(os.listdir(dir_path), ["a_spain.json", "a_spain.h5", "b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    def test_clean_data_rki_del_dir(self):

        self.set_dirs_and_files()

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

        cd.clean_data(False, True, False, False, False, False, self.path)

        dir_list = ['Spain', 'France', 'Italy', 'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure
        # Germany should be deleted, because no files where left after deletion
        self.assertEqual(len(os.listdir(self.path)), 10)

        self.assertEqual(os.listdir(self.path),
                         dir_list + ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json', 'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 4)
                self.assertEqual(os.listdir(dir_path), ["a_spain.json", "a_spain.h5", "b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    def test_clean_data_population(self):

        # test if writte fct works as expected

        self.set_dirs_and_files()

        cd.clean_data(False, False, False, False, True, False, self.path)

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
                self.assertEqual(len(os.listdir(dir_path)), 4)
                self.assertEqual(os.listdir(dir_path), ["a_spain.json", "a_spain.h5", "b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    def test_clean_data_population_hdf5(self):

        # test if writte fct works as expected

        self.set_dirs_and_files()

        cd.clean_data(False, False, False, False, True, True, self.path)

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
                self.assertEqual(len(os.listdir(dir_path)), 4)
                self.assertEqual(os.listdir(dir_path), ["a_spain.json", "a_spain.h5", "b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])


    def test_clean_data_population_del_dir(self):

        # test if writte fct works as expected

        self.set_dirs_and_files()

        dir_path = os.path.join(self.path, "Germany")
        files = os.listdir(dir_path)

        # delete all files except which will be deleted
        for item in files:
            if item == "PopulData.json" or item == "FullDataB.json" or item == "FullDataL.json":
                continue
            else:
                os.remove(os.path.join(dir_path, item))

        self.assertEqual(len(os.listdir(dir_path)), 3)

        cd.clean_data(False, False, False, False, True, False, self.path)

        dir_list = ['Spain', 'France', 'Italy', 'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 10)

        self.assertEqual(os.listdir(self.path),
                         dir_list + ['all_jh.json', 'all_jh.h5', 'FullJohnHopkins.json', 'FullJohnHopkins.h5'])

        for dir in dir_list:
            dir_path = os.path.join(self.path, dir)

            if dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 4)
                self.assertEqual(os.listdir(dir_path), ["a_spain.json", "a_spain.h5", "b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    def test_all_false(self):

        cd.clean_data(False, False, False, False, False, False, self.path)

        # test if writte fct works as expected

        self.set_dirs_and_files()

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
                self.assertEqual(len(os.listdir(dir_path)), 4)
                self.assertEqual(os.listdir(dir_path), ["a_spain.json", "a_spain.h5", "b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    def test_wrong_path(self):

        cd.clean_data(True, False, False, False, False, False, "/home/y")

        # TODO add some test: but what? - nothing is happening in this case


    def test_clean_data_jh(self):

        # test if writte fct works as expected

        self.set_dirs_and_files()

        cd.clean_data(False, False, True, False, False, False, self.path)

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
                self.assertEqual(len(os.listdir(dir_path)), 3)
                self.assertEqual(os.listdir(dir_path), ["a_spain.json", "a_spain.h5", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 1)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.h5"])

    def test_clean_data_jh_hdf5(self):

        # test if writte fct works as expected

        self.set_dirs_and_files()

        cd.clean_data(False, False, True, False, False, True, self.path)

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
                self.assertEqual(len(os.listdir(dir_path)), 3)
                self.assertEqual(os.listdir(dir_path), ["a_spain.json", "a_spain.h5", "b_jh.json"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 1)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json"])

    def test_clean_data_jh_both_endings(self):

        # test if writte fct works as expected

        self.set_dirs_and_files()

        cd.clean_data(False, False, True, False, False, False, self.path)
        cd.clean_data(False, False, True, False, False, True, self.path)

        dir_list = ['Germany', 'Spain']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 2)

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

            elif dir == "Spain":
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(os.listdir(dir_path), ["a_spain.json", "a_spain.h5"])

    def test_clean_data_spain(self):

        # test if writte fct works as expected

        self.set_dirs_and_files()

        cd.clean_data(False, False, False, True, False, False, self.path)

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
                self.assertEqual(len(os.listdir(dir_path)), 3)
                self.assertEqual(os.listdir(dir_path), ["a_spain.h5", "b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    def test_clean_data_spain_hdf5(self):

        # test if writte fct works as expected

        self.set_dirs_and_files()

        cd.clean_data(False, False, False, True, False, True, self.path)

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
                self.assertEqual(len(os.listdir(dir_path)), 3)
                self.assertEqual(os.listdir(dir_path), ["a_spain.json", "b_jh.json", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    def test_clean_data_spain_del_dir(self):

        # test if writte fct works as expected

        self.set_dirs_and_files()

        dir_path = os.path.join(self.path, "Spain")
        files = os.listdir(dir_path)

        # delete all files except which will be deleted
        for item in files:
            if item == "a_spain.json":
                continue
            else:
                os.remove(os.path.join(dir_path, item))

        cd.clean_data(False, False, False, True, False, False, self.path)

        dir_list = ['Germany', 'France', 'Italy', 'US', 'SouthKorea', 'China']

        # Test wanted folder and file structure

        self.assertEqual(len(os.listdir(self.path)), 10)

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

            else:
                self.assertEqual(len(os.listdir(dir_path)), 2)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.json", "c_jh.h5"])

    def test_clean_data_rki_john_hokins(self):

        # test if writte fct works as expected

        self.set_dirs_and_files()

        cd.clean_data(False, True, True, False, False, False, self.path)

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
                self.assertEqual(len(os.listdir(dir_path)), 3)
                self.assertEqual(os.listdir(dir_path), ["a_spain.json", "a_spain.h5", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 1)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.h5"])

    def test_clean_data_rki_john_hokins_spain(self):

        # test if writte fct works as expected

        self.set_dirs_and_files()

        cd.clean_data(False, True, True, True, False, False, self.path)

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
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(os.listdir(dir_path), ["a_spain.h5", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 1)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.h5"])

    def test_clean_data_rki_john_hokins_spain_population(self):

        # test if writte fct works as expected

        self.set_dirs_and_files()

        cd.clean_data(False, True, True, True, True, False, self.path)

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
                self.assertEqual(len(os.listdir(dir_path)), 2)
                self.assertEqual(os.listdir(dir_path), ["a_spain.h5", "b_jh.h5"])
            else:
                self.assertEqual(len(os.listdir(dir_path)), 1)
                if dir == "France":
                    self.assertEqual(os.listdir(dir_path), ["c_jh.h5"])


if __name__ == '__main__':
    unittest.main()