import unittest
from pyfakefs import fake_filesystem_unittest
import os
import sys
from io import StringIO
from unittest.mock import patch
from datetime import date, datetime
import pandas as pd

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd

class Test_getDataIntoPandasDataFrame(fake_filesystem_unittest.TestCase):

    path = '/home/x/'

    def setUp(self):
        self.setUpPyfakefs()

    def test_cli_correct_default(self):
        
        out_path_default = dd.defaultDict['out_folder']
        out_path_default = os.path.join(out_path_default, 'pydata')
        
        [read_data, out_form, out_folder] = gd.cli("spain")

        assert read_data == dd.defaultDict['read_data']
        assert out_form == dd.defaultDict['out_form']
        assert out_folder == out_path_default

        [read_data, out_form, out_folder] = gd.cli("population")

        assert read_data == dd.defaultDict['read_data']
        assert out_form == dd.defaultDict['out_form']
        assert out_folder == out_path_default

        [read_data, out_form, out_folder] = gd.cli("jh")

        assert read_data == dd.defaultDict['read_data']
        assert out_form == dd.defaultDict['out_form']
        assert out_folder == out_path_default

        [read_data, out_form, out_folder, make_plot] = gd.cli("rki")

        assert read_data == dd.defaultDict['read_data']
        assert out_form == dd.defaultDict['out_form']
        assert out_folder == out_path_default
        assert make_plot == dd.defaultDict['make_plot']

        [read_data, out_form, out_folder, make_plot] = gd.cli("rkiest")

        assert read_data == dd.defaultDict['read_data']
        assert out_form == dd.defaultDict['out_form']
        assert out_folder == out_path_default
        assert make_plot == dd.defaultDict['make_plot']

        [read_data, out_form, out_folder, end_date, start_date, update] = gd.cli("divi")

        assert read_data == dd.defaultDict['read_data']
        assert out_form == dd.defaultDict['out_form']
        assert out_folder == out_path_default
        assert end_date == dd.defaultDict['end_date']
        assert start_date == dd.defaultDict['start_date']
        assert update == dd.defaultDict['update_data']

        [read_data, out_form, out_folder, end_date, make_plot, start_date, update] = gd.cli("all")

        assert read_data == dd.defaultDict['read_data']
        assert out_form == dd.defaultDict['out_form']
        assert out_folder == out_path_default
        assert end_date == dd.defaultDict['end_date']
        assert start_date == dd.defaultDict['start_date']
        assert update == dd.defaultDict['update_data']
        assert make_plot == dd.defaultDict['make_plot']

    @patch('sys.stderr', new_callable=StringIO)
    def test_cli_correct_raise_exit(self, mock_stderr):

        with self.assertRaises(SystemExit) as cm:
            gd.cli("wrong_key")

        the_exception = cm.exception
        self.assertEqual(the_exception.code, "Wrong key or cli_dict.")

        sys.argv[1:] = ['-ff', 'wrong_format']

        with self.assertRaises(SystemExit) as cm:
            gd.cli("rki")
        self.assertRegexpMatches(mock_stderr.getvalue(), r"invalid choice")

        # update and read_data together does not work
        sys.argv[1:] = ['-u', '-r']
        with self.assertRaises(SystemExit) as cm:
            gd.cli("divi")

        the_exception = cm.exception
        self.assertEqual(the_exception.code, "You called the program with '--read-from-disk' and '--update'." +
                         "Please choose just one. Both together is not possible.")

        sys.argv[1:] = ['-u', '-r']
        with self.assertRaises(SystemExit) as cm:
            gd.cli("all")

        the_exception = cm.exception
        self.assertEqual(the_exception.code, "You called the program with '--read-from-disk' and '--update'." +
                         "Please choose just one. Both together is not possible.")

        sys.argv[1:] = ['--update',]
        with self.assertRaises(SystemExit) as cm:
            gd.cli("rki")
        self.assertRegexpMatches(mock_stderr.getvalue(), r"unrecognized arguments")

        sys.argv[1:] = ['--start_date']
        with self.assertRaises(SystemExit) as cm:
            gd.cli("rki")
        self.assertRegexpMatches(mock_stderr.getvalue(), r"unrecognized arguments")

        sys.argv[1:] = ['--end_date']
        with self.assertRaises(SystemExit) as cm:
            gd.cli("rki")
        self.assertRegexpMatches(mock_stderr.getvalue(), r"unrecognized arguments")

        sys.argv[1:] = ['--plot']
        with self.assertRaises(SystemExit) as cm:
            gd.cli("divi")
        self.assertRegexpMatches(mock_stderr.getvalue(), r"unrecognized arguments")

        sys.argv[1:] = ['--start-date', '2020,3,24']
        with self.assertRaises(SystemExit) as cm:
            gd.cli("divi")

        sys.argv[1:] = ['--end-date', '24-3-2020']
        with self.assertRaises(SystemExit) as cm:
            gd.cli("divi")


    def test_cli_set_different_values(self):

        folder = "some_folder"

        sys.argv[1:] = ['--read-from-disk', '--out-path', folder, '--file-format', 'hdf5']

        [read_data, out_form, out_folder] = gd.cli("spain")

        assert read_data == True
        assert out_form == 'hdf5'
        assert out_folder == "some_folder"

        [read_data, out_form, out_folder] = gd.cli("population")

        assert read_data == True
        assert out_form == 'hdf5'
        assert out_folder == "some_folder"

        [read_data, out_form, out_folder] = gd.cli("jh")

        assert read_data == True
        assert out_form == 'hdf5'
        assert out_folder == "some_folder"

        sys.argv[1:] = ['--read-from-disk', '--out-path', folder, '--file-format', 'hdf5', '--plot']

        [read_data, out_form, out_folder, make_plot] = gd.cli("rki")

        assert read_data == True
        assert out_form == 'hdf5'
        assert out_folder == "some_folder"
        assert make_plot == True

        sys.argv[1:] = ['--read-from-disk', '--out-path', folder, '--file-format', 'json', '--plot']

        [read_data, out_form, out_folder, make_plot] = gd.cli("rkiest")

        assert read_data == True
        assert out_form == 'json'
        assert out_folder == "some_folder"
        assert make_plot == True

        sys.argv[1:] = ['--out-path', folder, '--file-format', 'json', '--update',
                        '--start-date', '2020-11-24', '--end-date', '2020-11-26']

        [read_data, out_form, out_folder, end_date, start_date, update] = gd.cli("divi")

        assert read_data == dd.defaultDict['read_data']
        assert out_form == 'json'
        assert out_folder == "some_folder"
        assert end_date == date(2020,11,26)
        assert start_date == date(2020,11,24)
        assert update == True

        sys.argv[1:] = ['--out-path', folder, '--file-format', 'json', '--update', '--plot',
                        '--start-date', '2020-11-24', '--end-date', '2020-11-26']

        [read_data, out_form, out_folder, end_date, make_plot, start_date, update] = gd.cli("all")

        assert read_data == dd.defaultDict['read_data']
        assert out_form == 'json'
        assert out_folder == "some_folder"
        assert end_date == date(2020, 11, 26)
        assert start_date == date(2020, 11, 24)
        assert update == True
        assert make_plot == True

    def test_check_dir(self):

        self.assertFalse(os.path.exists(self.path))

        gd.check_dir(self.path)
        self.assertTrue(os.path.exists(self.path))

    test_string_json = ("""[{"Date":1606348800000,"col2":"d1"},{"Date":1595721600000,"col2":"d2"}]""")
    test_string_json_timeasstring = ("""[{"Date":"2020-11-26","col2":"d1"},{"Date":"2020-07-26","col2":"d2"}]""")

    def test_write_dataframe(self):

        gd.check_dir(self.path)

        d1 = date(2020,11,26)
        d1 = datetime.combine(d1, datetime.min.time())
        d2 = date(2020,7,26)
        d2 = datetime.combine(d2, datetime.min.time())

        d = {'Date': [d1, d2], 'col2': ["d1", "d2"]}
        df = pd.DataFrame(data=d)

        gd.write_dataframe(df, self.path, "test_json", 'json')

        file = "test_json.json"

        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(os.listdir(self.path), [file])

        file_with_path = os.path.join(self.path, file)
        f = open(file_with_path, "r")
        fread = f.read()
        self.assertEqual(fread, self.test_string_json)

        d1_in_sec = (d1- datetime(1970, 1, 1)).total_seconds()
        # date in milliseconds
        assert d1_in_sec*1000 == 1606348800000

        d2_in_sec = (d2 - datetime(1970, 1, 1)).total_seconds()
        assert d2_in_sec*1000 == 1595721600000

        gd.write_dataframe(df, self.path, "test_json_timeasstring", 'json_timeasstring')

        file2 = "test_json_timeasstring.json"

        self.assertEqual(len(os.listdir(self.path)), 2)
        self.assertEqual(os.listdir(self.path).sort(), [file, file2].sort())

        file_with_path = os.path.join(self.path, file2)
        f = open(file_with_path, "r")
        fread = f.read()
        self.assertEqual(fread, self.test_string_json_timeasstring)

        # TODO: Why does check of hdf5 not work?
        #gd.write_dataframe(df, self.path, "test_hdf", 'hdf5')
        # file = "test_hdf.h5"
        # self.assertEqual(len(os.listdir(self.path)), 3)
        # self.assertEqual(os.listdir(self.path), [file])


if __name__ == '__main__':
    unittest.main()