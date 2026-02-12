#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Henrik Zunker
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
import io
import os
import tempfile
import unittest
import zipfile
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd

from memilio.epidata.getContactData import (AGE_GROUP_LABELS,
                                            AGE_GROUP_LABELS_RKI,
                                            CONTACT_WORKBOOK_NAME,
                                            list_available_contact_countries,
                                            load_contact_matrix)


class TestGetContactData(unittest.TestCase):
    """Tests for loading contact matrices."""

    def setUp(self):
        """Track temporary files to clean up after each test."""
        self._temp_files = []

    def tearDown(self):
        for path in self._temp_files:
            try:
                os.remove(path)
            except FileNotFoundError:
                pass

    def _create_workbook(self, sheets: dict):
        """Create a temporary xlsx file with provided sheet_name."""
        fd, path = tempfile.mkstemp(suffix=".xlsx")
        os.close(fd)
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            for sheet_name, df in sheets.items():
                df.to_excel(writer, sheet_name=sheet_name, index=False)
        self._temp_files.append(path)
        return path

    def _create_zip_with_workbook(self, sheets: dict):
        """Create a ZIP in memory that contains MUestimates_all_locations_1.xlsx."""
        with io.BytesIO() as buf_zip:
            with zipfile.ZipFile(buf_zip, mode="w") as zf:
                with io.BytesIO() as buf_xlsx:
                    with pd.ExcelWriter(buf_xlsx, engine="openpyxl") as writer:
                        for sheet_name, df in sheets.items():
                            df.to_excel(
                                writer, sheet_name=sheet_name, index=False)
                    zf.writestr(CONTACT_WORKBOOK_NAME, buf_xlsx.getvalue())
            return buf_zip.getvalue()

    def test_list_available_explicit_path(self):
        """Ensure listing works with a custom workbook path."""
        data = pd.DataFrame(np.ones((16, 16)))
        contact_path = self._create_workbook({"Germany": data, "Spain": data})

        countries = list_available_contact_countries(contact_path=contact_path)
        self.assertEqual(set(countries), {"Germany", "Spain"})

    @patch('memilio.epidata.getContactData.requests.get')
    def test_list_available_downloads_zip(self, mock_get):
        """Listing without a path downloads the ZIP and reads workbook."""
        data = pd.DataFrame(np.ones((16, 16)))
        zip_bytes = self._create_zip_with_workbook(
            {"Germany": data, "Spain": data})
        mock_resp = Mock()
        mock_resp.content = zip_bytes
        mock_resp.raise_for_status = Mock()
        mock_get.return_value = mock_resp

        countries = list_available_contact_countries(contact_path=None)
        self.assertEqual(set(countries), {"Germany", "Spain"})
        mock_resp.raise_for_status.assert_called_once()
        mock_get.assert_called_once()

    def test_load_contact_matrix(self):
        """Contact matrix loads with full 16x16 matrix."""
        matrix_values = np.arange(16 * 16).reshape(16, 16)
        df = pd.DataFrame(matrix_values)
        contact_path = self._create_workbook({"Germany": df})

        matrix = load_contact_matrix(
            "Germany", contact_path=contact_path,
            reduce_to_rki_groups=False)

        self.assertEqual(matrix.shape, (16, 16))
        self.assertEqual(list(matrix.columns), AGE_GROUP_LABELS)
        self.assertEqual(list(matrix.index), AGE_GROUP_LABELS)
        self.assertEqual(matrix.iloc[0, 0], 0)
        self.assertEqual(matrix.iloc[-1, -1], 255)

    def test_load_contact_matrix_rki_groups(self):
        """Default is to aggregate to the 6-group RKI age groups."""
        matrix_values = np.arange(16 * 16).reshape(16, 16)
        df = pd.DataFrame(matrix_values)
        contact_path = self._create_workbook({"Germany": df})

        matrix = load_contact_matrix("Germany", contact_path=contact_path)

        self.assertEqual(matrix.shape, (6, 6))
        self.assertEqual(list(matrix.columns), AGE_GROUP_LABELS_RKI)
        self.assertEqual(list(matrix.index), AGE_GROUP_LABELS_RKI)
        # spot-check a few block means
        self.assertEqual(matrix.iloc[0, 0], 0)  # block with only (0,0)
        # rows 1,2 and col 0 -> values 16 and 32 => mean 24
        self.assertEqual(matrix.iloc[1, 0], 24)
        # rows 3,4,5,6 and cols 3,4,5,6 (15-34 vs 15-34)
        sub = matrix_values[3:7, 3:7].mean()
        self.assertEqual(matrix.iloc[2, 2], sub)

    @patch('memilio.epidata.getContactData.requests.get')
    def test_load_contact_matrix_downloads_no_path(self, mock_get):
        """When no path is given the workbook is downloaded from the DOI ZIP."""
        matrix_values = np.arange(16 * 16).reshape(16, 16)
        df = pd.DataFrame(matrix_values)
        zip_bytes = self._create_zip_with_workbook({"Germany": df})
        mock_resp = Mock()
        mock_resp.content = zip_bytes
        mock_resp.raise_for_status = Mock()
        mock_get.return_value = mock_resp

        matrix = load_contact_matrix("Germany", contact_path=None)
        self.assertEqual(matrix.shape, (6, 6))
        self.assertEqual(list(matrix.columns), AGE_GROUP_LABELS_RKI)
        self.assertEqual(list(matrix.index), AGE_GROUP_LABELS_RKI)
        self.assertEqual(matrix.iloc[0, 0], 0)
        self.assertEqual(matrix.iloc[1, 0], 24)
        mock_resp.raise_for_status.assert_called_once()
        mock_get.assert_called_once()

    def test_load_contact_matrix_unknown_country_raises(self):
        """Loading a non-existing sheet raises ValueError."""
        data = pd.DataFrame(np.zeros((16, 16)))
        contact_path = self._create_workbook({"Germany": data})

        with self.assertRaises(ValueError):
            load_contact_matrix("France", contact_path=contact_path)


if __name__ == '__main__':
    unittest.main()
