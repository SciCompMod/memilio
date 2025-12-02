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
"""
:strong:`getContactData.py`

Load an age-structured contact matrix for a chosen country based on
Prem et al., 2017. The module can download the supporting ZIP from
https://doi.org/10.1371/journal.pcbi.1005697.s002 (contains the
``MUestimates_all_locations_1.xlsx`` workbook) or read a defined local
workbook path. By default, downloads are done in memory and no
files are written.
"""

import io
import os
import zipfile
from typing import Iterable, List, Optional

import pandas as pd
import requests

CONTACT_ZIP_URL = (
    "https://journals.plos.org/ploscompbiol/article/file"
    "?id=10.1371/journal.pcbi.1005697.s002&type=supplementary"
)
CONTACT_WORKBOOK_NAME = "MUestimates_all_locations_1.xlsx"

# Only kept for tests or explicit overrides; regular calls should pass
# contact_path=None to trigger download.
DEFAULT_CONTACT_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                 CONTACT_WORKBOOK_NAME))

AGE_GROUP_LABELS = [
    "0-4",
    "5-9",
    "10-14",
    "15-19",
    "20-24",
    "25-29",
    "30-34",
    "35-39",
    "40-44",
    "45-49",
    "50-54",
    "55-59",
    "60-64",
    "65-69",
    "70-74",
    "75+",
]


def _normalize_country_name(country: str) -> str:
    """Return a case-insensitive key without whitespace or punctuation."""
    return "".join(ch for ch in country.casefold() if ch.isalnum())


def _download_contact_workbook(
        url: str = CONTACT_ZIP_URL, target_filename: str = CONTACT_WORKBOOK_NAME) -> bytes:
    """Download the ZIP from the DOI and return the workbook bytes in memory."""
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(response.content)) as zf:
        candidates = [name for name in zf.namelist()
                      if name.endswith(target_filename)]
        if not candidates:
            raise FileNotFoundError(
                f"'{target_filename}' not found in downloaded archive.")
        with zf.open(candidates[0]) as f:
            return f.read()


def _load_workbook_bytes(
        contact_path: Optional[str],
        url: str = CONTACT_ZIP_URL,
        target_filename: str = CONTACT_WORKBOOK_NAME) -> bytes:
    """Return workbook bytes either from a user path or by downloading the ZIP."""
    if contact_path:
        if not os.path.exists(contact_path):
            raise FileNotFoundError(
                f"Contact matrix file not found at {contact_path}")
        with open(contact_path, "rb") as f:
            return f.read()
    return _download_contact_workbook(url=url, target_filename=target_filename)


def list_available_contact_countries(
        contact_path: Optional[str] = None) -> List[str]:
    """List all country names available in the contact matrix workbook."""
    xls_bytes = _load_workbook_bytes(contact_path)
    xls = pd.ExcelFile(io.BytesIO(xls_bytes))
    return xls.sheet_names


def _select_sheet_name(country: str, sheet_names: Iterable[str]) -> str:
    lookup = {_normalize_country_name(name): name for name in sheet_names}
    key = _normalize_country_name(country)
    if key not in lookup:
        available = ", ".join(sheet_names)
        raise ValueError(
            f"Country '{country}' not found in contact matrices. "
            f"Available sheets: {available}")
    return lookup[key]


def load_contact_matrix(
        country: str, contact_path: Optional[str] = None) -> pd.DataFrame:
    """
    Load the all-locations contact matrix for the given country. If
    ``contact_path`` is not provided, the function downloads the
    ``MUestimates_all_locations_1.xlsx`` workbook from the DOI ZIP.

    :param country: Country name as listed in the workbook (case-insensitive).
    :param contact_path: Optional path to ``MUestimates_all_locations_1.xlsx``.
    :returns: DataFrame indexed by age group with floats.
    """
    xls_bytes = _load_workbook_bytes(contact_path)
    xls = pd.ExcelFile(io.BytesIO(xls_bytes))
    sheet_names = xls.sheet_names
    sheet = _select_sheet_name(country, sheet_names)
    df = pd.read_excel(xls, sheet_name=sheet, engine="openpyxl")

    # Ensure numeric values and trim potential trailing rows/cols.
    matrix = df.apply(pd.to_numeric, errors="coerce")
    matrix = matrix.iloc[:len(AGE_GROUP_LABELS), :len(AGE_GROUP_LABELS)]
    matrix.columns = AGE_GROUP_LABELS[:matrix.shape[1]]
    matrix.index = AGE_GROUP_LABELS[:matrix.shape[0]]

    if matrix.isnull().any().any():
        raise ValueError(
            f"Contact matrix for '{country}' contains non-numeric entries.")

    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(
            f"Contact matrix for '{country}' is not square: {matrix.shape}")

    return matrix
