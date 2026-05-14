#############################################################################
# Copyright (C) 2020-2026 MEmilio
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
Prem et al., 2017 (DOI: https://doi.org/10.1371/journal.pcbi.1005697).
The module can download the supporting ZIP from
https://doi.org/10.1371/journal.pcbi.1005697.s002 (contains the
``MUestimates_all_locations_1.xlsx`` and
``MUestimates_all_locations_2.xlsx`` workbooks) or read a defined local
workbook path. By default, downloads are done in memory and no files are
written.
"""

import io
import os
import zipfile
from typing import Optional
from collections.abc import Iterable

import numpy as np
import pandas as pd
import requests

CONTACT_ZIP_URL = (
    "https://journals.plos.org/ploscompbiol/article/file"
    "?id=10.1371/journal.pcbi.1005697.s002&type=supplementary"
)
CONTACT_WORKBOOK_NAMES = (
    "MUestimates_all_locations_1.xlsx",
    "MUestimates_all_locations_2.xlsx",
)

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

AGE_GROUP_LABELS_RKI = [
    "0-4",
    "5-14",
    "15-34",
    "35-59",
    "60-79",
    "80-99",
]


def _normalize_country_name(country: str):
    """
    Return a case-insensitive key without whitespace or punctuation.

    :param country: The country name to normalize.
    :returns: Normalized country name.
    """
    return "".join(ch for ch in country.casefold() if ch.isalnum())


def _download_contact_workbooks():
    """
    Download the ZIP and return all contact workbooks.

    :returns: List of workbook contents.
    """
    response = requests.get(CONTACT_ZIP_URL, timeout=30)
    response.raise_for_status()

    workbooks = []
    with zipfile.ZipFile(io.BytesIO(response.content)) as zf:
        for target_filename in CONTACT_WORKBOOK_NAMES:
            candidates = [name for name in zf.namelist()
                          if name.endswith(target_filename)]
            if not candidates:
                raise FileNotFoundError(
                    f"'{target_filename}' not found in downloaded workbook.")
            with zf.open(candidates[0]) as f:
                workbooks.append(f.read())

    return workbooks


def _load_workbooks_bytes(
        contact_path: Optional[str]):
    """
    Return one explicit workbook or, by default, all downloaded workbooks.

    :param contact_path: Optional local path to a single workbook.
    :returns: List of workbook contents.
    """
    if contact_path:
        if not os.path.exists(contact_path):
            raise FileNotFoundError(
                f"Contact matrix file not found at {contact_path}")
        with open(contact_path, "rb") as f:
            return [f.read()]
    return _download_contact_workbooks()


def list_available_contact_countries(
        contact_path: Optional[str] = None):
    """
    List all country names available in the contact matrix workbooks.

    :param contact_path: Optional local path to the workbook.
    :returns: List of all country names.
    """
    countries = []
    seen = set()
    for xls_bytes in _load_workbooks_bytes(contact_path):
        xls = pd.ExcelFile(io.BytesIO(xls_bytes))
        for sheet_name in xls.sheet_names:
            key = _normalize_country_name(sheet_name)
            if key not in seen:
                countries.append(sheet_name)
                seen.add(key)
    return countries


def _select_sheet_name(country: str, sheet_names: Iterable[str]):
    """
    Select the appropriate sheet name from the workbook for a given country.

    :param country: Country name as listed in the workbook (case-insensitive).
    :param sheet_names: Iterable of available sheet names.
    :returns: The exact sheet name as found in the workbook.
    """
    lookup = {_normalize_country_name(name): name for name in sheet_names}
    key = _normalize_country_name(country)
    if key not in lookup:
        available = ", ".join(sheet_names)
        raise ValueError(
            f"Country '{country}' not found in contact matrices. "
            f"Available sheets: {available}")
    return lookup[key]


def _read_contact_sheet(xls: pd.ExcelFile, sheet_name: str, country: str):
    """
    Read a contact sheet and extract its numeric 16x16 matrix.

    Some source workbooks contain an explicit header row, others contain only
    the matrix. Reading without headers and selecting the numeric block handles
    both formats.

    :param xls: Opened Excel workbook containing contact matrix sheets.
    :param sheet_name: Name of the sheet to read.
    :param country: Country name used for error messages.
    :returns: DataFrame containing the extracted 16x16 contact matrix.
    """
    df = pd.read_excel(
        xls,
        sheet_name=sheet_name,
        engine="openpyxl",
        header=None)

    return _extract_contact_matrix(df, country)


def _extract_contact_matrix(df: pd.DataFrame, country: str):
    """
    Extract the 16x16 numeric contact matrix from a raw Excel sheet.

    :param df: Raw sheet data read from the workbook without headers.
    :param country: Country name used for error messages.
    :returns: DataFrame with age-group labels as index and columns.
    """
    matrix_size = len(AGE_GROUP_LABELS)
    numeric = df.apply(pd.to_numeric, errors="coerce")
    max_row_start = numeric.shape[0] - matrix_size
    max_col_start = numeric.shape[1] - matrix_size

    for row_start in range(max_row_start, -1, -1):
        for col_start in range(max_col_start, -1, -1):
            matrix = numeric.iloc[
                row_start:row_start + matrix_size,
                col_start:col_start + matrix_size]
            if matrix.shape == (matrix_size, matrix_size):
                if not matrix.isnull().any().any():
                    matrix = matrix.copy()
                    matrix.columns = AGE_GROUP_LABELS
                    matrix.index = AGE_GROUP_LABELS
                    return matrix

    raise ValueError(
        f"Contact matrix for '{country}' does not contain a numeric "
        f"{matrix_size}x{matrix_size} block. Raw shape: {df.shape}")


def load_contact_matrix(
        country: str,
        contact_path: Optional[str] = None,
        reduce_to_rki_groups: bool = True,
        population: Optional[Iterable[float]] = None):
    """
    Load the all-locations contact matrix for the given country. If
    ``contact_path`` is not provided, the function downloads the
    ``MUestimates_all_locations_*.xlsx`` workbooks from Prem et al., 2017.

    :param country: Country name as listed in the workbook (case-insensitive).
    :param contact_path: Optional path to one ``MUestimates_all_locations``
      workbook.
    :param reduce_to_rki_groups: If True, aggregate to the six RKI age groups 
      (0-4, 5-14, 15-34, 35-59, 60-79, 80+ years). Default True.
    :param population: An iterable of 16 float values representing the
      population size for each original age group. Required if
      reduce_to_rki_groups is True.
    :returns: DataFrame indexed by age group with floats.
    """
    all_sheet_names = []
    for xls_bytes in _load_workbooks_bytes(contact_path):
        xls = pd.ExcelFile(io.BytesIO(xls_bytes))
        sheet_names = xls.sheet_names
        all_sheet_names.extend(sheet_names)
        try:
            sheet = _select_sheet_name(country, sheet_names)
        except ValueError:
            continue
        matrix = _read_contact_sheet(xls, sheet, country)
        break
    else:
        _select_sheet_name(country, all_sheet_names)

    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(
            f"Contact matrix for '{country}' is not square: {matrix.shape}")

    if reduce_to_rki_groups:
        if population is None:
            raise ValueError(
                "To reduce to RKI groups, the population distribution "
                "for the 16 original age groups must be provided via the "
                "'population' argument."
            )
        matrix = _aggregate_to_rki_age_groups(matrix, population)

    return matrix


def _aggregate_to_rki_age_groups(
        matrix: pd.DataFrame, population: Iterable[float]):
    """
    Aggregate an age-structured 16x16 contact matrix to the 6-group RKI
    scheme using population-weighted averages.
    Assumes the original columns/rows follow AGE_GROUP_LABELS order.
    Note: The source only provides data up to 70-74 and a 75+ group.
    We map 60-74 to the 60-79 RKI group and 75+ to the 80-99 RKI group.

    :param matrix: The 16x16 contact matrix to aggregate.
    :param population: Population size for each age group, in the same order
      as AGE_GROUP_LABELS.
    :returns: Aggregated 6x6 contact matrix.
    """
    pop_array = np.array(list(population), dtype=float)
    if len(pop_array) != len(AGE_GROUP_LABELS):
        raise ValueError(
            f"Population array must have {len(AGE_GROUP_LABELS)} elements.")

    groups = [
        [0],  # 0-4
        [1, 2],  # 5-14
        [3, 4, 5, 6],  # 15-34
        [7, 8, 9, 10, 11],  # 35-59
        [12, 13, 14],  # 60-79
        [15],  # 80-99 (source has 75+ only)
    ]

    aggregated = pd.DataFrame(
        index=AGE_GROUP_LABELS_RKI, columns=AGE_GROUP_LABELS_RKI, dtype=float)

    mat_values = matrix.values

    for i, row_indices in enumerate(groups):
        for j, col_indices in enumerate(groups):

            # Summing over the columns to get total contacts from group i to
            #  the new group j.
            contacts_summed_over_cols = mat_values[row_indices][
                :, col_indices].sum(
                axis=1)

            # Weight with the population of the original age groups in the new
            #  group j.
            pop_weights = pop_array[row_indices]
            weighted_contacts = contacts_summed_over_cols * pop_weights

            # Calculate the average: Sum of weighted contacts divided by the
            #  new total population per group.
            aggregated.iat[i, j] = weighted_contacts.sum() / pop_weights.sum()

    return aggregated
