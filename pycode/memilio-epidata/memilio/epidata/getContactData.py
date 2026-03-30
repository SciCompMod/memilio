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
``MUestimates_all_locations_1.xlsx`` workbook) or read a defined local
workbook path. By default, downloads are done in memory and no
files are written.
"""

import io
import os
import zipfile
from typing import Iterable, Optional

import numpy as np
import pandas as pd
import requests

CONTACT_ZIP_URL = (
    "https://journals.plos.org/ploscompbiol/article/file"
    "?id=10.1371/journal.pcbi.1005697.s002&type=supplementary"
)
CONTACT_WORKBOOK_NAME = "MUestimates_all_locations_1.xlsx"

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


def _download_contact_workbook(
        url: str = CONTACT_ZIP_URL, target_filename: str = CONTACT_WORKBOOK_NAME):
    """
    Download the ZIP from the url and return the workbook.

    :param url: URL to download the ZIP from.
    :param target_filename: Name of the workbook file within the ZIP.
    :returns: Content of the workbook.
    """
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(response.content)) as zf:
        candidates = [name for name in zf.namelist()
                      if name.endswith(target_filename)]
        if not candidates:
            raise FileNotFoundError(
                f"'{target_filename}' not found in downloaded workbook.")
        with zf.open(candidates[0]) as f:
            return f.read()


def _load_workbook_bytes(
        contact_path: Optional[str],
        url: str = CONTACT_ZIP_URL,
        target_filename: str = CONTACT_WORKBOOK_NAME):
    """
    Return workbook either from a user path or by downloading the ZIP.

    :param contact_path: Optional local path to the workbook.
    :param url: Url to download the ZIP from if no path is provided.
    :param target_filename: Name of the workbook file within the ZIP.
    :returns: Content of the workbook.
    """
    if contact_path:
        if not os.path.exists(contact_path):
            raise FileNotFoundError(
                f"Contact matrix file not found at {contact_path}")
        with open(contact_path, "rb") as f:
            return f.read()
    return _download_contact_workbook(url=url, target_filename=target_filename)


def list_available_contact_countries(
        contact_path: Optional[str] = None):
    """
    List all country names available in the contact matrix workbook.

    :param contact_path: Optional local path to the workbook.
    :returns: List of all country names.
    """
    xls_bytes = _load_workbook_bytes(contact_path)
    xls = pd.ExcelFile(io.BytesIO(xls_bytes))
    return xls.sheet_names


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


def load_contact_matrix(
        country: str,
        contact_path: Optional[str] = None,
        reduce_to_rki_groups: bool = True,
        population: Optional[Iterable[float]] = None):
    """
    Load the all-locations contact matrix for the given country. If
    ``contact_path`` is not provided, the function downloads the
    ``MUestimates_all_locations_1.xlsx`` workbook from Prem et al., 2017.

    :param country: Country name as listed in the workbook (case-insensitive).
    :param contact_path: Optional path to ``MUestimates_all_locations_1.xlsx``.
    :param reduce_to_rki_groups: If True, aggregate to the six RKI age groups 
      (0-4, 5-14, 15-34, 35-59, 60-79, 80+ years). Default True.
    :param population: An iterable of 16 float values representing the population 
      size for each original age group. Required if reduce_to_rki_groups is True.
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
    Aggregate an age-structured 16x16 contact matrix to the 6-group RKI scheme using 
    population-weighted averages.
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
