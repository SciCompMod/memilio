#############################################################################
# Copyright (C) 2020-2024 MEmilio
#
# Authors: Patrick Lenz
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
@file getVariantsData.py

@brief Data on SARS-CoV-2 variants is downloaded.

"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from memilio.epidata import defaultDict as dd
from memilio.epidata import getDataIntoPandasDataFrame as gd

# activate CoW for more predictable behaviour of pandas DataFrames
pd.options.mode.copy_on_write = True
# mpl.use('TKAgg')


def get_variants_data(read_data=dd.defaultDict['read_data'],
                      file_format=dd.defaultDict['file_format'],
                      out_folder=dd.defaultDict['out_folder'],
                      no_raw=dd.defaultDict['no_raw'],
                      make_plot=False,
                      transform_to_daily=True
                      ):

    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)

    filename = "FullData_variants"
    url = "https://opendata.ecdc.europa.eu/covid19/virusvariant/csv/data.csv"
    path = os.path.join(directory + filename + ".json")
    df_raw = gd.get_file(path, url, read_data, param_dict={}, interactive=True)

    if not df_raw.empty:
        if not no_raw:
            gd.write_dataframe(df_raw, directory, filename, file_format)
    else:
        raise gd.DataError("Something went wrong, dataframe is empty.")

    df_ger = df_raw[df_raw.country == "Germany"]
    # for now, only consider GISAID source
    df_ger = df_ger[df_ger.source == "GISAID"]
    variants = df_ger.variant.unique().tolist()
    # remove "UNK" from variants
    variants.remove("UNK")
    print(str(len(variants)) + " Variants")
    df_ger = df_ger[df_ger.percent_variant.notna()]

    all_dates = np.array([], dtype='datetime64')
    weekdays = [[] for _ in range(7)]
    for weekday in range(1, 8):
        all_dates = np.concatenate((all_dates, pd.to_datetime(
            df_ger.year_week + '-' + str(weekday), format="%G-%V-%u").unique()))
        weekdays[weekday-1] = pd.to_datetime(
            df_ger.year_week + '-' + str(weekday), format="%G-%V-%u").unique()

    df_out = pd.DataFrame()
    df_out['Date'] = all_dates

    for v in variants:
        # create column for eacht variant
        df_out[v] = 0.
        # read percentage for each variant and apply those values to df out
        df_sub = df_ger[df_ger.variant == v]
        df_sub_values = df_sub.percent_variant.values[:]
        for weekday in range(1, 8):
            dates = pd.to_datetime(df_sub.year_week.unique(
            ) + '-' + str(weekday), format="%G-%V-%u").unique()
            df_out.loc[df_out.Date.isin(dates), v] += df_sub_values

    if transform_to_daily:
        # transform from weekly data to daily data
        df_out = df_out.set_index('Date').resample('D').ffill()
        df_out.reset_index(inplace=True)

    if make_plot:
        fig, axs = plt.subplots(4, 6)
        for v in range(len(variants)):
            axs[int((v)/6)][(v) % 6].fill_between(df_out.Date,
                                                  df_out[variants[v]], 0., label=variants[v])
            axs[int((v)/6)][(v) % 6].legend()

        plt.savefig('variants_test.png')
        plt.show()

    return df_out


def sanitize_data(df_variants):
    pass


def plot_variants_data(df_variants, variants='wildtype_alpha_delta'):

    fig, ax = plt.subplots()

    if variants == 'all':
        ax.plot(df_variants.iloc[:, 0], df_variants.iloc[:,
                1:-1], label=df_variants.columns[1:-1])
        filename = 'all_variants'

    elif variants == 'wildtype_alpha_delta':
        variants_of_interest = ['Other', 'B.1.617.2', 'B.1.1.7']
        ax.plot(df_variants.iloc[:, 0],
                df_variants.loc[:, variants_of_interest], label=variants_of_interest)
        filename = 'wildtype_alpha_delta'

    if not os.path.isdir(f'plots/variants'):
        os.makedirs(f'plots/variants')
    plt.tight_layout()
    plt.legend()
    plt.savefig(f'plots/variants/{filename}.png', format='png',
                dpi=500)

    # df_variants.plot(
    #     x='Date', y=df_variants.columns[1:-1], kind='line')
    # plt.show()
    plt.close()


def main():
    """ Main program entry."""

    df_variants = get_variants_data(make_plot=False, transform_to_daily=True)
    plot_variants_data(df_variants, variants='wildtype_alpha_delta')


if __name__ == "__main__":
    main()
