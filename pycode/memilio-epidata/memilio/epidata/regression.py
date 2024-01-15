#############################################################################
# Copyright (C) 2020-2024 MEmilio
#
# Authors: Anna Wendler, PATRICK LENZ !
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

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

from datetime import date

from memilio.epidata import defaultDict as dd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import getCaseData as gcd
from memilio.epidata import getNPIData as gnd
from memilio.epidata import modifyDataframeSeries as mdfs

pd.options.mode.copy_on_write = True


def compute_R_eff(out_folder=dd.defaultDict['out_folder']):

    directory = out_folder
    directory = os.path.join(directory, 'Germany/')
    gd.check_dir(directory)
    filepath = os.path.join(
        directory, 'cases_all_county_ma7.json')

    if not os.path.exists(filepath):
        # use moving average to avoid influence of reporting delay
        # we use (by default) rep_date = False, so we consider cases by reference date
        # TODO: discuss if ma should be used here (we use ref_date and not rep_date)
        arg_dict = {'moving_average': 7, 'files': 'all_county'}
        gcd.get_case_data(**arg_dict)

    # get case data per county
    # TODO: discuss if we want to include age groups

    df_cases = pd.read_json(filepath)
    df_cases = df_cases[df_cases.ID_County == 1001].reset_index(drop=True)

    # create df for effective reproduction number
    # use Date and ID_County from case data df, add column for R_eff where initially all values are 0
    df_r_eff = df_cases.iloc[:, :2]
    df_r_eff.insert(len(df_r_eff.columns), 'R_eff', 0)

    counties = df_cases['ID_County'].unique()
    counties = [1001]
    # to check number of dates and counties with low cases
    incidence_counter = 0
    for county in counties:  # counties or [counties[0]]

        start_date = df_cases['Date'][0]
        # start on 5th day of available cases so that we can get case data from 4 days before
        date_counter = 5
        # for every date from then on compute incidence and incidence from 4 days ago
        # TODO: can we compute incidence more effectively by shifting date by 1 day in another dataframe and then subtract dataframes from each other?
        while date_counter < df_cases.loc[df_cases['ID_County'] == county].shape[0]:
            incidence_today = df_cases.loc[(df_cases['ID_County'] == county) & (df_cases['Date'] == start_date +
                                                                                pd.Timedelta(date_counter, 'd')), 'Confirmed'].item() - df_cases.loc[(df_cases['ID_County'] == county) & (df_cases['Date'] == start_date +
                                                                                                                                                                                          pd.Timedelta(date_counter-1, 'd')), 'Confirmed'].item()

            incidence_4_days_ago = df_cases.loc[(df_cases['ID_County'] == county) & (df_cases['Date'] == start_date +
                                                                                     pd.Timedelta(date_counter-4, 'd')), 'Confirmed'].item() - df_cases.loc[(df_cases['ID_County'] == county) & (df_cases['Date'] == start_date +
                                                                                                                                                                                                 pd.Timedelta(date_counter-5, 'd')), 'Confirmed'].item()

            # check if incidence_today and incidence_4_days_ago are >= 1, only then compute R_eff
            # TODO: discuss if we should also check values in between
            # TODO: 1 seems arbitrary to me, maybe consider population in respective county
            if (incidence_today >= 1 and incidence_4_days_ago >= 1):
                # compute R_eff and store in df_r_eff
                df_r_eff.loc[(df_r_eff['ID_County'] == county) & (df_r_eff['Date'] == start_date +
                                                                  pd.Timedelta(date_counter, 'd')), 'R_eff'] = incidence_today / incidence_4_days_ago

            if incidence_today < 1.0:
                incidence_counter += 1

            date_counter += 1

    # drop all rows where R_eff = 0
    # TODO: to dicuss if this is what we want
    df_r_eff.drop(df_r_eff[df_r_eff['R_eff'] == 0.0].index, inplace=True)
    df_r_eff.reset_index(inplace=True, drop=True)
    df_r_eff = mdfs.extract_subframe_based_on_dates(df_r_eff, date(2020,3,1), date(2022,2,15))

    # get number of days and counties where incidence < 1.0
    # only useful results if we compute R for all counties
    num_low_cases = incidence_counter
    num_total = df_cases.shape[0]-1
    print('Ratio low cases: ', num_low_cases/num_total)

    if True:
        gd.write_dataframe(df_r_eff, directory, "r_eff_county1001", "json")

    return df_r_eff


def regression_model(out_folder=dd.defaultDict['out_folder']):
    
    directory = out_folder
    directory = os.path.join(directory, 'Germany/')
    gd.check_dir(directory)
    filepath = os.path.join(
        directory, 'germany_counties_npi_maincat.csv')
    
    if not os.path.exists(filepath):
        df_npis = gnd.get_npi_data(start_date=date(2020, 1, 1),
                                   fine_resolution=0, file_format='csv')
    else:
        df_npis = pd.read_csv(filepath)
    df_npis = df_npis[df_npis.ID_County == 1001]
        
    filepath = os.path.join(
        directory, "r_eff_county1001.json")
    
    if not os.path.exists(filepath):
        df_r = compute_R_eff()
    else:
        df_r = pd.read_json(filepath)

    # remove dates from df_npis which are not in df_r
    df_npis = df_npis[df_npis['Date'].astype(str).isin(df_r.Date.astype(str))]
    
    Y = df_r['R_eff']
    columns = ['M01a', 'M01b', 'M02a', 'M02b',
       'M03', 'M04', 'M05', 'M06', 'M07', 'M08', 'M09', 'M10', 'M11', 'M12',
       'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20', 'M21']
    X = np.array([df_npis[column] for column in columns]).T
    X = sm.add_constant(X)
    plt.plot(df_r.Date, df_r.R_eff, marker='o')

    model = sm.GLM(Y, X, family=sm.families.Gamma(
        sm.families.links.Log()))
    
    results = model.fit()
    print(results.summary())
    fig = sm.graphics.plot_partregress_grid(results)
    plt.show()


    # run simple regression model where log(R) depends on NPIs only
    # use fine_resolution=0 in NPI df for simplicity for now

    return


def main():
    regression_model()


if __name__ == "__main__":
    main()
