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

from memilio.epidata import modifyDataframeSeries as mdfs
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getNPIData as gnd
from memilio.epidata import getCaseData as gcd
from memilio.epidata import getPopulationData as gpd
from memilio.epidata import getVaccinationData as gvd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import defaultDict as dd
from memilio.epidata import getVariantsData as gvsd
from memilio.epidata import progress_indicator
from datetime import date, datetime, timedelta
import os
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statsmodels.api as sm

pd.options.mode.copy_on_write = True
# mpl.use('TKAgg')


def compute_R_eff(iter_ids, state_level, out_folder=dd.defaultDict['out_folder']):
    with progress_indicator.Spinner(message="Computing r-value"):
        dict_entry = dd.EngEng['idState'] if state_level else dd.EngEng['idCounty']
        # TODO: discuss how we want to compute R value

        start_date_npis = date(2020, 3, 1)
        end_date_npis = date(2022, 2, 15)

        directory = out_folder
        directory = os.path.join(directory, 'Germany/')
        gd.check_dir(directory)
        if not state_level:
            filepath = os.path.join(
                directory, 'cases_all_county.json')
        else:
            filepath = os.path.join(
                directory, 'cases_all_state.json')

        if not os.path.exists(filepath):
            # use moving average to avoid influence of reporting delay
            # we use (by default) rep_date = False, so we consider cases by reference date
            # TODO: discuss if ma should be used here (we use ref_date and not rep_date)
            if not state_level:
                arg_dict = {'files': 'all_county'}
            else:
                arg_dict = {'files': 'all_state'}
            gcd.get_case_data(**arg_dict)

        # get case data per county or state
        # TODO: discuss if we want to include age groups
        df_cases = pd.read_json(filepath)
        df_cases = df_cases[df_cases[dict_entry].isin(
            iter_ids)].reset_index(drop=True)
        df_cases.drop(['Deaths', 'Recovered'], inplace=True, axis=1)
        df_cases = mdfs.impute_and_reduce_df(
            df_cases, group_by_cols={
                dict_entry: df_cases[dict_entry].unique()},
            mod_cols=["Confirmed"])

        # create df for effective reproduction number
        # use Date and ID_County from case data df, add column for R_eff where initially all values are 0
        df_r_eff = df_cases.iloc[:, :2]
        df_r_eff.insert(len(df_r_eff.columns), 'R_eff', 0)

        _, axs = plt.subplots(2)
        for idx in iter_ids:
            # This is a copy of a subframe
            df_cases_county = df_cases[df_cases[dict_entry] == idx]
            # get incidence (This is not actually the incidence, population and "per 100000" cancels out)
            df_cases_county["Incidence"] = df_cases_county["Confirmed"].diff(
                periods=7)
            # dont use values where incidence is 0 ((or close to 0)?), to prevent high errors
            # TODO: find out what value seems reasonable here (maybe percentage of county population?)
            df_cases_county.loc[df_cases_county.Incidence <
                                100, "Incidence"] = np.nan
            # R_t = #neue Fälle(t-6,t)/#neue Fälle(t-10,t-4) (See RKI paper)
            df_r_eff.loc[df_r_eff[dict_entry] == idx, "R_eff"] += (df_cases_county["Incidence"]/(
                df_cases_county["Incidence"].shift(4))).replace([np.nan, np.inf], 0)
            # df_r_eff.drop(df_r_eff[df_r_eff['R_eff'] == 0.0].index, inplace=True)
            if True:
                axs[0].plot(df_cases_county["Date"],
                            df_cases_county['Incidence'])

                axs[1].scatter(df_r_eff[df_r_eff[dict_entry] == idx]["Date"], np.log(
                    df_r_eff[df_r_eff[dict_entry] == idx]['R_eff']), marker='.')

        # plt.show()
        # drop all rows where R_eff = 0
        # TODO: to dicuss if this is what we want
        df_r_eff.drop(df_r_eff[df_r_eff['R_eff'] == 0.0].index, inplace=True)
        df_r_eff.reset_index(inplace=True, drop=True)
        df_r_eff = mdfs.extract_subframe_based_on_dates(
            df_r_eff, start_date_npis, end_date_npis)

        # get number of days and counties where incidence < 1.0
        # only useful results if we compute R for all counties
        num_low_cases = len(np.where(df_r_eff.R_eff < 1)[0])
        num_total = df_cases.shape[0]-1
        print('Ratio low cases: ', num_low_cases/num_total)

        if True:
            gd.write_dataframe(df_r_eff, "", filepath, "json")

    return df_r_eff


def compute_R_eff_old_method(counties, out_folder=dd.defaultDict['out_folder']):
    """!!!This function is not used!!!"""
    # TODO: discuss how we want to compute R value

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
    df_cases = df_cases[df_cases[self.dict_entry].isin(
        counties)].reset_index(drop=True)

    # create df for effective reproduction number
    # use Date and ID_County from case data df, add column for R_eff where initially all values are 0
    df_r_eff = df_cases.iloc[:, :2]
    df_r_eff.insert(len(df_r_eff.columns), 'R_eff', 0)

    counties = df_cases['ID_County'].unique()
    # to check number of dates and counties with low cases
    incidence_counter = 0
    for county in counties:  # counties or [counties[0]]

        start_date = df_cases['Date'][0]
        # start on 5th day of available cases so that we can get case data from 4 days before
        date_counter = 5
        # for every date from then on compute incidence and incidence from 4 days ago
        # TODO: can we compute incidence more effectively by shifting date by 1 day in another dataframe and then subtract dataframes from each other?

        # TODO: check if we can use the below code from getNPIData here
        # from NPI data:
        # df_infec_local['Incidence'] = (pd.Series(
        #         df_infec_local_repeat_first_entry).diff(periods=7) /
        #         pop_local * 100000)[7:].values

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
    df_r_eff = mdfs.extract_subframe_based_on_dates(
        df_r_eff, date(2020, 3, 1), date(2022, 2, 15))

    # get number of days and counties where incidence < 1.0
    # only useful results if we compute R for all counties
    num_low_cases = incidence_counter
    num_total = df_cases.shape[0]-1
    print('Ratio low cases: ', num_low_cases/num_total)

    if True:
        gd.write_dataframe(df_r_eff, directory,
                           "r_eff_county_multiple_c", "json")

    return df_r_eff


class NPIRegression():

    def __init__(self, counties, min_date='2020-03-01', max_date='2022-03-01', fine_resolution=0, delay=0, fixed_effects=False, clustering=None, rki=False, state_level=False):
        self.counties = counties
        self.min_date = min_date
        self.max_date = max_date
        self.fine_resolution = fine_resolution
        self.delay = delay
        self.fixed_effects = fixed_effects
        if clustering != None:
            self.clustering = clustering
        else:
            self.clustering = 'clustering_xxx'
        self.rki = rki
        self.state_level = state_level
        self.states = geoger.get_state_ids()
        self.iter_ids = self.states if state_level else self.counties
        self.dict_entry = dd.EngEng['idState'] if state_level else dd.EngEng['idCounty']

    # read data that is relevant for regression and store in dataframes

    def read_data(self, out_folder=dd.defaultDict['out_folder']):
        with progress_indicator.Spinner(message="Read in data"):
            directory = out_folder
            directory = os.path.join(directory, 'Germany/')
            gd.check_dir(directory)

            if self.rki:
                filepath = os.path.join(directory + 'cluster_rki.json')
                self.clustering = 'cluster_rki'

                if not os.path.exists(filepath):
                    print('Clustering analogous to RKI not found.')
                    return

            else:
                filepath = os.path.join(directory + self.clustering + '.json')

            if not os.path.exists(filepath):

                print('Clustered NPI data not found. Using NPI data.')

                if self.fine_resolution == 0:
                    filepath = os.path.join(
                        directory, 'germany_counties_npi_maincat.csv')

                if self.fine_resolution == 2:
                    filepath = os.path.join(
                        directory, 'germany_counties_npi_subcat.csv')

                if not os.path.exists(filepath):
                    print('NPI data not found. Running download script.')
                    self.df_npis = gnd.get_npi_data(start_date=date(2020, 1, 1),
                                                    fine_resolution=self.fine_resolution, file_format='csv')

                self.df_npis = pd.read_csv(filepath)
                self.used_npis = self.df_npis.columns.to_list()
                self.used_npis.remove('Unnamed: 0')
                self.used_npis.remove('Date')
                self.used_npis.remove('ID_County')

            else:
                self.df_npis = pd.read_json(filepath)

                self.used_npis = self.df_npis.columns.to_list()
                self.used_npis.remove('Date')
                self.used_npis.remove('ID_County')

            # read population data
            filepath = os.path.join(
                directory, 'county_current_population.json')

            if not os.path.exists(filepath):
                df_population = gpd.get_population_data(file_format='json')
            else:
                df_population = pd.read_json(filepath)
            if not self.state_level:
                df_population = df_population[df_population[self.dict_entry].isin(
                    self.counties)]
            else:
                county_to_state_map = geoger.get_countyid_to_stateid_map()
                # aggregate npi data to state level
                # npis are county population weighted mean of all npis in state
                df_population[self.dict_entry] = df_population[dd.EngEng['idState']
                                                               if not self.state_level else dd.EngEng['idCounty']].map(county_to_state_map)
                state_pop_dict = pd.Series(df_population.groupby(self.dict_entry).sum(
                )['Population'].values, index=df_population['ID_State'].unique()).to_dict()
                df_population['Ratio per State'] = df_population['Population'] / \
                    (df_population['ID_State'].map(state_pop_dict))
                ratio_dict = pd.Series(
                    df_population['Ratio per State'].values, index=df_population['ID_County']).to_dict()
                self.df_npis.iloc[:, 2:] = self.df_npis.iloc[:, 2:].mul(
                    self.df_npis['ID_County'].astype(float).map(ratio_dict).values, axis=0)
                self.df_npis['ID_State'] = self.df_npis.ID_County.map(
                    county_to_state_map)
                self.df_npis = self.df_npis.groupby(['Date', 'ID_State']).sum(
                ).reset_index(drop=False).drop('ID_County', axis=1)

                df_population.drop(
                    dd.EngEng['idState'] if not self.state_level else dd.EngEng['idCounty'], axis=1, inplace=True)
                df_population = df_population.groupby(
                    self.dict_entry).sum().reset_index(drop=False)

            # read vaccination data
            filepath = os.path.join(
                directory, 'vacc_county_all_dates.json')

            if not os.path.exists(filepath):
                self.df_vaccinations = gvd.get_vaccination_data(
                    start_date=date(2020, 1, 1), file_format='json')
            else:
                self.df_vaccinations = pd.read_json(filepath)
            if not self.state_level:
                self.df_vaccinations = self.df_vaccinations[self.df_vaccinations[self.dict_entry].isin(
                    self.counties)]
            else:
                self.df_vaccinations = self.df_vaccinations.groupby(
                    ['Date', 'ID_State']).sum().reset_index()
                self.df_vaccinations.drop('ID_County', axis=1, inplace=True)

            # read inkar data
            if not self.state_level:
                # TODO: Which year do we want to use for INKAR data?
                filepath = os.path.join(
                    directory, 'inkar_2022.xls')
                # TODO: aggregate on state level
                self.df_inkar = pd.read_excel(filepath)

            # read values for effective reproduction number
            if self.state_level == False:
                filepath = os.path.join(
                    directory, "r_eff_county_multiple_c.json")
            else:
                filepath = os.path.join(
                    directory, "r_eff_state_multiple_c.json")

            if not os.path.exists(filepath):
                self.df_r = compute_R_eff(
                    iter_ids=self.iter_ids, state_level=self.state_level)
            else:
                # TODO: not working
                self.df_r = pd.read_json(filepath=self.states)

            # shift column 'Date' for a certain number of days given by delay
            # a positive delay means that the effect of NPIs occurs after they are being implemented
            # a negative delay means that the effect of NPIs occurs before they are being implemented
            # this is why we shift the 'Date'-column of df_r by delay days
            self.df_r['Date'] -= pd.Timedelta(days=self.delay)

        with progress_indicator.Spinner(message="Preparing vaccination data"):
            # computing proportion of vaccinated individuals by dividing by respective population of county
            # TODO: discuss if this is what we want (in contrast to absolute values)
            # TODO: discuss if we should save this in some file as well
            unique_ids = self.df_vaccinations[self.dict_entry].unique()
            self.all_vacc_states = self.df_vaccinations.keys()[3:]
            for idx in unique_ids:
                for vacc_state in self.all_vacc_states:
                    self.df_vaccinations.loc[(self.df_vaccinations[self.dict_entry] == idx), vacc_state] = self.df_vaccinations.loc[self.df_vaccinations[self.dict_entry]
                                                                                                                                    == idx, vacc_state]/df_population.loc[df_population[self.dict_entry] == idx, 'Population'].iloc[0]
            # add rows for missing dates if self.min_date is smaller than first date in self.df_vaccinations
            min_date_vacc = min(self.df_vaccinations.Date).strftime("%Y-%m-%d")
            if self.min_date < min_date_vacc:
                # list of all df we want to concatenate
                df_list = []
                # create df for each date that needs to be added
                for date in pd.date_range(start=pd.to_datetime(self.min_date), end=pd.to_datetime(min_date_vacc)-pd.Timedelta(days=1)):
                    # get correct structure for subdataframe
                    df_add_date = self.df_vaccinations.loc[self.df_vaccinations['Date']
                                                           == min_date_vacc]
                    # set correct date
                    df_add_date['Date'] = date
                    # set all values vor vaccinations to 0
                    df_add_date.iloc[:, 3:] = 0
                    df_list.append(df_add_date)
                df_list.append(self.df_vaccinations)
                self.df_vaccinations = pd.concat(df_list)

        with progress_indicator.Spinner(message="Preparing seasonality data"):
            # variable for seasonality
            self.df_seasonality = self.df_vaccinations.loc[:, [
                self.dict_entry, dd.EngEng['date']]]
            self.df_seasonality.insert(2, "sin", np.sin(
                2*np.pi*self.df_seasonality['Date'].dt.dayofyear/365))
            self.df_seasonality.insert(3, "cos", np.cos(
                2*np.pi*self.df_seasonality['Date'].dt.dayofyear/365))

        # variable for virus variant
        with progress_indicator.Spinner(message="Preparing variant data"):
            df_var = gvsd.get_variants_data(transform_to_daily=True)
            # use only alpha and delta (because we want to assess effect of alpha and delta wrt to wildtype)
            # TODO: discuss if this is the right approach (if we want to estimate effects of variants at all with our model)
            self.variants = ['B.1.617.2', 'B.1.1.7']
            # self.variants = df_var.iloc[:, 1:].columns.to_list()
            # add column for counties
            df_var[self.dict_entry] = 0
            # initialize dataframe
            self.df_variants = pd.DataFrame(columns=df_var.columns)
            for idx in self.iter_ids:
                # append same data for each county
                df_var[self.dict_entry] = idx
                self.df_variants = pd.concat([self.df_variants, df_var])

            # fig, ax = plt.subplots()
            # ax.plot(self.df_variants.iloc[:, 0],
            #         self.df_variants.iloc[:, 1:-1], label=self.df_variants.columns[1:-1])

            # if not os.path.isdir(f'plots'):
            #     os.makedirs(f'plots')
            # plt.tight_layout()
            # plt.legend()
            # plt.savefig(f'plots/variants.png', format='png',
            #             dpi=500)

            # plt.close()

            # self.df_variants.plot(
            #     x='Date', y=self.df_variants.columns[1:-1], kind='line')
            # plt.show()

        if not self.state_level:
            with progress_indicator.Spinner(message="Preparing age structure data"):
                # variables for age structure
                self.df_agestructure = self.df_vaccinations.loc[:, [
                    dd.EngEng['idCounty'], dd.EngEng['date']]]

                self.age_categories = ['Unter 18',
                                       'Über 65', 'Altersdurchschnitt']
                for county in self.counties:
                    self.df_agestructure.loc[self.df_agestructure['ID_County'] == county,
                                             self.age_categories[0]] = self.df_inkar.loc[self.df_inkar['Kennziffer'] == county]['Einwohner unter 6 Jahre'].values[0] + self.df_inkar.loc[self.df_inkar['Kennziffer'] ==
                                                                                                                                                                                         county]['Einwohner von 6 bis unter 18 Jahren'].values[0]
                    self.df_agestructure.loc[self.df_agestructure['ID_County'] == county,
                                             self.age_categories[1]] = self.df_inkar.loc[self.df_inkar['Kennziffer'] == county]['Einwohner 65 Jahre und älter'].values[0]
                    self.df_agestructure.loc[self.df_agestructure['ID_County'] == county,
                                             self.age_categories[2]] = self.df_inkar.loc[self.df_inkar['Kennziffer'] == county]['Durchschnittsalter der Bevölkerung'].values[0]

        if not self.state_level:
            with progress_indicator.Spinner(message="Preparing region types"):
                # variables for region types
                # create dataframe with len of df_vaccination to assign region type for every date and countyID
                self.df_regions = self.df_vaccinations.loc[:, [
                    self.dict_entry, dd.EngEng['date']]]
                # create column for every region type
                # self.region_types = ['Stadtregion - Metropole', 'Stadtregion - Regiopole und Großstadt', 'Stadtregion - Mittelstadt, städtischer Raum',
                #                      'Stadtregion - Kleinstädtischer, dörflicher Raum', 'Ländliche Region - Zentrale Stadt', 'Ländliche Region - Städtischer Raum',
                #                      'Ländliche Region - Kleinstädtischer, dörflicher Raum']
                self.region_types = ['71', '72', '73', '74', '75', '76', '77']
                # create columns for region types
                for r_id in range(7):
                    self.df_regions[self.region_types[r_id]] = 0
                # get region type of each county
                for county in self.counties:
                    for key, val in dd.RegioStaR7ToCountyID.items():
                        if county in val:
                            break
                    else:
                        print("WARNING: County " + str(county) +
                              " not found in RegioStaR7 dictionary")
                    # add 1 or each county in the respective column
                    self.df_regions.loc[self.df_regions[self.dict_entry] ==
                                        county, self.region_types[key-1]] += 1

        if self.rki:
            filepath = os.path.join(directory, 'holidays.json')
            df_holidays_temp = pd.read_json(filepath)

            dict_states = {1: 'SH', 2: 'HH', 3: 'NI', 4: 'HB', 5: 'NW', 6: 'HE', 7: 'RP', 8: 'BW',
                           9: 'BY', 10: 'SL', 11: 'BE', 12: 'BB', 13: 'MV', 14: 'SN', 15: 'ST', 16: 'TH'}

            self.df_holidays = self.df_vaccinations.loc[:, [
                self.dict_entry, dd.EngEng['date']]]

            # variable for school holidays
            self.df_holidays['schoolholidays'] = 0
            # variable for after the school holidays
            self.df_holidays['after_holiday_effect'] = 0
            # variable for second half of school holidays
            self.df_holidays['second_half_effect'] = 0

            for state in self.states:
                # get lists for start and end dates of holidays for particular state
                start_dates = list(df_holidays_temp.loc[(df_holidays_temp['stateCode'] == dict_states[state]) & (
                    (df_holidays_temp['year'] == 2020) | (df_holidays_temp['year'] == 2021)), 'start'])
                end_dates = list(df_holidays_temp.loc[(df_holidays_temp['stateCode'] == dict_states[state]) & (
                    (df_holidays_temp['year'] == 2020) | (df_holidays_temp['year'] == 2021)), 'end'])
                for holiday in range(len(start_dates)):
                    start_date = start_dates[holiday]
                    end_date = end_dates[holiday]

                    length_of_holidays = datetime.strptime(end_date, '%Y-%m-%d').date(
                    )-datetime.strptime(start_date, '%Y-%m-%d').date()

                    self.df_holidays.loc[((self.df_holidays['ID_State'] == state) & (
                        self.df_holidays['Date'] >= start_date) & (
                        self.df_holidays['Date'] <= end_date)), 'schoolholidays'] = 1

                    # only add an after effect to holidays if holidays are at least 5 days long
                    if length_of_holidays.days >= 5:
                        self.df_holidays.loc[((self.df_holidays['ID_State'] == state) & (
                            self.df_holidays['Date'] > end_date) & (
                            self.df_holidays['Date'] <= np.datetime64(datetime.strptime(end_date, '%Y-%m-%d').date() + timedelta(days=5)))), 'after_holiday_effect'] = 1

                    # only add an effect for the second half of school holidays if they are at least 12 days long
                    print('length of holidays: ', length_of_holidays.days)
                    if length_of_holidays.days >= 12:
                        self.df_holidays.loc[((self.df_holidays['ID_State'] == state) & (
                            self.df_holidays['Date'] >= np.datetime64(datetime.strptime(start_date, '%Y-%m-%d').date() + timedelta(days=int(0.5*length_of_holidays.days)))) & (
                            self.df_holidays['Date'] <= end_date)), 'second_half_effect'] = 1

            # variable for Easter and Christmas
            dates_easter_christmas = [["2020-04-10", "2020-04-13"], ["2020-12-24", "2021-01-01"], [
                "2021-04-02", "2021-04-05"], ["2021-12-24", "2022-01-01"]]

            self.df_holidays['easter_christmas'] = 0

            for holiday in range(len(dates_easter_christmas)):
                start_date = dates_easter_christmas[holiday][0]
                end_date = dates_easter_christmas[holiday][1]
                self.df_holidays.loc[((self.df_holidays['ID_State'] == state) & (
                    self.df_holidays['Date'] >= start_date) & (
                    self.df_holidays['Date'] <= end_date)), 'easter_christmas'] = 1

            self.holiday_variables = [
                'schoolholidays', 'after_holiday_effect', 'second_half_effect', 'easter_christmas']

        with progress_indicator.Spinner(message="Handle datetimes"):
            # TODO: check if this works as expected
            # make dates consistent, use df_vaccinations as reference

            # first remove all dates which are not in df_npis
            min_date_npis = min(self.df_npis.Date)
            max_date_npis = max(self.df_npis.Date)
            if not isinstance(self.df_npis.Date[0], str):
                min_date_npis = datetime.strftime(min_date_npis, '%Y-%m-%d')
                max_date_npis = datetime.strftime(max_date_npis, '%Y-%m-%d')
            self.min_date = max(min_date_npis, self.min_date)
            self.max_date = min(max_date_npis, self.max_date)

            self.df_r = mdfs.extract_subframe_based_on_dates(
                self.df_r, self.min_date, self.max_date)
            self.df_vaccinations = mdfs.extract_subframe_based_on_dates(
                self.df_vaccinations, self.min_date, self.max_date)
            self.df_npis = mdfs.extract_subframe_based_on_dates(
                self.df_npis, self.min_date, self.max_date)
            self.df_seasonality = mdfs.extract_subframe_based_on_dates(
                self.df_seasonality, self.min_date, self.max_date)
            self.df_variants = mdfs.extract_subframe_based_on_dates(
                self.df_variants, self.min_date, self.max_date)
            if (not self.state_level) or (not self.rki):
                self.df_agestructure = mdfs.extract_subframe_based_on_dates(
                    self.df_agestructure, self.min_date, self.max_date)
                self.df_regions = mdfs.extract_subframe_based_on_dates(
                    self.df_regions, self.min_date, self.max_date)
            if self.rki:
                self.df_holidays = mdfs.extract_subframe_based_on_dates(
                    self.df_holidays, self.min_date, self.max_date)

            self.df_npis['Date'] = pd.to_datetime(self.df_npis['Date'])
            self.df_vaccinations['Date'] = pd.to_datetime(
                self.df_vaccinations['Date'])
            self.df_seasonality['Date'] = pd.to_datetime(
                self.df_seasonality['Date'])
            self.df_variants['Date'] = pd.to_datetime(
                self.df_variants['Date'])
            if (not self.state_level) or (not self.rki):
                self.df_agestructure['Date'] = pd.to_datetime(
                    self.df_agestructure['Date'])
                self.df_regions['Date'] = pd.to_datetime(
                    self.df_regions['Date'])
            if self.rki:
                self.df_holidays['Date'] = pd.to_datetime(
                    self.df_holidays['Date'])

            # remove all dates where the r-value could not be computed
            # check if whole counties are removed due to missing R value
            removed_counties = []
            for idx in self.iter_ids:
                df_r_c = self.df_r[self.df_r[self.dict_entry] == idx]
                # use region dataframe for reference (smallest of other dfs)
                df_seasonality_c = self.df_seasonality[self.df_seasonality[self.dict_entry] == idx]
                missing_dates = set(df_seasonality_c.Date)-set(df_r_c.Date)
                if len(df_r_c.Date) == 0:
                    removed_counties.append(idx)
                self.df_npis = self.df_npis.drop(self.df_npis[(self.df_npis.Date.isin(
                    missing_dates)) & (self.df_npis[self.dict_entry] == idx)].index)
                self.df_vaccinations = self.df_vaccinations.drop(self.df_vaccinations[(
                    self.df_vaccinations.Date.isin(missing_dates)) & (self.df_vaccinations[self.dict_entry] == idx)].index)
                self.df_seasonality = self.df_seasonality.drop(self.df_seasonality[(self.df_seasonality.Date.isin(
                    missing_dates)) & (self.df_seasonality[self.dict_entry] == idx)].index)
                self.df_variants = self.df_variants.drop(self.df_variants[(self.df_variants.Date.isin(
                    missing_dates)) & (self.df_variants[self.dict_entry] == idx)].index)
                if (not self.state_level) or (not self.rki):
                    self.df_agestructure = self.df_agestructure.drop(self.df_agestructure[(self.df_agestructure.Date.isin(
                        missing_dates)) & (self.df_agestructure[self.dict_entry] == idx)].index)
                    self.df_regions = self.df_regions.drop(self.df_regions[(self.df_regions.Date.isin(
                        missing_dates)) & (self.df_regions[self.dict_entry] == idx)].index)
                if self.rki:
                    self.df_holidays = self.df_holidays.drop(self.df_holidays[(self.df_holidays.Date.isin(
                        missing_dates)) & (self.df_holidays[self.dict_entry] == idx)].index)

            # drop columns where no NPIs are assigned
            num_dropped = 0
            null_idx = np.where(self.df_npis[self.used_npis].sum() == 0)[0]
            for null_item in null_idx:
                dropped = self.used_npis.pop(null_item-num_dropped)
                print('No Data for NPI: ' + str(dropped))
                self.df_npis = self.df_npis.drop(dropped, axis=1)
                num_dropped += 1

            # sort all dataframes first by county and then by date so that all dataframes are aligned
            self.df_r.sort_values([self.dict_entry, 'Date'])
            self.df_npis.sort_values(
                [self.dict_entry, 'Date'])
            self.df_vaccinations.sort_values(
                [self.dict_entry, 'Date'])
            self.df_seasonality.sort_values(
                [self.dict_entry, 'Date'])
            self.df_variants.sort_values(
                [self.dict_entry, 'Date'])
            if (not self.state_level) or (not self.rki):
                self.df_agestructure.sort_values(
                    [self.dict_entry, 'Date'])
                self.df_regions.sort_values(
                    [self.dict_entry, 'Date'])
            if self.rki:
                self.df_holidays.sort_values(
                    [self.dict_entry, 'Date'])

            # reset index so that we can concatenate dataframes without problems
            self.df_r.reset_index(inplace=True)
            self.df_npis.reset_index(inplace=True)
            self.df_vaccinations.reset_index(inplace=True)
            self.df_seasonality.reset_index(inplace=True)
            self.df_variants.reset_index(inplace=True)
            if (not self.state_level) or (not self.rki):
                self.df_agestructure.reset_index(inplace=True)
                self.df_regions.reset_index(inplace=True)
            if self.rki:
                self.df_holidays.reset_index(inplace=True)

    # define variables that will be used in every regression

    def set_up_model(self, delay=0):

        self.read_data()

        # set up endogenous variable
        self.df_variants.sort_values([self.dict_entry, 'Date'], inplace=True)
        self.df_r.sort_values([self.dict_entry, 'Date'], inplace=True)
        self.Y = self.df_r['R_eff']
        # if fixed effects is True, the R-value is adjusted by setting the effects of the variants and seasonality
        # as below
        if self.fixed_effects or self.rki:
            variants_considered = ['Other', 'B.1.617.2', 'B.1.1.7']
            # consider (known) variant data to r-value, so the effect does not have to be estimated in regression
            variants = self.df_variants.loc[:, variants_considered]/100
            # add all variants not considered to 'Other' column (with an effect of 0% to wildtype)
            variants['Other'] += self.df_variants.loc[:, ~
                                                      self.df_variants.columns.isin(variants_considered)].iloc[:, 2:-1].sum(axis=1)/100
            variants['Other'] *= 1.
            variants['B.1.1.7'] *= 1.3
            variants['B.1.617.2'] *= 1.6
            self.Y /= (variants.sum(axis=1))

            if not self.rki:
                # add seasonality as a multiplicative factor (TODO: find values / other formula)
                # for now use simple cos
                beta0 = 1
                beta1 = 0.5
                self.Y *= (beta0*(1+beta1*np.cos(2*np.pi *
                                                 self.df_r.Date.dt.day_of_year.values/365)))

        # TODO: discuss which vaccination states we want to include
        self.used_vacc_states = list(self.all_vacc_states[0:1])

        # TODO: discuss which region type we want to use as reference
        # for now use Metropole
        # self.reference_region = 'Stadtregion - Metropole'
        if not self.state_level:
            self.reference_region = '71'
            self.region_types.remove(self.reference_region)

        # TODO: references for other variables?

        # define df with all variables that go into the model
        if self.fixed_effects:
            self.df_allvariables = pd.DataFrame([self.df_vaccinations[vacc_state]
                                                for vacc_state in self.used_vacc_states] +
                                                [self.df_regions[region_type]
                                                for region_type in self.region_types] +
                                                # [self.df_seasonality['sin'], self.df_seasonality['cos']] +
                                                [self.df_agestructure[age_category]
                                                 for age_category in self.age_categories] +
                                                [self.df_npis[npi] for npi in self.used_npis]).transpose()  # + [self.df_variants[variant] for variant in self.variants]).transpose()
        elif self.rki or self.state_level:
            self.df_allvariables = pd.DataFrame([self.df_vaccinations[vacc_state]
                                                 for vacc_state in self.used_vacc_states] +
                                                [self.df_seasonality['sin'], self.df_seasonality['cos']] +
                                                [self.df_npis[npi] for npi in self.used_npis] +
                                                [self.df_holidays[holiday] for holiday in self.holiday_variables]).transpose()
        else:
            self.df_allvariables = pd.DataFrame([self.df_vaccinations[vacc_state]
                                                 for vacc_state in self.used_vacc_states] +
                                                [self.df_regions[region_type]
                                                 for region_type in self.region_types] +
                                                [self.df_seasonality['sin'], self.df_seasonality['cos']] +
                                                [self.df_agestructure[age_category]
                                                 for age_category in self.age_categories] +
                                                [self.df_npis[npi] for npi in self.used_npis] +
                                                [self.df_variants[variant] for variant in self.variants]).transpose()

    # define variables for regression according to input and fit model

    def do_regression(self, variables):

        X = np.array([self.df_allvariables[variable]
                     for variable in variables]).T

        # TODO: check why there is not always a constant added automatically, do we have linear dependence somewhere?
        # with has_constant = 'add' we force it to add a constant but can this lead to issues?
        X = sm.add_constant(X, has_constant='add')

        # plt.plot(df_r.Date, df_r.R_eff, marker='o')

        # do regression
        model = sm.GLM(self.Y, X, family=sm.families.Gamma(
            sm.families.links.Log()))

        results = model.fit()

        return results

    def regression_with_all_variables(self, plot=False):
        # set up regression model
        self.set_up_model()

        # define variables that will be used in backward selection
        if self.fixed_effects:
            regression_variables = self.used_vacc_states + \
                self.region_types + self.age_categories + self.used_npis
        elif self.rki:
            regression_variables = self.used_vacc_states + \
                ['sin', 'cos'] + \
                self.used_npis + \
                self.holiday_variables
        else:
            regression_variables = self.used_vacc_states + \
                self.region_types + ['sin', 'cos'] + \
                self.age_categories + self.used_npis + self.variants

        # do regression with all NPIs
        results = self.do_regression(regression_variables)
        # store pvalues in dataframe
        self.df_pvalues = pd.DataFrame({"pvalues": results.pvalues})
        # always keep variable 'const' in regression model
        self.fixed_variables = ['const']
        # add column with column names to df
        self.df_pvalues.insert(
            1, "columns", self.fixed_variables + regression_variables)
        # append coefficients and lower and upper boundary of confidence intervals to df_pvalues
        self.df_pvalues.insert(2, "coeffs", list(results.params))
        self.df_pvalues.insert(
            3, "conf_int_min", list(results.conf_int()[0]))
        self.df_pvalues.insert(
            4, "conf_int_max", list(results.conf_int()[1]))

        # this shouldn't be needed anymore because we check before if there is an NPI that is never active
        # drop rows with pvalue that is NaN
        # TODO: check if there are NaNs from other variables such as variants
        self.df_pvalues.dropna(inplace=True)

        if plot:
            self.plot_confidence_intervals(iteration=0)

        # compute AIC and BIC as reference for later
        aic = results.aic
        bic = results.bic_llf

        return self.df_pvalues, results, aic, bic

    # do backward selection to find out which NPIs have significant impact on R value

    def backward_selection(self, plot=False):

        self.df_pvalues, results, aic_initial, bic_initial = self.regression_with_all_variables()

        aic_min = aic_initial
        bic_min = bic_initial
        print('AIC init: ', aic_min)
        print('BIC init: ', bic_min)

        # counter for iterations in backward selection
        iteration = 0

        # count how often an NPI was selected by highest p value but not removed due to unclear AIC/BIC in a row
        counter_not_removed = 0

        # list with NPIs that were removed
        removed_list = []

        # TODO: think about how to decide when backwards selection is "done"
        while (counter_not_removed < 3) and (len(self.df_pvalues) > 5+len(self.fixed_variables)):
            iteration += 1

            # choose NPI of interest which is chosen according to the n-th highest pvalue
            # n is determined by the counter_not_removed which is set accordingly if a NPI was removed or not in the previous iteration, see below
            variable_of_interest = self.df_pvalues.sort_values(
                'pvalues', ascending=False).iloc[counter_not_removed].name

            # if variable_of_interest is 'const' use variable with next highest pvalue
            if self.df_pvalues['columns'][variable_of_interest] == 'const':
                variable_of_interest = self.df_pvalues.sort_values(
                    'pvalues', ascending=False).iloc[counter_not_removed+1].name

            # create view of df_pvalues where we remove variable_of_interest and that will be used for regression_model
            df_view = self.df_pvalues[~self.df_pvalues.index.isin(
                [variable_of_interest])]
            print("NPI of interest: ",
                  self.df_pvalues['columns'][variable_of_interest])

            # do new regression and compute AIC and BIC
            # [num_fixed_variables:] because we do want to keep 'const' in regression model
            num_fixed_variables = len(self.fixed_variables)
            results = self.do_regression(list(
                df_view['columns'][num_fixed_variables:]))

            # compute AIC and BIC
            aic = results.aic
            bic = results.bic_llf
            print('AIC: ', aic)
            print('BIC: ', bic)

            # check if AIC and BIC have decreased compared to before
            if (aic < aic_min) and (bic < bic_min):
                if plot:
                    # plot pvalues
                    self.plot_pvalues(iteration,
                                      variable_of_interest, removed=True)

                # set new reference values for AIC and BIC
                aic_min = aic
                bic_min = bic

                # add variable_of_interest to removed_list
                removed_list.append(
                    self.df_pvalues['columns'][variable_of_interest])
                print('Removed ',
                      self.df_pvalues['columns'][variable_of_interest])
                # change df_pvalues to df_view because we actually want to remove variable_of_interest
                self.df_pvalues = df_view[:]
                self.df_pvalues["coeffs"] = list(results.params)
                self.df_pvalues["conf_int_min"] = list(results.conf_int()[0])
                self.df_pvalues["conf_int_max"] = list(
                    results.conf_int()[1])

                # set counter_not_removed  = 0 because we want to count how many times in a row a selected NPI was not removed
                # due to unclear AIC and BIC
                # also, in this case we want to select variable_of_interest by taking the NPI with the highest pvalue of remaining NPIs
                counter_not_removed = 0

                self.plot_confidence_intervals(iteration)

            else:
                if plot:
                    # plot pvalues
                    self.plot_pvalues(iteration,
                                      variable_of_interest, removed=False)

                if aic < aic_min:
                    print("BIC didn't decrease, don't remove {}".format(
                        self.df_pvalues['columns'][variable_of_interest]))
                elif bic < bic_min:
                    print("AIC didn't decrease, don't remove {}".format(
                        self.df_pvalues['columns'][variable_of_interest]))
                else:
                    print("AIC and BIC didn't decrease, don't remove {}".format(
                        self.df_pvalues['columns'][variable_of_interest]))

                # increase counter_not_removed , we select variable_of_interest by taking the NPI with the next highest pvalue
                counter_not_removed += 1

        print('Removed NPIs: ', removed_list)

        # do one last regression here to make sure that df_pvalues and results are matching
        # (i.e. also if in last loop no variable was removed)
        results = self.do_regression(
            self.df_pvalues['columns'][num_fixed_variables:])

        self.plot_confidence_intervals(iteration='final')
        self.plot_active_countydays(iteration='final')

        return self.df_pvalues, results, aic_initial, aic_min

    # plot coefficients and confidence intervals per independent variable
    def plot_confidence_intervals(self, iteration):
        num_plots = len(self.df_pvalues)//25
        split_variables = np.array_split(
            range(len(self.df_pvalues)), num_plots)
        for plot_number in range(num_plots):

            # determine begin and end for plotted variables
            plotted_variables_begin = sum([len(x)
                                          for x in split_variables][:plot_number])
            plotted_variables_end = sum([len(x)
                                        for x in split_variables][:plot_number+1])

            fig, ax = plt.subplots()

            # vertical line at x=0
            ax.axvline(x=1, color='gray')

            for i in range(plotted_variables_begin, plotted_variables_end):
                # plot cofidence interval
                ax.plot((np.exp(self.df_pvalues['conf_int_min'].iloc[i]),
                         np.exp(self.df_pvalues['conf_int_max'].iloc[i])), (i, i), '-o', color='teal', markersize=3)
                # plot coefficient
                ax.scatter(np.exp(self.df_pvalues['coeffs'].iloc[i]),
                           i, color='teal', marker='x')

                # different colors for background
                if i % 2 == 0:
                    ax.axhspan(i-0.5, i+0.5, facecolor='gray', alpha=0.1)

            ax.set_yticks(range(plotted_variables_begin, plotted_variables_end),
                          list(self.df_pvalues['columns'])[plotted_variables_begin: plotted_variables_end], fontsize=5)
            ax.invert_yaxis()

            ax.set_xlabel('Effect on R-value')
            ax.set_ylabel('Variables')

            if not os.path.isdir(f'plots/{self.min_date}to{self.max_date}/{self.clustering}/regression_results'):
                os.makedirs(
                    f'plots/{self.min_date}to{self.max_date}/{self.clustering}/regression_results')
            plt.tight_layout()
            plt.savefig(f'plots/{self.min_date}to{self.max_date}/{self.clustering}/regression_results/regression_results_iteration_{iteration}_plot{plot_number}.png', format='png',
                        dpi=500)

            plt.close()

    # plot pvalues and variable_of_interest for iteration in backward selection
    def plot_pvalues(self, iteration, variable_of_interest, removed):
        num_plots = 3
        # plot pvalues
        fig, (ax1, ax2, ax3) = plt.subplots(1, num_plots, sharex=True)
        axes = [ax1, ax2, ax3]

        split_variables = np.array_split(
            range(len(self.df_pvalues)), num_plots)

        for plot_number in range(num_plots):

            # determine begin and end for plotted variables
            plotted_variables_begin = sum([len(x)
                                          for x in split_variables][:plot_number])
            plotted_variables_end = sum([len(x)
                                        for x in split_variables][:plot_number+1])
            axes[plot_number].barh(range(plotted_variables_begin, plotted_variables_end),
                                   self.df_pvalues['pvalues'][plotted_variables_begin: plotted_variables_end], zorder=1)

            # add vertical lines for orientation
            axes[plot_number].axvline(0.25, color='gray', alpha=0.2, zorder=0)
            axes[plot_number].axvline(0.5, color='gray', alpha=0.2, zorder=0)
            axes[plot_number].axvline(0.75, color='gray', alpha=0.2, zorder=0)

            # set y axis
            axes[plot_number].set_yticks(range(plotted_variables_begin, plotted_variables_end), list(
                self.df_pvalues['columns'][plotted_variables_begin: plotted_variables_end]), fontsize=5)
            axes[plot_number].invert_yaxis()

            axes[plot_number].set_xlim(0, self.df_pvalues['pvalues'].max()*1.4)

        # get index of npi_of interest and change color of that bar
        index = self.df_pvalues.index.get_loc(
            variable_of_interest)
        # get plot where variable_of_interest is displayed
        plot_of_interest = [i for i in range(
            num_plots) if index in split_variables[i]][0]
        # adjust index with respect to plot_of_interest
        index = index - sum([len(x)
                             for x in split_variables][:plot_of_interest])

        # if variable_of_interest was removed change color to green
        if removed:
            axes[plot_of_interest].get_children()[index].set_color('g')
            labels = ['Variable of interest was removed']
            handles = [plt.Rectangle((0, 0), 1, 1, color='g')]

        # if variable_of_interest was not removed change color to red
        else:
            axes[plot_of_interest].get_children()[index].set_color('r')
            labels = ['Variable of interest was not removed']
            handles = [plt.Rectangle((0, 0), 1, 1, color='r')]

        axes[1].set_xlabel('P-values')
        axes[0].set_ylabel('Variables')

        plt.legend(handles, labels, bbox_to_anchor=(1, 0), loc="lower right",
                   bbox_transform=fig.transFigure, ncol=1, fontsize=8)

        if not os.path.isdir(f'plots/{self.min_date}to{self.max_date}/{self.clustering}/pvalues'):
            os.makedirs(
                f'plots/{self.min_date}to{self.max_date}/{self.clustering}/pvalues')
        plt.tight_layout()
        plt.savefig(f'plots/{self.min_date}to{self.max_date}/{self.clustering}/pvalues/pvalues_iteration_{iteration}.png', format='png',
                    dpi=500)

        plt.close()

    def plot_active_countydays(self, iteration):
        num_plots = len(self.df_pvalues)//25
        split_variables = np.array_split(
            range(len(self.df_pvalues)), num_plots)
        for plot_number in range(num_plots):

            # determine begin and end for plotted variables
            plotted_variables_begin = sum([len(x)
                                          for x in split_variables][:plot_number])
            plotted_variables_end = sum([len(x)
                                        for x in split_variables][:plot_number+1])

            fig, ax = plt.subplots()

            # vertical line at x=0
            ax.axvline(color='gray')

            npi_names_all = []
            county_days_all = []
            county_days_total = self.df_npis.count().iloc[0]
            for i in range(plotted_variables_begin, plotted_variables_end):
                # get number of active county day for all NPIs/ clusters of NPIs
                npi_name = self.df_pvalues['columns'].iloc[i]

                if npi_name in list(self.df_npis.columns[3:]):
                    npi_names_all.append(npi_name)
                    county_days_all.append(self.df_npis[npi_name].sum())

            for j in range(len(self.df_npis.columns[3:])):

                ax.barh(npi_names_all, county_days_all/county_days_total)

                ax.set_yticks(range(len(npi_names_all)),
                              npi_names_all, fontsize=5)
                ax.invert_yaxis()

                ax.set_xlabel('Percentage of active county days')
                ax.set_ylabel('NPIs')

            if not os.path.isdir(f'plots/{self.min_date}to{self.max_date}/{self.clustering}/regression_results'):
                os.makedirs(
                    f'plots/{self.min_date}to{self.max_date}/{self.clustering}/regression_results')
            plt.tight_layout()
            plt.savefig(f'plots/{self.min_date}to{self.max_date}/{self.clustering}/regression_results/active_countydays_iteration_{iteration}_plot{plot_number}.png', format='png',
                        dpi=500)

            plt.close()


def main():
    counties = geoger.get_county_ids(merge_eisenach=True, merge_berlin=True)

    min_date = '2020-03-01'  # as in rki
    max_date = '2021-08-31'  # as in rki

    fine_resolution = 2

    delay = 0

    fixed_effects = False

    clustering = 'clustered_npis.json'

    rki = True

    if rki:
        npi_regression = NPIRegression(
            counties, min_date, max_date, fine_resolution, delay, fixed_effects, clustering=None, rki=True, state_level=True)

        df_pvalues, results, aic_initial, aic_final = npi_regression.regression_with_all_variables(
            plot=True)

    else:
        npi_regression = NPIRegression(
            counties, min_date, max_date, fine_resolution, delay, fixed_effects, clustering=None, rki=False)

        df_pvalues, results, aic_initial, aic_final = npi_regression.backward_selection(
            plot=True)


if __name__ == "__main__":
    main()
