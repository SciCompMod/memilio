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
from datetime import date, datetime
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

pd.options.mode.copy_on_write = True

def compute_R_eff(counties, out_folder=dd.defaultDict['out_folder']):
    # TODO: discuss how we want to compute R value

    start_date_npis = date(2020, 3, 1)
    end_date_npis = date(2022, 2, 15)

    directory = out_folder
    directory = os.path.join(directory, 'Germany/')
    gd.check_dir(directory)
    filepath = os.path.join(
        directory, 'cases_all_county.json')

    if not os.path.exists(filepath):
        # use moving average to avoid influence of reporting delay
        # we use (by default) rep_date = False, so we consider cases by reference date
        # TODO: discuss if ma should be used here (we use ref_date and not rep_date)
        arg_dict = {'files': 'all_county'}
        gcd.get_case_data(**arg_dict)

    # get case data per county
    # TODO: discuss if we want to include age groups
    df_cases = pd.read_json(filepath)
    df_cases = df_cases[df_cases.ID_County.isin(
        counties)].reset_index(drop=True)
    df_cases.drop(['Deaths', 'Recovered'], inplace = True, axis=1)
    df_cases = mdfs.impute_and_reduce_df(
                    df_cases, group_by_cols={dd.EngEng['idCounty']: df_cases[dd.EngEng['idCounty']].unique()},
                    mod_cols=["Confirmed"])
    

    filepath = os.path.join(directory, 'county_current_population.json')
    if not os.path.exists(filepath):
        df_pop = gpd.get_population_data()
    else:
        df_pop = pd.read_json(filepath)



    # create df for effective reproduction number
    # use Date and ID_County from case data df, add column for R_eff where initially all values are 0
    df_r_eff = df_cases.iloc[:, :2]
    df_r_eff.insert(len(df_r_eff.columns), 'R_eff', 0)

    
    for county in counties:
        # This is a copy of a subframe
        df_cases_county = df_cases[df_cases[dd.EngEng['idCounty']]==county]
        pop_county = df_pop[df_pop[dd.EngEng['idCounty']]==county]["Population"].values[0]
        # get confirmed infections of last 7 days
        df_cases_county["Incidence"] = df_cases_county["Confirmed"].diff(periods = 7)/pop_county*100000
        r_local =  (df_cases_county["Incidence"]/(df_cases_county["Incidence"].shift(4))).replace([np.nan, np.inf], 0)
        df_r_eff.loc[df_r_eff.ID_County==county, "R_eff"] += r_local
        
        
    # drop all rows where R_eff = 0
    # TODO: to dicuss if this is what we want
    df_r_eff.drop(df_r_eff[df_r_eff['R_eff'] == 0.0].index, inplace=True)
    df_r_eff.reset_index(inplace=True, drop=True)
    df_r_eff = mdfs.extract_subframe_based_on_dates(
        df_r_eff, start_date_npis, end_date_npis)

    # get number of days and counties where incidence < 1.0
    # only useful results if we compute R for all counties
    num_low_cases = len(np.where(df_r_eff.R_eff<1)[0])
    num_total = df_cases.shape[0]-1
    print('Ratio low cases: ', num_low_cases/num_total)

    if True:
        gd.write_dataframe(df_r_eff, directory,
                           "r_eff_county_multiple_c", "json")

    return df_r_eff


def compute_R_eff_old(counties, out_folder=dd.defaultDict['out_folder']):
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
    df_cases = df_cases[df_cases.ID_County.isin(
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

    def __init__(self, counties):
        self.counties = counties

    # read data that is relevant for regression and store in dataframes
    def read_data(self, out_folder=dd.defaultDict['out_folder']):
        # for now, all data is only read for county 1001 to make testing easier
        directory = out_folder
        directory = os.path.join(directory, 'Germany/')
        gd.check_dir(directory)

        # read NPI data
        filepath = os.path.join(
            directory, 'germany_counties_npi_maincat.csv')

        if not os.path.exists(filepath):
            self.df_npis = gnd.get_npi_data(start_date=date(2020, 1, 1),
                                            fine_resolution=0, file_format='csv')
        else:
            self.df_npis = pd.read_csv(filepath)

        self.df_npis = self.df_npis[self.df_npis.ID_County.isin(self.counties)]

        # read population data
        filepath = os.path.join(
            directory, 'county_current_population.json')

        if not os.path.exists(filepath):
            df_population = gpd.get_population_data(
                start_date=date(2020, 1, 1), file_format='json')
        else:
            df_population = pd.read_json(filepath)
        df_population = df_population[df_population.ID_County.isin(
            self.counties)]

        # read vaccination data
        filepath = os.path.join(
            directory, 'vacc_county_all_dates.json')

        if not os.path.exists(filepath):
            self.df_vaccinations = gvd.get_vaccination_data(
                start_date=date(2020, 1, 1), file_format='json')
        else:
            self.df_vaccinations = pd.read_json(filepath)
        self.df_vaccinations = self.df_vaccinations[self.df_vaccinations.ID_County.isin(
            self.counties)]

        # computing proportion of vaccinated individuals by dividing by respective population of county
        # TODO: discuss if this is what we want (in contrast to absolute values)
        # TODO: discuss if we should save this in some file as well
        counties = self.df_vaccinations['ID_County'].unique()
        self.all_vacc_states = self.df_vaccinations.keys()[3:]
        for county in counties:
            for vacc_state in self.all_vacc_states:
                self.df_vaccinations.loc[(self.df_vaccinations['ID_County'] == county), vacc_state] = self.df_vaccinations.loc[self.df_vaccinations['ID_County']
                                                                                                                               == county, vacc_state]/df_population.loc[df_population.ID_County == county, 'Population'].iloc[0]

        # variable for seasonality
        self.df_seasonality = self.df_vaccinations.loc[:, [
            dd.EngEng['idCounty'], dd.EngEng['date']]]
        self.df_seasonality.insert(2, "sin", np.sin(
            2*np.pi*self.df_seasonality['Date'].dt.dayofyear/365))
        self.df_seasonality.insert(3, "cos", np.cos(
            2*np.pi*self.df_seasonality['Date'].dt.dayofyear/365))

        # variable for virus variant
        # tbd

        # variables for age structure
        # TODO: Which year do we want to use for INKAR data?
        filepath = os.path.join(
            directory, 'inkar_2022.xls')

        self.df_inkar = pd.read_excel(filepath)
        self.df_agestructure = self.df_vaccinations.loc[:, [
            dd.EngEng['idCounty'], dd.EngEng['date']]]

        self.age_categories = ['Unter 18', 'Über 65', 'Altersdurchschnitt']
        for county in self.counties:
            self.df_agestructure.loc[self.df_agestructure['ID_County'] == county,
                                     self.age_categories[0]] = self.df_inkar.loc[self.df_inkar['Kennziffer'] == county]['Einwohner unter 6 Jahre'].values[0] + self.df_inkar.loc[self.df_inkar['Kennziffer'] ==
                                                                                                                                                                                 county]['Einwohner von 6 bis unter 18 Jahren'].values[0]
            self.df_agestructure.loc[self.df_agestructure['ID_County'] == county,
                                     self.age_categories[1]] = self.df_inkar.loc[self.df_inkar['Kennziffer'] == county]['Einwohner 65 Jahre und älter'].values[0]
            self.df_agestructure.loc[self.df_agestructure['ID_County'] == county,
                                     self.age_categories[2]] = self.df_inkar.loc[self.df_inkar['Kennziffer'] == county]['Durchschnittsalter der Bevölkerung'].values[0]

        # variables for region types
        # create dataframe with len of df_vaccination to assign region type for every date and countyID
        self.df_regions = self.df_vaccinations.loc[:, [
            dd.EngEng['idCounty'], dd.EngEng['date']]]
        # create column for every region type
        self.region_types = ['Stadtregion - Metropole', 'Stadtregion - Regiopole und Großstadt', 'Stadtregion - Mittelstadt, städtischer Raum',
                             'Stadtregion - Kleinstädtischer, dörflicher Raum', 'Ländliche Region - Zentrale Stadt', 'Ländliche Region - Städtischer Raum',
                             'Ländliche Region - Kleinstädtischer, dörflicher Raum']
        for r_id in range(7):
            self.df_regions[self.region_types[r_id]] = [1 if self.df_regions.iloc[_id, 0]
                                                        in dd.RegioStaR7ToCountyID[r_id + 1] else 0 for _id in range(len(self.df_regions))]

        # read values for effective reproduction number
        filepath = os.path.join(
            directory, "r_eff_county_multiple_c.json")

        if not os.path.exists(filepath):
            self.df_r = compute_R_eff(counties=self.counties)
        else:
            self.df_r = pd.read_json(filepath)

        # make dates consistent, use df_vaccinations as reference

        # TODO: define timeframe we want to investigate
        # TODO: make dataframes consistent wrt timeframe, need to insert 0 entries before vaccination timeframe
        # for now, use df_vaccinations as reference because this has the least entries

        # remove dates from df_r which are not in df_vaccinations and vice versa
        # first remove all dates which are not in df_vaccination or df_npis
        min_date_vacc = min(self.df_vaccinations.Date)
        min_date_npis = min(self.df_npis.Date)
        min_date = max(min_date_npis, min_date_vacc.strftime("%Y-%m-%d"))
        max_date_vacc = max(self.df_vaccinations.Date)
        max_date_npis = max(self.df_npis.Date)
        max_date = min(max_date_npis, max_date_vacc.strftime("%Y-%m-%d"))

        self.df_r = mdfs.extract_subframe_based_on_dates(
            self.df_r, min_date, max_date)
        self.df_vaccinations = mdfs.extract_subframe_based_on_dates(
            self.df_vaccinations, min_date, max_date)
        self.df_npis = mdfs.extract_subframe_based_on_dates(
            self.df_npis, min_date, max_date)
        self.df_regions = mdfs.extract_subframe_based_on_dates(
            self.df_regions, min_date, max_date)
        self.df_seasonality = mdfs.extract_subframe_based_on_dates(
            self.df_seasonality, min_date, max_date)
        self.df_agestructure = mdfs.extract_subframe_based_on_dates(
            self.df_agestructure, min_date, max_date)

        self.df_npis['Date'] = pd.to_datetime(self.df_npis['Date'])
        self.df_vaccinations['Date'] = pd.to_datetime(
            self.df_vaccinations['Date'])
        self.df_regions['Date'] = pd.to_datetime(self.df_regions['Date'])
        self.df_seasonality['Date'] = pd.to_datetime(
            self.df_seasonality['Date'])
        self.df_agestructure['Date'] = pd.to_datetime(
            self.df_agestructure['Date'])

        # remove all dates where the r-value could not be computed
        for county in self.counties:
            df_r_c = self.df_r[self.df_r.ID_County == county]
            # use region dataframe for reference (smallest of other dfs)
            df_region_c = self.df_regions[self.df_regions.ID_County == county]
            missing_dates = set(df_region_c.Date)-set(df_r_c.Date)
            self.df_npis = self.df_npis.drop(self.df_npis[(self.df_npis.Date.isin(
                missing_dates)) & (self.df_npis.ID_County == county)].index)
            self.df_vaccinations = self.df_vaccinations.drop(self.df_vaccinations[(
                self.df_vaccinations.Date.isin(missing_dates)) & (self.df_vaccinations.ID_County == county)].index)
            self.df_regions = self.df_regions.drop(self.df_regions[(self.df_regions.Date.isin(
                missing_dates)) & (self.df_regions.ID_County == county)].index)
            self.df_seasonality = self.df_seasonality.drop(self.df_seasonality[(self.df_seasonality.Date.isin(
                missing_dates)) & (self.df_seasonality.ID_County == county)].index)
            self.df_agestructure = self.df_agestructure.drop(self.df_agestructure[(self.df_agestructure.Date.isin(
                missing_dates)) & (self.df_agestructure.ID_County == county)].index)

        # drop columns where no NPIs are assigned
        num_dropped = 0
        null_idx = np.where(self.df_npis[self.used_npis].sum() == 0)[0]
        for null_item in null_idx:
            dropped = self.used_npis.pop(null_item-num_dropped)
            print('No Data for NPI: ' + str(dropped))
            self.df_npis = self.df_npis.drop(dropped, axis=1)
            num_dropped += 1

        # sort all dataframes first by county and then by date so that all dataframes are aligned
        self.df_r.sort_values(['ID_County', 'Date']).reset_index(inplace=True)
        self.df_npis.sort_values(
            ['ID_County', 'Date']).reset_index(inplace=True)
        self.df_vaccinations.sort_values(
            ['ID_County', 'Date']).reset_index(inplace=True)
        self.df_regions.sort_values(
            ['ID_County', 'Date']).reset_index(inplace=True)
        self.df_seasonality.sort_values(
            ['ID_County', 'Date']).reset_index(inplace=True)
        self.df_agestructure.sort_values(
            ['ID_County', 'Date']).reset_index(inplace=True)

        # reset index so that we can concatenate dataframes without problems
        self.df_r.reset_index(inplace=True)
        self.df_npis.reset_index(inplace=True)
        self.df_vaccinations.reset_index(inplace=True)
        self.df_regions.reset_index(inplace=True)
        self.df_seasonality.reset_index(inplace=True)
        self.df_agestructure.reset_index(inplace=True)

    # define variables that will be used in every regression

    def set_up_model(self):

        self.read_data()

        # set up endogenous variable
        self.Y = self.df_r['R_eff']

        # TODO: discuss which vaccination states we want to include
        self.used_vacc_states = list(self.all_vacc_states[0:3])

        # TODO: discuss which region type we want to use as reference
        # for now use Metropole
        self.reference_region = 'Stadtregion - Metropole'
        self.region_types.remove(self.reference_region)

        # define df with all variables that go into the model
        self.df_allvariables = pd.DataFrame([self.df_vaccinations[vacc_state]
                                             for vacc_state in self.used_vacc_states] +
                                            [self.df_regions[region_type]
                                             for region_type in self.region_types] +
                                            [self.df_seasonality['sin'], self.df_seasonality['cos']] +
                                            [self.df_agestructure[age_category]
                                             for age_category in self.age_categories] +
                                            [self.df_npis[npi] for npi in self.used_npis]).transpose()

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

    # do backward selection to find out which NPIs have significant impact on R value
    def backward_selection(self, plot=False):

        # initial set of NPIs
        # use fine_resolution=0 for now for simplicity
        # TODO: move this to constructor?
        self.used_npis = ['M01a', 'M01b', 'M02a', 'M02b',
                          'M03', 'M04', 'M05', 'M06', 'M07', 'M08', 'M09', 'M10', 'M11', 'M12',
                          'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20', 'M21']

        # set up regression model
        self.set_up_model()

        # define variables that will be used in backward selection
        regression_variables = self.used_vacc_states + self.region_types + \
            ['sin', 'cos'] + \
            self.age_categories + self.used_npis

        # do regression with all NPIs
        results = self.do_regression(regression_variables)
        # store pvalues in dataframe
        self.df_pvalues = pd.DataFrame({"pvalues": results.pvalues})
        # always keep variable 'const' in regression model
        fixed_variables = ['const']
        # add column with column names to df
        self.df_pvalues.insert(
            1, "columns", fixed_variables + regression_variables)
        # this shouldn't be needed anymore because we check before if there is an NPI that is never active
        # drop rows with pvalue that is NaN
        # self.df_pvalues.dropna(inplace=True)

        # compute AIC and BIC as reference for later
        aic_min = results.aic
        bic_min = results.bic_llf
        print('AIC init: ', aic_min)
        print('BIC init: ', bic_min)

        # count how often an NPI was selected by highest p value but not removed due to unclear AIC/BIC in a row
        counter_not_removed = 0

        # list with NPIs that were removed
        removed_list = []
        iteration = 0
        # TODO: think about how to decide when backwards selection is "done"
        while (counter_not_removed < 7) and (len(self.df_pvalues) > 5+len(fixed_variables)):
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
            num_fixed_variables = len(fixed_variables)
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

                # set counter_not_removed  = 0 because we want to count how many times in a row a selected NPI was not removed
                # due to unclear AIC and BIC
                # also, in this case we want to select variable_of_interest by taking the NPI with the highest pvalue of remaining NPIs
                counter_not_removed = 0

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

        # append coefficients and lower and upper boundary of confidence intervals to df_pvalues
        self.df_pvalues.insert(2, "coeffs", list(results.params))
        self.df_pvalues.insert(3, "conf_int_min", list(results.conf_int()[0]))
        self.df_pvalues.insert(4, "conf_int_max", list(results.conf_int()[1]))

        return self.df_pvalues, results

    # plot coefficients and confidence intervals per independent variable
    def plot_confidence_intervals(self):

        fig, ax = plt.subplots()
        for i in range(len(self.df_pvalues)):
            ax.plot((self.df_pvalues['conf_int_min'].iloc[i],
                    self.df_pvalues['conf_int_max'].iloc[i]), (i, i), '-o', color='teal', markersize=3)
            ax.scatter(self.df_pvalues['coeffs'].iloc[i],
                       i, color='teal', marker='x')

        ax.axvline(color='gray')

        ax.set_yticks(range(0, len(self.df_pvalues)),
                      list(self.df_pvalues['columns']), fontsize=5)
        ax.invert_yaxis()

        ax.set_xlabel('Values of coefficients')
        ax.set_ylabel('Variables')

        if not os.path.isdir('plots'):
            os.makedirs('plots')
        plt.tight_layout()
        plt.savefig('plots/regression_results.png', format='png',
                    dpi=500)

        plt.close()

    # plot pvalues and variable_of_interest for iteration in backward selection
    def plot_pvalues(self, iteration, variable_of_interest, removed):
        # plot pvalues
        fig, ax = plt.subplots()
        ax.barh(range(len(self.df_pvalues)), self.df_pvalues['pvalues'])
        # get index of npi_of interest and change color of that bar
        index = self.df_pvalues.index.get_loc(variable_of_interest)

        # if variable_of_interest was removed change color to green
        if removed:
            ax.get_children()[index].set_color('g')
            labels = ['Variable of interest was removed']
            handles = [plt.Rectangle((0, 0), 1, 1, color='g')]

        # if variable_of_interest was not removed change color to red
        else:
            ax.get_children()[index].set_color('r')
            labels = ['Variable of interest was not removed']
            handles = [plt.Rectangle((0, 0), 1, 1, color='r')]

        plt.legend(handles, labels, loc='lower right')

        ax.set_yticks(range(0, len(self.df_pvalues)), list(
            self.df_pvalues['columns']))
        ax.invert_yaxis()

        ax.set_xlabel('P-values')
        ax.set_ylabel('Variables')

        if not os.path.isdir('plots'):
            os.makedirs('plots')
        plt.tight_layout()
        plt.savefig(f'plots/pvalues_iteration{iteration}.png', format='png',
                    dpi=500)

        plt.close()


def main():
    # 4 randomly chosen counties of each regioStar7 type
    counties = geoger.get_county_ids(merge_eisenach=True, merge_berlin=True)

    npi_regression = NPIRegression(counties)

    df_pvalues, results = npi_regression.backward_selection(plot=True)
    npi_regression.plot_confidence_intervals()


if __name__ == "__main__":
    main()
