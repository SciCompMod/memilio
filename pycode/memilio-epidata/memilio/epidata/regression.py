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
from memilio.epidata import getNPIData as gnd
from memilio.epidata import getCaseData as gcd
from memilio.epidata import getPopulationData as gpd
from memilio.epidata import getVaccinationData as gvd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import defaultDict as dd
from datetime import date
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm


pd.options.mode.copy_on_write = True


def compute_R_eff(out_folder=dd.defaultDict['out_folder']):
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
    df_r_eff = mdfs.extract_subframe_based_on_dates(
        df_r_eff, date(2020, 3, 1), date(2022, 2, 15))

    # get number of days and counties where incidence < 1.0
    # only useful results if we compute R for all counties
    num_low_cases = incidence_counter
    num_total = df_cases.shape[0]-1
    print('Ratio low cases: ', num_low_cases/num_total)

    if True:
        gd.write_dataframe(df_r_eff, directory, "r_eff_county1001", "json")

    return df_r_eff


# TODO: Refactor this function into one where set up the model with all variables that are fixed and will not be evaluated
# in the backward selection and one function where we solve the model with an input of the variables that go into
# the backward selection, i.e. the NPIs? At the moment we are reading the data every time we solve the model, not efficient
# TODO: Discuss which variables should go in backward selection, also vaccination states?
def regression_model(columns, out_folder=dd.defaultDict['out_folder']):

    # for now, all data is only read for county 1001 to make testing easier
    directory = out_folder
    directory = os.path.join(directory, 'Germany/')
    gd.check_dir(directory)

    # read NPI data
    filepath = os.path.join(
        directory, 'germany_counties_npi_maincat.csv')

    if not os.path.exists(filepath):
        df_npis = gnd.get_npi_data(start_date=date(2020, 1, 1),
                                   fine_resolution=0, file_format='csv')
    else:
        df_npis = pd.read_csv(filepath)
    df_npis = df_npis[df_npis.ID_County == 1001]

    # read population data
    filepath = os.path.join(
        directory, 'county_current_population.json')

    if not os.path.exists(filepath):
        df_population = gpd.get_population_data(
            start_date=date(2020, 1, 1), file_format='json')
    else:
        df_population = pd.read_json(filepath)
    df_population = df_population[df_population.ID_County == 1001]

    # read vaccination data
    filepath = os.path.join(
        directory, 'vacc_county_all_dates.json')

    if not os.path.exists(filepath):
        df_vaccinations = gvd.get_vaccination_data(
            start_date=date(2020, 1, 1), file_format='json')
    else:
        df_vaccinations = pd.read_json(filepath)
    df_vaccinations = df_vaccinations[df_vaccinations.ID_County == 1001]

    # TODO: discuss if this is what we want (in contrast to absolute values)
    # computing proportion of vaccinated individuals by dividing by respective population of county
    counties = df_vaccinations['ID_County'].unique()
    vacc_states = df_vaccinations.keys()[3:]
    for county in counties:
        for vacc_state in vacc_states:
            df_vaccinations.loc[(df_vaccinations['ID_County'] == county), vacc_state] = df_vaccinations.loc[df_vaccinations['ID_County']
                                                                                                            == county, vacc_state]/df_population.loc[df_population.ID_County == county, 'Population'].iloc[0]

    # variable for seasonality

    # variables for region types
    # tbd

    # read values for effective reproduction number
    filepath = os.path.join(
        directory, "r_eff_county1001.json")

    if not os.path.exists(filepath):
        df_r = compute_R_eff()
    else:
        df_r = pd.read_json(filepath)

    # make dates consistent, use df_vaccinations as reference

    # TODO: define timeframe we want to investigate
    # TODO: make dataframes consistent wrt timeframe, need to insert 0 entries before vaccination timeframe
    # for now, use df_vaccinations as reference because this has the least entries

    # remove dates from df_r which are not in df_vaccinations and vice versa
    df_r = df_r[df_r['Date'].astype(
        str).isin(df_vaccinations.Date.astype(str))]
    df_vaccinations = df_vaccinations[df_vaccinations['Date'].astype(
        str).isin(df_r.Date.astype(str))]
    # remove dates from df_npis which are not in df_vaccinations
    df_npis = df_npis[df_npis['Date'].astype(
        str).isin(df_vaccinations.Date.astype(str))]

    # set up regression model
    Y = df_r['R_eff']

    # TODO: discuss which vaccination states we want to include
    # for now use Vacc_completed and Vacc_refreshed and use Vacc_partially as reference
    vacc_states = vacc_states[1:3]
    X_vaccinations = np.array([df_vaccinations[vacc_state]
                              for vacc_state in vacc_states]).T
    X_npis = np.array([df_npis[column] for column in columns]).T

    X = np.concatenate((X_vaccinations, X_npis), axis=1)
    # TODO: check why there is not always a constant added automatically, do we have linear dependence somewhere?
    # with has_constant = 'add' we force it to add a constant but can this lead to issues?
    X = sm.add_constant(X, has_constant='add')

    plt.plot(df_r.Date, df_r.R_eff, marker='o')

    # do regression
    model = sm.GLM(Y, X, family=sm.families.Gamma(
        sm.families.links.Log()))

    results = model.fit()

    # TODO: returning vacc_states is not nice, refactor this function?
    return results, list(vacc_states)


def backward_selection(plot=False):

    # initial set of NPIs
    # use fine_resolution=0 for now for simplicity
    # include 'const' as this is an additional variable that we want to evaluate according to pvalue (at least for now)
    column_names = ['M01a', 'M01b', 'M02a', 'M02b',
                    'M03', 'M04', 'M05', 'M06', 'M07', 'M08', 'M09', 'M10', 'M11', 'M12',
                    'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20', 'M21']

    # do regression with all NPIs
    results, vacc_states = regression_model(column_names)
    # store pvalues in dataframe
    df_pvalues = pd.DataFrame({"pvalues": results.pvalues})
    # add column with column names to df
    non_npi_variables = ['const'] + vacc_states
    df_pvalues.insert(1, "columns", non_npi_variables + column_names)
    # drop rows with pvalue that is NaN
    # TODO: check why we get NaNs here in the first place
    df_pvalues.dropna(inplace=True)

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
    while (counter_not_removed < 7) and (len(df_pvalues) > 5):
        iteration += 1

        # choose NPI of interest which is chosen according to the n-th highest pvalue
        # n is determined by the counter_not_removed which is set accordingly if a NPI was removed or not in the previous iteration, see below
        npi_of_interest = df_pvalues.sort_values(
            'pvalues', ascending=False).iloc[counter_not_removed].name

        # if npi_of_interest is in non_npi_variables, take variable with next higher pvalue
        counter_non_npi_var = 0
        while df_pvalues['columns'][npi_of_interest] in non_npi_variables:
            counter_non_npi_var += 1
            # if npi_of_interest in non_npi_variables:
            npi_of_interest = df_pvalues.sort_values(
                'pvalues', ascending=False).iloc[counter_not_removed + counter_non_npi_var].name
            # adjust counter_not_removed in case the new npi_of_interest doesn't get removed
        counter_not_removed += counter_non_npi_var

        # plot_pvalues(df_pvalues, iteration, npi_of_interest)

        # create view of df_pvalues where we remove npi_of_interest and that will be used for regression_model
        df_view = df_pvalues[~df_pvalues.index.isin([npi_of_interest])]
        print("NPI of interest: ", df_pvalues['columns'][npi_of_interest])

        # do new regression and compute AIC and BIC
        # [1:] because we do only want NPIs as input for regression model, not 'const' or vacc_states
        num_non_npi_variables = len(non_npi_variables)
        results, vacc_states = regression_model(
            df_view['columns'][num_non_npi_variables:])
        aic = results.aic
        bic = results.bic_llf
        print('AIC: ', aic)
        print('BIC: ', bic)

        # check if AIC and BIC have decreased compared to before
        if (aic < aic_min) and (bic < bic_min):
            if plot:
                # plot pvalues
                plot_pvalues(df_pvalues, iteration,
                             npi_of_interest, removed=True)

            # set new reference values for AIC and BIC
            aic_min = aic
            bic_min = bic

            # add npi_of_interest to removed_list
            removed_list.append(df_pvalues['columns'][npi_of_interest])
            print('Removed ', df_pvalues['columns'][npi_of_interest])
            # change df_pvalues to df_view because we actually want to remove npi_of_interest
            df_pvalues = df_view[:]

            # set counter_not_removed  = 0 because we want to count how many times in a row a selected NPI was not removed
            # due to unclear AIC and BIC
            # also, in this case we want to select npi_of_interest by taking the NPI with the highest pvalue of remaining NPIs
            counter_not_removed = 0

        else:
            if plot:
                # plot pvalues
                plot_pvalues(df_pvalues, iteration,
                             npi_of_interest, removed=False)

            if aic < aic_min:
                print("BIC didn't decrease, don't remove {}".format(
                    df_pvalues['columns'][npi_of_interest]))
            elif bic < bic_min:
                print("AIC didn't decrease, don't remove {}".format(
                    df_pvalues['columns'][npi_of_interest]))
            else:
                print("AIC and BIC didn't decrease, don't remove {}".format(
                    df_pvalues['columns'][npi_of_interest]))

            # increase counter_not_removed , we select npi_of_interest by taking the NPI with the next highest pvalue
            counter_not_removed += 1

    print(removed_list)

    # do one last regression here to make sure that df_pvalues and results are matching
    # (i.e. also if in last loop no NPI was removed)
    results, vacc_states = regression_model(
        df_pvalues['columns'][num_non_npi_variables:])

    # append coefficients and lower and upper boundary of confidence intervals to df_pvalues
    df_pvalues.insert(2, "coeffs", list(results.params))
    df_pvalues.insert(3, "conf_int_min", list(results.conf_int()[0]))
    df_pvalues.insert(4, "conf_int_max", list(results.conf_int()[1]))

    return df_pvalues, results


def plot_confidence_intervals(df_pvalues):

    # plot coefficients and confidence intervals per NPI
    fig, ax = plt.subplots()
    for i in range(len(df_pvalues)):
        ax.plot((df_pvalues['conf_int_min'][i],
                 df_pvalues['conf_int_max'][i]), (i, i), '-o', color='teal', markersize=3)
        ax.scatter(df_pvalues['coeffs'][i], i, color='teal', marker='x')

    ax.set_yticks(range(0, len(df_pvalues)), list(df_pvalues['columns']))
    ax.invert_yaxis()

    ax.set_xlabel('Values of coefficients')
    ax.set_ylabel('Variables')

    if not os.path.isdir('plots'):
        os.makedirs('plots')
    plt.savefig('plots/regression_results.png', format='png',
                dpi=500)

    plt.close()


def plot_pvalues(df_pvalues, iteration, npi_of_interest, removed):
    # plot pvalues
    fig, ax = plt.subplots()
    ax.barh(range(len(df_pvalues)), df_pvalues['pvalues'])
    # get index of npi_of interest and change color of that bar
    index = df_pvalues.index.get_loc(npi_of_interest)

    # if npi_of_interest was removed change color to green
    if removed:
        ax.get_children()[index].set_color('g')
        labels = ['NPI of interest was removed']
        handles = [plt.Rectangle((0, 0), 1, 1, color='g')]

    # if npi_of_interest was not removed change color to red
    else:
        ax.get_children()[index].set_color('r')
        labels = ['NPI of interest was not removed']
        handles = [plt.Rectangle((0, 0), 1, 1, color='r')]

    plt.legend(handles, labels, loc='lower right')

    ax.set_yticks(range(0, len(df_pvalues)), list(
        df_pvalues['columns']))
    ax.invert_yaxis()

    ax.set_xlabel('P-values')
    ax.set_ylabel('Variables')

    if not os.path.isdir('plots'):
        os.makedirs('plots')
    plt.savefig(f'plots/pvalues_iteration{iteration}.png', format='png',
                dpi=500)

    plt.close()


def main():
    df_pvalues, results = backward_selection(plot=True)
    plot_confidence_intervals(df_pvalues)


if __name__ == "__main__":
    main()
