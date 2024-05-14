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
from regression import NPIRegression
from datetime import date, datetime
import os
import math
import scipy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statsmodels.api as sm

pd.options.mode.copy_on_write = True


def investigate_influence_of_delay(counties, min_date, max_date, fine_resolution, delay_list, plot=True):

    aic_initial_all = []
    aic_final_all = []

    for delay in delay_list:

        npi_regression = NPIRegression(
            counties, min_date, max_date, fine_resolution, delay)

        df_pvalues, results, aic_initial, aic_final = npi_regression.backward_selection(
            plot=False)

        aic_initial_all.append(aic_initial)
        aic_final_all.append(aic_final)

    plot_aic_for_delay(delay_list, aic_initial_all, aic_final_all)

    if plot:
        plot_aic_for_delay(delay_list, aic_initial_all, aic_final_all)


def plot_aic_for_delay(delay_list, aic_initial_all, aic_final_all):

    fig, ax = plt.subplots()

    ax.plot(delay_list, aic_initial_all, label='Initial model')
    ax.plot(delay_list, aic_final_all, label='Final model')

    ax.set_xticks(delay_list)

    ax.set_xlabel('Delay')
    ax.set_ylabel('AIC')

    plt.legend()

    if not os.path.isdir(f'plots/delay'):
        os.makedirs(
            f'plots/delay')
    plt.tight_layout()
    plt.savefig(f'plots/delay/aic.png', format='png',
                dpi=500)

    plt.close()


def compute_chow_statistic(results_ur1, results_ur2, results_r):
    # TODO: is pseudo_rsquared a good choice?
    pseudo_rsquared_ur1 = results_ur1.pseudo_rsquared()
    pseudo_rsquared_ur2 = results_ur2.pseudo_rsquared()
    pseudo_rsquared_r = results_r.pseudo_rsquared()

    df_ur1 = results_ur1.df_model
    df_ur2 = results_ur2.df_model
    df_r = results_r.df_model

    chow = ((pseudo_rsquared_r-pseudo_rsquared_ur1-pseudo_rsquared_ur2)/(df_r -
            df_ur1-df_ur2)) / ((pseudo_rsquared_ur1+pseudo_rsquared_ur2)/(df_ur1+df_ur2))

    pvalue = scipy.stats.f.cdf(chow, df_r-df_ur1-df_ur2, df_ur1+df_ur2)


def apply_chow_test(counties, min_date, max_date, variant_change_date):
    # model containing only first variant
    unrestricted_model_1 = NPIRegression(
        counties, min_date, variant_change_date, fine_resolution=2, delay=0, fixed_effects=False)
    df_pvalues_ur1, results_ur1, aic_initial_ur1, aic_final_ur1 = unrestricted_model_1.regression_with_all_variables(
        plot=True)

    # model containing only second variant
    unrestricted_model_2 = NPIRegression(
        counties, variant_change_date, max_date, fine_resolution=2, delay=0, fixed_effects=False)
    df_pvalues_ur2, results_ur2, aic_initial_ur2, aic_final_ur2 = unrestricted_model_2.regression_with_all_variables(
        plot=True)

    # model containing both variants
    restricted_model = NPIRegression(
        counties, min_date, max_date, fine_resolution=2, delay=0, fixed_effects=False)
    df_pvalues_r, results_r, aic_initial_r, aic_final_r = restricted_model.regression_with_all_variables(
        plot=True)

    pvalue = compute_chow_statistic(results_ur1, results_ur2, results_r)

    return pvalue


def main():
    counties = geoger.get_county_ids(merge_eisenach=True, merge_berlin=True)

    fine_resolution = 2

    investigate_delay = False
    if investigate_delay:
        min_date = '2020-03-01'
        max_date = '2020-05-01'

        delay_list = [delay for delay in range(-10, 11)]
        investigate_influence_of_delay(
            counties, min_date, max_date, fine_resolution, delay_list, plot=True)

    do_chow_test = True
    if do_chow_test:
        # TODO: find out exact dates, these are rough estimates from plot
        min_date = '2020-03-01'
        max_date = '2021-07-01'
        change_wildtype_alpha = '2021-02-01'
        pvalue = apply_chow_test(
            counties, min_date, max_date, change_wildtype_alpha)

        print(f'P value for Chow test is {pvalue}.')
        print(0)


if __name__ == "__main__":
    main()
