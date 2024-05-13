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


def chow_test():
    pass


def main():
    counties = geoger.get_county_ids(merge_eisenach=True, merge_berlin=True)

    min_date = '2020-03-01'
    max_date = '2020-05-01'

    fine_resolution = 2

    delay = 0

    fixed_effects = False

    npi_regression = NPIRegression(
        counties, min_date, max_date, fine_resolution, delay, fixed_effects)

    df_pvalues, results, aic_initial, aic_final = npi_regression.backward_selection(
        plot=True)

    investigate_delay = False
    if investigate_delay:
        delay_list = [delay for delay in range(-10, 11)]
        investigate_influence_of_delay(
            counties, min_date, max_date, fine_resolution, delay_list, plot=True)


if __name__ == "__main__":
    main()
