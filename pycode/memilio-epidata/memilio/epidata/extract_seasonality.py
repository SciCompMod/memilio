import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import os
from memilio.epidata import defaultDict as dd
from memilio.epidata import regression as r
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import getPopulationData as gpd
from statsmodels.tsa.seasonal import seasonal_decompose
from statsmodels.tsa.seasonal import MSTL, STL
from scipy import optimize
from datetime import datetime
pd.options.mode.copy_on_write = True
mpl.use('QtAgg')


def naive_model(df_r):
    df_r.set_index('Date', inplace=True)
    data = df_r[['R_eff']].copy()
    decompose_result_mult = seasonal_decompose(
        data, model="multiplicative", period=365)
    decompose_result_mult.plot()
    plt.savefig('plots/seasonality_decomposition.png', dpi=500)
    x = df_r.index
    y = decompose_result_mult.seasonal.values
    plot_seasonality(df_r, x, y)
    return dict(zip(x, y))


def plot_seasonality(df_r, x, y, ):
    p, cov = optimize.curve_fit(seasonality_func, np.array(
        df_r.index.dayofyear.to_list()), y)
    plt.figure(figsize=(6, 4))
    plt.scatter(x, y, label='Data', c='b', alpha=0.5)
    plt.plot(x, seasonality_func(np.array(df_r.index.dayofyear.to_list()),
             p[0], p[1], p[2], p[3], p[4]), label='Fitted function', c='r')
    plt.axhline(y.mean(), c='b')
    print(p)
    y_new = y-seasonality_func(np.array(df_r.index.dayofyear.to_list()),
                               p[0], p[1], p[2], p[3], p[4])+1
    plt.scatter(x, y_new, c='y', alpha=0.3, label='Resolved data')
    plt.fill_between(x, y.mean()+y.std(), y.mean()-y.std(),
                     color='aqua', alpha=0.3, label='Scattering with seasonality')
    plt.axhline(y_new.mean(), c='g')
    plt.fill_between(x, y_new.mean()+y_new.std(), y_new.mean()-y_new.std(),
                     color='lime', alpha=0.2, label='Scattering without Seasonality')
    plt.legend(loc='best')

    plt.savefig('plots/seasonality.png', dpi=500)

    plt.show()


def LOESS_model(df_r):
    # locally estimated scatterplot smoothing
    data = df_r[['R_eff']].copy()
    stl = STL(data, period=365)
    res = stl.fit()
    fig = res.plot()
    x = df_r.index
    y = res.seasonal.values
    plot_seasonality(df_r, x, y)


def seasonality_func(x, a, b, c, d, e):
    return a*np.sin(2*np.pi*x/365-b)+c*np.cos(2*np.pi*x/365-d)+e


def adjust_r(df_r, filename):
    df_pop = gpd.get_population_data()
    pop_total = df_pop.Population.sum()
    if 'state' in filename:
        # get population per State
        df_pop['ID_State'] = df_pop['ID_County'].map(
            geoger.get_countyid_to_stateid_map())
        df_pop = df_pop[['ID_State', 'Population']].groupby('ID_State').sum()
        df_pop.reset_index(drop=False, inplace=True)
        ratio_by_state = dict(
            zip(df_pop.ID_State, df_pop.Population/pop_total))
        # take mean
        df_r['R_eff'] *= df_r['ID_State'].map(ratio_by_state)
        df_r.drop('ID_State', axis=1, inplace=True)
    elif 'county' in filename:
        # To be tested
        # TODO: population weighted mean is not good. Gives different results on state and county level.
        # High values in the devision can cause Problems even in small counties. Maybe find another method.
        # State level could be a good approx. since the errors on our calculation are smaller.
        ratio_by_county = dict(zip(df_pop['ID_County'], df_pop['Population']))
        df_r['R_eff'] *= df_r['ID_County'].map(ratio_by_county)
        df_r.drop('ID_County', axis=1, inplace=True)
    df_r = df_r.groupby('Date').sum()
    df_r.reset_index(drop=False, inplace=True)
    return df_r


if __name__ == '__main__':
    # mpl.use('TKagg')
    out_folder = dd.defaultDict['out_folder']
    directory = os.path.join(out_folder, 'Germany/')
    filename = "r_eff_county_multiple_c"
    filepath = os.path.join(directory, filename + ".json")
    if False:
        # if os.path.exists(filepath):
        df_r = pd.read_json(filepath)
    else:
        df_r = r.compute_R_eff(
            iter_ids=geoger.get_county_ids(), state_level=False)
    r_orig = df_r[:]
    # adjust r_value so it can be used later. Take (population weighted) mean, remove county column etc.
    df_r = adjust_r(df_r, filename)
    data = naive_model(df_r)
    plt.figure()
    plt.scatter(r_orig.Date, np.log(r_orig.R_eff), c='maroon',
                alpha=0.3, marker='.', label='R_value with seasonality')
    r_orig['R_eff'] /= r_orig['Date'].map(data)
    plt.scatter(r_orig.Date, np.log(r_orig.R_eff), c='dodgerblue',
                alpha=0.3, marker='.', label='R_value without seasonality')
    plt.legend()
    plt.grid(True)
    gd.write_dataframe(r_orig, directory, filename+'_season', "json")
