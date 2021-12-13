#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Martin J. Kuehn
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
@file customPlot.py
@brief Plots time series data as provided by the scripts.
"""
import matplotlib.pyplot as plt 
import numpy as np
from epidemiology.epidata import getDataIntoPandasDataFrame as gd
import pandas as pd


def plotList(x, y, legend, title='', xlabel='', ylabel='', fig_name=''):
    """! Plots time series data.

    @param x x values to plot.
    @param y List of curve(s) to plot.
    @param legend List of legend(s) for y value curve(s).
    @param title title of the plot.
    @param xlabel labels for x values.
    @param ylabel labels for y values.
    @param fig_name name of the figure to save.
    @return dataframe with imputed dates (and moving average if requested)
    """
    if len(y) > 3:
        fig, ax = plt.subplots(figsize=(10,6))
    else:
        fig, ax = plt.subplots(figsize=(8,6))
    
    for i in range(len(y)):
        ax.plot(x, y[i], label=legend[i])
    if title != '': 
        ax.set_title(title, fontsize=18)
    tick_range = (np.arange(int(len(x) / 30) + 1) * 30) # only plot every 30th day
    tick_range[-1] -= 1 # adapt for last element
    if isinstance(x, pd.Series):
        x = x.values
    ax.set_xticks(x[tick_range])
    ax.set_xticklabels(x[tick_range], rotation=90, fontsize=8)
    if xlabel != '':
        ax.set_xlabel(xlabel, fontsize=12)
    if ylabel != '':
        ax.set_ylabel(ylabel, fontsize=12)
    
    if len(y) > 3:
        ax.legend(fontsize=12, bbox_to_anchor=(1,1), loc="upper left")
    else:
        ax.legend(fontsize=12)
    
    fig.subplots_adjust(bottom=0.2) # adjust to show full xlabel information

    gd.check_dir('figures/')
    plt.savefig('figures/' +  fig_name + '.jpg', bbox_inches='tight')
