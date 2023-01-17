#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import numpy as np
import pandas as pd

from memilio.epidata import getDataIntoPandasDataFrame as gd


def plotList(
        x, y, legend, title='', xlabel='', ylabel='', linewidth=1,
        xticks_idx='default', loc_legend='upper_left', fig_size=(10,6), fig_name='customPlot', dpi=300, outercolor='white', innercolor='white'):
    """! Plots time series data.

    @param x x values to plot.
    @param y List of curve(s) to plot.
    @param legend List of legend(s) for y value curve(s).
    @param title title of the plot.
    @param xlabel labels for x values.
    @param ylabel labels for y values.
    @param linewidth [Default: 1] Width of plotted curves.
    @param xticks_idx [Default: 'default'] List of xticks indices to be used. 
        If 'default' is used, then every 30th x-value is printed.
    @param dpi Dots per inch of figure.
    @param outercolor Outer color of figure.
    @oaram innercolor Inner color of figure.
    @param fig_name name of the figure to save.
    """
    fig, ax = plt.subplots(figsize=fig_size, facecolor=outercolor)
    
    for i in range(len(y)):
        ax.plot(x, y[i], linewidth=linewidth, label=legend[i])
    if title != '': 
        ax.set_title(title, fontsize=18)

    ax.set_facecolor(innercolor)

    if isinstance(x, pd.Series):
        x = x.values        
    if xticks_idx == 'default':
        tick_range = (np.arange(int(len(x) / 30) + 1) * 30) # only plot every 30th day
        tick_range[-1] -= 1 # adapt for last element
        ax.set_xticks(x[tick_range])
        ax.set_xticklabels(x[tick_range], rotation=90, fontsize=8)
    else:
        ax.set_xticks([x[i] for i in xticks_idx])
        ax.set_xticklabels([x[i] for i in xticks_idx], rotation=90, fontsize=10)

    if xlabel != '':
        ax.set_xlabel(xlabel, fontsize=12)
    if ylabel != '':
        ax.set_ylabel(ylabel, fontsize=12)
    
    if len(y) > 3:
        ax.legend(fontsize=12, bbox_to_anchor=(1,1), loc=loc_legend)
    else:
        ax.legend(fontsize=12, loc=loc_legend)
    
    fig.subplots_adjust(bottom=0.2) # adjust to show full xlabel information

    gd.check_dir('figures/')
    plt.savefig('figures/' +  fig_name + '.jpg', bbox_inches='tight', dpi=dpi)
