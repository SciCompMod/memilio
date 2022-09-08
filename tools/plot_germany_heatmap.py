#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
#
# Authors: Daniel Abele, Martin J. Kuehn
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
# import matplotlib.pyplot as plt
from gc import get_count
import matplotlib.colors as pclrs
import pandas as pd
import datetime as dt
import numpy as np
import os
import h5py
# in case of a necessary manual installation of GDAL and Fiona on Windows, see
# https://stackoverflow.com/questions/69521550/importerror-the-read-file-function-requires-the-fiona-package-but-it-is-no
import geopandas as gpd

def print_manual_download(filename, url):
    """! Prints message to ask the user to manually download a file.

    @param filename Filename of file needed.
    @param url URL to file.
    """ 
    print(
        'This script needs manual downloading of files. Please download '
         + filename + ' from ' + url + 'and move it the extracted folder'
        'to the current working directory under tools/.')


def merge_eisenach(map_data : gpd.GeoDataFrame):
    """! Merges geometries for Eisenach with Wartburgkreis of Geopandas 
    dataframe.

    @param map_data GeoPandasDataFrame from county shape file
    """ 
    wartburg = map_data.ARS == '16063'
    eisenach = map_data.ARS == '16056'
    map_data.loc[wartburg, 'geometry'] = [map_data[wartburg].geometry.values[0].union(map_data[eisenach].geometry.values[0])]
    # remove Eisenach and return
    return map_data.drop(map_data[eisenach].index.values[0])

def get_county_data(file, date, column, filters=None, output='sum', file_format='json'):
    """ Reads county data from a json or hdf5 file. County data can be time
    series data where a particular date is extracted or single time data which
    is used rightaway. If the file contains multiple features per county, a 
    particular columns need to be specified.

    @param file Path and filename of file to be read in, relative from current
        directory.
    @param data Date to be extracted from data frame or empty if single or 
        no date information is in the input table.
    @param column Column with values that will be plotted.
    @param filters Dictionary with columns and values for filtering rows.
        If a column or all filters are None then all values are taken.
    @param output [Either 'sum' or 'matrix'] If 'sum' is chosen all selected 
        values for one county will be summed up and only 'column' is returned
        for each county. If 'matrix' is chosen, then also the filter columns
        will be returned with the corresponding entries for each selected 
        criterion or value.
    @param file_format File format; either json or h5.
    """
    input_file = os.path.join(os.getcwd(), file[1])
    if file_format == 'json':
        df = pd.read_json(input_file + '.' + file_format)
    elif file_format == 'h5':
        # todo
        df = h5py.File(input_file + '.' + file_format, 'r')
    else:
        print('Error')
        # raise gd.DataError("Download of Vaccination Data failed. File is empty.")


    date_str = date.strftime('%Y-%m-%d')

    if 'Date' in df.columns:
        if df['Date'].nunique() > 1:
            df = df[df['Date'] == date_str]
    dffilter = pd.Series([True for i in range(len(df))], index = df.index)
    if filters != None:
        for col, vals in filters.items():
            if vals == None:
                vals = df[col].unique()
            dffilter = dffilter & (df[col].isin(vals))

    if output == 'sum':
        return df[dffilter].groupby(['ID_County']).agg({column : sum}).reset_index()
    elif output == 'matrix':
        if filters != None:
            return df[dffilter].loc[:,['ID_County'] +  list(filters.keys()) + [column]]
        else:
            return df
    else:
        print('Error')
        # raise gd.DataError("Chose output form.")

    return pd.DataFrame()

def plot(
    files : list,
    date : dt.date,
    column : str = '',
    filters = None,
    relative = True
):


    df = []
    for file in files.items():
        df.append(get_county_data(file, date=date, column=column, filters=filters))
        df[-1].rename(columns={column : 'Count'}, inplace=True)

    # read and filter population data
    if relative:
        pop_data = pd.read_json(os.path.join(os.getcwd(), 'tools/data/county_current_population_dim401.json'))
        if 'Age_RKI' in filters.keys() and filters['Age_RKI'] != None:
            print('Error. Functionality needs to be implemented.')
        # TODO: 
        # Use create_intervals_mapping() (extract first) and 
        # functionality from get_vaccination_data()

        # merge Eisenach
        # TODO: change back to .values addition for more age groups!
        pop_data = pop_data[['ID_County', 'Population']]
        pop_data.loc[pop_data['ID_County']==16063, 'Population'] += pop_data.loc[pop_data.ID_County==16056, 'Population'] 
        pop_data = pop_data[pop_data.ID_County != 16056]
    
    #read and filter map data
    try:
        map_data = gpd.read_file(os.path.join(os.getcwd(), 'tools/shapes/vg2500_01-01.utm32s.shape/vg2500/vg2500_krs.shp'))
    except FileNotFoundError:
        print_manual_download('Georeferenzierung: UTM32s, Format: shape (ZIP, 3 MB)', 'https://gdz.bkg.bund.de/index.php/default/digitale-geodaten/verwaltungsgebiete/verwaltungsgebiete-1-2-500-000-stand-01-01-vg2500.html')
    
    map_data = merge_eisenach(map_data)
    
    map_filter = pd.Series([False for i in range(len(map_data))], index = map_data.index)
    if 'ID_State' in filters.keys() and filters['ID_State'] != None:
        for state in filters['ID_State']:
            map_filter = map_filter | (map_data.ARS.str.startswith('{:02}'.format(state)))
    map_data = map_data[map_filter]

    def get_plot_data(e):
        if relative:
            pop = pop_data.loc[pop_data['ID_County']==int(e.ARS), 'Population'].values[0]
        local_out = []
        for i in range(len(files.items())):
            local_abs = df[i].loc[df[i]['ID_County']==int(e.ARS), 'Count'].sum()
            if not relative:
                local_out.append(local_abs)
            else:
                local_out.append(local_abs / pop)

        return local_out
    map_data[list(files.keys())] = map_data.apply(get_plot_data, axis=1, result_type='expand')

    # calc min and max and map to non-uniform scale
    map_data_columns = map_data[list(files.keys())]
    min_value = min(map_data_columns.min().values)
    max_value = max(map_data_columns.max().values)
    s = 0.7
    
    def normalize_lin(x, src, dest):
        return (x - src[0]) / (src[1] - src[0]) * (dest[1] - dest[0]) + dest[0]
    def normalize_sine(x):
        if (x < s):
            i = [min_value, s]
            j = [-np.pi / 2, 0]
        else:
            i = [s, max_value]
            j = [0, np.pi / 2]
        map_x = normalize_lin(x, i, j)
        sin_map_x = np.sin(map_x)
        y = normalize_lin(sin_map_x, [-1, 1], [0, 1])
        return y
    def normalize_sine_inv(y):
        map_y = normalize_lin(y, [0, 1], [-1, 1])
        asin_map_y = np.arcsin(map_y)
        if (asin_map_y < 0):
            i = [min_value, s]
            j = [-np.pi / 2, 0]
        else:
            i = [s, max_value]
            j = [0, np.pi / 2]
        x = normalize_lin(asin_map_y, j, i)
        return x
    e = 2
    def normalize_skewed_sine(x):
        x = normalize_lin(x, (min_value, max_value), (-np.pi/2, np.pi/2))
        x = np.sin(x)
        def p(x):
            return (-1)**(e-1) * (0.5**e) * (x - 1)**e + 1
        x = p(x)
        return x
    def normalize_skewed_sine_inv(y):
        def p_inv(y):
            return -((-1)**e * (1-y) / (0.5**e))**(1/e) + 1
        y = p_inv(y)
        y = np.arcsin(y)
        y = normalize_lin(y, (-np.pi/2, np.pi/2), (min_value, max_value))
        return y

    color_norm = pclrs.FuncNorm((np.vectorize(normalize_sine), np.vectorize(normalize_sine_inv)), vmin = min_value, vmax = max_value)

    # file name extensions for filters
    part_filename_age = 'allages'
    if not filter_agegroup is None:
        part_filename_age = 'age{0}'.format(filter_agegroup)    
    part_filename_state = 'allstates'
    if not filter_state is None:
        part_filename_state = 'state{0}'.format(filter_state)
    part_filename_vacc = 'vacc{0}'.format(filter_vacc_state)
    part_filename_date = date.strftime('%Y%m%d')
    part_filename = '{}_{}_{}_{}'.format(part_filename_date, part_filename_vacc, part_filename_age, part_filename_state)

    #interactive html files
    def save_interactive(col, filename):
        map_data.explore(col, legend = True, vmin = min_value, vmax = max_value).save(filename) # TODO: norm

    save_interactive('VaccRateRKI', 'vaccrate_{}_rki.html'.format(part_filename))
    save_interactive('VaccRateNew', 'vaccrate_{}_new.html'.format(part_filename))

    from matplotlib.gridspec import GridSpec
    fig = plt.figure(constrained_layout=True)
    gs = GridSpec(3, 3, figure=fig, wspace = 0.1, width_ratios = [0.05, 1, 1], height_ratios = [0.15, 1, 0.15])
    cax = fig.add_subplot(gs[1,0])
    ax1 = fig.add_subplot(gs[:,1])
    ax2 = fig.add_subplot(gs[:,2])
    map_data.plot('VaccRateRKI', ax = ax1, cax = cax, legend = True, vmin = min_value, vmax = max_value, norm = color_norm)
    ax1.set_axis_off()
    map_data.plot('VaccRateNew', ax = ax2, legend = False, vmin = min_value, vmax = max_value, norm = color_norm)
    ax2.set_axis_off()

    plt.savefig('vaccrate_' + part_filename + "_compare.svg")

    # plt.show()

if __name__ == '__main__':

    files_vacc = {'Reported data' : 'tools/data/raw/all_county_agevacc_vacc_all_dates', 'Sanitized data' : 'tools/data/sanitized/all_county_agevacc_vacc_all_dates'}
    age_groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80+']

    for vacc in ['Vacc_partially']:
        for state in [None]: #[None, 9]:
            # plot(files_vacc, dt.date(2021, 11, 18), column = vacc, filters={'ID_State' : [1], 'Age_RKI' : ['05-11','12-17']})
            plot(files_vacc, dt.date(2021, 11, 18), column = vacc, filters={'ID_State' : [1], 'Age_RKI' : None})