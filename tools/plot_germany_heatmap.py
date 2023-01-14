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
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getDataIntoPandasDataFrame as gd
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
    @return Frame with Wartburgkreis and Eisenach merged. 
    """
    wartburg = map_data.ARS == '16063'
    eisenach = map_data.ARS == '16056'
    map_data.loc[wartburg, 'geometry'] = [
        map_data[wartburg].geometry.values[0].union(
            map_data[eisenach].geometry.values[0])]
    # remove Eisenach and return
    return map_data.drop(map_data[eisenach].index.values[0])


def read_data(
        file, region_spec, column, date, filters=None, output='sum',
        file_format='json'):
    """ Reads data from a general json or specific hdf5 file as output by the 
    MEmilio simulation framework and extracts parts of the data if filters are
     applied. The input regional data can be a time series data where a 
    particular date is extracted or single time data which is used rightaway.
    If the file contains multiple features per region,
    particular columns need to be specified.

    @param file Path and filename of file to be read in, relative from current
        directory.
    @param region_spec Specificier for region, for json, e.g., column name
        of county IDs for hdf5 level of nesting X where to retrieve
        fileX.keys(), X could either by empty, such that file.keys() is taken
        or a list of nested keys [X0, X1, ...] such that 
        file[X0][X1][...].keys()
    @param column Column with values that will be plotted.
    @param data Date to be extracted from data frame if it contains a time 
        series or if single or no date information is in the input table, 
        provide date=None. For json, pd.date Objects need to be used, for
        h5, provide an integer from the beginning of the time series.
    @param filters Dictionary with columns and values for filtering rows.
        None for a column or all filters means that all values are taken. 
        For MEmilio h5 Files only 'Group' (AgeGroups from 'Group1' to 'Total') 
        and 'InfectionState' (from 0 to InfetionState::Count-1) can be filtered.
    @param output [Either 'sum' or 'matrix'] If 'sum' is chosen all selected 
        values for one county will be summed up and only 'column' is returned
        for each county. If 'matrix' is chosen, then also the filter columns
        will be returned with the corresponding entries for each selected 
        criterion or value.
    @param file_format File format; either json or h5.
    """
    input_file = os.path.join(os.getcwd(), str(file))
    if file_format == 'json':
        df = pd.read_json(input_file + '.' + file_format)

        date_str = date.strftime('%Y-%m-%d')

        # If data contains a time series, extract the date, 
        if date is not None:
            if 'Date' in df.columns:
                if df['Date'].nunique() > 1:
                    df = df[df['Date'] == date_str]
            else:
                raise gd.DataError('Please provide date column.')   


        # Filter data frame according to filter criterion provided 
        dffilter = pd.Series([True for i in range(len(df))], index = df.index)
        if filters != None:
            for col, vals in filters.items():
                if vals == None:
                    vals = df[col].unique()
                dffilter = dffilter & (df[col].isin(vals))
                

        # Aggregated or matrix output
        df.rename(columns={column : 'Count'}, inplace=True)
        if output == 'sum':
            return df[dffilter].groupby(region_spec).agg({'Count' : sum}).reset_index()
        elif output == 'matrix':
            if filters != None:
                return df[dffilter].loc[:,[region_spec] +  list(filters.keys()) + ['Count']]
            else:
                return df
        else:
            raise gd.DataError("Chose output form.")

    elif file_format == 'h5':
        h5file = h5py.File(input_file + '.' + file_format, 'r')
        # for nested groups, extract the level where the regions to be read in
        # are stored
        if region_spec is not None:
            for elem in region_spec:
                h5file = h5file[elem]

        regions = list(h5file.keys())

        # set no filtering if filters were set to None
        if filters == None:
            filters['Group'] = list(h5file[regions[i]].keys())[:-2] # remove 'Time' and 'Total'
            filters['InfectionState'] = list(range(h5file[regions[i]]['Group1'].shape[1]))

        InfectionStateList = [j for j in filters['InfectionState']]

        # Create data frame to store results to plot
        df = pd.DataFrame(columns=['Time', 'Region', 'Group', 'InfectionState', 'Count'])
        df['Time'] = np.array([date for j in 
                               range(len(regions) * len(filters['Group']) * len(InfectionStateList))])
        
        for i in range(len(regions)):
            # set region identifier in data frame
            region_start_idx = i * len(filters['Group']) * len(InfectionStateList)
            region_end_idx = (i+1) * len(filters['Group']) * len(InfectionStateList)- 1
            df.loc[region_start_idx:region_end_idx, 'Region'] = regions[i]

            # get date index (allows different time points in regions as long 
            # the particular point in time is available)
            date_idx = np.where(np.abs(h5file[regions[i]]['Time'][:] - date)<1e-13)[0]
            if len(date_idx) == 1:
                k = 0
                for group in filters['Group']:
                    # Add (Age)Group identifier
                    df.loc[region_start_idx + k*len(InfectionStateList):
                            region_start_idx+(k+1)*len(InfectionStateList)-1,
                            'Group'] = group

                    # Add InfectionState identifier
                    df.loc[region_start_idx + k*len(InfectionStateList):
                            region_start_idx+(k+1)*len(InfectionStateList)-1,
                            'InfectionState'] = InfectionStateList

                    # Add counts
                    df.loc[region_start_idx + k*len(InfectionStateList):
                            region_start_idx+(k+1)*len(InfectionStateList)-1,
                            'Count'] = np.array(h5file[regions[i]][group][
                        date_idx][0][filters['InfectionState']])     
                    
                    k += 1               
            else:
                raise gd.ValueError("Time point not found for region " + str(regions[i]) + ".")                

        
        # Aggregated or matrix output
        if output == 'sum':
            return df.groupby('Region').agg({'Count' : sum}).reset_index()
        elif output == 'matrix':
            return df
        else:
            raise gd.DataError("Chose output form.")
            
    else:
        raise gd.DataError("Data could not be read in.")

    raise gd.DataError('No data frame can be returned.')
    return pd.DataFrame()

def scale_data_relative(df, path_population):
    """! Scales a population-related data frame relative to the size of the  
    local populations or subpopulations (e.g., if not all age groups are
    considered).

    The first column in the input data frame and the population data frame to
    be read need to be named according to the region identifiers. All regions
    of the data frame to be scaled need to be available in the population
    data set.

    @param df Data frame with population-related information. 
    @param path_population Relative path to population data set.
    @param 
    @return Scaled data set.
    """

    pop_data = pd.read_json(os.path.join(os.getcwd(), path_population))

    # Merge population data of Eisenach (if counted separately) with Wartburgkreis
    if 16063 in pop_data[df.columns[0]].values:
        for i in range(1, len(pop_data.columns)):
            pop_data.loc[pop_data[df.columns[0]] == 16063, pop_data.columns[i]
                         ] += pop_data.loc[pop_data.ID_County == 16056, pop_data.columns[i]]
        pop_data = pop_data[pop_data.ID_County != 16056]   


    return df 

def plot(
    data : pd.DataFrame,
    xlabel : str,
    ylabel : str
):
    # read and filter population data
    # if relative:
        # TODO: 
        # Use create_intervals_mapping() (extract first) and 
        # functionality from get_vaccination_data()

    
    #read and filter map data
    if data.iloc[:, 0].isin(geoger.get_county_ids()).all():
        try:
            map_data = gpd.read_file(os.path.join(os.getcwd(), 'tools/shapes/vg2500_01-01.utm32s.shape/vg2500/vg2500_krs.shp'))
        except FileNotFoundError:
            print_manual_download('Georeferenzierung: UTM32s, Format: shape (ZIP, 3 MB)', 'https://gdz.bkg.bund.de/index.php/default/digitale-geodaten/verwaltungsgebiete/verwaltungsgebiete-1-2-500-000-stand-01-01-vg2500.html')
    else:
        gd.raiseDataError('Provide shape files regions to be plotted.')
    
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

    files_input = {'Reported data' : 'tools/raw/test', 'Sanitized data' : 'tools/raw/test_copy'}
    age_groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80+']
    relative = True

    i = 0
    for file in files_input.values():
        # df.append(read_data(file, region_spec=None, column=None, date=1, filters={
        #           'Group': ['Group1', 'Group2'], 'InfectionState': [3, 4]}, file_format='h5'))
        df = read_data(
            file, region_spec='ID_County',
            column='Vacc_partially',
            date=dt.date(2021, 11, 18),
            filters={'ID_State': [1],
                     'Age_RKI': None},
                file_format='json')

        if relative:
            df = scale_data_relative(df, 'tools/raw/county_current_population.json')

        if i == 0:
            dfs_all = pd.DataFrame(df.iloc[:, 0])

        dfs_all['Count ' + str(i)] = df['Count']


    x=15
    # for vacc in ['Vacc_partially']:
    #     for state in [None]: #[None, 9]:
    #         # plot(files_vacc, dt.date(2021, 11, 18), column = vacc, filters={'ID_State' : [1], 'Age_RKI' : ['05-11','12-17']})
    #         plot(files_vacc, region_spec=None, column=None, date=1, filters={'Group': ['Group1' , 'Group2'], 'InfectionState': [3,4]}, file_format='h5')
    #         plot(files_vacc, region_spec='ID_County', column=vacc, date=dt.date(2021, 11, 18), filters={'ID_State': [1], 'Age_RKI': None}, file_format='h5')

