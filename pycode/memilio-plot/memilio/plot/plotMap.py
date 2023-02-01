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
import os

# in case of a necessary manual installation of GDAL and Fiona on Windows, see
# https://stackoverflow.com/questions/69521550/importerror-the-read-file-function-requires-the-fiona-package-but-it-is-no
import geopandas
import h5py
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import modifyDataframeSeries as mdfs


def print_manual_download(filename, url):
    """! Prints message to ask the user to manually download a file.

    @param[in] filename Filename of file needed.
    @param[in] url URL to file.
    """
    print(
        'This script needs manual downloading of files. Please download '
        + filename + ' from ' + url + 'and move it the extracted folder'
        'to the current working directory under tools/.')


def merge_eisenach(map_data: geopandas.GeoDataFrame):
    """! Merges geometries for Eisenach with Wartburgkreis of Geopandas 
    dataframe.

    @param[in,out] map_data GeoPandasDataFrame from county shape file.
    @return Frame with Wartburgkreis and Eisenach merged. 
    """
    wartburg = map_data.ARS == '16063'
    eisenach = map_data.ARS == '16056'
    map_data.loc[wartburg, 'geometry'] = [
        map_data[wartburg].geometry.values[0].union(
            map_data[eisenach].geometry.values[0])]
    # Remove Eisenach and return.
    return map_data.drop(map_data[eisenach].index.values[0])


def extract_data(
        file, region_spec, column, date, filters=None, output='sum',
        file_format='json'):
    """ Reads data from a general json or specific hdf5 file as output by the 
    MEmilio simulation framework and extracts parts of the data if filters are
     applied. The input regional data can be a time series data where a 
    particular date is extracted or single time data which is used rightaway.
    If the file contains multiple features per region,
    particular columns need to be specified.

    @param[in] file Path and filename of file to be read in, relative from current
        directory.
    @param[in] region_spec Specificier for region, for json, e.g., column name
        of county IDs for hdf5 level of nesting X where to retrieve
        fileX.keys(), X could either by empty, such that file.keys() is taken
        or a list of nested keys [X0, X1, ...] such that 
        file[X0][X1][...].keys()
    @param[in] column Column with values that will be plotted.
    @param[in] data Date to be extracted from data frame if it contains a time 
        series or if single or no date information is in the input table, 
        provide date=None. For json, pd.date Objects need to be used, for
        h5, provide an integer from the beginning of the time series.
    @param[in] filters Dictionary with columns and values for filtering rows.
        None for a column or all filters means that all values are taken. 
        For MEmilio h5 Files only 'Group' (AgeGroups from 'Group1' to 'Total') 
        and 'InfectionState' (from 0 to InfetionState::Count-1) can be filtered.
    @param[in] output [Either 'sum' or 'matrix'] If 'sum' is chosen all selected 
        values for one county will be summed up and only 'column' is returned
        for each county. If 'matrix' is chosen, then also the filter columns
        will be returned with the corresponding entries for each selected 
        criterion or value.
    @param[in] file_format File format; either json or h5.
    @return Extracted data set.
    """
    input_file = os.path.join(os.getcwd(), str(file))
    if file_format == 'json':
        df = pd.read_json(input_file + '.' + file_format)

        date_str = date.strftime('%Y-%m-%d')

        # If data contains a time series, extract the date.
        if date is not None:
            if 'Date' in df.columns:
                if df['Date'].nunique() > 1:
                    df = df[df['Date'] == date_str]
            else:
                raise gd.DataError('Please provide date column.')

        # Filter data frame according to filter criterion provided.
        dffilter = pd.Series([True for i in range(len(df))], index=df.index)
        if filters != None:
            for col, vals in filters.items():
                if vals == None:
                    vals = df[col].unique()
                dffilter = dffilter & (df[col].isin(vals))

        # Aggregated or matrix output.
        df.rename(columns={column: 'Count'}, inplace=True)
        if output == 'sum':
            return df[dffilter].groupby(region_spec).agg(
                {'Count': sum}).reset_index()
        elif output == 'matrix':
            if filters != None:
                return df[dffilter].loc[:, [region_spec] +
                                        list(filters.keys()) + ['Count']]
            else:
                return df
        else:
            raise gd.DataError("Chose output form.")

    elif file_format == 'h5':
        h5file = h5py.File(input_file + '.' + file_format, 'r')
        # For nested groups, extract the level where the regions to be read in
        # are stored.
        if region_spec is not None:
            for elem in region_spec:
                h5file = h5file[elem]

        regions = list(h5file.keys())

        # Set no filtering if filters were set to None.
        if filters == None:
            filters['Group'] = list(h5file[regions[i]].keys())[
                :-2]  # Remove 'Time' and 'Total'.
            filters['InfectionState'] = list(
                range(h5file[regions[i]]['Group1'].shape[1]))

        InfectionStateList = [j for j in filters['InfectionState']]

        # Create data frame to store results to plot.
        df = pd.DataFrame(
            columns=['Time', 'Region', 'Group', 'InfectionState', 'Count'])
        df['Time'] = np.array(
            [date
             for j in range(
                 len(regions) * len(filters['Group'])
                 * len(InfectionStateList))])

        for i in range(len(regions)):
            # Set region identifier in data frame.
            region_start_idx = i * \
                len(filters['Group']) * len(InfectionStateList)
            region_end_idx = (
                i+1) * len(filters['Group']) * len(InfectionStateList) - 1
            df.loc[region_start_idx:region_end_idx, 'Region'] = regions[i]

            # Get date index (allows different time points in regions as long
            # the particular point in time is available).
            date_idx = np.where(
                np.abs(h5file[regions[i]]['Time'][:] - date) < 1e-13)[0]
            if len(date_idx) == 1:
                k = 0
                for group in filters['Group']:
                    # Add (Age)Group identifier.
                    df.loc[region_start_idx + k*len(InfectionStateList):
                           region_start_idx+(k+1)*len(InfectionStateList)-1,
                           'Group'] = group

                    # Add InfectionState identifier.
                    df.loc[region_start_idx + k*len(InfectionStateList):
                           region_start_idx+(k+1)*len(InfectionStateList)-1,
                           'InfectionState'] = InfectionStateList

                    # Add counts.
                    df.loc[region_start_idx + k*len(InfectionStateList):
                           region_start_idx+(k+1)*len(InfectionStateList)-1,
                           'Count'] = np.array(h5file[regions[i]][group][
                               date_idx][0][filters['InfectionState']])

                    k += 1
            else:
                raise gd.ValueError(
                    "Time point not found for region " + str(regions[i]) + ".")

        # Aggregated or matrix output.
        if output == 'sum':
            return df.groupby('Region').agg({'Count': sum}).reset_index()
        elif output == 'matrix':
            return df
        else:
            raise gd.DataError("Chose output form.")

    else:
        raise gd.DataError("Data could not be read in.")


def scale_dataframe_relative(df, age_groups, df_population):
    """! Scales a population-related data frame relative to the size of the  
    local populations or subpopulations (e.g., if not all age groups are
    considered).

    The first column in the input data frame and the population data frame to
    be read need to be named according to the region identifiers. All regions
    of the data frame to be scaled need to be available in the population
    data set.

    @param[in,out] df Data frame with population-related information. 
    @param[in] age_groups Considered age groups of population data taken into 
        summation.
    @param[in] df_population Data frame with population data set.
    @return Scaled data set.
    """

    # Merge population data of Eisenach (if counted separately) with Wartburgkreis.
    if 16056 in df_population[df.columns[0]].values:
        for i in range(1, len(df_population.columns)):
            df_population.loc[df_population[df.columns[0]] == 16063, df_population.columns[i]
                              ] += df_population.loc[df_population.ID_County == 16056, df_population.columns[i]]
        df_population = df_population[df_population.ID_County != 16056]

    df_population_agegroups = pd.DataFrame(
        columns=[df_population.columns[0]] + age_groups)
    # Extrapolate on oldest age group with maximumg age 100.
    for region_id in df.iloc[:, 0]:
        df_population_agegroups.loc[len(df_population_agegroups.index), :] = [region_id] + list(
            mdfs.fit_age_group_intervals(df_population[df_population.iloc[:, 0] == region_id].iloc[:, 2:], age_groups))

    def scale_row(elem):
        population_local_sum = df_population_agegroups[
            df_population_agegroups[df.columns[0]] == elem[0]].iloc[
                :, 1:].sum(axis=1)
        return elem['Count'] / population_local_sum.values[0]

    df_population_agegroups.loc[:, age_groups].sum(axis=1)

    df['Count (rel)'] = df.apply(scale_row, axis=1)

    return df


# Save interactive html files.
def save_interactive(col, filename, map_data, scale_colors):
    """! Plots region-specific information in an interactive html map.

    @param[in] col The column that will be plotted.
    @param[in] filename Filename with path that determines the output directory.
    @param[in] map_data Geopandas file with plot data.
    @param[in] scale_colors Array of min-max-values to scale colorbar.
    """
    map_data.explore(col, legend=True, vmin=scale_colors[0],
                     vmax=scale_colors[1]).save(filename)


def plot_map(data: pd.DataFrame,
             scale_colors: np.array([0, 1]),
             legend: list = [],
             title: str = '',
             plot_colorbar: bool = True,
             output_path: str = '',
             fig_name: str = 'customPlot',
             dpi: int = 300,
             outercolor='white'):
    """! Plots region-specific information onto a interactive html map and
    returning svg and png image. Allows the comparisons of a variable list of
    data sets.

    @param[in] data Data to be plotted. First column must contain regional 
        specifier, following columns will be plotted for comparison.
    @param[in] scale_colors Array of min-max-values to scale colorbar.
    @param[in] legend List subtitles for different columns. Can be list of 
        empty strings.
    @param[in] title Title of the plot.
    @param[in] plot_colorbar Defines if a colorbar will be plotted.
    @param[in] output_path Output path for the figure.   
    @param[in] fig_name Name of the figure created.
    @param[in] dpi Dots-per-inch value for the exported figure.
    @param[in] outercolor Background color of the plot image.
    """
    region_classifier = data.columns[0]

    data_columns = data.columns[1:]
    # Read and filter map data.
    if data[region_classifier].isin(geoger.get_county_ids()).all():
        try:
            map_data = geopandas.read_file(
                os.path.join(
                    os.getcwd(),
                    'tools/vg2500_12-31.utm32s.shape/vg2500/VG2500_KRS.shp'))
            if '16056' in map_data.ARS.values:
                map_data = merge_eisenach(map_data)
            # Remove information for plot.
            map_data = map_data[['ARS', 'GEN', 'NUTS', 'geometry']]
            # Use string values as in shape data file.
            data[region_classifier] = data[region_classifier].astype(
                'str').str.zfill(5)
        except FileNotFoundError:
            print_manual_download(
                'Georeferenzierung: UTM32s, Format: shape (ZIP, 5 MB)',
                'https://gdz.bkg.bund.de/index.php/default/verwaltungsgebiete-1-2-500-000-stand-31-12-vg2500-12-31.html')
    else:
        raise gd.DataError('Provide shape files regions to be plotted.')

    # Remove regions that are not input data table.
    map_data = map_data[map_data.ARS.isin(data[region_classifier])]

    map_data[data_columns] = data.loc[:, data_columns]

    for i in range(len(data_columns)):
        if legend[i] == '':
            fname = 'data_column_' + str(i)
        else:
            fname = str(legend[i].replace(' ', '_'))
        save_interactive(data[data_columns[i]], os.path.join(
            output_path, fname) + '.html', map_data, scale_colors)

    fig = plt.figure(figsize=(4 * len(data_columns), 6), facecolor=outercolor)
    # Use n+2 many columns (1: legend + 2: empty space + 3-n: data sets) and
    # n+2 rows where the top row is used for a potential title, the second row
    # for the content and all other rows have height zero.
    height_ratios = [0.05, 1, 0]
    if len(data_columns) > 1:
        height_ratios = height_ratios + [
            0.0 for i in range(len(data_columns)-1)]
    gs = GridSpec(
        len(data_columns) + 2, len(data_columns) + 2, figure=fig, wspace=0.1,
        width_ratios=[0.05, 0.05] + [1 for i in range(len(data_columns))],
        height_ratios=height_ratios)

    # Use top row for title.
    tax = fig.add_subplot(gs[0, :])
    tax.set_axis_off()
    tax.set_title(title, fontsize=18)
    if plot_colorbar:
        # Prepare colorbar.
        cax = fig.add_subplot(gs[1, 0])
    else:
        cax = None

    for i in range(len(data_columns)):

        ax = fig.add_subplot(gs[:, i+2])
        if cax is not None:
            map_data.plot(data_columns[i], ax=ax, cax=cax, legend=True,
                          vmin=scale_colors[0], vmax=scale_colors[1])
        else:
            # Do not plot colorbar.
            map_data.plot(data_columns[i], ax=ax, legend=False,
                          vmin=scale_colors[0], vmax=scale_colors[1])

        ax.set_title(legend[i], fontsize=12)
        ax.set_axis_off()

    plt.subplots_adjust(bottom=0.1)

    plt.savefig(os.path.join(output_path, fig_name + '.png'), dpi=dpi)
    plt.savefig(os.path.join(output_path, fig_name + '.svg'), dpi=dpi)

    plt.show()
