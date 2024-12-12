
import datetime as dt
import os.path

import numpy as np
import pandas as pd

import geopandas
from matplotlib.gridspec import GridSpec

from memilio.epidata import geoModificationGermany as geoger

import memilio.epidata.getPopulationData as gpd
from memilio.epidata import getDataIntoPandasDataFrame as gd
import memilio.plot.plotMap as pm

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def plot_map_nrw(data: pd.DataFrame,
             scale_colors: np.array([0, 1]),
             legend: list = [],
             title: str = '',
             plot_colorbar: bool = True,
             output_path: str = '',
             fig_name: str = 'customPlot',
             dpi: int = 300,
             outercolor='white',
             log_scale=False):
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
    @param[in] log_scale Defines if the colorbar is plotted in log scale.
    """
    region_classifier = data.columns[0]
    region_data = data[region_classifier].to_numpy().astype(int)

    data_columns = data.columns[1:]
    # Read and filter map data.
    if np.isin(region_data, geoger.get_county_ids()).all():
        try:
            map_data = geopandas.read_file(
                os.path.join(
                    os.getcwd(),
                    'tools/vg2500_12-31.utm32s.shape/vg2500/VG2500_KRS.shp'))
            if '16056' in map_data.ARS.values:
                map_data = pm.merge_eisenach(map_data)
            # Remove information for plot.
            map_data = map_data[['ARS', 'GEN', 'NUTS', 'geometry']]
            # Use string values as in shape data file.
            data[region_classifier] = data[region_classifier].astype(
                'str').str.zfill(5)
        except FileNotFoundError:
            pm.print_manual_download(
                'Georeferenzierung: UTM32s, Format: shape (ZIP, 5 MB)',
                'https://gdz.bkg.bund.de/index.php/default/verwaltungsgebiete-1-2-500-000-stand-31-12-vg2500-12-31.html')
    else:
        raise gd.DataError('Provide shape files regions to be plotted.')

    # Remove regions that are not input data table.
    map_data = map_data[map_data.ARS.isin(data[region_classifier])]

    data['new_index'] = map_data.index.array
    data = data.set_index('new_index')

    map_data[data_columns] = data.loc[:, data_columns] 

    for i in range(len(data_columns)):
        if legend[i] == '':
            fname = 'data_column_' + str(i)
        else:
            fname = str(legend[i].replace(' ', '_'))
        pm.save_interactive(data[data_columns[i]], os.path.join(
            output_path, fname) + '.html', map_data, scale_colors)

    fig = plt.figure(figsize=(3.5 * len(data_columns), 3), facecolor=outercolor)
    # Use n+2 many columns (1: legend + 2: empty space + 3-n: data sets) and
    # n+2 rows where the top row is used for a potential title, the second row
    # for the content and all other rows have height zero.
    height_ratios = [0.05, 1, 0]
    if len(data_columns) > 1:
        height_ratios = height_ratios + [
            0.0 for i in range(len(data_columns)-1)]
    gs = GridSpec(
        len(data_columns) + 2, len(data_columns) + 2, figure=fig,
        width_ratios=[1 for i in range(len(data_columns))] + [0.1, 0.2],
        height_ratios=height_ratios)

    # Use top row for title.
    tax = fig.add_subplot(gs[0, :])
    tax.set_axis_off()
    tax.set_title(title, fontsize=16)
    if plot_colorbar:
        # Prepare colorbar.
        cax = fig.add_subplot(gs[1, -2])

    else:
        cax = None

    if log_scale:
        norm = mcolors.LogNorm(vmin=scale_colors[0], vmax=scale_colors[1])

    for i in range(len(data_columns)):

        cmap = 'viridis'
        ax = fig.add_subplot(gs[1, i])
        if log_scale:
            map_data.plot(data_columns[i], ax=ax, legend=False,
                          norm=norm, cmap=cmap)

        elif cax is not None:
            map_data.plot(data_columns[i], ax=ax, cax=cax, legend=True,
                          vmin=scale_colors[0], vmax=scale_colors[1])
        else:
            # Do not plot colorbar.
            map_data.plot(data_columns[i], ax=ax, legend=False,
                          vmin=scale_colors[0], vmax=scale_colors[1])

        ax.set_title(legend[i], fontsize=12)
        ax.set_axis_off()

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_ticks([scale_colors[0], scale_colors[1]])
    cbar.set_ticklabels([f'{scale_colors[0]:.4e}', f'{scale_colors[1]:.4e}'])

    plt.subplots_adjust(bottom=0.1)
    plt.savefig(os.path.join(output_path, fig_name + '.png'), dpi=dpi)


if __name__ == '__main__':

    files_input = {'Data set 1': 'cpp/build/ode_result_nrw', 
                   'Data set 2': 'cpp/build/graph_result_nrw'} # Result file of equation-based model has to be first
    file_format = 'h5'
    # Define age groups which will be considered through filtering
    # Keep keys and values as well as its assignment constant, remove entries
    # if only part of the population should be plotted or considered, e.g., by
    # setting:
    # age_groups = {1: '5-14', 2: '15-34'}
    age_groups = {0: '0-4', 1: '5-14', 2: '15-34',
                  3: '35-59', 4: '60-79', 5: '80+'}
    if len(age_groups) == 6:
        filter_age = None
    else:
        if file_format == 'json':
            filter_age = [val for val in age_groups.values()]
        else:
            filter_age = ['Group' + str(key) for key in age_groups.keys()]

    relative = True

    date = 14

    i = 0
    for file in files_input.values():
        # MEmilio backend hdf5 example
        if(i == 0): # result file of equation-based model has to be first
            df = pm.extract_data(
                file, region_spec=None, column=None, date=date,
                filters={'Group': filter_age, 'InfectionState': [2]},
                output='matrix',
                file_format=file_format)
            df['Group'] = df.Group.str.extract('(\d+)')
            df['Group'] = df['Group'].apply(pd.to_numeric, errors='coerce')
            # df['Region'] = df['Group']
            df['Region'] = (df['Group']-1) // len(age_groups)
            df = df.groupby(['Region'], as_index=False).agg({'Count': "sum"})

            ids = geoger.get_county_ids()
            ids = [id for id in ids if str(id).startswith('5')]
            # ids = [5111, 5112, 5113]

            if len(ids) != len(df):
                raise gd.DataError("Data is not compatible with number of NRW counties.")

            df['Region'] = ids
        else:
            df = pm.extract_data(
                file, region_spec=None, column=None, date=date,
                filters={'Group': filter_age, 'InfectionState': [2]},
                file_format=file_format)

        df = df.apply(pd.to_numeric, errors='coerce')

        if relative:

            try:
                population = pd.read_json(
                    'data/pydata/Germany/county_current_population.json')
            # pandas>1.5 raise FileNotFoundError instead of ValueError
            except (ValueError, FileNotFoundError):
                print("Population data was not found. Download it from the internet.")
                population = gpd.get_population_data(
                    read_data=False, file_format=file_format,
                    out_folder='data/pydata/Germany/', no_raw=True,
                    merge_eisenach=True)

            # For fitting of different age groups we need format ">X".
            age_group_values = list(age_groups.values())
            age_group_values[-1] = age_group_values[-1].replace('80+', '>79')
            # scale data
            df = pm.scale_dataframe_relative(df, age_group_values, population)

        if i == 0:
            dfs_all = pd.DataFrame(df.iloc[:, 0])

        dfs_all[df.columns[-1] + ' ' + str(i)] = df[df.columns[-1]]
        i += 1

    min_val = dfs_all[dfs_all.columns[1:]].min().min()
    max_val = dfs_all[dfs_all.columns[1:]].max().max()

    plot_map_nrw(
        dfs_all, scale_colors=[min_val, max_val],
        legend=['', ''],
        title='NRW - Simulation Day 10', plot_colorbar=True,
        output_path=os.path.dirname(__file__),
        fig_name='NRWPlot', dpi=300,
        outercolor='white', 
        log_scale=True)
