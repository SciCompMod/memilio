import datetime as dt
import os.path
import h5py

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
from matplotlib.ticker import FormatStrFormatter

compartments = {'Susceptible': 0,
                'Exposed': 1,
                'Infected': 2,
                'Recovered': 3}


def plot_map_nrw(data: pd.DataFrame,
                 scale_colors: np.array([0, 1]),
                 legend: list = [],
                 title: str = '',
                 plot_colorbar: bool = True,
                 output_path: str = '',
                 fig_name: str = 'customPlot',
                 dpi: int = 300,
                 outercolor='white',
                 log_scale=False,
                 cmap='viridis',
                 fontsize=10):
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

    fig = plt.figure(figsize=(3.5 * len(data_columns), 3),
                     facecolor=outercolor)
    # Use n+2 many columns (1: legend + 2: empty space + 3-n: data sets) and
    # n+2 rows where the top row is used for a potential title, the second row
    # for the content and all other rows have height zero.
    height_ratios = [0.05, 1]
    # if len(data_columns) > 1:
    #     height_ratios = height_ratios + [
    #         0.0 for i in range(len(data_columns)-1)]
    if plot_colorbar:
        gs = GridSpec(
            2, len(data_columns)+2, figure=fig,
            width_ratios=[1 for i in range(len(data_columns))]+[0.1, 0.2],
            height_ratios=height_ratios)
    else:
        gs = GridSpec(
            2, len(data_columns), figure=fig,
            width_ratios=[1 for i in range(len(data_columns))],
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
    else:
        norm = mcolors.TwoSlopeNorm(
            vmin=scale_colors[0],
            vmax=scale_colors[1],
            vcenter=0)

    for i in range(len(data_columns)):

        ax = fig.add_subplot(gs[1, i])
        if log_scale:
            map_data.plot(
                data_columns[i],
                ax=ax, legend=False, norm=norm, cmap=cmap, edgecolor='black',
                linewidth=0.1)

        elif cax is not None:
            map_data.plot(
                data_columns[i],
                ax=ax, cax=cax, legend=True, vmin=scale_colors[0],
                vmax=scale_colors[1],
                cmap=cmap, edgecolor='black', linewidth=0.1)
        else:
            # Do not plot colorbar.
            map_data.plot(
                data_columns[i],
                ax=ax, legend=False, vmin=scale_colors[0],
                vmax=scale_colors[1],
                cmap=cmap, edgecolor='black', linewidth=0.1)

        ax.set_title(legend[i], fontsize=fontsize)
        ax.set_axis_off()

    if plot_colorbar:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_ticks([scale_colors[0], scale_colors[1]])
        cbar.set_ticklabels(
            [f'{scale_colors[0]:.3e}', f'{scale_colors[1]:.3e}'],
            fontsize=7)

    plt.subplots_adjust(bottom=0.1, left=0.1)
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, fig_name + '.png'), dpi=dpi)
    plt.close()


def plot_maps(files, output_dir, legend, name=''):

    dfs = {}
    min_vals = []
    max_vals = []

    for date in range(10, 101, 20):
        dfs[date] = extract_nrw_data_and_combine(
            files=files, date=date, relative=True)
        min_vals.append(dfs[date].drop(columns='Region').min().min())
        max_vals.append(dfs[date].drop(columns='Region').max().max())

    min_val = min(min_vals)
    max_val = max(max_vals)

    cmap = 'viridis'
    norm = mcolors.LogNorm(vmin=min_val, vmax=max_val)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_fig, cax = plt.subplots(figsize=(12, 1))
    cbar = plt.colorbar(sm, orientation='horizontal', cax=cax)
    cbar.ax.tick_params(labelsize=12)
    # cbar.set_ticks([min_val, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, max_val])
    # cbar.ax.xaxis.set_major_formatter(FormatStrFormatter('%.5f'))
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.3)
    plt.savefig(os.path.join(output_dir, 'colorbar.png'), dpi=300)
    plt.close()

    for date in range(10, 101, 20):
        plot_map_nrw(
            dfs[date], scale_colors=[min_val, max_val],
            legend=legend,
            title='NRW - Simulation Day '+str(date), plot_colorbar=False,
            output_path=output_dir,
            fig_name=name+str(date), dpi=900,
            outercolor='white',
            log_scale=True,
            fontsize=13)


def plot_difference(files, output_dir, name='difference2D'):
    fig = plt.figure()

    df_dif = pd.DataFrame(columns=['Time', 'difference', 'absolute value'])
    dates = [i for i in range(100)]
    df_dif['Time'] = dates
    df_dif.set_index('Time', inplace=True)

    total_population = 18190422.

    for date in range(100):
        dfs_all = extract_nrw_data_and_combine(
            files=files, date=date, relative=False)
        df_dif.loc[date, 'difference'] = (
            dfs_all[dfs_all.columns[1]] - dfs_all[dfs_all.columns[2]]).sum() / total_population
        df_dif.loc[date, 'absolute value'] = (
            dfs_all[dfs_all.columns[1]] - dfs_all[dfs_all.columns[2]]).abs().sum() / total_population

    df_dif['difference'].plot(label='Relative difference')
    df_dif['absolute value'].plot(
        label='Relative difference in absolute value')
    plt.xlim(left=0., right=101.)
    plt.tight_layout()
    plt.legend()
    plt.grid(linestyle='--')
    plt.savefig(os.path.join(output_dir, name))
    plt.close()


def plot_difference_maps(files, output_dir):

    df_dif1 = pd.DataFrame(columns=['Region'])
    df_dif2 = pd.DataFrame(columns=['Region'])

    for date in range(10, 51, 10):
        dfs_all = extract_nrw_data_and_combine(
            files=files, date=date, relative=True)
        df_dif1['Region'] = dfs_all['Region']
        df_dif1['Count (rel)' + str(date)] = dfs_all[dfs_all.columns[1]
                                                     ] - dfs_all[dfs_all.columns[2]]

    for date in range(60, 101, 10):
        dfs_all = extract_nrw_data_and_combine(
            files=files, date=date, relative=True)
        df_dif2['Region'] = dfs_all['Region']

        df_dif2['Count (rel)' + str(date)] = dfs_all[dfs_all.columns[1]
                                                     ] - dfs_all[dfs_all.columns[2]]

    min_val1 = df_dif1.drop(columns=['Region']).min().min()
    max_val1 = df_dif1.drop(columns=['Region']).max().max()
    min_val2 = df_dif2.drop(columns=['Region']).min().min()
    max_val2 = df_dif2.drop(columns=['Region']).max().max()
    min_val = min(min_val1, min_val2)
    max_val = max(max_val1, max_val2)
    maximum_abs = abs(max([min_val, max_val], key=abs))

    cmap = 'seismic'
    norm = mcolors.TwoSlopeNorm(vmin=-maximum_abs, vmax=maximum_abs, vcenter=0)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_fig, cax = plt.subplots(figsize=(12, 1))
    cbar = plt.colorbar(sm, orientation='horizontal', cax=cax)
    ticks = np.linspace(-maximum_abs, maximum_abs, 10)
    cbar.ax.tick_params(labelsize=13)
    cbar.set_ticks(ticks)
    cbar.ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.3)
    plt.savefig(os.path.join(output_dir, 'colorbar_difference.png'), dpi=900)
    plt.close()

    plot_map_nrw(
        df_dif1, scale_colors=[-maximum_abs, maximum_abs],
        legend=['Day ' + str(date) for date in range(10, 51, 10)],
        title='', plot_colorbar=False,
        output_path=output_dir,
        fig_name="difference10-50", dpi=900,
        outercolor='white',
        log_scale=False,
        cmap='seismic',
        fontsize=17)

    plot_map_nrw(
        df_dif2, scale_colors=[-maximum_abs, maximum_abs],
        legend=['Day ' + str(date) for date in range(60, 101, 10)],
        title='', plot_colorbar=False,
        output_path=output_dir,
        fig_name="difference60-100", dpi=900,
        outercolor='white',
        log_scale=False,
        cmap='seismic')


def extract_nrw_data_and_combine(files, date, relative=True):
    age_groups = {0: '0-4', 1: '5-14', 2: '15-34',
                  3: '35-59', 4: '60-79', 5: '80+'}
    filter_age = None

    i = 0
    for file in files.values():
        model_type = os.path.basename(file).split('_')[0]
        if model_type == 'ode':
            df = pm.extract_data(
                file, region_spec=None, column=None, date=date,
                filters={'Group': filter_age, 'InfectionState': [2]},
                output='matrix',
                file_format=file_format)
            df['Group'] = df.Group.str.extract(r'(\d+)')
            df['Group'] = df['Group'].apply(pd.to_numeric, errors='coerce')
            df['Region'] = (df['Group']-1) // len(age_groups)
            df = df.groupby(['Region'], as_index=False).agg({'Count': "sum"})

            ids = geoger.get_county_ids()
            ids = [id for id in ids if str(id).startswith('5')]

            if len(ids) != len(df):
                raise gd.DataError(
                    "Data is not compatible with number of NRW counties.")

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

    return dfs_all


def plot_total_compartment(
        files, output_dir, legend, compartment='Infected', name='', ax=None,
        print_legend=True):

    colors = ['#2ca02c', '#ff7f0e', '#9C180D']
    file_idx = 0
    if ax is None:
        fig, ax = plt.subplots()
    ax.grid(True, linestyle='--')
    for file in files.values():
        model_type = os.path.basename(file).split('_')[0]
        # Load data.
        h5file = h5py.File(file + '.h5', 'r')
        if model_type == 'ode':
            dates = h5file['1']['Time'][:]
            data = h5file['1']['Total'][:, compartments[compartment]]
            ax.plot(
                dates, data, label=legend[file_idx],
                linewidth=2, color=colors[file_idx])
            ax.set_title(compartment, fontsize=8)
            # ax.set_ylim(bottom=0.)
            ax.set_xlim(left=0., right=dates.max()+1)
        else:
            df = pd.DataFrame()
            regions = list(h5file.keys())
            population = 0
            for i in range(len(regions)):
                for comp in compartments.keys():
                    population += h5file[regions[i]
                                         ]['Total'][:, compartments[comp]]
                df['Region'+str(i)] = h5file[regions[i]
                                             ]['Total'][:, compartments[compartment]]
            df['Total'] = df.sum(axis=1)
            df['Time'] = h5file[regions[0]]['Time'][:]  # hardcoded
            ax.plot(
                df['Time'],
                df['Total'],
                label=legend[file_idx],
                linewidth=1.5, color=colors[file_idx],
                linestyle='--')
            ax.set_title(compartment, fontsize=8)
            ax.set_ylim(bottom=0.)
            ax.set_xlim(left=0., right=df['Time'].max()+1)

        file_idx = file_idx+1
    ax.tick_params(labelsize=7, )
    plt.tight_layout()
    if print_legend:
        plt.legend()
    plt.savefig(os.path.join(output_dir, name + '.png'), dpi=300)

    return ax


def compare_compartments(files, output_dir, legend):
    fig, axs = plt.subplots(
        2, 2, sharex='all')
    axs = axs.flatten()
    for i, compartment in enumerate(compartments.keys()):
        plot_total_compartment(
            files=files, output_dir=output_dir, legend=legend,
            compartment=compartment, ax=axs[i],
            print_legend=False)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)
    lines, labels = axs[0].get_legend_handles_labels()
    fig.legend(lines, labels, ncol=len(models), loc='center',
               fontsize=10, bbox_to_anchor=(0.5, 0.05))
    plt.savefig(
        os.path.join(output_dir, 'compare_all_compartments.png'),
        dpi=300)
    plt.close()


if __name__ == '__main__':

    results = {'Model B': 'results/ode_result_wang_nrw',
               'Model C': 'results/ode_result_nrw',
               'Model D': 'results/graph_result_nrw'}

    file_format = 'h5'

    models = ['Model B (ODE)',
              'Model C (ODE)',
              'Model D (Graph-ODE)']

    plot_dir = os.path.join(os.path.dirname(__file__), '../Plots')

    # plot_maps(files=results, output_dir=plot_dir,
    #           legend=models, name='NRWAdaptiveDay')
    # plot_difference_maps(
    #     files={key: value for key, value in results.items()
    #            if key in {'Model C', 'Model D'}},
    #     output_dir=plot_dir)
    # plot_difference(
    #     files={key: value for key, value in results.items()
    #            if key in {'Model C', 'Model D'}},
    #     output_dir=plot_dir)
    # compare_compartments(files=results, output_dir=plot_dir,  legend=models)
    plot_difference(files={'Old result': 'cpp/build/ode_result_nrw_test', 'New result': 'cpp/build/ode_result_nrw'}, output_dir=plot_dir, name='difference_ode_results_test')
    # plot_difference(files={'Old result': 'results/graph_result_nrw', 'New result': 'cpp/build/graph_result_nrw'}, output_dir=plot_dir, name='difference_graph_results')