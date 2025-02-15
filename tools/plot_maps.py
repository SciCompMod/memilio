import datetime as dt
import os
import imageio
import warnings

import geopandas
import numpy as np
import pandas as pd

from matplotlib.gridspec import GridSpec
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import SymLogNorm, LinearSegmentedColormap, Normalize
from mpl_toolkits.axes_grid1.inset_locator import inset_axes  # Local import
import matplotlib as mpl

from tqdm.auto import tqdm


# Custom modules (make sure these are in your PYTHONPATH)
import memilio.epidata.getPopulationData as gpd_data
import memilio.plot.plotMap as pm
from memilio.epidata import geoModificationGermany as geoger

# Ignore FutureWarnings
warnings.filterwarnings("ignore", category=FutureWarning)


def print_manual_download(filename: str, url: str) -> None:
    """
    Inform the user that a file must be manually downloaded.

    Parameters:
        filename (str): Name of the required file.
        url (str): URL where the file can be downloaded.
    """
    message = (
        "This script requires manual downloading of files. Please download "
        f"{filename} from {url} and move the extracted folder to the current working "
        "directory under tools/."
    )
    print(message)


def merge_eisenach(map_data: geopandas.GeoDataFrame) -> geopandas.GeoDataFrame:
    """
    Merge Eisenach with Wartburgkreis in the GeoDataFrame.

    Parameters:
        map_data (GeoDataFrame): GeoDataFrame from the county shape file.

    Returns:
        GeoDataFrame: Modified GeoDataFrame with Eisenach merged into Wartburgkreis.
    """
    wartburg = map_data.ARS == '16063'
    eisenach = map_data.ARS == '16056'
    merged_geom = map_data[wartburg].geometry.values[0].union(
        map_data[eisenach].geometry.values[0]
    )
    map_data.loc[wartburg, 'geometry'] = [merged_geom]
    return map_data.drop(map_data[eisenach].index.values[0])


def plot_risk_map(path_results: str, path_plots: str, days: list, percentile: str) -> None:
    """
    Plot risk maps using the provided results.

    Parameters:
        path_results (str): Path to the results.
        path_plots (str): Directory where plots will be saved.
        days (list): List of days for which to plot.
        percentile (str): Percentile key for results.
    """
    if not os.path.exists(path_plots):
        os.makedirs(path_plots)
    path_risk_results = os.path.join(path_results, "risk")
    plot_maps(
        path_results=path_risk_results,
        path_plots=path_plots,
        compartments=[0],
        percentile=percentile,
        days=days,
        min_val=0,
        max_val=1,
        filename="risk_map",
        relative=False,
        age_groups={0: '0-4'}
    )


def plot_icu_map(path_results: str, path_plots: str, days: list, percentile: str, max_val: float, icu_cap: float) -> None:
    """
    Plot ICU maps based on the provided results.

    Parameters:
        path_results (str): Path to the results.
        path_plots (str): Directory where plots will be saved.
        days (list): List of days for which to plot.
        percentile (str): Percentile key for results.
        max_val (float): Maximum value for color normalization.
        icu_cap (float): ICU capacity value.
    """
    if not os.path.exists(path_plots):
        os.makedirs(path_plots)
    plot_maps(
        path_results=path_results,
        path_plots=path_plots,
        compartments=[7],
        percentile=percentile,
        days=days,
        min_val=0,
        max_val=max_val,
        filename="icu_map",
        relative=True,
        icu_cap=icu_cap
    )


def create_colorbar(path_plots: str, norm: Normalize, title: str) -> None:
    """
    Create and save a colorbar figure.

    Parameters:
        path_plots (str): Directory where the colorbar image will be saved.
        norm (Normalize): Normalization for the colormap.
        title (str): Title for the colorbar image file.
    """
    colors = ["green", "yellow", "red", "purple"]
    cmap = LinearSegmentedColormap.from_list("my_colormap", colors)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # Create a horizontal colorbar
    cbar_fig, ax = plt.subplots(figsize=(8, 1))
    plt.colorbar(sm, orientation='horizontal', cax=ax)
    ax.tick_params(labelsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(path_plots, f"{title}_colorbar.png"), dpi=300)
    plt.clf()


def plot_map(norm: SymLogNorm,
             data: pd.DataFrame,
             scale_colors: np.array,
             legend: list = [],
             title: str = '',
             plot_colorbar: bool = True,
             output_path: str = '',
             fig_name: str = 'customPlot',
             dpi: int = 300,
             outercolor: str = 'white') -> None:
    """
    Plot region-specific data on a map and save the figure as a PNG image.
    The first column in `data` must contain regional identifiers, and subsequent
    columns will be visualized.

    Parameters:
        norm (SymLogNorm): Normalization for the colormap.
        data (pd.DataFrame): DataFrame with region IDs in the first column and data columns thereafter.
        scale_colors (np.array): Array with minimum and maximum values for color scaling.
        legend (list, optional): List of subtitles for each data column.
        title (str, optional): Title of the plot.
        plot_colorbar (bool, optional): Whether to plot a colorbar.
        output_path (str, optional): Path where the figure will be saved.
        fig_name (str, optional): Base name for the saved figure file.
        dpi (int, optional): Resolution in dots per inch.
        outercolor (str, optional): Background color for the figure.
    """
    # Ensure a legend entry exists for each data column.
    data_columns = data.columns[1:]
    if legend is None or len(legend) < len(data_columns):
        legend = (legend or []) + [''] * \
            (len(data_columns) - len(legend or []))

    # Convert region identifiers to the required format.
    region_col = data.columns[0]
    region_ids = data[region_col].to_numpy().astype(int)

    # Create custom colormap.
    colors = ["green", "yellow", "red", "purple"]
    cmap = LinearSegmentedColormap.from_list("custom_colormap", colors)

    # Obtain mapping from county IDs to state IDs.
    county_to_state = geoger.get_countyid_to_stateid_map(merge_eisenach=True)

    # Check if all provided region IDs are valid.
    if np.isin(region_ids, geoger.get_county_ids()).all():
        try:
            shape_path = os.path.join(
                os.getcwd(),
                'tools/vg2500_12-31.utm32s.shape',
                'vg2500',
                'VG2500_KRS.shp'
            )
            map_data = geopandas.read_file(shape_path)
            if '16056' in map_data.ARS.values:
                map_data = merge_eisenach(map_data)
            map_data = map_data[['ARS', 'GEN', 'NUTS', 'geometry']]
            data[region_col] = data[region_col].astype(str).str.zfill(5)
        except FileNotFoundError:
            print_manual_download(
                'Georeferenzierung: UTM32s, Format: shape (ZIP, 5 MB)',
                'https://gdz.bkg.bund.de/index.php/default/verwaltungsgebiete-1-2-500-000-stand-31-12-vg2500-12-31.html'
            )
            return
    else:
        raise ValueError(
            "Provided regional identifiers do not match the expected shape file regions.")

    # Merge input data with the map geometries.
    map_data = map_data[map_data.ARS.isin(data[region_col])].copy()
    map_data = map_data.merge(data, left_on='ARS', right_on=region_col)

    # Add state IDs based on county mapping.
    map_data['state_id'] = map_data['ARS'].astype(
        int).map(county_to_state).fillna(0).astype(int)

    # Dissolve counties into states for plotting state borders.
    states = map_data.dissolve(by='state_id')

    # Create the figure and set the title.
    num_columns = len(data_columns)
    fig = plt.figure(figsize=(4 * num_columns, 5), facecolor=outercolor)
    fig.suptitle(title, fontsize=16)

    # Create grid layout: first column for colorbar (if used), then one per data column.
    gs = GridSpec(
        1, num_columns + 1, figure=fig, wspace=0.1,
        width_ratios=[0.05] + [1] * num_columns, height_ratios=[1]
    )

    if plot_colorbar:
        cax = fig.add_subplot(gs[0, 0])
        sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, cax=cax)
    else:
        cax = None

    # Plot each data column.
    for i, col in enumerate(data_columns):
        ax = fig.add_subplot(gs[0, i + 1])
        map_data.plot(
            column=col,
            ax=ax,
            cax=cax,
            legend=False,
            norm=norm,
            cmap=cmap,
            edgecolor='gray',
            linewidth=0.5
        )
        # Overlay state borders.
        states.boundary.plot(ax=ax, edgecolor='darkgray', linewidth=1.0)
        if legend[i]:
            ax.set_title(legend[i], fontsize=12)
        ax.set_axis_off()

    plt.subplots_adjust(left=0, right=1, top=0.9, bottom=0)
    output_file = os.path.join(output_path, f"{fig_name}.png")
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', pad_inches=0)
    plt.close(fig)


def plot_maps(path_results: str,
              path_plots: str,
              compartments: list,
              percentile: str,
              days: list,
              min_val: float,
              max_val: float,
              filename: str = "data",
              relative: bool = True,
              age_groups: dict = {0: '0-4', 1: '5-14',
                                  2: '15-34', 3: '35-59', 4: '60-79', 5: '80+'},
              flows: bool = False,
              path_results2: str = "",
              icu_cap: float = 1) -> None:
    """
    Create maps for the specified days using the extracted data.

    Parameters:
        path_results (str): Path to the primary results.
        path_plots (str): Directory where plots will be saved.
        compartments (list): List of compartments to extract.
        percentile (str): Percentile key for the results.
        days (list): List of days for which to generate maps.
        min_val (float): Minimum value for color normalization.
        max_val (float): Maximum value for color normalization.
        filename (str, optional): Base filename for the maps.
        relative (bool, optional): Whether to scale data relative to the population.
        age_groups (dict, optional): Mapping of age group indices to age range labels.
        flows (bool, optional): If True, compute daily flows from cumulative data.
        path_results2 (str, optional): Optional second results path for comparison.
        icu_cap (float, optional): ICU capacity for relative scaling.
    """
    progress_bar = tqdm(total=len(days))
    path = os.path.join(path_results, percentile, "Results")
    get_max_val = 0

    # Define the input files.
    files_input = {'Data set 1': path}
    if path_results2:
        files_input['Data set 2'] = os.path.join(
            path_results2, percentile, "Results")

    file_format = 'h5'
    # Determine the filter for age groups.
    if len(age_groups) == 6:
        filter_age = None
    else:
        filter_age = (
            [val for val in age_groups.values()]
            if file_format == 'json'
            else ['Group' + str(key + 1) for key in age_groups.keys()]
        )

    norm = SymLogNorm(linthresh=1, linscale=0.7, vmin=min_val, vmax=max_val)
    create_colorbar(path_plots, norm, filename)

    for day in days:
        i = 0
        for file in files_input.values():
            # Extract data for the current day.
            df = pm.extract_data(
                file,
                region_spec=None,
                column=None,
                date=day,
                filters={'Group': filter_age, 'InfectionState': compartments},
                file_format=file_format
            )
            df = df.apply(pd.to_numeric, errors='coerce')

            if flows:
                if day > 0:
                    df_previous = pm.extract_data(
                        file,
                        region_spec=None,
                        column=None,
                        date=day - 1,
                        filters={'Group': filter_age,
                                 'InfectionState': compartments},
                        file_format=file_format
                    )
                    df['Count'] = df['Count'] - df_previous['Count']
                    if df['Count'].min() < 0:
                        print("Negative values in data.")
                        df['Count'] = df['Count'].clip(lower=0)

            if relative:
                try:
                    population = pd.read_json(
                        'data/pydata/Germany/county_current_population.json')
                except (ValueError, FileNotFoundError):
                    print(
                        "Population data was not found. Downloading from the internet...")
                    population = gpd_data.get_population_data(
                        read_data=False,
                        file_format=file_format,
                        out_folder='data/pydata/Germany/',
                        no_raw=True,
                        split_gender=False,
                        merge_eisenach=True
                    )
                # Adjust the age group format if necessary.
                age_group_values = list(age_groups.values())
                age_group_values[-1] = age_group_values[-1].replace(
                    '80+', '>79')
                df = pm.scale_dataframe_relative(
                    df, age_group_values, population)
                df['Count'] = df['Count (rel)'] * 100_000 / icu_cap
                df = df.drop(columns=['Count (rel)'])

            if i == 0:
                dfs_all = pd.DataFrame(df.iloc[:, 0])
            dfs_all[df.columns[-1] + ' ' + str(i)] = df[df.columns[-1]]
            i += 1

        fn = f"{filename}_day_{day}"
        dfs_all = dfs_all.apply(pd.to_numeric, errors='coerce')
        dfs_all_sorted = dfs_all.sort_values(
            by='Region').reset_index(drop=True)

        # Update maximum value if necessary.
        if dfs_all_sorted['Count 0'].max() > get_max_val:
            if not path_results2:
                get_max_val = dfs_all_sorted['Count 0'].max()
            else:
                get_max_val = max(
                    dfs_all_sorted['Count 0'].max(), dfs_all_sorted['Count 1'].max())

        plot_map(
            norm=norm,
            data=dfs_all_sorted,
            scale_colors=[min_val, max_val],
            legend=['', ''],
            title='',
            plot_colorbar=False,
            output_path=path_plots,
            fig_name=fn,
            dpi=300,
            outercolor='white'
        )
        progress_bar.update()
    progress_bar.close()
    print("max value:", get_max_val)


def create_combined_plots(kmax: str = '0.80', days: list = None, base_dir: str = None) -> None:
    """
    Creates combined plots by overlaying ICU images with a risk inset for each configuration.
    Searches in the folder 'plots/ICUCap_9.000000' for subdirectories corresponding to different
    configurations and creates a combined plot per configuration.

    Parameters:
        kmax (str): The kmax value as a string (e.g., '0.60').
        days (list, optional): List of days to include in the combined plot. Defaults to [80, 100, 120, 160].
        base_dir (str, optional): Base directory. Defaults to the current working directory.
    """
    if base_dir is None:
        base_dir = os.getcwd()
    if days is None:
        days = [80, 100, 120, 160]

    paths = [
        f'fs_bf0_kmin_0.000000_kmax_{kmax}0000',
        f'fs_kmin_0.000000_kmax_{kmax}0000',
        f'rp_kmin_0.000000_kmax_{kmax}0000',
        f'rp_735_kmin_0.000000_kmax_{kmax}0000'
    ]

    for path in paths:
        path_plots = os.path.join(base_dir, 'plots/ICUCap_9.000000', path)

        # List only files (exclude directories)
        list_all_files = [file for file in os.listdir(path_plots)
                          if os.path.isfile(os.path.join(path_plots, file))]

        # Filter files that correspond to the specified days.
        filtered_files = []
        for file in list_all_files:
            day_in_file = ''.join(filter(str.isdigit, file))
            if day_in_file.isdigit() and int(day_in_file) in days:
                filtered_files.append(file)

        # Create subplots: one subplot per day.
        fig, axs = plt.subplots(1, len(days), figsize=(16, 4))

        for i, day in enumerate(days):
            # Look for ICU and risk files containing the day as substring.
            icu_files = [
                file for file in filtered_files if 'icu' in file and str(day) in file]
            risk_files = [
                file for file in filtered_files if 'risk' in file and str(day) in file]
            if not icu_files or not risk_files:
                print(
                    f"Missing ICU or risk file for day {day} in {path_plots}")
                continue
            icu_file = icu_files[0]
            risk_file = risk_files[0]

            icu_img = plt.imread(os.path.join(path_plots, icu_file))
            risk_img = plt.imread(os.path.join(path_plots, risk_file))

            title = f'Day {day}' if day != 198 else 'Day 200'
            axs[i].imshow(icu_img)
            axs[i].set_title(title, fontsize=18)
            axs[i].axis('off')

            # Create an inset axis for the risk plot.
            axins = inset_axes(
                axs[i],
                width="40%",   # Relative width of inset
                height="40%",  # Relative height of inset
                loc='upper right',
                bbox_to_anchor=(0, 0, 1.4, 1),
                bbox_transform=axs[i].transAxes,
                borderpad=0
            )
            axins.imshow(risk_img)
            axins.axis('off')

        plt.tight_layout()
        save_dir = os.path.join(base_dir, 'plots/ICUCap_9.000000', 'new')
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(
            save_dir, f'{path}_combined_plot_adjusted.png')
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        plt.close(fig)
    plt.show()


if __name__ == '__main__':
    # Configuration and paths
    modes = ["FeedbackDamping"]
    path_cwd = os.getcwd()
    icu_cap = [6, 9, 12, 15]
    cap_indx = 1

    kmax = '0.80'
    path_results = os.path.join(
        path_cwd, "results", "ICUCap_9.000000", "rho_1.000000",
        "BlendingFactorRegional_0.000000", f"kmin_0.000000_kmax_{kmax}0000", "FeedbackDamping"
    )
    path_plots = os.path.join(
        path_cwd, "plots", f"ICUCap_{icu_cap[cap_indx]}.000000",
        f"fs_bf0_kmin_0.000000_kmax_{kmax}0000"
    )
    path_data = os.path.join(path_cwd, "data")
    num_days = 198
    percentile = "p50"
    days = list(range(0, num_days + 1, 20))
    days.append(num_days)

    # Generate maps for risk and ICU
    plot_risk_map(path_results, path_plots, days, percentile)
    # plot_icu_map(path_results, path_plots, days,
    #              percentile, 1, icu_cap[cap_indx])

    # Create combined plots
    # create_combined_plots(kmax, [80, 100, 120, 160], path_cwd)
