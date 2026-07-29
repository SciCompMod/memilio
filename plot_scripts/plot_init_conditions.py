#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Anna Wendler
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
import h5py
import os
import numpy as np
import matplotlib.pyplot as plt


STYLE = {
    "groundtruth":    {"label": "Groundtruth",    "color": "#0072B2", "linestyle": "-",  "linewidth": 5., "alpha":0.3  },
    "detailed":       {"label": "Detailed long",  "color": "#E69F00", "linestyle": "--", "linewidth": 2., "alpha":1.},
    "detailed_short": {"label": "Detailed short", "color": "#009E73", "linestyle": "-.", "linewidth": 1.5, "alpha":1.},
    "simple":         {"label": "Simple",         "color": "#D55E00", "linestyle": ":",  "linewidth": 2., "alpha":1.},
}
# order matches input file lists
FILE_KEYS = [ "detailed", "detailed_short", "simple"]
FILE_KEYS_GROUNDTRUTH = ["groundtruth", "detailed", "detailed_short", "simple"]

SCATTERSIZE = 12
TITLE_FONTSIZE = 11
T0_COLOR = "gray"
T0_LABEL = r"$t_0$"


PANEL_WIDTH = 3.4
PANEL_HEIGHT = 4


def get_t0_from_dir_name(dir_name):
    t0_string = [x for x in dir_name.split("_") if "t0ide" in x  or "t0" in x]
    t0 = float(t0_string[0].split("=")[-1])
    return t0


def load_h5_total(filepath):
    """Load the single dataset's 'Total' array and 'Time' array from an h5 file."""
    with h5py.File(str(filepath) + ".h5", "r") as h5file:
        if len(list(h5file.keys())) > 1:
            print("File should contain one dataset.")
            return
        group_key = list(h5file.keys())[0]
        if len(list(h5file[group_key].keys())) > 3:
            print("Expected only one group.")
            return
        data = h5file[group_key]
        total = data["Total"][:, :]
        dates = data["Time"][:]
    return dates, total


def _style_axis(ax, title):
    ax.set_title(title, fontsize=TITLE_FONTSIZE, pad=12)
    ax.set_axisbelow(True)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))


def _add_shared_legend(fig, used_keys, include_t0=True):

    styles = [STYLE[k] for k in used_keys]
    handles = [plt.Line2D([0], [0], **style)
               for style in styles]
    if include_t0:
        handles.append(plt.Line2D(
            [0], [0], color=T0_COLOR, alpha=0.5, label=T0_LABEL))

    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.08),
               ncol=len(handles), fancybox=False, shadow=False, frameon=False)


def _save(fig, save_dir, filename):
    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        fig.savefig(os.path.join(save_dir, filename),
                    bbox_inches="tight", dpi=500)
    plt.close(fig)


def plot_compartments(files, fileending, save_dir="", zoom = False):
    """
    Plots simulation results (S, I, R) comparing groundtruth/detailed/simple.

    @param[in] files List of three files [groundtruth, detailed, simple] with
        compartment simulation results, in that order.
    @param[in] fileending Suffix for the saved plot filename.
    @param[in] save_dir Directory to save the plot in; not saved if empty.
    """
    compartment_names = {0: "Susceptible", 1: "Infected", 2: "Recovered"}
    num_plots = len(compartment_names)

    fig, axs = plt.subplots(1, num_plots, figsize=(PANEL_WIDTH * num_plots, PANEL_HEIGHT),
                            sharex="all", num="Compare files")

    plotted_keys = []
    plot_min = 1e7*np.ones(num_plots)
    plot_max = np.zeros(num_plots)
    all_dates = []
    all_totals = []

    for key, filepath in zip(FILE_KEYS_GROUNDTRUTH, files):
        dates, total = load_h5_total(filepath)
        all_dates.append(dates)
        all_totals.append(total)
        # if key == "groundtruth":
        #     # groundtruth_total = total
        #     continue
        plotted_keys.append(key)
        style = STYLE[key]
        for i in range(num_plots):
            _style_axis(axs[i], compartment_names[i])
            if np.min(total[:, i]) < plot_min[i]:
                plot_min[i] = np.min(total[:, i])
            if np.max(total[:, i]) > plot_max[i]:
                plot_max[i] = np.max(total[:, i])
            axs[i].plot(dates, total[:, i], **style)

    t0 = get_t0_from_dir_name(files[0])
    for i in range(num_plots):
        axs[i].vlines(t0, plot_min[i], plot_max[i],
                      color=T0_COLOR, alpha=0.5)

    if zoom:
        zoom_min = t0 - 7
        zoom_max = t0 + 7
        for i in range(num_plots):
            axs[i].set_xlim(zoom_min, zoom_max)
            
            # Calculate y-limits based on data within zoom range
            zoom_plot_min = 1e7
            zoom_plot_max = 0
            for dates, total in zip(all_dates, all_totals):
                mask = (dates >= zoom_min) & (dates <= zoom_max)
                if np.any(mask):
                    zoom_plot_min = min(zoom_plot_min, np.min(total[mask, i]))
                    zoom_plot_max = max(zoom_plot_max, np.max(total[mask, i]))
            
            # Add some padding to the y-limits
            y_padding = (zoom_plot_max - zoom_plot_min) * 0.1
            axs[i].set_ylim(zoom_plot_min - y_padding, zoom_plot_max + y_padding)

    _add_shared_legend(fig, plotted_keys)

    fig.supxlabel("Simulation time [days]", y=0.04)
    fig.supylabel("Number of individuals", y=0.54)
    plt.tight_layout()

    if zoom:
        _save(fig, save_dir, f"compare_compartments_{fileending}_zoom.png")
    else:
        _save(fig, save_dir, f"compare_compartments_{fileending}.png")


def plot_flow_S_to_I(files, fileending, save_dir=""):
    """
    Plots the Susceptible -> Infected flow, comparing groundtruth/detailed/simple.

    @param[in] files List of three files [groundtruth, detailed, simple].
    @param[in] fileending Suffix for the saved plot filename.
    @param[in] save_dir Directory to save the plot in; not saved if empty.
    """
    fig, ax = plt.subplots(1, 1, figsize=(
        1.5 * PANEL_WIDTH, PANEL_HEIGHT), num="Compare files")

    plotted_keys = []
    plot_min = 1e7
    plot_max = 0

    t0 = get_t0_from_dir_name(files[0])

    _style_axis(ax, "Susceptible to Infected flow")

    for key, filepath in zip(FILE_KEYS, files):
        dates, total = load_h5_total(filepath)
        if key == "detailed":
            tinit_detailed = dates[0]
        plotted_keys.append(key)
        style = STYLE[key]

        t0_index = np.where(np.isclose(dates, t0))[0][0]
        ax.scatter(dates[:t0_index+1], total[:, 0][:t0_index+1], color=style["color"],
                   marker="o",
                   s=SCATTERSIZE,
                   linewidths=0,
                   edgecolors="none")

        # ax.plot(dates[:t0_index+1], total[:, 0][:t0_index+1], color=style["color"],
        #         marker="o",
        #         linewidth=12)

    # tinit =
    ax.set_xlim(tinit_detailed-2, t0+2)
    # ax.vlines(t0, plot_min, plot_max,
    #           color=T0_COLOR, alpha=0.5)

    handles = [plt.Line2D([0], [0], color=STYLE[k]["color"],
                          linestyle="-", linewidth=LINEWIDTH,
                          label=STYLE[k]["label"])
               for k in plotted_keys]

    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.08),
               ncol=len(handles), fancybox=False, shadow=False, frameon=False)

    fig.supxlabel("Simulation time [days]",  y=0.04)
    fig.supylabel("Number of individuals", y=0.65)
    plt.tight_layout()

    _save(fig, save_dir, f"compare_flow_{fileending}.png")


def plot_infectionage_distribution(files, fileending, save_dir=""):
    """
    Plots the infected-per-infection-age distribution, comparing
    groundtruth/detailed/simple.

    @param[in] files List of three files [groundtruth, detailed, simple].
    @param[in] fileending Suffix for the saved plot filename.
    @param[in] save_dir Directory to save the plot in; not saved if empty.
    """
    fig, ax = plt.subplots(1, 1, figsize=(
        1.5*PANEL_WIDTH, PANEL_HEIGHT), num="Compare files")

    plotted_keys = []

    _style_axis(ax, r"Infected per infection age at $t_0$")

    for key, filepath in zip(FILE_KEYS, files):
        dates, total = load_h5_total(filepath)
        if key == "groundtruth":
            continue
        plotted_keys.append(key)
        style = STYLE[key]
        ax.scatter(dates, total, color=style["color"],
                   linestyle=style["linestyle"], marker="o",
                   s=SCATTERSIZE,
                   linewidths=0,
                   edgecolors="none")

    handles = [plt.Line2D([0], [0], color=STYLE[k]["color"],
                          linestyle="-", linewidth=LINEWIDTH,
                          label=STYLE[k]["label"])
               for k in plotted_keys]

    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.08),
               ncol=len(handles), fancybox=False, shadow=False, frameon=False)

    fig.supxlabel("Infection age [days]",  y=0.04)
    fig.supylabel("Number of individuals")
    plt.tight_layout()

    _save(fig, save_dir, f"infected_per_infection_age_{fileending}.png")


def plot_init_conditions(files_flows, files_infage_dist, fileending, save_dir=""):
    """
    Combines the Susceptible -> Infected flow plot and the infected-per-
    infection-age distribution plot into one figure with a shared legend.

    @param[in] files_flows List of flow files in the order
        [detailed, detailed_short, simple].
    @param[in] files_infage_dist List of infection-age distribution files in
        the order [detailed, detailed_short, simple].
    @param[in] fileending Suffix for the saved plot filename.
    @param[in] save_dir Directory to save the plot in; not saved if empty.
    """
    if len(files_flows) != len(FILE_KEYS):
        raise ValueError(
            f"Expected {len(FILE_KEYS)} flow files, got {len(files_flows)}."
        )

    if len(files_infage_dist) != len(FILE_KEYS):
        raise ValueError(
            "Expected "
            f"{len(FILE_KEYS)} infection-age distribution files, "
            f"got {len(files_infage_dist)}."
        )

    fig, axs = plt.subplots(
        1,
        2,
        figsize=(3 * PANEL_WIDTH, PANEL_HEIGHT),
        num="Compare initial conditions",
    )

    _style_axis(axs[0], "Susceptible to Infected flow")
    _style_axis(axs[1], r"Infected per infection age at $t_0$")

    plotted_keys = []
    t0 = get_t0_from_dir_name(files_flows[0])
    tinit_detailed = None

    # Plot the flow data.
    for key, filepath in zip(FILE_KEYS, files_flows):
        dates_flow, total_flow = load_h5_total(filepath)

        if key == "detailed":
            tinit_detailed = dates_flow[0]

        style = STYLE[key]
        plotted_keys.append(key)

        t0_indices = np.where(np.isclose(dates_flow, t0))[0]
        if len(t0_indices) == 0:
            raise ValueError(
                f"Could not find t0={t0} in flow file '{filepath}'."
            )

        t0_index = t0_indices[0]

        axs[0].scatter(
            dates_flow[:t0_index + 1],
            total_flow[:t0_index + 1, 0],
            color=style["color"],
            marker="o",
            s=SCATTERSIZE,
            linewidths=0,
            edgecolors="none",
            zorder=2,
        )

    # Plot the infection-age distributions.
    for key, filepath in zip(FILE_KEYS, files_infage_dist):
        dates_age, total_age = load_h5_total(filepath)
        style = STYLE[key]

        axs[1].scatter(
            dates_age,
            total_age,
            color=style["color"],
            marker="o",
            s=SCATTERSIZE,
            linewidths=0,
            edgecolors="none",
            zorder=2,
        )

    if tinit_detailed is not None:
        axs[0].set_xlim(tinit_detailed - 2, t0 + 2)

    axs[0].set_xlabel("Simulation time [days]")
    axs[1].set_xlabel("Infection age [days]")
    fig.supylabel("Number of individuals", y=0.57)

    styles = [STYLE[k] for k in plotted_keys]

    handles = [
        plt.Line2D(
            [0],
            [0],
            **style,
        )
        for style in styles
    ]

    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.02),
        ncol=len(handles),
        fancybox=False,
        shadow=False,
        frameon=False,
    )

    plt.tight_layout(rect=(0, 0.10, 1, 1))

    _save(fig, save_dir, f"init_conditions_{fileending}.png")


def subfolders_scandir(path):
    with os.scandir(path) as it:
        return [entry.name for entry in it if entry.is_dir()]


if __name__ == "__main__":

    root_dir = os.path.join(os.path.dirname(__file__), "../simulation_results")
    main_dir = "2026-07-29/compare_different_inits_exp" 

    relevant_dir = os.path.join(root_dir, main_dir)
    sub_dirs = subfolders_scandir(relevant_dir)

    gregory_order = 3
    ide_exponent = 2

    for sub_dir in sub_dirs:
        print(main_dir + "/" + sub_dir)

        result_dir = os.path.join(os.path.dirname(__file__),
                                  f"../simulation_results/{main_dir}/{sub_dir}")
        plot_dir = os.path.join(os.path.dirname(__file__),
                                f"../plots/{main_dir}/{sub_dir}/")

        files = os.listdir(result_dir)

        if f"detailed_dt=1e-{ide_exponent}_gregoryorder={gregory_order}.h5" not in files:
            continue

        # print(ide_exponent)

        base = f"dt=1e-{ide_exponent}"

        plot_compartments(
            [os.path.join(result_dir, f"groundtruth_{base}"),
             os.path.join(result_dir, f"detailed_{base}_gregoryorder={gregory_order}"),
             os.path.join(result_dir, f"detailed_short_{base}_gregoryorder={gregory_order}"),
             os.path.join(result_dir, f"simple_{base}_gregoryorder={gregory_order}")],
            fileending=f"dt=1e-{ide_exponent}", save_dir=plot_dir)

        plot_compartments(
                    [os.path.join(result_dir, f"groundtruth_{base}"),
                     os.path.join(result_dir, f"detailed_{base}_gregoryorder={gregory_order}"),
                     os.path.join(result_dir, f"detailed_short_{base}_gregoryorder={gregory_order}"),
                     os.path.join(result_dir, f"simple_{base}_gregoryorder={gregory_order}")],
                    fileending=f"dt=1e-{ide_exponent}", save_dir=plot_dir, zoom=True)

        # plot_flow_S_to_I(
        #     [os.path.join(result_dir, f"detailed_flows_{base}"),
        #      os.path.join(result_dir, f"detailed_short_flows_{base}"),
        #      os.path.join(result_dir, f"simple_flows_{base}")],
        #     fileending=f"dt=1e-{ide_exponent}", save_dir=plot_dir)

        # plot_infectionage_distribution(
        #     [os.path.join(
        #         result_dir, f"detailed_infectionagedistribution_{base}"),
        #         os.path.join(
        #             result_dir, f"detailed_short_infectionagedistribution_{base}"),
        #         os.path.join(result_dir, f"simple_infectionagedistribution_{base}")],
        #     fileending=f"dt=1e-{ide_exponent}", save_dir=plot_dir)

        plot_init_conditions(
            [
             os.path.join(result_dir, f"detailed_flows_{base}_gregoryorder={gregory_order}"),
             os.path.join(result_dir, f"detailed_short_flows_{base}_gregoryorder={gregory_order}"),
             os.path.join(result_dir, f"simple_flows_{base}_gregoryorder={gregory_order}")], 
            [
             os.path.join(
                 result_dir, f"detailed_infectionagedistribution_{base}_gregoryorder={gregory_order}"),
             os.path.join(
                 result_dir, f"detailed_short_infectionagedistribution_{base}_gregoryorder={gregory_order}"),
             os.path.join(result_dir, f"simple_infectionagedistribution_{base}_gregoryorder={gregory_order}")],
            fileending=f"dt=1e-{ide_exponent}",
            save_dir=plot_dir,
        )
