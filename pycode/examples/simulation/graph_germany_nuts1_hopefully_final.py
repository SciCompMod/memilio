#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Carlotta Gerstein
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
os.environ["KERAS_BACKEND"] = "tensorflow"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import pickle
from scipy.stats import truncnorm

from matplotlib.patches import Patch

import bayesflow as bf
import keras

import memilio.simulation as mio
import memilio.simulation.osecir as osecir
from memilio.simulation.osecir import Model, interpolate_simulation_result
from memilio.epidata import defaultDict as dd

import geopandas as gpd

name = "nuts1_hf"

region_ids = [region_id for region_id in dd.State.keys()]
inference_params = ['scaling_factor', 'damping_values', 't_E', 't_ISy', 't_ISev',
                    't_Cr', 'mu_CR', 'mu_IH', 'mu_HU', 'mu_UD', 'transmission_prob']
summary_vars = [f'fed_state{i}' for i in range(16)] + ["state"]

bounds = {
    'scaling_factor': (1.0, 10.0),
    't_E': (1.0, 5.2),
    't_ISy': (4.0, 10.0),
    't_ISev': (5.0, 10.0),
    't_Cr': (9.0, 17.0),
    'mu_CR': (0.0, 0.4),
    'mu_IH': (0.0, 0.2),
    'mu_HU': (0.0, 0.4),
    'mu_UD': (0.0, 0.4),
    'transmission_prob': (0.0, 0.2)
}
SPIKE_SCALE = 0.4
SLAB_SCALE = 0.2
DATE_TIME = datetime.date(year=2020, month=10, day=1)
NUM_DAMPING_POINTS = 3


def set_fontsize(base_fontsize=17):
    fontsize = base_fontsize
    plt.rcParams.update({
        'font.size': fontsize,
        'axes.titlesize': fontsize * 1,
        'axes.labelsize': fontsize,
        'xtick.labelsize': fontsize * 0.8,
        'ytick.labelsize': fontsize * 0.8,
        'legend.fontsize': fontsize * 0.8,
        'font.family': "Arial"
    })


plt.style.use('default')

dpi = 300

colors = {"Blue": "#155489",
          "Medium blue": "#64A7DD",
          "Light blue": "#B4DCF6",
          "Lilac blue": "#AECCFF",
          "Turquoise": "#76DCEC",
          "Light green": "#B6E6B1",
          "Medium green": "#54B48C",
          "Green": "#5D8A2B",
          "Teal": "#20A398",
          "Yellow": "#FBD263",
          "Orange": "#E89A63",
          "Rose": "#CF7768",
          "Red": "#A34427",
          "Purple": "#741194",
          "Grey": "#C0BFBF",
          "Dark grey": "#616060",
          "Light grey": "#F1F1F1"}

def plot_region_fit(
    data: np.ndarray,
    region: int,
    true_data=None,
    ax=None,
    label=None,
    color="red",
    only_80q=False
):
    if data.ndim != 3:
        raise ValueError("Array not of shape (samples, time_points, regions)")
    if true_data is not None:
        if true_data.shape != data.shape[1:]:
            raise ValueError("True data shape does not match data shape")
    n_samples, n_time, n_regions = data.shape
    if not (0 <= region < n_regions):
        raise IndexError

    x = np.arange(n_time)
    vals = data[:, :, region]  # (samples, time_points)

    qs_80 = np.quantile(vals, q=[0.1, 0.9], axis=0)
    qs_90 = np.quantile(vals, q=[0.05, 0.95], axis=0)
    qs_95 = np.quantile(vals, q=[0.025, 0.975], axis=0)

    med = np.median(vals, axis=0)

    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(
        x, med, lw=2, label=label or f"{dd.State[region_ids[region]]}", color=color)
    ax.fill_between(x, qs_90[0], qs_90[1], linewidth=0, alpha=0.3,
                    color=color)
    ax.fill_between(x, qs_95[0], qs_95[1], alpha=0.1,
                    color=color)

    if true_data is not None:
        true_vals = true_data[:, region]  # (time_points,)
        ax.plot(x, true_vals, lw=2, color="black", label="Reported data")

    ax.set_xlabel("Time")
    ax.set_ylabel("ICU cases [#]")
    ax.set_title(f"{dd.State[region_ids[region]]}")


def plot_aggregated_over_regions(
    data: np.ndarray,
    region_agg=np.sum,
    true_data=None,
    ax=None,
    label=None,
    color='red',
    only_80q=False
):
    if data.ndim != 3:
        raise ValueError("Array not of shape (samples, time_points, regions)")
    if true_data is not None:
        if true_data.shape != data.shape[1:]:
            raise ValueError("True data shape does not match data shape")

    # Aggregate over regions
    agg_over_regions = region_agg(data, axis=-1)  # (samples, time_points)

    qs_80 = np.quantile(agg_over_regions, q=[0.1, 0.9], axis=0)
    qs_90 = np.quantile(agg_over_regions, q=[0.05, 0.95], axis=0)
    qs_95 = np.quantile(agg_over_regions, q=[0.025, 0.975], axis=0)

    # Aggregate over samples
    agg_median = np.median(agg_over_regions, axis=0)        # (time_points, )

    x = np.arange(agg_median.shape[0])
    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(x, agg_median,
            label=label or "Aggregated simulation", color=color)
    # ax.fill_between(x, qs_80[0], qs_80[1], alpha=0.5,
    #                 color=color, label="80% CI")
    # if not only_80q:
    ax.fill_between(x, qs_90[0], qs_90[1], linewidth=0, alpha=0.4,
                    color=color)
    ax.fill_between(x, qs_95[0], qs_95[1], alpha=0.25,
                    color=color)
    if true_data is not None:
        true_vals = region_agg(true_data, axis=-1)  # (time_points,)
        ax.scatter(x, true_vals, color="black", label="Reported data", marker='x')

    # Update x-axis to show dates in YYYY-MM format
    start_date = datetime.date(2020, 10, 1)
    date_labels = [(start_date + datetime.timedelta(days=int(day))).strftime('%Y-%m-%d') for day in x]
    ax.set_xticks(x[::7])  # Set ticks every 7 days
    ax.set_xticklabels(date_labels[::7], rotation=45)

    ax.set_ylabel("ICU cases [#]")

def plot_icu_on_germany(simulations, synthetic, with_aug):
    med = np.median(simulations, axis=0)

    population = pd.read_json('data/Germany/pydata/county_current_population_states.json')
    values = med / population['Population'].to_numpy()[None, :] * 100000

    map_data = gpd.read_file(os.path.join(os.getcwd(), 'tools/vg2500_12-31.utm32s.shape/vg2500/VG2500_KRS.shp'))
    fedstate_data = gpd.read_file(os.path.join(os.getcwd(), 'tools/vg2500_12-31.utm32s.shape/vg2500/VG2500_LAN.shp'))

    fig, axes = plt.subplots(1, 2, figsize=(12, 8), layout="constrained")
    vmin = 0
    vmax = 15.5

    plot_map(values[0], map_data, fedstate_data, axes[0], "Median ICU", vmin, vmax)
    plot_map(values[-1], map_data, fedstate_data, axes[1], "Median ICU", vmin, vmax)

    plt.savefig(f"{name}/median_icu_germany_{name}{synthetic}{with_aug}.png", bbox_inches='tight', dpi=dpi)


def plot_map(values, map_data, fedstate_data, ax, label, vmin, vmax):
    map_data[label] = map_data['ARS'].map({f"{region_id:05d}": values[region_id // 1000 - 1]  for region_id in dd.County.keys()})
    # map_data[label] = map_data['ARS'].map(dict(zip([f"{region_id:02d}" for region_id in region_ids], values)))

    map_data['state_id'] = map_data['ARS'].astype(
        int).map({county: county // 1000 for county in dd.County.keys()}).fillna(0).astype(int)

    states = map_data.dissolve(by='state_id')

    map_data.plot(
        column=f"{label}",
        cmap='Reds',
        linewidth=0.,
        ax=ax,
        legend=False,
        vmin=vmin,
        vmax=vmax,
    )
    states.boundary.plot(ax=ax, edgecolor='darkgray', linewidth=0.5)

    ax.axis('off')

def plot_all_regions(simulations, divi_data, synthetic, with_aug):
    n_regions = simulations.shape[-1]
    fig, ax = plt.subplots(nrows=4, ncols=4, figsize=(25, 25), layout="constrained")
    ax = ax.flatten()
    for i in range(n_regions):
        plot_region_fit(
            simulations, region=i, true_data=divi_data, label="Median", ax=ax[i], color=colors["Red"]
        )
    plt.savefig(f'{name}/federal_states_{name}{synthetic}{with_aug}.png', dpi=dpi)
    plt.close()

def calibration_curves_per_region(
    data: np.ndarray,
    true_data: np.ndarray,
    levels=np.linspace(0.01, 0.99, 20),
    ax=None,
    max_regions=None,
    cmap=plt.cm.Blues,
    linewidth=1.5,
    legend=True,
    with_ideal=True,
):
    """
    Per-region calibration curves, each region in a different shade of a colormap.

    data: (samples, time, regions)
    true_data: (time, regions)
    max_regions: limit number of regions shown
    cmap: matplotlib colormap for line shades
    """
    if data.ndim != 3:
        raise ValueError("Array not of shape (samples, time_points, regions)")
    if true_data.shape != data.shape[1:]:
        raise ValueError("True data shape does not match data shape")

    n_samples, n_time, n_regions = data.shape
    if max_regions is None:
        R = n_regions
    else:
        R = min(max_regions, n_regions)

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4))

    colors = [cmap(i / (R + 1)) for i in range(1, R + 1)]

    x = np.asarray(levels)
    for r, col in zip(range(R), colors):
        emp = []
        for nominal in levels:
            q_low = (1.0 - nominal) / 2.0
            q_high = 1.0 - q_low
            lo = np.quantile(data[:, :, r], q_low, axis=0)
            hi = np.quantile(data[:, :, r], q_high, axis=0)
            hits = (true_data[:, r] >= lo) & (true_data[:, r] <= hi)
            emp.append(hits.mean())
        emp = np.asarray(emp)
        # , label=f"Region {r+1}")
        ax.plot(x, emp, lw=linewidth, color=col, alpha=0.5)

    if with_ideal:
        ideal_line = ax.plot([0, 1], [0, 1], linestyle="--",
                             lw=1.2, color="black", label="Ideal")[0]
    else:
        ideal_line = None

    if legend:
        # Custom legend: one patch for regions, one line for ideal
        region_patch = Patch(color=colors[-1], label="Regions")
        ax.legend(handles=[region_patch, ideal_line],
                  frameon=True, ncol=1)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Nominal level")
    ax.set_ylabel("Empirical coverage")
    ax.set_title("Calibration per region")
    return ax


def calibration_median_mad_over_regions(
    data: np.ndarray,
    true_data: np.ndarray,
    levels=np.linspace(0.01, 0.99, 20),
    ax=None,
    color="tab:blue",
    alpha=0.25,
    linewidth=2.0,
    with_ideal=True,
):
    """
    Compute per region empirical coverage at each nominal level,
    then summarize across regions with median and MAD.

    data: (samples, time, regions)
    true_data: (time, regions)
    """
    if data.ndim != 3:
        raise ValueError("Array not of shape (samples, time_points, regions)")
    if true_data.shape != data.shape[1:]:
        raise ValueError("True data shape does not match data shape")

    n_samples, n_time, n_regions = data.shape
    L = len(levels)
    per_region = np.empty((n_regions, L), dtype=float)

    # per level coverage per region
    for j, nominal in enumerate(levels):
        q_low = (1.0 - nominal) / 2.0
        q_high = 1.0 - q_low
        lo = np.quantile(data, q_low, axis=0)   # (time, regions)
        hi = np.quantile(data, q_high, axis=0)  # (time, regions)
        hits = (true_data >= lo) & (true_data <= hi)  # (time, regions)
        # mean over time for each region
        per_region[:, j] = hits.mean(axis=0)

    med = np.median(per_region, axis=0)                # (levels,)
    mad = np.median(np.abs(per_region - med[None, :]), axis=0)

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4))

    x = np.asarray(levels)
    ax.fill_between(x, med - mad, med + mad,
                    alpha=alpha, color=color, label=None)
    ax.plot(x, med, lw=linewidth, color=color, label="Median across regions")

    if with_ideal:
        ax.plot([0, 1], [0, 1], linestyle="--",
                lw=1.2, color="black", label="Ideal")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Nominal level")
    ax.set_ylabel("Empirical coverage")
    ax.set_title("Calibration median and MAD across regions")
    ax.legend()
    return ax, {"levels": x, "median": med, "mad": mad}


def plot_damping_values(damping_values, synthetic):

    med = np.median(damping_values, axis=0)
    mad = np.median(np.abs(damping_values - med), axis=0)

    # Extend for step plotting
    med_extended = np.hstack([med, med[:, -1][:, None]])
    mad_extended = np.hstack([mad, mad[:, -1][:, None]])

    # Plot damping values per region
    fig, axes = plt.subplots(4, 4, figsize=(15, 10), constrained_layout=True)
    axes = axes.flatten()
    x = np.arange(15, 61, 15)  # Time steps from 15 to 60

    for i, ax in enumerate(axes):
        if i < 16:
            ax.stairs(med[i], edges=x, lw=2, color='red', baseline=None)
            ax.fill_between(
                x, med_extended[i] - mad_extended[i], med_extended[i] + mad_extended[i],
                alpha=0.25, color='red', step='post'
            )
            ax.set_title(f"{dd.State[i+1]}")
            ax.set_xlabel("Time")
            ax.set_ylabel("Damping Value")
        else:
            ax.axis('off')  # Hide unused subplots

    plt.suptitle("Damping Values per Region")
    plt.savefig(f"{name}/damping_values{name}{synthetic}.png", dpi=dpi)

    # Combined plot for all regions
    fig, ax = plt.subplots(figsize=(10, 6))
    cmap = plt.cm.get_cmap("viridis", 16)  # Colormap with 16 distinct colors

    for i in range(16):
        ax.stairs(
            med[i], edges=x, lw=2, label=f"{dd.State[i+1]}",
            color=cmap(i), baseline=None
        )

    ax.set_title("Damping Values per Region (Combined)")
    ax.set_xlabel("Time")
    ax.set_ylabel("Damping Value")
    ax.legend(loc="upper right", ncol=2)
    plt.savefig(f"{name}/damping_values_combined{name}{synthetic}.png", dpi=dpi)


def best_fit_sim(results, true_data, ax=None):
    rmse_per_region = np.sqrt(np.mean((results - true_data[None, :, :])**2, axis=1))
    rmse_agg = np.sum(rmse_per_region, axis=-1)
    rmse_agg_idx = np.argmin(rmse_agg)

    agg_over_regions = np.sum(results, axis=-1)  # (samples, time_points)
    true_vals = np.sum(true_data, axis=-1)  # (time_points,)

    rmse = np.sqrt(np.mean((agg_over_regions - true_vals[None, :])**2, axis=1))
    best_fit_index = np.argmin(rmse)

    x = np.arange(agg_over_regions.shape[1])
    ax.plot(x, agg_over_regions[rmse_agg_idx], lw=2, label="First rmse then aggregated", color="red")
    ax.plot(x, agg_over_regions[best_fit_index], lw=2,
            label="First aggregated then rmse", color="blue", linestyle="--")
    ax.plot(x, true_vals, lw=2, color="black", label="Reported data")

    ax.set_xlabel("Time")
    ax.set_ylabel("ICU")
    ax.legend()

    return best_fit_index

class Simulation:
    """ """

    def __init__(self, data_dir, start_date, results_dir):
        self.num_groups = 1
        self.data_dir = data_dir
        self.start_date = start_date
        self.results_dir = results_dir
        if not os.path.exists(self.results_dir):
            os.makedirs(self.results_dir)

    def set_covid_parameters(self, model, t_E, t_ISy, t_ISev, t_Cr, mu_CR, mu_IH, mu_HU, mu_UD, transmission_prob):
        model.parameters.TimeExposed[mio.AgeGroup(0)] = t_E
        model.parameters.TimeInfectedNoSymptoms[mio.AgeGroup(0)] = 5.2 - t_E
        model.parameters.TimeInfectedSymptoms[mio.AgeGroup(0)] = t_ISy
        model.parameters.TimeInfectedSevere[mio.AgeGroup(0)] = t_ISev
        model.parameters.TimeInfectedCritical[mio.AgeGroup(0)] = t_Cr

        # probabilities
        model.parameters.TransmissionProbabilityOnContact[mio.AgeGroup(
            0)] = transmission_prob
        model.parameters.RelativeTransmissionNoSymptoms[mio.AgeGroup(0)] = 1

        model.parameters.RecoveredPerInfectedNoSymptoms[mio.AgeGroup(
            0)] = mu_CR
        model.parameters.SeverePerInfectedSymptoms[mio.AgeGroup(0)] = mu_IH
        model.parameters.CriticalPerSevere[mio.AgeGroup(0)] = mu_HU
        model.parameters.DeathsPerCritical[mio.AgeGroup(0)] = mu_UD

        # start day is set to the n-th day of the year
        model.parameters.StartDay = self.start_date.timetuple().tm_yday

        model.parameters.Seasonality = mio.UncertainValue(0.2)

    def set_contact_matrices(self, model):
        contact_matrices = mio.ContactMatrixGroup(1, self.num_groups)

        baseline = np.ones((self.num_groups, self.num_groups)) * 7.95
        minimum = np.zeros((self.num_groups, self.num_groups))
        contact_matrices[0] = mio.ContactMatrix(baseline, minimum)
        model.parameters.ContactPatterns.cont_freq_mat = contact_matrices

    def set_npis(self, params, end_date, damping_values):
        start_damping_1 = DATE_TIME + datetime.timedelta(days=15)
        start_damping_2 = DATE_TIME + datetime.timedelta(days=30)
        start_damping_3 = DATE_TIME + datetime.timedelta(days=45)

        if start_damping_1 < end_date:
            start_date = (start_damping_1 - self.start_date).days
            params.ContactPatterns.cont_freq_mat[0].add_damping(
                mio.Damping(np.r_[damping_values[0]], t=start_date))

        if start_damping_2 < end_date:
            start_date = (start_damping_2 - self.start_date).days
            params.ContactPatterns.cont_freq_mat[0].add_damping(
                mio.Damping(np.r_[damping_values[1]], t=start_date))

        if start_damping_3 < end_date:
            start_date = (start_damping_3 - self.start_date).days
            params.ContactPatterns.cont_freq_mat[0].add_damping(
                mio.Damping(np.r_[damping_values[2]], t=start_date))

    def get_graph(self, end_date, scaling_factor, t_E, t_ISy, t_ISev, t_Cr, mu_CR, mu_IH, mu_HU, mu_UD, transmission_prob):
        print("Initializing model...")
        model = Model(self.num_groups)
        self.set_covid_parameters(
            model, t_E, t_ISy, t_ISev, t_Cr, mu_CR, mu_IH, mu_HU, mu_UD, transmission_prob)
        self.set_contact_matrices(model)
        print("Model initialized.")

        graph = osecir.ModelGraph()

        scaling_factor_infected = [scaling_factor]
        scaling_factor_icu = 1.0

        data_dir_Germany = os.path.join(self.data_dir, "Germany")
        mobility_data_file = os.path.join(
            data_dir_Germany, "mobility", "commuter_mobility_2022_states.txt")
        pydata_dir = os.path.join(data_dir_Germany, "pydata")

        path_population_data = os.path.join(
            pydata_dir, "county_current_population_states.json")

        print("Setting nodes...")
        mio.osecir.set_nodes_states(
            model.parameters,
            mio.Date(self.start_date.year,
                     self.start_date.month, self.start_date.day),
            mio.Date(end_date.year,
                     end_date.month, end_date.day), pydata_dir,
            path_population_data, False, graph, scaling_factor_infected,
            scaling_factor_icu, 1.0, 0, False)

        print("Setting edges...")
        mio.osecir.set_edges(mobility_data_file, graph, 1)

        print("Graph created.")

        return graph

    def run(self, num_days_sim, scaling_factor, damping_values, t_E, t_ISy, t_ISev, t_Cr, mu_CR, mu_IH, mu_HU, mu_UD, transmission_prob, export_timeseries):
        mio.set_log_level(mio.LogLevel.Warning)
        end_date = self.start_date + datetime.timedelta(days=num_days_sim)

        graph = self.get_graph(end_date, scaling_factor, t_E, t_ISy, t_ISev,
                               t_Cr, mu_CR, mu_IH, mu_HU, mu_UD, transmission_prob)

        mobility_graph = osecir.MobilityGraph()
        for node_idx in range(graph.num_nodes):
            node = graph.get_node(node_idx)

            self.set_npis(
                node.property.parameters, 
                end_date, 
                damping_values[node_idx])
            mobility_graph.add_node(node.id, node.property)
        for edge_idx in range(graph.num_edges):
            mobility_graph.add_edge(
                graph.get_edge(edge_idx).start_node_idx,
                graph.get_edge(edge_idx).end_node_idx,
                graph.get_edge(edge_idx).property)
        mobility_sim = osecir.MobilitySimulation(mobility_graph, t0=0, dt=0.5)
        mobility_sim.advance(num_days_sim)

        if export_timeseries:
            graph = mobility_sim.graph
            node_ids = [graph.get_node(i).id for i in range(graph.num_nodes)]

            osecir.save_results(
                [osecir.interpolate_simulation_result(graph)], [[graph.get_node(node_indx).property.model
                    for node_indx in range(graph.num_nodes)]], node_ids, self.results_dir, True, False
            )

        results = {}
        for node_idx in range(mobility_sim.graph.num_nodes):
            results[f'fed_state{node_idx}'] = osecir.interpolate_simulation_result(
                mobility_sim.graph.get_node(node_idx).property.result)

        return results


def run_germany_nuts1_simulation(scaling_factor, damping_values, t_E, t_ISy, t_ISev, t_Cr, mu_CR, mu_IH, mu_HU, mu_UD, transmission_prob, export_timeseries=False):
    mio.set_log_level(mio.LogLevel.Warning)
    file_path = os.path.dirname(os.path.abspath(__file__))

    sim = Simulation(
        data_dir=os.path.join(file_path, "../../../data"),
        start_date=DATE_TIME,
        results_dir=os.path.join(file_path, "../../../results_osecir"))
    num_days_sim = 60

    results = sim.run(num_days_sim, scaling_factor, damping_values, t_E, t_ISy,
                      t_ISev, t_Cr, mu_CR, mu_IH, mu_HU, mu_UD, transmission_prob, export_timeseries)

    return results


def prior():
    damping_values = np.zeros((NUM_DAMPING_POINTS, 16))
    for i in range(NUM_DAMPING_POINTS):
        mean = np.random.uniform(0, 1)
        scale = 0.1
        a, b = (0 - mean) / scale, (1 - mean) / scale
        damping_values[i] = truncnorm.rvs(
            a=a, b=b, loc=mean, scale=scale, size=16
        )
    return {
        'scaling_factor': np.random.uniform(*bounds['scaling_factor']),
        'damping_values': np.transpose(damping_values),
        't_E': np.random.uniform(*bounds['t_E']),
        't_ISy': np.random.uniform(*bounds['t_ISy']),
        't_ISev': np.random.uniform(*bounds['t_ISev']),
        't_Cr': np.random.uniform(*bounds['t_Cr']),
        'mu_CR': np.random.uniform(*bounds['mu_CR']),
        'mu_IH': np.random.uniform(*bounds['mu_IH']),
        'mu_HU': np.random.uniform(*bounds['mu_HU']),
        'mu_UD': np.random.uniform(*bounds['mu_UD']),
        'transmission_prob': np.random.uniform(*bounds['transmission_prob'])
    }


def load_divi_data():
    file_path = os.path.dirname(os.path.abspath(__file__))
    divi_path = os.path.join(file_path, "../../../data/Germany/pydata")

    data = pd.read_json(os.path.join(divi_path, "state_divi_ma7.json"))
    data = data[data['Date'] >= np.datetime64(DATE_TIME)]
    data = data[data['Date'] <= np.datetime64(DATE_TIME + datetime.timedelta(days=60))]
    data = data.sort_values(by=['ID_State', 'Date'])
    divi_data = data.pivot(index='Date', columns='ID_State', values='ICU')
    divi_dict = {}
    for i in range(len(region_ids)):
        divi_dict[f"fed_state{i}"] = divi_data[region_ids[i]].to_numpy()[None, :, None]
    return divi_dict



def extract_observables(simulation_results, observable_index=7):
    for key in simulation_results.keys():
        if key not in inference_params:
            simulation_results[key] = simulation_results[key][:, :, observable_index][..., np.newaxis]
    return simulation_results


def create_train_data(filename, number_samples=1000):

    simulator = bf.simulators.make_simulator(
        [prior, run_germany_nuts1_simulation]
    )
    trainings_data = simulator.sample(number_samples)
    trainings_data = extract_observables(trainings_data)
    with open(filename, 'wb') as f:
        pickle.dump(trainings_data, f, pickle.HIGHEST_PROTOCOL)


def load_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)

def is_state_key(k: str) -> bool:
    return 'fed_state' in k

def apply_aug(d: dict, aug) -> dict:
    return {k: np.clip(aug(v), 0, None) if is_state_key(k) else v for k, v in d.items()}

def concat_dicts(base: dict, new: dict) -> dict:
    missing = set(base) - set(new)
    if missing:
        raise KeyError(f"new dict missing keys: {sorted(missing)}")
    for k in base:
        base[k] = np.concatenate([base[k], new[k]])
    return base


def combine_results(dict_list):
    combined = {}
    for d in dict_list:
        combined = concat_dicts(combined, d) if combined else d
    return combined


def skip_2weeks(d: dict) -> dict:
    return {k: v[:, 14:, :] if is_state_key(k) else v for k, v in d.items()}

def aggregate_states(d: dict) -> None:
    d["state"] = np.sum([d[f"fed_state{r}"] for r in range(16)], axis=0)


def get_workflow():

    simulator = bf.make_simulator(
        [prior, run_germany_nuts1_simulation]
    )
    adapter = (
        bf.Adapter()
        .to_array()
        .convert_dtype("float64", "float32")
        .constrain("scaling_factor", lower=bounds["scaling_factor"][0], upper=bounds["scaling_factor"][1])
        .constrain("damping_values", lower=0.0, upper=1.0)
        .constrain("t_E", lower=bounds["t_E"][0], upper=bounds["t_E"][1])
        .constrain("t_ISy", lower=bounds["t_ISy"][0], upper=bounds["t_ISy"][1])
        .constrain("t_ISev", lower=bounds["t_ISev"][0], upper=bounds["t_ISev"][1])
        .constrain("t_Cr", lower=bounds["t_Cr"][0], upper=bounds["t_Cr"][1])
        .constrain("mu_CR", lower=bounds["mu_CR"][0], upper=bounds["mu_CR"][1])
        .constrain("mu_IH", lower=bounds["mu_IH"][0], upper=bounds["mu_IH"][1])
        .constrain("mu_HU", lower=bounds["mu_HU"][0], upper=bounds["mu_HU"][1])
        .constrain("mu_UD", lower=bounds["mu_UD"][0], upper=bounds["mu_UD"][1])
        .constrain("transmission_prob", lower=bounds["transmission_prob"][0], upper=bounds["transmission_prob"][1])
        .concatenate(
            ["scaling_factor", "damping_values", "t_E", "t_ISy", "t_ISev", "t_Cr",
             "mu_CR", "mu_IH", "mu_HU", "mu_UD", "transmission_prob"],
            into="inference_variables",
            axis=-1
        )
        .concatenate(summary_vars, into="summary_variables", axis=-1)
    )

    summary_network = bf.networks.TimeSeriesNetwork(
        summary_dim=(len(bounds)+16*NUM_DAMPING_POINTS)*2, dropout=0.1
    )
    inference_network = bf.networks.FlowMatching(subnet_kwargs={'widths': (512, 512, 512, 512, 512)})

    # aug = bf.augmentations.NNPE(spike_scale=SPIKE_SCALE, slab_scale=SLAB_SCALE, per_dimension=False)
    workflow = bf.BasicWorkflow(
        simulator=simulator,
        adapter=adapter,
        summary_network=summary_network,
        inference_network=inference_network,
        standardize='all'
        # augmentations={f'fed_state{i}': aug for i in range(len(region_ids))}
    )

    return workflow


def run_training(num_training_files=20):
    train_template = name+"/trainings_data{i}_"+name+".pickle"
    val_path = f"{name}/validation_data_{name}.pickle"

    aug = bf.augmentations.NNPE(
        spike_scale=SPIKE_SCALE, slab_scale=SLAB_SCALE, per_dimension=False
    )

    # training data
    train_files = [train_template.format(i=i) for i in range(1, 1+num_training_files)]
    trainings_data = None
    for p in train_files:
        d = load_pickle(p)
        d = apply_aug(d, aug=aug)  # only on region keys
        d['damping_values'] = d['damping_values'].reshape((d['damping_values'].shape[0], -1))
        if trainings_data is None:
            trainings_data = d
        else:
            trainings_data = concat_dicts(trainings_data, d)
    aggregate_states(trainings_data)

    # validation data
    validation_data = apply_aug(load_pickle(val_path), aug=aug)
    aggregate_states(validation_data)
    validation_data['damping_values'] = validation_data['damping_values'].reshape((validation_data['damping_values'].shape[0], -1))

    # check data
    workflow = get_workflow()
    print("summary_variables shape:", workflow.adapter(trainings_data)["summary_variables"].shape)
    print("inference_variables shape:", workflow.adapter(trainings_data)["inference_variables"].shape)

    history = workflow.fit_offline(
        data=trainings_data, epochs=300, batch_size=64, validation_data=validation_data
    )

    workflow.approximator.save(
        filepath=os.path.join(f"{name}/model_{name}.keras")
    )

    plots = workflow.plot_default_diagnostics(
        test_data=validation_data, calibration_ecdf_kwargs={'difference': True, 'stacked': True}
    )
    plots['losses'].savefig(f'{name}/losses_{name}.png', dpi=dpi)
    plots['recovery'].savefig(f'{name}/recovery_{name}.png', dpi=dpi)
    plots['calibration_ecdf'].savefig(f'{name}/calibration_ecdf_{name}.png', dpi=dpi)
    #plots['z_score_contraction'].savefig(f'{name}/z_score_contraction_{name}.png', dpi=dpi)


def run_inference(num_samples=1000, on_synthetic_data=False):
    val_path = f"{name}/validation_data_{name}.pickle"
    synthetic = "_synthetic" if on_synthetic_data else ""

    aug = bf.augmentations.NNPE(
        spike_scale=SPIKE_SCALE, slab_scale=SLAB_SCALE, per_dimension=False
    )

    # validation data
    validation_data = load_pickle(val_path)  # synthetic data

    if on_synthetic_data:
        # validation data
        validation_data = apply_aug(validation_data, aug=aug)
        validation_data['damping_values'] = validation_data['damping_values'].reshape((validation_data['damping_values'].shape[0], -1))
        aggregate_states(validation_data)
        validation_data = {k: v[0][np.newaxis, ...] if is_state_key(k) else v for k, v in validation_data.items()}  # only one dataset
        validation_data = {k: v[0][np.newaxis, ...] if k== "state" else v for k, v in validation_data.items()}  # only one dataset
        divi_dict = validation_data

        divi_data = np.concatenate(
            [divi_dict[f'fed_state{i}'] for i in range(len(region_ids))], axis=-1
        )[0]  # only one dataset
    else:
        divi_dict = load_divi_data()
        aggregate_states(divi_dict)
        divi_data = np.concatenate(
            [divi_dict[f'fed_state{i}'] for i in range(len(region_ids))], axis=-1
        )[0]

    with open(f'{name}/divi_data_{name}{synthetic}.pickle', 'wb') as f:
        pickle.dump(divi_data, f, pickle.HIGHEST_PROTOCOL)

    workflow = get_workflow()
    workflow.approximator = keras.models.load_model(
        filepath=os.path.join(f"{name}/model_{name}.keras")
    )

    if os.path.exists(f'{name}/sims_{name}{synthetic}_with_aug.pickle') and os.path.exists(f'{name}/sims_{name}{synthetic}.pickle') and os.path.exists(f'{name}/samples_{name}{synthetic}.pickle'):
        simulations = load_pickle(f'{name}/sims_{name}{synthetic}.pickle')
        simulations_aug = load_pickle(
            f'{name}/sims_{name}{synthetic}_with_aug.pickle')
        samples = load_pickle(f'{name}/samples_{name}{synthetic}.pickle')
        print("loaded simulations from file")
    else:
        samples = workflow.sample(conditions=divi_dict, num_samples=num_samples)
        with open(f'{name}/samples_{name}{synthetic}.pickle', 'wb') as f:
            pickle.dump(samples, f, pickle.HIGHEST_PROTOCOL)
        samples['damping_values'] = samples['damping_values'].reshape((samples['damping_values'].shape[0], num_samples, 16, NUM_DAMPING_POINTS))
        results = []
        for i in range(num_samples):  # we only have one dataset for inference here
            result = run_germany_nuts1_simulation(
                scaling_factor=samples['scaling_factor'][0, i],
                damping_values=samples['damping_values'][0, i],
                t_E=samples['t_E'][0, i], t_ISy=samples['t_ISy'][0, i],
                t_ISev=samples['t_ISev'][0, i], t_Cr=samples['t_Cr'][0, i],
                mu_CR=samples['mu_CR'][0, i], mu_IH=samples['mu_IH'][0, i],
                mu_HU=samples['mu_HU'][0, i], mu_UD=samples['mu_UD'][0, i],
                transmission_prob=samples['transmission_prob'][0, i]
            )
            for key in result.keys():
                result[key] = np.array(result[key])[None, ...]  # add sample axis
            results.append(result)
        results = combine_results(results)
        results = extract_observables(results)
        results_aug = apply_aug(results, aug=aug)

        # get sims in shape (samples, time, regions)
        simulations = np.zeros((num_samples, divi_data.shape[0], divi_data.shape[1]))
        simulations_aug = np.zeros((num_samples, divi_data.shape[0], divi_data.shape[1]))
        for i in range(num_samples):
            simulations[i] = np.concatenate([results[key][i] for key in results.keys()], axis=-1)
            simulations_aug[i] = np.concatenate([results_aug[key][i] for key in results.keys()], axis=-1)

        fig, ax = plt.subplots(1, 1, figsize=(8, 6), constrained_layout=True)
        sim_idx = best_fit_sim(simulations_aug, true_data=divi_data, ax=ax)
        plt.savefig(f'{name}/best_sims_{name}{synthetic}.png', dpi=dpi)
        plt.close()

        best_sim = run_germany_nuts1_simulation(
            scaling_factor=samples['scaling_factor'][0, sim_idx],
            damping_values=samples['damping_values'][0, sim_idx],
            t_E=samples['t_E'][0, sim_idx], t_ISy=samples['t_ISy'][0, sim_idx],
            t_ISev=samples['t_ISev'][0, sim_idx], t_Cr=samples['t_Cr'][0, sim_idx],
            mu_CR=samples['mu_CR'][0, sim_idx], mu_IH=samples['mu_IH'][0, sim_idx],
            mu_HU=samples['mu_HU'][0, sim_idx], mu_UD=samples['mu_UD'][0, sim_idx],
            transmission_prob=samples['transmission_prob'][0, sim_idx], export_timeseries = True
        )
        
        best_samples = {}
        for key in inference_params:
            best_samples[key] = samples[key][0, sim_idx]
        with open(f'{name}/best_samples_{name}.pickle', 'wb') as f:
            pickle.dump(best_samples, f, pickle.HIGHEST_PROTOCOL)

        # save sims
        with open(f'{name}/sims_{name}{synthetic}.pickle', 'wb') as f:
            pickle.dump(simulations, f, pickle.HIGHEST_PROTOCOL)
        with open(f'{name}/sims_{name}{synthetic}_with_aug.pickle', 'wb') as f:
            pickle.dump(simulations_aug, f, pickle.HIGHEST_PROTOCOL)

        samples['damping_values'] = samples['damping_values'].reshape(
        (samples['damping_values'].shape[0], samples['damping_values'].shape[1], -1))
        validation_data['damping_values'] = validation_data['damping_values'].reshape(
        (validation_data['damping_values'].shape[0], -1))

        plot = bf.diagnostics.pairs_posterior(
            samples, priors=validation_data, dataset_id=0)
        plot.savefig(f'{name}/pairs_posterior_{name}{synthetic}.png', dpi=dpi)

    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    # Plot with augmentation
    plot_aggregated_over_regions(
        simulations_aug, true_data=divi_data, label="Aggregated Simulation", ax=ax, color=colors["Red"]
    )
    lines, labels = ax.get_legend_handles_labels()
    fig.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.135, 0.99), ncol=1)
    plt.savefig(f'{name}/test_region_aggregated_{name}{synthetic}.png', dpi=dpi)
    plt.close()

    fig, axis = plt.subplots(1, 2, figsize=(10, 4), sharex=True, layout="constrained")
    ax = calibration_curves_per_region(simulations, divi_data, ax=axis[0])
    ax, stats = calibration_median_mad_over_regions(simulations, divi_data, ax=axis[1])
    plt.savefig(f'{name}/calibration_per_region_{name}{synthetic}.png', dpi=dpi)
    plt.close()
    fig, axis = plt.subplots(1, 2, figsize=(10, 4), sharex=True, layout="constrained")
    ax = calibration_curves_per_region(simulations_aug, divi_data, ax=axis[0])
    ax, stats = calibration_median_mad_over_regions(simulations_aug, divi_data, ax=axis[1])
    plt.savefig(f'{name}/calibration_per_region_{name}{synthetic}_with_aug.png', dpi=dpi)
    plt.close()

    plot_icu_on_germany(simulations, synthetic, with_aug="")
    plot_icu_on_germany(simulations_aug, synthetic, with_aug="_with_aug")

    simulation_agg = np.sum(simulations, axis=-1, keepdims=True)  # sum over regions
    simulation_aug_agg = np.sum(simulations_aug, axis=-1, keepdims=True)

    rmse = bf.diagnostics.metrics.root_mean_squared_error(np.swapaxes(simulation_agg, 0,1), np.sum(divi_data, axis=-1, keepdims=True), normalize=False)
    rmse_aug = bf.diagnostics.metrics.root_mean_squared_error(np.swapaxes(simulation_aug_agg, 0,1), np.sum(divi_data, axis=-1, keepdims=True), normalize=False)
    print("Mean RMSE over regions:", rmse["values"].mean())
    print("Mean RMSE over regions (with aug):", rmse_aug["values"].mean())

    cal_error = bf.diagnostics.metrics.calibration_error(np.swapaxes(simulation_agg, 0,1), np.sum(divi_data, axis=-1, keepdims=True))
    cal_error_aug = bf.diagnostics.metrics.calibration_error(np.swapaxes(simulation_aug_agg, 0,1), np.sum(divi_data, axis=-1, keepdims=True))
    print("Mean Calibration Error over regions:", cal_error["values"].mean())
    print("Mean Calibration Error over regions (with aug):", cal_error_aug["values"].mean())


if __name__ == "__main__":

    set_fontsize()

    if not os.path.exists(name):
        os.makedirs(name)
    # create_train_data(filename=f'{name}/validation_data_{name}.pickle', number_samples=100)
    # run_training(num_training_files=20)
    run_inference(on_synthetic_data=True, num_samples=100)
    run_inference(on_synthetic_data=False)