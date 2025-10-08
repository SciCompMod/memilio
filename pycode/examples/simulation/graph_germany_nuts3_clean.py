#######################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Henrik Zunker
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
from memilio.simulation.osecir import Model, Simulation, interpolate_simulation_result
from memilio.epidata import defaultDict as dd



excluded_ids = [11001, 11002, 11003, 11004, 11005, 11006,
                11007, 11008, 11009, 11010, 11011, 11012, 16056]
no_icu_ids = [7338, 9374, 9473, 9573]
region_ids = [region_id for region_id in dd.County.keys()
              if region_id not in excluded_ids]

inference_params = ['damping_values', 't_E', 't_ISy', 't_ISev',
                    't_Cr', 'mu_CR', 'mu_IH', 'mu_HU', 'mu_UD', 'transmission_prob']
summary_vars = ['state'] + [f'fed_state{i}' for i in range(16)] + [
    f'region{i}' for i in range(len(region_ids)) if region_ids[i] not in no_icu_ids]

bounds = {
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
SPIKE_SCALE = 0.85
SLAB_SCALE = 0.4
DATE_TIME = datetime.date(year=2020, month=10, day=1)

# %%
def plot_region_median_mad(
    data: np.ndarray,
    region: int,
    true_data=None,
    ax=None,
    label=None,
    color="red"
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
    med = np.median(vals, axis=0)
    mad = np.median(np.abs(vals - med), axis=0)

    if ax is None:
        fig, ax = plt.subplots()

    line, = ax.plot(
        x, med, lw=2, label=label or f"Region {region}", color=color)
    band = ax.fill_between(x, med - mad, med + mad, alpha=0.25, color=color)
    if true_data is not None:
        true_vals = true_data[:, region]  # (time_points,)
        ax.plot(x, true_vals, lw=2, color="black", label="True data")

    ax.set_xlabel("Time", fontsize=12)
    ax.set_ylabel("ICU", fontsize=12)
    ax.set_title(f"Region {region}", fontsize=12)
    if label is not None:
        ax.legend(fontsize=11, loc="upper right")
    return line, band


def plot_aggregated_over_regions(
    data: np.ndarray,
    region_agg=np.sum,
    true_data=None,
    ax=None,
    label=None,
    color='red'
):
    if data.ndim != 3:
        raise ValueError("Array not of shape (samples, time_points, regions)")
    if true_data is not None:
        if true_data.shape != data.shape[1:]:
            raise ValueError("True data shape does not match data shape")

    # Aggregate over regions
    agg_over_regions = region_agg(data, axis=-1)  # (samples, time_points)

    # Aggregate over samples
    agg_median = np.median(agg_over_regions, axis=0)        # (time_points, )
    agg_mad = np.median(
        np.abs(agg_over_regions - agg_median[None]),
        axis=0
    )

    x = np.arange(agg_median.shape[0])
    if ax is None:
        fig, ax = plt.subplots()

    line, = ax.plot(x, agg_median, lw=2,
                    label=label or "Aggregated over regions", color=color)
    band = ax.fill_between(
        x,
        agg_median - agg_mad,
        agg_median + agg_mad,
        alpha=0.25,
        color=color
    )
    if true_data is not None:
        true_vals = region_agg(true_data, axis=-1)  # (time_points,)
        ax.plot(x, true_vals, lw=2, color="black", label="True data")

    ax.set_xlabel("Time", fontsize=12)
    ax.set_ylabel("ICU", fontsize=12)
    if label is not None:
        ax.legend(fontsize=11)
    return line, band

# %%


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
                  frameon=True, ncol=1, fontsize=12)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Nominal level", fontsize=12)
    ax.set_ylabel("Empirical coverage", fontsize=12)
    ax.set_title("Calibration per region", fontsize=12)
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
    ax.set_xlabel("Nominal level", fontsize=12)
    ax.set_ylabel("Empirical coverage", fontsize=12)
    ax.set_title("Calibration median and MAD across regions", fontsize=12)
    ax.legend(fontsize=12)
    return ax, {"levels": x, "median": med, "mad": mad}


class Simulation:  # todo: correct class?
    """ """

    def __init__(self, data_dir, start_date, results_dir):
        self.num_groups = 1
        self.data_dir = data_dir
        self.start_date = start_date
        self.results_dir = results_dir
        if not os.path.exists(self.results_dir):
            os.makedirs(self.results_dir)

    def set_covid_parameters(self, model, t_E, t_ISy, t_ISev, t_Cr, mu_CR, mu_IH, mu_HU, mu_UD, transmission_prob):
        """

        :param model: 

        """
        model.parameters.TimeExposed[mio.AgeGroup(0)] = t_E
        model.parameters.TimeInfectedNoSymptoms[mio.AgeGroup(0)] = 5.2 - t_E  # todo: correct?
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
        """

        :param model: 

        """
        contact_matrices = mio.ContactMatrixGroup(1, self.num_groups)

        baseline = np.ones((self.num_groups, self.num_groups)) * 7.95
        minimum = np.zeros((self.num_groups, self.num_groups))
        contact_matrices[0] = mio.ContactMatrix(baseline, minimum)
        model.parameters.ContactPatterns.cont_freq_mat = contact_matrices

    def set_npis(self, params, end_date, damping_value):
        """

        :param params: 
        :param end_date: 

        """

        start_damping = DATE_TIME + datetime.timedelta(days=7)

        if start_damping < end_date:
            start_date = (start_damping - self.start_date).days
            params.ContactPatterns.cont_freq_mat[0].add_damping(
                mio.Damping(np.r_[damping_value], t=start_date))

    def get_graph(self, end_date, t_E, t_ISy, t_ISev, t_Cr, mu_CR, mu_IH, mu_HU, mu_UD, transmission_prob):
        """

        :param end_date: 

        """
        print("Initializing model...")
        model = Model(self.num_groups)
        self.set_covid_parameters(
            model, t_E, t_ISy, t_ISev, t_Cr, mu_CR, mu_IH, mu_HU, mu_UD, transmission_prob)
        self.set_contact_matrices(model)
        print("Model initialized.")

        graph = osecir.ModelGraph()

        scaling_factor_infected = [2.5]
        scaling_factor_icu = 1.0

        data_dir_Germany = os.path.join(self.data_dir, "Germany")
        mobility_data_file = os.path.join(
            data_dir_Germany, "mobility", "commuter_mobility_2022.txt")
        pydata_dir = os.path.join(data_dir_Germany, "pydata")

        path_population_data = os.path.join(pydata_dir,
                                            "county_current_population.json")

        print("Setting nodes...")
        mio.osecir.set_nodes(
            model.parameters,
            mio.Date(self.start_date.year,
                     self.start_date.month, self.start_date.day),
            mio.Date(end_date.year,
                     end_date.month, end_date.day), pydata_dir,
            path_population_data, True, graph, scaling_factor_infected,
            scaling_factor_icu, 0, 0, False)

        print("Setting edges...")
        mio.osecir.set_edges(mobility_data_file, graph, 1)

        print("Graph created.")

        return graph

    def run(self, num_days_sim, damping_values, t_E, t_ISy, t_ISev, t_Cr, mu_CR, mu_IH, mu_HU, mu_UD, transmission_prob, save_graph=True):
        """

        :param num_days_sim: 
        :param num_runs:  (Default value = 10)
        :param save_graph:  (Default value = True)
        :param create_gif:  (Default value = True)

        """
        mio.set_log_level(mio.LogLevel.Warning)
        end_date = self.start_date + datetime.timedelta(days=num_days_sim)

        graph = self.get_graph(end_date, t_E, t_ISy, t_ISev,
                               t_Cr, mu_CR, mu_IH, mu_HU, mu_UD, transmission_prob)

        mobility_graph = osecir.MobilityGraph()
        for node_idx in range(graph.num_nodes):
            node = graph.get_node(node_idx)
            self.set_npis(node.property.parameters, end_date,
                          damping_values[node.id // 1000 - 1]) 
            mobility_graph.add_node(node.id, node.property)
        for edge_idx in range(graph.num_edges):
            mobility_graph.add_edge(
                graph.get_edge(edge_idx).start_node_idx,
                graph.get_edge(edge_idx).end_node_idx,
                graph.get_edge(edge_idx).property)
        mobility_sim = osecir.MobilitySimulation(mobility_graph, t0=0, dt=0.5)
        mobility_sim.advance(num_days_sim)

        results = {}
        for node_idx in range(mobility_sim.graph.num_nodes):
            node = mobility_sim.graph.get_node(node_idx)
            if node.id in no_icu_ids:
                results[f'no_icu_region{node_idx}'] = osecir.interpolate_simulation_result(
                    node.property.result)
            else:
                results[f'region{node_idx}'] = osecir.interpolate_simulation_result(
                    node.property.result)

        return results


def run_germany_nuts3_simulation(damping_values, t_E, t_ISy, t_ISev, t_Cr, mu_CR, mu_IH, mu_HU, mu_UD, transmission_prob):
    mio.set_log_level(mio.LogLevel.Warning)
    file_path = os.path.dirname(os.path.abspath(__file__))

    sim = Simulation(
        data_dir=os.path.join(file_path, "../../../data"),
        start_date=DATE_TIME,
        results_dir=os.path.join(file_path, "../../../results_osecir"))
    num_days_sim = 60

    results = sim.run(num_days_sim, damping_values, t_E, t_ISy,
                      t_ISev, t_Cr, mu_CR, mu_IH, mu_HU, mu_UD, transmission_prob)

    return results

def prior():
    mean = np.random.uniform(0, 1)
    scale = 0.1
    a, b = (0 - mean) / scale, (1 - mean) / scale
    damping_values = truncnorm.rvs(
        a=a, b=b, loc=mean, scale=scale, size=16
    )
    return {
        'damping_values': damping_values,
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

    data = pd.read_json(os.path.join(divi_path, "county_divi_ma7.json"))
    data = data[data['Date'] >= np.datetime64(DATE_TIME)]
    data = data[data['Date'] <= np.datetime64(DATE_TIME + datetime.timedelta(days=60))]
    data = data.sort_values(by=['ID_County', 'Date'])
    divi_data = data.pivot(index='Date', columns='ID_County', values='ICU')
    divi_dict = {}
    for i, region_id in enumerate(region_ids):
        if region_id not in no_icu_ids:
            divi_dict[f"region{i}"] = divi_data[region_id].to_numpy()[None, :, None]
        else:
            divi_dict[f"no_icu_region{i}"] = np.zeros((1, divi_data.shape[0], 1))
    return divi_dict


def extract_observables(simulation_results, observable_index=7):
    for key in simulation_results.keys():
        if key not in inference_params:
            simulation_results[key] = simulation_results[key][:, :, observable_index][..., np.newaxis]
    return simulation_results


def create_train_data(filename, number_samples=1000):

    simulator = bf.simulators.make_simulator(
        [prior, run_germany_nuts3_simulation]
    )
    trainings_data = simulator.sample(number_samples)
    trainings_data = extract_observables(trainings_data)
    with open(filename, 'wb') as f:
        pickle.dump(trainings_data, f, pickle.HIGHEST_PROTOCOL)


def load_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)

def is_region_key(k: str) -> bool:
    return 'region' in k

def apply_aug(d: dict, aug) -> dict:
    return {k: np.clip(aug(v), 0, None) if is_region_key(k) else v for k, v in d.items()}

def concat_dicts(base: dict, new: dict) -> dict:
    missing = set(base) - set(new)
    if missing:
        raise KeyError(f"new dict missing keys: {sorted(missing)}")
    for k in base:
        base[k] = np.concatenate([base[k], new[k]])
    return base

def region_keys_sorted(d: dict):
    def idx(k):
        # handles "regionN" and "no_icu_regionN"
        return int(k.split("region")[-1])
    return sorted([k for k in d if is_region_key(k)], key=idx)

def aggregate_states(d: dict) -> None:
    n_regions = len(region_ids)
    # per state
    for state in range(16):
        idxs = [
            r for r in range(n_regions)
            if region_ids[r] // 1000 == state + 1
        ]
        d[f"fed_state{state}"] = np.sum([d[f"region{r}"] if region_ids[r] not in no_icu_ids else d[f"no_icu_region{r}"] for r in idxs], axis=0)
    # all allowed regions
    d["state"] = np.sum([d[f"fed_state{r}"] for r in range(16)], axis=0)

def combine_results(dict_list):
    combined = {}
    for d in dict_list:
        combined = concat_dicts(combined, d) if combined else d
    return combined

def skip_2weeks(d:dict) -> dict:
    return {k: v[:, 14:, :] if is_region_key(k) else v for k, v in d.items()}

def get_workflow():

    simulator = bf.make_simulator(
        [prior, run_germany_nuts3_simulation]
    )
    adapter = (
        bf.Adapter()
        .to_array()
        .convert_dtype("float64", "float32")
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
            ["damping_values", "t_E", "t_ISy", "t_ISev", "t_Cr",
             "mu_CR", "mu_IH", "mu_HU", "mu_UD", "transmission_prob"],
            into="inference_variables",
            axis=-1
        )
        .concatenate(summary_vars, into="summary_variables", axis=-1)
    )

    summary_network = bf.networks.FusionTransformer(
        summary_dim=(len(bounds)+16)*2, dropout=0.1
    )
    inference_network = bf.networks.FlowMatching()

    # aug = bf.augmentations.NNPE(spike_scale=SPIKE_SCALE, slab_scale=SLAB_SCALE, per_dimension=False)
    workflow = bf.BasicWorkflow(
        simulator=simulator,
        adapter=adapter,
        summary_network=summary_network,
        inference_network=inference_network,
        standardize='all'
        # augmentations={f'region{i}': aug for i in range(len(region_ids))}
    )

    return workflow


def run_training(name, num_training_files=20):
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
        d = skip_2weeks(d) 
        if trainings_data is None:
            trainings_data = d
        else:
            trainings_data = concat_dicts(trainings_data, d)
    aggregate_states(trainings_data)

    # validation data
    validation_data = apply_aug(load_pickle(val_path), aug=aug)
    validation_data = skip_2weeks(validation_data)
    aggregate_states(validation_data)

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
    plots['losses'].savefig(f'{name}/losses_{name}.png')
    plots['recovery'].savefig(f'{name}/recovery_{name}.png')
    plots['calibration_ecdf'].savefig(f'{name}/calibration_ecdf_{name}.png')
    #plots['z_score_contraction'].savefig(f'{name}/z_score_contraction_{name}.png')


def run_inference(name, num_samples=100, on_synthetic_data=False, apply_augmentation=True):
    val_path = f"{name}/validation_data_{name}.pickle"
    synthetic = "_synthetic" if on_synthetic_data else ""
    with_aug = "_with_aug" if apply_augmentation else ""

    aug = bf.augmentations.NNPE(
        spike_scale=SPIKE_SCALE, slab_scale=SLAB_SCALE, per_dimension=False
    )

    if on_synthetic_data:
        # validation data
        validation_data = load_pickle(val_path)
        validation_data = apply_aug(validation_data, aug=aug)
        validation_data_skip2w = skip_2weeks(validation_data)
        aggregate_states(validation_data_skip2w)
        divi_dict = validation_data
        divi_region_keys = region_keys_sorted(divi_dict)

        divi_data = np.concatenate(
            [divi_dict[key] for key in divi_region_keys], axis=-1
        )[0]  # only one dataset
    else:
        divi_dict = load_divi_data()
        validation_data_skip2w = skip_2weeks(divi_dict)
        aggregate_states(validation_data_skip2w)
        divi_region_keys = region_keys_sorted(divi_dict)
        divi_data = np.concatenate(
            [divi_dict[key] for key in divi_region_keys], axis=-1
        )[0]

    workflow = get_workflow()
    workflow.approximator = keras.models.load_model(
        filepath=os.path.join(f"{name}/model_{name}.keras")
    )

    if False: #os.path.exists(f'{name}/sims_{name}{synthetic}{with_aug}.pickle'):
        simulations = load_pickle(f'{name}/sims_{name}{synthetic}{with_aug}.pickle')
        print("loaded simulations from file")
    else:
        samples = workflow.sample(conditions=validation_data_skip2w, num_samples=num_samples)
        results = []
        for i in range(num_samples):  # we only have one dataset for inference here
            result = run_germany_nuts3_simulation(
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
        if apply_augmentation:
            results = apply_aug(results, aug=aug)

        # get sims in shape (samples, time, regions)
        simulations = np.zeros((num_samples, divi_data.shape[0], divi_data.shape[1]))
        for i in range(num_samples):
            simulations[i] = np.concatenate([results[key][i] for key in divi_region_keys], axis=-1)

        # save sims
        with open(f'{name}/sims_{name}{synthetic}{with_aug}.pickle', 'wb') as f:
            pickle.dump(simulations, f, pickle.HIGHEST_PROTOCOL)

    # plot simulations
    fig, ax = plt.subplots(nrows=2, ncols=5, figsize=(12, 5), layout="constrained")
    ax = ax.flatten()
    rand_index = np.random.choice(simulations.shape[-1], replace=False, size=len(ax))
    for i, a in enumerate(ax):
        plot_region_median_mad(
            simulations, region=rand_index[i], true_data=divi_data, label=r"Median $\pm$ Mad", ax=a
        )
    plt.savefig(f'{name}/random_regions_{name}{synthetic}{with_aug}.png')
    plt.close()

    plot_aggregated_over_regions(simulations, true_data=divi_data, label="Region Aggregated Median $\pm$ Mad")
    plt.savefig(f'{name}/region_aggregated_{name}{synthetic}{with_aug}.png')
    plt.close()


    fig, axis = plt.subplots(1, 2, figsize=(10, 4), sharex=True, layout="constrained")
    ax = calibration_curves_per_region(simulations, divi_data, ax=axis[0])
    ax, stats = calibration_median_mad_over_regions(simulations, divi_data, ax=axis[1])
    plt.savefig(f'{name}/calibration_per_region_{name}{synthetic}{with_aug}.png')
    plt.close()

    # plot = bf.diagnostics.pairs_posterior(simulations, priors=validation_data, dataset_id=0)
    # plot.savefig(f'pairs_posterior_wcovidparams_oct{synthetic}_ma7_noise.png')


if __name__ == "__main__":
    name = "skip2w"

    if not os.path.exists(name):
        os.makedirs(name)
    # create_train_data(filename=f'{name}/trainings_data1_{name}.pickle', number_samples=10)
    # run_training(name=name, num_training_files=20)
    run_inference(name=name, on_synthetic_data=True)
    run_inference(name=name, on_synthetic_data=True, apply_augmentation=False)
    run_inference(name=name, on_synthetic_data=False)
    run_inference(name=name, on_synthetic_data=False, apply_augmentation=False)