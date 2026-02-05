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

name = "nuts0"

region_ids = [0]
inference_params = ['scaling_factor', 'damping_values', 't_E', 't_ISy', 't_ISev',
                    't_Cr', 'mu_CR', 'mu_IH', 'mu_HU', 'mu_UD', 'transmission_prob']
summary_vars = ['state0']

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

    qs_50 = np.quantile(vals, q=[0.25, 0.75], axis=0)
    qs_90 = np.quantile(vals, q=[0.05, 0.95], axis=0)
    qs_95 = np.quantile(vals, q=[0.025, 0.975], axis=0)

    med = np.median(vals, axis=0)

    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(
        x, med, lw=3, label="Simulation", color=color, zorder=1)
    ax.fill_between(x, qs_50[0], qs_50[1], alpha=0.6, linewidth=0,
                    color=color)
    ax.fill_between(x, qs_90[0], qs_90[1], linewidth=0, alpha=0.4,
                    color=color)
    ax.fill_between(x, qs_95[0], qs_95[1], linewidth=0, alpha=0.2,
                    color=color)
    if true_data is not None:
        true_vals = true_data[:, region]  # (time_points,)
        ax.scatter(x, true_vals, color="black", label="Reported data", marker='x', zorder=2)

    # Update x-axis to show dates in YYYY-MM format
    start_date = datetime.date(2020, 10, 1)
    date_labels = [(start_date + datetime.timedelta(days=int(day))).strftime('%Y-%m-%d') for day in x]
    ax.set_xticks(x[::14])  # Set ticks every 7 days
    ax.set_xticklabels([])
    # ax.set_xticklabels(date_labels[::14], rotation=45)

    ax.set_ylabel("ICU cases [#]")


def plot_icu_on_germany(simulations):
    med = np.median(simulations, axis=0)

    population = pd.read_json('data/Germany/pydata/county_current_population_germany.json')
    values = med / population['Population'].to_numpy()[None, :] * 100000

    map_data = gpd.read_file(os.path.join(os.getcwd(), 'tools/vg2500_12-31.utm32s.shape/vg2500/VG2500_KRS.shp'))
    fedstate_data = gpd.read_file(os.path.join(os.getcwd(), 'tools/vg2500_12-31.utm32s.shape/vg2500/VG2500_LAN.shp'))
    state_data = gpd.read_file(os.path.join(os.getcwd(), 'tools/vg2500_12-31.utm32s.shape/vg2500/VG2500_STA.shp'))

    fig, axes = plt.subplots(1, 2, figsize=(4.5, 2.5), layout="constrained")
    vmin = 0
    vmax = 8
    
    plot_map(values[0], state_data, map_data, fedstate_data, axes[0], "Median ICU", vmin, vmax)
    plot_map(values[-1], state_data, map_data, fedstate_data, axes[1], "Median ICU", vmin, vmax)

    plt.savefig(f"{name}/median_icu_germany_nuts0.png", bbox_inches='tight', dpi=dpi)


def plot_map(values, map_data, ax, label, vmin, vmax):
    map_data[label] = map_data['ARS'].map({f"{region_id:05d}": values[0]  for region_id in dd.County.keys()})

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
        start_damping_2 = DATE_TIME + datetime.timedelta(days=25)
        start_damping_3 = DATE_TIME + datetime.timedelta(days=35)

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
        pydata_dir = os.path.join(data_dir_Germany, "pydata")

        path_population_data = os.path.join(
            pydata_dir, "county_current_population_germany.json")

        print("Setting nodes...")
        mio.osecir.set_node_germany(
            model.parameters,
            mio.Date(self.start_date.year,
                     self.start_date.month, self.start_date.day),
            mio.Date(end_date.year,
                     end_date.month, end_date.day), pydata_dir,
            path_population_data, False, graph, scaling_factor_infected,
            scaling_factor_icu, 1.0, 0, False)

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
            results[f'state{node_idx}'] = osecir.interpolate_simulation_result(
                mobility_sim.graph.get_node(node_idx).property.result)

        return results


def run_germany_nuts0_simulation(scaling_factor, damping_values, t_E, t_ISy, t_ISev, t_Cr, mu_CR, mu_IH, mu_HU, mu_UD, transmission_prob, export_timeseries=False):
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
    damping_values = np.zeros((NUM_DAMPING_POINTS, 1))
    for i in range(NUM_DAMPING_POINTS):
        mean = np.random.uniform(0, 1)
        scale = 0.1
        a, b = (0 - mean) / scale, (1 - mean) / scale
        damping_values[i] = truncnorm.rvs(
            a=a, b=b, loc=mean, scale=scale, size=1
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

    data = pd.read_json(os.path.join(divi_path, "germany_divi_ma7.json"))
    data = data[data['Date'] >= np.datetime64(DATE_TIME)]
    data = data[data['Date'] <= np.datetime64(DATE_TIME + datetime.timedelta(days=60))]
    data = data.sort_values(by=['Date'])
    divi_dict = {}
    divi_dict[f"state0"] = data['ICU'].to_numpy()[None, :, None]
    return divi_dict



def extract_observables(simulation_results, observable_index=7):
    for key in simulation_results.keys():
        if key not in inference_params:
            simulation_results[key] = simulation_results[key][:, :, observable_index][..., np.newaxis]
    return simulation_results


def create_train_data(filename, number_samples=1000):

    simulator = bf.simulators.make_simulator(
        [prior, run_germany_nuts0_simulation]
    )
    trainings_data = simulator.sample(number_samples)
    trainings_data = extract_observables(trainings_data)
    with open(filename, 'wb') as f:
        pickle.dump(trainings_data, f, pickle.HIGHEST_PROTOCOL)


def load_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)

def is_state_key(k: str) -> bool:
    return 'state' in k

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


def get_workflow():

    simulator = bf.make_simulator(
        [prior, run_germany_nuts0_simulation]
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
        summary_dim=(len(bounds)+1*NUM_DAMPING_POINTS)*2, dropout=0.1
    )
    inference_network = bf.networks.FlowMatching()

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

    # validation data
    validation_data = apply_aug(load_pickle(val_path), aug=aug)
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

def run_inference(num_samples=1000):

    aug = bf.augmentations.NNPE(
        spike_scale=SPIKE_SCALE, slab_scale=SLAB_SCALE, per_dimension=False
    )

    divi_dict = load_divi_data()
    divi_data = np.concatenate(
        [divi_dict[f'state{i}'] for i in range(len(region_ids))], axis=-1
    )[0]

    workflow = get_workflow()
    workflow.approximator = keras.models.load_model(
        filepath=os.path.join(f"{name}/model_{name}.keras")
    )

    samples = workflow.sample(conditions=divi_dict, num_samples=num_samples)
    samples['damping_values'] = samples['damping_values'].reshape((samples['damping_values'].shape[0], num_samples, 1, NUM_DAMPING_POINTS))
    results = []
    for i in range(num_samples):  # we only have one dataset for inference here
        result = run_germany_nuts0_simulation(
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
    results = apply_aug(results, aug=aug)

    # get sims in shape (samples, time, regions)
    simulations = np.zeros((num_samples, divi_data.shape[0], divi_data.shape[1]))
    for i in range(num_samples):
        simulations[i] = np.concatenate([results[key][i] for key in results.keys()], axis=-1)


    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    # Plot with augmentation
    plot_region_fit(
        simulations, region=0, true_data=divi_data, label="Median", ax=ax, color=colors["Red"]
    )    
    lines, labels = ax.get_legend_handles_labels()
    fig.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.115, 0.99), ncol=1)
    plt.savefig(f'{name}/region_aggregated_{name}.png', dpi=dpi)
    plt.close()

    plot_icu_on_germany(simulations)

if __name__ == "__main__":

    set_fontsize()

    if not os.path.exists(name):
        os.makedirs(name)
    create_train_data(filename=f'{name}/validation_data_{name}.pickle', number_samples=100)
    create_train_data(filename=f'{name}/trainings_data1_{name}.pickle', number_samples=20000)
    run_training(num_training_files=1)
    run_inference(num_samples=1000)