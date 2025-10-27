#############################################################################
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
import numpy as np
import datetime
import os
import memilio.simulation as mio
import memilio.simulation.osecir as osecir
import matplotlib.pyplot as plt

import pandas as pd
from memilio.epidata import defaultDict as dd

excluded_ids = [11001, 11002, 11003, 11004, 11005, 11006, 11007, 11008, 11009, 11010, 11011, 11012, 16056]
no_icu_ids = [7338, 9374, 9473, 9573]
region_ids = [region_id for region_id in dd.County.keys() if region_id not in excluded_ids]

inference_params = ['damping_values', 't_E', 't_ISy', 't_ISev', 'transmission_prob']

def plot_region_median_mad(
    data: np.ndarray,
    region: int,
    true_data = None,
    ax = None,
    label = None,
    color = "red"
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

    line, = ax.plot(x, med, lw=2, label=label or f"Region {region}", color=color)
    band = ax.fill_between(x, med - mad, med + mad, alpha=0.25, color=color)
    if true_data is not None:
        true_vals = true_data[:, region]  # (time_points,)
        ax.plot(x, true_vals, lw=2, color="black", label="True data")

    ax.set_xlabel("Time")
    ax.set_ylabel("ICU")
    ax.set_title(f"Region {region}")
    if label is not None:
        ax.legend()
    return line, band


def plot_aggregated_over_regions(
    data: np.ndarray,
    region_agg = np.sum,
    true_data = None,
    ax = None,
    label = None,
    color = 'red'
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

    line, = ax.plot(x, agg_median, lw=2, label=label or "Aggregated over regions", color=color)
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

    ax.set_xlabel("Time")
    ax.set_ylabel("ICU")
    if label is not None:
        ax.legend()

    return line, band

class Simulation:
    """ """

    def __init__(self, data_dir, start_date, results_dir):
        self.num_groups = 1
        self.data_dir = data_dir
        self.start_date = start_date
        self.results_dir = results_dir
        if not os.path.exists(self.results_dir):
            os.makedirs(self.results_dir)

    def set_covid_parameters(self, model):
        """
        :param model: 
        """
        model.parameters.TimeExposed[mio.AgeGroup(0)] =  3.3
        model.parameters.TimeInfectedNoSymptoms[mio.AgeGroup(0)] = 1.9
        model.parameters.TimeInfectedSymptoms[mio.AgeGroup(0)] = 6.9
        model.parameters.TimeInfectedSevere[mio.AgeGroup(0)] = 8
        model.parameters.TimeInfectedCritical[mio.AgeGroup(0)] = 10.5

        # probabilities
        model.parameters.TransmissionProbabilityOnContact[mio.AgeGroup(0)] = 0.08
        model.parameters.RelativeTransmissionNoSymptoms[mio.AgeGroup(0)] = 1

        model.parameters.RecoveredPerInfectedNoSymptoms[mio.AgeGroup(0)] = 0.2069
        model.parameters.SeverePerInfectedSymptoms[mio.AgeGroup(0)] = 0.07864
        model.parameters.CriticalPerSevere[mio.AgeGroup(0)] = 0.17318
        model.parameters.DeathsPerCritical[mio.AgeGroup(0)] = 0.21718

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

        start_damping = datetime.date(
            year=2020, month=10, day=8)

        if start_damping < end_date:
            start_date = (start_damping - self.start_date).days
            params.ContactPatterns.cont_freq_mat[0].add_damping(mio.Damping(np.r_[damping_value], t=start_date))

    def get_graph(self, end_date):
        """
        :param end_date: 
        """
        print("Initializing model...")
        model = osecir.Model(self.num_groups)
        self.set_covid_parameters(model)
        self.set_contact_matrices(model)
        print("Model initialized.")

        graph = osecir.ModelGraph()

        scaling_factor_infected = [2.5]
        scaling_factor_icu = 1.0
        tnt_capacity_factor = 7.5 / 100000.

        data_dir_Germany = os.path.join(self.data_dir, "Germany")
        mobility_data_file = os.path.join(
            data_dir_Germany, "mobility", "commuter_migration_scaled_dim401.txt")
        pydata_dir = os.path.join(data_dir_Germany, "pydata")

        path_population_data = os.path.join(pydata_dir,
                                            "county_current_population.json")

        print("Setting nodes...")
        osecir.set_nodes(
            model.parameters,
            mio.Date(self.start_date.year,
                     self.start_date.month, self.start_date.day),
            mio.Date(end_date.year,
                     end_date.month, end_date.day), pydata_dir,
            path_population_data, True, graph, scaling_factor_infected,
            scaling_factor_icu, 0, 0, False)

        print("Setting edges...")
        osecir.set_edges(
            mobility_data_file, graph, 1)

        print("Graph created.")

        return graph

    def run(self, num_days_sim, save_graph=True):
        """
        :param num_days_sim: 
        :param num_runs:  (Default value = 10)
        :param save_graph:  (Default value = True)
        :param create_gif:  (Default value = True)
        """
        mio.set_log_level(mio.LogLevel.Warning)
        end_date = self.start_date + datetime.timedelta(days=num_days_sim)

        graph = self.get_graph(end_date)

        mobility_graph = osecir.MobilityGraph()
        for node_idx in range(graph.num_nodes):
            node = graph.get_node(node_idx)
            self.set_npis(node.property.parameters, end_date, 0.3)
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
                results[f'no_icu_region{node_idx}'] = osecir.interpolate_simulation_result(node.property.result)
            else:
                results[f'region{node_idx}'] = osecir.interpolate_simulation_result(node.property.result)

        return results

def run_germany_nuts3_simulation():
    mio.set_log_level(mio.LogLevel.Warning)
    file_path = os.path.dirname(os.path.abspath(__file__))

    sim = Simulation(
        data_dir=os.path.join(file_path, "../../../data"),
        start_date=datetime.date(year=2020, month=10, day=1),
        results_dir=os.path.join(file_path, "../../../results_osecir"))
    num_days_sim = 60

    results = sim.run(num_days_sim)

    return results

def load_divi_data():
    file_path = os.path.dirname(os.path.abspath(__file__))
    divi_path = os.path.join(file_path, "../../../data/Germany/pydata")

    data = pd.read_json(os.path.join(divi_path, "county_divi_all_dates.json"))
    data = data[data['Date']>= np.datetime64(datetime.date(2020, 10, 1))]
    data = data[data['Date'] <= np.datetime64(datetime.date(2020, 10, 1) + datetime.timedelta(days=50))]
    data = data.sort_values(by=['ID_County', 'Date'])
    divi_data = data.pivot(index='Date', columns='ID_County', values='ICU')
    divi_dict = {f"region{i}": divi_data[region_id].to_numpy()[None, :, None] for i, region_id in enumerate(region_ids) if region_id not in no_icu_ids}

    return divi_data.to_numpy(), divi_dict

if __name__ == "__main__":
    run_germany_nuts3_simulation()