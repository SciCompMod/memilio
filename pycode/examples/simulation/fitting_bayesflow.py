#############################################################################
# Copyright (C) 2020-2026 MEmilio
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
import memilio.simulation.halle_abm as halle_abm
import memilio.simulation as mio
import keras
import bayesflow as bf
import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
os.environ["KERAS_BACKEND"] = "tensorflow"


name = "fitting_abm"

inference_params = ['viral_shedding_rate', 'dark_figure']
summary_vars = ['deaths', 'severe', 'critical']

bounds = {
    'viral_shedding_rate': (1.4, 2.0),
    'dark_figure': (2.5, 5.5)
}
DATE_TIME = mio.Date(2022, 7, 1)


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


def plot_results(results, validation_data, name, synthetic=False):
    synthetic = "_synthetic" if synthetic else ""
    for key in summary_vars:
        num_days = results[key].shape[0]
        plt.figure(figsize=(6, 4))
        print(validation_data[key])
        plt.plot(np.arange(num_days), np.median(
            results[key], axis=1), color=colors["Blue"])
        plt.plot(np.arange(num_days),
                 validation_data[key], color=colors["Red"])
        plt.title(key)
        plt.xlabel("time [days]")
        plt.ylabel("number of people")
        plt.savefig(f"{name}/{key}{synthetic}.png", dpi=dpi)
        plt.close()


def run_abm(viral_shedding_rate, dark_figure):
    mio.set_log_level(mio.LogLevel.Warning)
    file_path = os.path.dirname(os.path.abspath(__file__))

    person_file = 'data/Germany/Halle/halle_population_data.csv'
    contact_dir = 'data/Germany/contacts'
    history_file = 'data/Germany/Halle/260304_infections_vaccines_halle.csv'
    history_lookback_days = 0
    allow_missing_history = True
    start_date = DATE_TIME
    num_days_sim = 30

    result = halle_abm.run_once(person_file, contact_dir, history_file, history_lookback_days, allow_missing_history, [
                                viral_shedding_rate, dark_figure], start_date, num_days_sim)
    result = np.vstack([np.asarray(result[0]), np.asarray(result[1]).T])
    return {'deaths': result[1],
            'severe': result[2],
            'critical': result[3]}


def prior():
    return {
        'viral_shedding_rate': np.random.uniform(*bounds['viral_shedding_rate']),
        'dark_figure': np.random.uniform(*bounds['dark_figure'])
    }


def load_real_data():

    return None


def create_train_data(filename, number_samples=1000):

    simulator = bf.simulators.make_simulator(
        [prior, run_abm]
    )
    trainings_data = simulator.sample(number_samples)
    with open(filename, 'wb') as f:
        pickle.dump(trainings_data, f, pickle.HIGHEST_PROTOCOL)


def load_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def is_state_key(k: str) -> bool:
    return 'state' in k


def concat_dicts(base: dict, new: dict) -> dict:
    missing = set(base) - set(new)
    if missing:
        raise KeyError(f"new dict missing keys: {sorted(missing)}")
    for k in base:
        base[k] = np.concatenate([base[k], new[k]], axis=-1)
    return base


def combine_results(dict_list):
    combined = {}
    for d in dict_list:
        combined = concat_dicts(combined, d) if combined else d
    return combined


def reshape_results(results):
    for key in summary_vars:
        results[key] = results[key][..., np.newaxis]
    return results


def get_workflow():

    simulator = bf.make_simulator(
        [prior, run_abm]
    )
    adapter = (
        bf.Adapter()
        .to_array()
        .convert_dtype("float64", "float32")
        .constrain("viral_shedding_rate", lower=bounds["viral_shedding_rate"][0], upper=bounds["viral_shedding_rate"][1])
        .constrain("dark_figure", lower=bounds["dark_figure"][0], upper=bounds["dark_figure"][1])
        .concatenate(
            ["viral_shedding_rate", "dark_figure"],
            into="inference_variables",
            axis=-1
        )
        .concatenate(summary_vars, into="summary_variables", axis=-1)
    )

    summary_network = bf.networks.TimeSeriesNetwork(
        summary_dim=len(bounds)*2, dropout=0.1
    )
    inference_network = bf.networks.FlowMatching()

    workflow = bf.BasicWorkflow(
        simulator=simulator,
        adapter=adapter,
        summary_network=summary_network,
        inference_network=inference_network,
        standardize='all'
    )

    return workflow


def run_training(num_training_files=20):
    train_template = name+"/trainings_data{i}_"+name+".pickle"
    val_path = f"{name}/validation_data_{name}.pickle"

    # training data
    train_files = [train_template.format(i=i)
                   for i in range(1, 1+num_training_files)]
    trainings_data = None
    for p in train_files:
        d = reshape_results(load_pickle(p))
        if trainings_data is None:
            trainings_data = d
        else:
            trainings_data = concat_dicts(trainings_data, d)

    # validation data
    validation_data = reshape_results(load_pickle(val_path))

    # check data
    workflow = get_workflow()
    print("summary_variables shape:", workflow.adapter(
        trainings_data)["summary_variables"].shape)
    print("inference_variables shape:", workflow.adapter(
        trainings_data)["inference_variables"].shape)

    history = workflow.fit_offline(
        data=trainings_data, epochs=3, batch_size=64, validation_data=validation_data
    )

    workflow.approximator.save(
        filepath=os.path.join(f"{name}/model_{name}.keras")
    )


def run_inference(name, num_samples=1000, on_synthetic_data=False):
    val_path = f"{name}/validation_data_{name}.pickle"
    synthetic = "_synthetic" if on_synthetic_data else ""

    validation_data = reshape_results(load_pickle(val_path))
    if on_synthetic_data:
        validation_dict = {
            key: validation_data[key][0, ...] for key in summary_vars}
        real_data = np.concatenate([validation_data[key][0]
                                   for key in summary_vars], axis=-1)

    else:
        real_data = load_real_data()  # to be determined

    workflow = get_workflow()
    workflow.approximator = keras.models.load_model(
        filepath=os.path.join(f"{name}/model_{name}.keras")
    )

    if os.path.exists(f'{name}/results_{name}{synthetic}.pickle'):
        results = load_pickle(f'{name}/results_{name}{synthetic}.pickle')
        print("loaded results from file")
    else:
        samples = workflow.sample(
            conditions=validation_dict, num_samples=num_samples)
        results = []
        for i in range(num_samples):  # we only have one dataset for inference here
            result = run_abm(
                viral_shedding_rate=samples['viral_shedding_rate'][0, i][0],
                dark_figure=samples['dark_figure'][0, i][0]
            )
            results.append(reshape_results(result))
        results = combine_results(results)

        with open(f'{name}/results_{name}{synthetic}.pickle', 'wb') as f:
            pickle.dump(results, f, pickle.HIGHEST_PROTOCOL)

    plot_results(results, validation_dict, name, synthetic=on_synthetic_data)


if __name__ == "__main__":

    set_fontsize()

    if not os.path.exists(name):
        os.makedirs(name)
    create_train_data(
        filename=f'{name}/validation_data_{name}.pickle', number_samples=100)
    create_train_data(
        filename=f'{name}/trainings_data1_{name}.pickle', number_samples=1000)
    run_training(num_training_files=1)
    run_inference(name=name, num_samples=100, on_synthetic_data=True)
