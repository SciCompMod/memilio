import os
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def get_mean_err(data, plot_compartments, mode):
    labels = data['labels']
    num_samples = labels.shape[0]

    # reverse log scale
    labels = np.expm1(labels)

    # split at 80%
    split_point = int(num_samples * 0.8)

    # Calculate mean across the first 80%
    mean_80_percent = np.mean(labels[:split_point, :, :], axis=0)

    # Extract the last 20% of data
    last_20_percent = labels[split_point:, :, :]

    aggregated_last_20 = np.sum(last_20_percent.reshape(last_20_percent.shape[0],
                                                        last_20_percent.shape[1],
                                                        6, 8), axis=2)

    # Calculate MAPE for each element in the last 20%
    # mape_values = np.abs(
    #     (last_20_percent - mean_80_percent) / last_20_percent) * 100
    mape_values = []
    for ref in last_20_percent:
        add_val = np.abs((ref - mean_80_percent) / ref) * 100
        if add_val.max() > 1e7 or add_val.min() < -1e7:
            test = 1
        if np.isnan(add_val).any():
            indx = np.where(np.isnan(add_val))
            num_nans = len(indx[0])
            test = 2
        mape_values.append(add_val)

    mape_values = np.array(mape_values)
    # replace nan values with 0
    # mape_values = np.nan_to_num(mape_values)

    # plot mape values in a histogram
    # plt.hist(mape_values.flatten(), bins=100)
    # plt.xlabel("MAPE")
    # plt.ylabel("Frequency")
    # plt.yscale('log')
    # plt.show()

    # calculate mean ( -> Sum up and divide by n)
    mape = np.mean(mape_values)
    median = np.median(mape_values)

    if plot_compartments:
        # Define compartment names
        compartments = [
            'Susceptible',
            'Exposed',
            'InfectedNoSymptoms',
            'InfectedSymptoms',
            'InfectedSevere',
            'InfectedCritical',
            'Recovered',
            'Dead'
        ]
        num_compartments = len(compartments)

        # Create a figure with a 2x4 grid of subplots (one subplot per compartment)
        fig, axes = plt.subplots(2, 4, figsize=(16, 8))
        axes = axes.flatten()

        # Choose 200 random runs (or all runs if fewer than 200)
        total_runs = aggregated_last_20.shape[0]
        n_to_plot = 200 if total_runs >= 200 else total_runs
        random_indices = np.random.choice(
            total_runs, size=n_to_plot, replace=False)

        # Plot the time series for each compartment over the 30 days
        for comp in range(num_compartments):
            for run in random_indices:
                # aggregated_last_20[run, :, comp] is the time series for the given run and compartment
                axes[comp].plot(aggregated_last_20[run, :, comp],
                                alpha=0.7,
                                label=None)
            axes[comp].set_xlabel("Days")
            axes[comp].set_ylabel("Value")
            axes[comp].set_title(compartments[comp])
            # Optionally, add a legend in the first subplot
            if comp == 0:
                axes[comp].legend()
        plt.tight_layout()

        # Save the plot in the desired directory
        save_dir = "/localdata1/gnn_paper_2024/images/without_spatial_res/plot_labels/secir_groups/"
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(
            save_dir, f"data_secir_agegroups_{mode}_compartments.png")
        plt.savefig(save_path)
        plt.clf()

    # 10 highest MAPE values
    # print("10 highest MAPE values:")
    # for i in range(10):
    #     max_idx = np.unravel_index(mape_values.argmax(), mape_values.shape)

    # print(f"MAPE for the last 20%: {median:.2f}%")
    return mape


path_datasets = "/localdata1/gnn_paper_2024/data/one_population/with_agegroups_Germany/"
modes = ["nodamp"]  # ,"with_damp"]  #

for mode in modes:
    path_datasets = os.path.join(path_datasets, mode)
    datasets = os.listdir(path_datasets)
    # remove every non .pickle file
    datasets = [dataset for dataset in datasets if dataset.endswith('.pickle')]
    for dataset in datasets:
        with open(os.path.join(path_datasets, dataset), 'rb') as f:
            data = pickle.load(f)

        # calculate MAPE for both datasets
        mape_old = get_mean_err(data, True, dataset)
        print(f"{dataset}_{mode}_MAPE for the last 20%: {mape_old:.2f}%")
