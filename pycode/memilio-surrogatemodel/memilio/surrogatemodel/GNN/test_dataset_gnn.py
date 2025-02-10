import pickle
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

path_datasets = "/localdata1/gnn_paper_2024/data/GNNs/"


# list all files in path_datasets
datasets = os.listdir(path_datasets)

days = 90

for dataset in datasets:
    # only specific datasets
    if dataset != f"GNN_data_{days}days_nodeswithvariance_1k.pickle" and dataset != f"GNN_data_{days}days_samenodes_1k.pickle":
        continue
    path_dataset = os.path.join(path_datasets, dataset)
    modes = ["non_log", "log"]

    # read the data
    with open(path_dataset, 'rb') as f:
        data = pickle.load(f)

    for mode in modes:
        labels = data['labels']
        num_samples = labels.shape[0]

        # reverse log scale
        if mode == "non_log":
            labels = np.expm1(labels)

        # split at 80%
        split_point = int(num_samples * 0.8)

        # Calculate mean across the first 80%
        mean_80_percent = np.mean(labels[:split_point, :, :], axis=0)

        # Extract the last 20% of data
        last_20_percent = labels[split_point:, :, :]

        if mode == "non_log":
            compartments = ['Susceptible',
                            'Exposed',
                            'InfectedNoSymptoms',
                            'InfectedSymptoms',
                            'InfectedSevere',
                            'InfectedCritical',
                            'Recovered',
                            'Dead']

            last_20_percent_summed = np.sum(last_20_percent, axis=-1)
            last_20_percent_sum_groups = np.split(
                last_20_percent_summed, 6, axis=1)

            fig, axes = plt.subplots(2, 4, figsize=(16, 8))
            axes = axes.flatten()

            for comp in range(8):
                for run in range(last_20_percent_summed.shape[0]):
                    comp_groups = [group[run, comp, :]
                                   for group in last_20_percent_sum_groups]
                    mean_total = np.sum(comp_groups, axis=0)
                    axes[comp].plot(mean_total)

                axes[comp].set_xlabel("Days")
                axes[comp].set_ylabel("Value")
                axes[comp].set_title(f"{compartments[comp]}")

            plt.tight_layout()
            plt.savefig(
                f"/localdata1/gnn_paper_2024/images/with_spatial_res/plot_labels_gnn//{dataset}_{mode}_compartments.png")
            plt.clf()

        # Calculate MAPE for each element in the last 20%
        mape_values = []
        for i in range(len(last_20_percent)):
            if np.any(last_20_percent[i] == 0):
                # zero_count = np.sum(last_20_percent[i] == 0)
                # print(
                #     f"There are {zero_count} zero values in last_20_percent[i].")
                continue
            mape_values.append(np.abs(
                (last_20_percent[i] - mean_80_percent) / last_20_percent[i]) * 100)

        # plot mape values in a histogram
        # plt.hist(mape_values.flatten(), bins=100)
        # plt.xlabel("MAPE")
        # plt.ylabel("Frequency")
        # plt.yscale('log')
        # plt.show()

        # calculate mean ( -> Sum up and divide by n)
        mape = np.mean(mape_values)
        median = np.median(mape_values)

        # make a histogram of the mape values
        # plt.hist(mape_values, bins=100)
        # plt.xlabel("MAPE")
        # plt.ylabel("Frequency")
        # plt.yscale('log')
        # plt.title(f"{dataset}_{mode}")
        # plt.show()

        print(f"{dataset}_{mode}_MAPE for the last 20%: {mape:.2f}%")
