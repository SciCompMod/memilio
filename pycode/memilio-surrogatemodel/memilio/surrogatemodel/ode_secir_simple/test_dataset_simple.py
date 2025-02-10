import pickle
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

path_datasets = "/localdata1/gnn_paper_2024/data/one_population/without_agegroups/"

# list all files in path_datasets
datasets = os.listdir(path_datasets)

for dataset in datasets:
    path_dataset = os.path.join(path_datasets, dataset)
    modes = ["non_log", "log"]
    # read the data
    with open(path_dataset, 'rb') as f:
        data = pickle.load(f)

    for mode in modes:
        labels = data['labels']
        num_samples = labels.shape[0]

        # labels are log scaled. Reverse this using np.expm1
        if mode == "non_log":
            labels = np.expm1(labels)

        # split at 80%
        split_point = int(num_samples * 0.8)

        # Calculate mean across the first 80%
        mean_80_percent = np.mean(labels[:split_point, :, :], axis=0)

        # Extract the last 20% of data
        last_20_percent = labels[split_point:, :, :]

        # Calculate MAPE for each element in the last 20%
        mape_values = np.abs(
            (last_20_percent - mean_80_percent) / last_20_percent) * 100

        # plot mape values in a histogram
        # plt.hist(mape_values.flatten(), bins=100)
        # plt.xlabel("MAPE")
        # plt.ylabel("Frequency")
        # plt.yscale('log')
        # plt.show()

        # calculate mean ( -> Sum up and divide by n)
        mape = np.mean(mape_values)
        median = np.median(mape_values)

        print(f"{dataset}_{mode}_MAPE for the last 20%: {mape:.2f}%")
