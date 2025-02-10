import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def get_mean_err(data):
    labels = data['labels']
    num_samples = labels.shape[0]

    # reverse log scale
    # labels = np.expm1(labels)

    # split at 80%
    split_point = int(num_samples * 0.8)

    # Calculate mean across the first 80%
    mean_80_percent = np.mean(labels[:split_point, :, :], axis=0)

    # Extract the last 20% of data
    last_20_percent = labels[split_point:, :, :]

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
    plt.hist(mape_values.flatten(), bins=100)
    plt.xlabel("MAPE")
    plt.ylabel("Frequency")
    plt.yscale('log')
    # plt.show()

    # calculate mean ( -> Sum up and divide by n)
    mape = np.mean(mape_values)
    median = np.median(mape_values)

    # 10 highest MAPE values
    # print("10 highest MAPE values:")
    # for i in range(10):
    #     max_idx = np.unravel_index(mape_values.argmax(), mape_values.shape)

    # print(f"MAPE for the last 20%: {median:.2f}%")
    return mape


path_dataset_old = "/localdata1/gnn_paper_2024/data/one_population/with_agegroups_Germany/nodamp//data_secir_groups_60days_I_based_Germany_100k_nodamp.pickle"
path_dataset_ibased = "/localdata1/gnn_paper_2024/data/one_population/with_agegroups_Germany/nodamp//data_secir_groups_60days_I_based_Germany_100k_nodamp.pickle"

# read the data
with open(path_dataset_old, 'rb') as f:
    data_old = pickle.load(f)

with open(path_dataset_ibased, 'rb') as f:
    data_ibased = pickle.load(f)

# calculate MAPE for both datasets
mape_old = get_mean_err(data_old)
mape_ibased = get_mean_err(data_ibased)

print(f"Mean MAPE for old dataset: {mape_old}%")
print(f"Mean MAPE for I-based dataset: {mape_ibased}%")
