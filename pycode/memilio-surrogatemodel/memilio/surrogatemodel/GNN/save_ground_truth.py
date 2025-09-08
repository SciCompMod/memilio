import os
import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
import h5py


def read_case_data():
    df = pd.read_hdf(
        "/localdata1/hege_mn/memilio/data/Germany/pydata/Results_rki.h5")
    return df


def read_results_h5(path, group_key='Group1'):
    with h5py.File(path, 'r') as f:
        keys = list(f.keys())
        res = {}
        for i, key in enumerate(keys):
            group = f[key]
            total = group[group_key][()]
            # remove confirmed compartments
            sum_inf_no_symp = np.sum(total[:, [2, 3]], axis=1)
            sum_inf_symp = np.sum(total[:, [4, 5]], axis=1)
            total[:, 2] = sum_inf_no_symp
            total[:, 4] = sum_inf_symp
            total = np.delete(total, [3, 5], axis=1)
            res[key] = total
    return res


def main():
    path_rki_h5 = "/localdata1/hege_mn/memilio/data/Germany/pydata/Results_rki.h5"
    num_age_groups = 6
    all_age_data_list = []
    group_keys = ['Group1', 'Group2', 'Group3', 'Group4', 'Group5', 'Group6']
    for age in range(num_age_groups):
        age_data_dict = read_results_h5(path_rki_h5, group_keys[age])
        age_data_np = np.array(list(age_data_dict.values()))
        all_age_data_list.append(age_data_np)
    # Combine age groups to get shape (num_nodes, timesteps, num_features=48)
    all_age_data = np.stack(all_age_data_list, axis=-1)
    print(f"Shape of all_age_data: {all_age_data.shape}")
    num_nodes, timesteps, _, _ = all_age_data.shape

    # Save the ground truth data
    output_path = "/localdata1/hege_mn/memilio/data/Germany/pydata/ground_truth_all_nodes.pickle"
    with open(output_path, 'wb') as f:
        pd.to_pickle(all_age_data, f)
    print(f"Ground truth data saved to {output_path}")


if __name__ == "__main__":
    main()
