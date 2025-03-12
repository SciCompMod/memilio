import os
import pickle
import numpy as np
import matplotlib.pyplot as plt


def load_data(file_path):
    """Load and return data from a pickle file."""
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
    return data


def compute_mape(labels, mode):
    """
    Compute the Mean Absolute Percentage Error (MAPE) for the last 20% of the data.

    If mode is 'non_log', reverse the log scaling using np.expm1.

    Returns:
        last_20_percent: The last 20% slice of the labels.
        mape: The overall mean absolute percentage error.
        median: The median absolute percentage error.
    """
    if mode == "non_log":
        labels = np.expm1(labels)

    num_samples = labels.shape[0]
    split_point = int(num_samples * 0.8)

    # Calculate the mean over the first 80% of the data
    mean_80_percent = np.mean(labels[:split_point, :, :], axis=0)
    # Extract the last 20% of the data
    last_20_percent = labels[split_point:, :, :]

    # Calculate MAPE for each element in the last 20%
    mape_values = np.abs(
        (last_20_percent - mean_80_percent) / last_20_percent) * 100
    mape = np.mean(mape_values)
    median = np.median(mape_values)

    return last_20_percent, mape, median


def plot_compartments(last_20_percent, dataset, mode, save_dir):
    """
    Plot the compartment time series for each run using data from the last 20%.
    The plot is saved to the specified directory.
    """
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

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    for comp_idx in range(num_compartments):
        for run in range(200):
            axes[comp_idx].plot(last_20_percent[run, :, comp_idx], alpha=0.7)
        axes[comp_idx].set_xlabel("Days")
        axes[comp_idx].set_ylabel("Value")
        axes[comp_idx].set_title(compartments[comp_idx])
        # Optionally add legend to the first subplot if needed
        if comp_idx == 0:
            axes[comp_idx].legend()

    plt.tight_layout()
    os.makedirs(save_dir, exist_ok=True)
    save_path = os.path.join(save_dir, f"{dataset}_{mode}_compartments.png")
    plt.savefig(save_path)
    plt.clf()


def plot_mean_std(last_20_percent, dataset, mode, save_dir):
    """
    Plot the mean and standard deviation of each compartment time series using data from the last 20%.
    The plot is saved to the specified directory.
    """
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
    runs = last_20_percent.shape[0]
    days = last_20_percent.shape[1]

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    for comp_idx in range(num_compartments):
        # Extract data for the current compartment (shape: runs x days)
        compartment_data = last_20_percent[:, :, comp_idx]
        # Compute mean and standard deviation across runs for each day
        mean_ts = np.mean(compartment_data, axis=0)
        std_ts = np.std(compartment_data, axis=0)

        days_range = np.arange(days)
        axes[comp_idx].plot(days_range, mean_ts, label='Mean')
        axes[comp_idx].fill_between(days_range, mean_ts - std_ts, mean_ts + std_ts,
                                    alpha=0.3, label='Std Dev')
        axes[comp_idx].set_xlabel("Days")
        axes[comp_idx].set_ylabel("Value")
        axes[comp_idx].set_title(compartments[comp_idx])
        axes[comp_idx].legend()

    plt.tight_layout()
    os.makedirs(save_dir, exist_ok=True)
    save_path = os.path.join(
        save_dir, f"{dataset}_{mode}_mean_std_compartments.png")
    plt.savefig(save_path)
    plt.clf()


def process_dataset(dataset_file, base_dataset_path, save_dir):
    """
    Process a single dataset:
      - Load the data.
      - Compute MAPE for both 'non_log' and 'log' modes.
      - Plot compartment time series for 'non_log' mode.
      - Plot mean and standard deviation for each compartment for 'non_log' mode.
    """
    dataset_path = os.path.join(base_dataset_path, dataset_file)
    data = load_data(dataset_path)

    modes = ["non_log", "log"]
    for mode in modes:
        labels = data['labels']
        last_20_percent, mape, median = compute_mape(labels, mode)
        print(f"{dataset_file}_{mode}_MAPE for the last 20%: {mape:.2f}%")

        if mode == "non_log":
            plot_compartments(last_20_percent, dataset_file, mode, save_dir)
            plot_mean_std(last_20_percent, dataset_file, mode, save_dir)


def main():
    base_dataset_path = "/localdata1/gnn_paper_2024/data/one_population/without_agegroups/"
    save_dir = "/localdata1/gnn_paper_2024/images/without_spatial_res/plot_labels/secir_simple/"

    dataset_files = os.listdir(base_dataset_path)
    # only use specific datasets 'data_secir_simple_30days_I_based_10k.pickle' and 'data_secir_simple_30days_10k.pickle'
    # dataset_files = ['data_secir_simple_30days_I_based_10k.pickle',
    #                  'data_secir_simple_30days_10k.pickle']
    for dataset_file in dataset_files:
        process_dataset(dataset_file, base_dataset_path, save_dir)


if __name__ == "__main__":
    main()
