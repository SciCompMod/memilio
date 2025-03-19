import os
import pickle
import geopandas as gpd
import numpy as np
import tensorflow as tf
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from keras.layers import Dense
from keras.models import Model
from matplotlib.cm import ScalarMappable
from spektral.layers import ARMAConv
from spektral.utils.convolution import normalized_laplacian, rescale_laplacian

# Parameters (adjust as needed)
DAY = 30  # Example: 30-day prediction
NUMBER_OF_NODES = 400
# In training, labels for 30 days with 6 age groups and 8 compartments were used:
OUTPUT_DIM = DAY * 6 * 8  # 1440
INPUT_WIDTH = 5          # 5 input days
INPUT_FEATURE_DIM = INPUT_WIDTH * 48  # 5 days * 48 features = 240

# Paths (adjust as needed)
CWD = os.getcwd()
COMMUTER_FILE_PATH = os.path.join(
    CWD, 'data', 'mobility', 'commuter_mobility.txt')
PATH_WEIGHTS = '/localdata1/gnn_paper_2024/data/results/saved_models/saved_models_GNN/'
WEIGHTS_FILENAME = f"GNN_{DAY}days_nodeswithvariance_1k_test.pickle"


def get_adjacency_matrix(number_of_nodes, commuter_file_path):
    """
    Reads the commuter data and creates a binary adjacency matrix.
    """
    commuter_data = pd.read_csv(commuter_file_path, sep=" ", header=None)
    sub_matrix = commuter_data.iloc[:number_of_nodes, :number_of_nodes]
    adjacency_matrix = sub_matrix.to_numpy()
    # Set all values > 0 to 1
    adjacency_matrix[adjacency_matrix > 0] = 1
    return adjacency_matrix


class Net(Model):
    """
    The GNN model with seven ARMAConv layers and a final Dense layer.
    """

    def __init__(self, output_dim, channels=512, number_of_layers=7, layer=ARMAConv):
        super().__init__()
        self.convs = []
        for _ in range(number_of_layers):
            self.convs.append(layer(channels, activation='relu'))
        self.dense = Dense(output_dim, activation='linear')

    def call(self, inputs):
        # Expecting inputs as a tuple: (x, a)
        x, a = inputs  # x: Node features, a: Adjacency matrix
        for conv in self.convs:
            x = conv([x, a])
        output = self.dense(x)
        return output


def load_model_and_weights(a):
    """
    Creates the model, builds it, and loads the saved weights.
    The parameter a is the (scaled) adjacency matrix.
    """
    model = Net(output_dim=OUTPUT_DIM)

    # Create the model variables by running a dummy forward pass
    dummy_x = tf.random.normal((1, NUMBER_OF_NODES, INPUT_FEATURE_DIM))
    _ = model([dummy_x, a])
    print("Model successfully built. Number of weight tensors:", len(model.weights))

    # Load the weights from the pickle file
    weight_path = os.path.join(PATH_WEIGHTS, WEIGHTS_FILENAME)
    if not os.path.exists(weight_path):
        raise FileNotFoundError(f"Weight file not found: {weight_path}")
    with open(weight_path, "rb") as fp:
        weights = pickle.load(fp)
    model.set_weights(weights)
    print("Model weights successfully loaded.")
    return model


def prepare_input_data():
    """
    Prepares sample input data.
    Expected shape: (batch, NUMBER_OF_NODES, INPUT_FEATURE_DIM)
    If a file 'sample_input.npy' exists in the current directory, it is loaded;
    otherwise, dummy data is generated.
    """
    sample_file = os.path.join(CWD, "sample_input.npy")
    if os.path.exists(sample_file):
        x = np.load(sample_file)
        print(f"Input data loaded from {sample_file}.")
    else:
        x = np.random.rand(
            NUMBER_OF_NODES, INPUT_FEATURE_DIM).astype('float32')
        print("Dummy input data generated.")
    # Add batch dimension
    x = np.expand_dims(x, axis=0)
    return x


def load_dataset(day):
    """
    Loads the dataset from a pickle file.
    The expected file is located in '/localdata1/gnn_paper_2024/data/GNNs/' and contains a dictionary
    with keys 'inputs' and 'labels'.
    """
    path_datasets = '/localdata1/gnn_paper_2024/data/GNNs/'
    dataset_file = os.path.join(
        path_datasets, f'GNN_data_{day}days_nodeswithvariance_1k.pickle')
    if not os.path.exists(dataset_file):
        raise FileNotFoundError(f"Dataset file not found: {dataset_file}")
    with open(dataset_file, 'rb') as f:
        data_secir = pickle.load(f)
    print(f"Dataset loaded from {dataset_file}.")
    return data_secir['inputs'], data_secir['labels']


def compute_metrics(model, test_inputs, test_labels, a):
    """
    Uses the provided model to predict on test_inputs and computes the Mean Absolute Percentage Error (MAPE)
    on both the scaled data and the reversed (original) scale.

    In addition, it computes:
      - mape_per_run: one error value per run (flattening regions and features),
      - mape_per_region: error per region for each run,
      - mean_mape_per_region: the mean MAPE per region averaged over all runs (shape: (n_regions,))
    """
    # If test_inputs has 4 dimensions (e.g., shape: (n_runs, 48, 5, 400)),
    # then transpose and reshape to the expected shape: (n_runs, NUMBER_OF_NODES, INPUT_FEATURE_DIM)
    if len(test_inputs.shape) == 4:
        # Assuming original shape: (n_runs, 48, 5, 400)
        # Transpose to (n_runs, 400, 5, 48) and then reshape to (n_runs, 400, 5*48)
        # Now (n_runs, 400, 5, 48)
        test_inputs = np.transpose(test_inputs, (0, 3, 2, 1))
        test_inputs = test_inputs.reshape(
            test_inputs.shape[0], test_inputs.shape[1], -1)  # (n_runs, 400, 240)
    if len(test_labels.shape) == 4:
        test_labels = np.transpose(test_labels, (0, 3, 2, 1))
        test_labels = test_labels.reshape(
            test_labels.shape[0], test_labels.shape[1], -1)

    # Ensure the adjacency matrix 'a' is repeated for each sample in the batch
    if len(a.shape) == 2:
        n_runs = test_inputs.shape[0]
        a_batch = np.repeat(np.expand_dims(a, axis=0), n_runs, axis=0)
    else:
        a_batch = a

    # Make predictions using the model
    pred = model.predict([test_inputs, a_batch])

    # Reverse the log scaling transformation using np.expm1
    pred_reversed = np.expm1(pred)
    labels_reversed = np.expm1(np.asarray(test_labels))

    # Compute MAPE per day (averaging across runs)
    mape_per_day = 100 * \
        np.mean(np.abs((test_labels - pred) / test_labels), axis=0)
    mape_reversed_per_day = 100 * \
        np.mean(np.abs((labels_reversed - pred_reversed) / labels_reversed), axis=0)

    # Get dimensions from test_labels
    n_runs = test_labels.shape[0]
    n_regions = test_labels.shape[1]  # expected to be 400
    n_features = test_labels.shape[2]

    # Compute MAPE per run (flattening regions and features)
    test_labels_flat = test_labels.reshape(n_runs, -1)
    pred_flat = pred.reshape(n_runs, -1)
    labels_reversed_flat = labels_reversed.reshape(n_runs, -1)
    pred_reversed_flat = pred_reversed.reshape(n_runs, -1)
    mape_per_run = 100 * \
        np.mean(np.abs((test_labels_flat - pred_flat) / test_labels_flat), axis=1)
    mape_reversed_per_run = 100 * \
        np.mean(np.abs((labels_reversed_flat - pred_reversed_flat) /
                labels_reversed_flat), axis=1)

    # Compute MAPE per region for each run: shape (n_runs, n_regions)
    mape_per_region = 100 * \
        np.mean(np.abs((test_labels - pred) / test_labels), axis=2)
    mape_reversed_per_region = 100 * \
        np.mean(np.abs((labels_reversed - pred_reversed) / labels_reversed), axis=2)

    # Compute mean MAPE per region over all runs, resulting in shape (n_regions,)
    mean_mape_per_region = np.mean(mape_per_region, axis=0)
    mean_mape_reversed_per_region = np.mean(mape_reversed_per_region, axis=0)

    return mean_mape_per_region, mean_mape_reversed_per_region


def plot_map(mean_mape, mean_mape_reversed, size_colorbar=40, size_ticks=30, size_edge=0.5):
    """
    Reads a shapefile and plots the mean MAPE per region with a continuous colorbar.
    Two separate plots are generated:
      - One for the scaled (log scale) data.
      - One for the reversed (original) data.

    The shapefile is expected to contain at least 400 regions; if it contains more,
    only the first 400 are used.

    Parameters:
        mean_mape (list or array): MAPE values for the log-scaled plot.
        mean_mape_reversed (list or array): MAPE values for the original-scale plot.
        size_title (int): Font size for plot titles.
        size_colorbar (int): Font size for colorbar labels.
        size_ticks (int): Font size for tick labels.
        size_edge (float): Line width for region borders.
    """

    # Adjust the path to your shapefile
    shapefile_path = os.path.join(
        CWD, 'tools', 'vg2500_12-31.utm32s.shape', 'vg2500', 'VG2500_KRS.shp')

    if not os.path.exists(shapefile_path):
        raise FileNotFoundError(f"Shapefile not found: {shapefile_path}")

    # Load shapefile
    map_data = gpd.read_file(shapefile_path)

    # Directory to save plots
    save_dir_plots = "/localdata1/gnn_paper_2024/images/with_spatial_res/mean_mape_spatially_resolved/"
    os.makedirs(save_dir_plots, exist_ok=True)

    # Use only first 400 regions if more are present
    if len(map_data) > 400:
        map_data = map_data.iloc[:400].copy()

    # Ensure mean_mape and mean_mape_reversed have correct length
    if len(mean_mape) != len(map_data) or len(mean_mape_reversed) != len(map_data):
        raise ValueError(
            "Length of 'mean_mape' or 'mean_mape_reversed' does not match the number of regions.")

    # Add MAPE values as new columns
    map_data['mean_mape'] = mean_mape
    map_data['mean_mape_reversed'] = mean_mape_reversed

    cmap = 'plasma'

    def plot_and_save(column, filename, colorbar_label):
        """ Helper function to plot and save maps with a colorbar. """
        fig, ax = plt.subplots(figsize=(10, 10))
        norm = Normalize(vmin=map_data[column].min(),
                         vmax=map_data[column].max())
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm._A = []  # Dummy array for colorbar compatibility

        # Plot the regions
        map_data.plot(column=column, ax=ax, cmap=cmap,
                      edgecolor='black', linewidth=size_edge, legend=False)

        # Add colorbar
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label(colorbar_label, fontsize=size_colorbar)
        cbar.ax.tick_params(labelsize=size_ticks)

        # Remove axis and save
        ax.set_axis_off()
        plt.tight_layout()
        plt.savefig(os.path.join(save_dir_plots, filename), dpi=300)
        plt.close(fig)

    # Plot for log-scaled data
    plot_and_save('mean_mape', "GNN_mean_mape_log_scaled.png",
                  "MAPE (log. scale)")

    # Plot for reversed (original scale) data
    plot_and_save('mean_mape_reversed',
                  "GNN_mean_mape_unscaled.png", "MAPE (orig. scale)")


def predict_and_plot(model, dataset_path, a):
    """
    Loads a real dataset from the specified pickle file, makes a prediction using the model,
    and plots the predictions versus the true labels for each compartment.

    Additionally, it plots the 5 historical input days as a dotted line (on the left)
    preceding the 30-day predictions/labels.

    Parameters:
      - model: the loaded Keras model.
      - dataset_path: full path to the pickle file containing the dataset.
      - a: the (scaled) adjacency matrix.

    The function assumes the dataset is a dictionary with keys 'inputs' and 'labels'.  
    If the data are 4D, they are transposed and reshaped to (n_samples, NUMBER_OF_NODES, INPUT_FEATURE_DIM).
    """
    # -------------------------
    # Load the dataset
    # -------------------------
    with open(dataset_path, 'rb') as f:
        data = pickle.load(f)
    inputs = data['inputs']
    labels = data['labels']

    # If inputs/labels have 4 dimensions (e.g., shape: (n_samples, 48, 5, 400)),
    # transpose and reshape them to (n_samples, NUMBER_OF_NODES, INPUT_FEATURE_DIM)
    if len(inputs.shape) == 4:
        # Original assumed shape: (n_samples, 48, 5, 400)
        # Transpose to (n_samples, 400, 5, 48) and then reshape to (n_samples, 400, 240)
        inputs = np.transpose(inputs, (0, 3, 2, 1))
        inputs = inputs.reshape(inputs.shape[0], inputs.shape[1], -1)
    if len(labels.shape) == 4:
        labels = np.transpose(labels, (0, 3, 2, 1))
        labels = labels.reshape(labels.shape[0], labels.shape[1], -1)

    # -------------------------
    # Make predictions
    # -------------------------
    # Ensure the adjacency matrix 'a' is repeated for each sample in the batch.
    # In this case, there is only 1 sample.
    a_batch = np.repeat(np.expand_dims(a, axis=0), inputs.shape[0], axis=0)
    # Make predictions using the model.
    pred = model.predict([inputs, a_batch])
    # Reverse the log scaling transformation for predictions and true labels.
    pred_reversed = np.expm1(pred)
    labels_reversed = np.expm1(labels)

    # -------------------------
    # Process predictions and labels:
    # Sum over regions (axis=1) so that shape becomes (1, 1440)
    pred_sum_regions = np.sum(pred_reversed, axis=1)
    labels_sum_regions = np.sum(labels_reversed, axis=1)
    # Remove the sample dimension (only one sample)
    pred_sum_regions = pred_sum_regions[0]  # shape: (1440,)
    labels_sum_regions = labels_sum_regions[0]  # shape: (1440,)
    # Reshape into (days, age_groups, compartments), i.e. (30, 6, 8)
    pred_reshaped = pred_sum_regions.reshape(30, 6, 8)
    labels_reshaped = labels_sum_regions.reshape(30, 6, 8)
    # Sum over the age groups (axis=1) to get values per compartment for each day → (30, 8)
    pred_compartments = np.sum(pred_reshaped, axis=1)
    labels_compartments = np.sum(labels_reshaped, axis=1)

    # -------------------------
    # Process inputs (historical data):
    # The inputs (after reshaping) have shape (1, 400, 240), where 240 = 5*48.
    # Sum over regions (axis=1) → shape becomes (1, 240)
    # reverse transformation as well
    inputs_sum_regions = np.sum(np.expm1(inputs), axis=1)
    inputs_sum_regions = inputs_sum_regions[0]  # shape: (240,)
    # Reshape into (days, age_groups, compartments) for the 5 input days → (5, 6, 8)
    inputs_reshaped = inputs_sum_regions.reshape(5, 6, 8)
    # Sum over age groups to get (5, 8)
    inputs_compartments = np.sum(inputs_reshaped, axis=1)

    # -------------------------
    # Plotting:
    # We want to plot:
    #   - The 5 historical input days as a dotted line (for x-values: -4, -3, -2, -1, 0),
    #   - The true labels (30 days) as a solid line (for x-values: 1 to 30),
    #   - The predictions (30 days) as a dashed line (for x-values: 1 to 30).
    # Create a figure with 8 subplots (one per compartment)
    # -------------------------
    compartments = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                    'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical',
                    'Recovered', 'Dead']

    # Create subplots (2 rows x 4 columns)
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()

    # Define x-axes:
    # Historical input days: we'll label them as days -4, -3, -2, -1, 0 (5 points)
    x_inputs = np.arange(-4, 1)  # i.e. -4, -3, -2, -1, 0
    # Predictions/labels days: days 1 to 30
    x_pred = np.arange(1, 31)

    for i in range(8):
        ax = axes[i]
        # Plot historical input data as dotted line
        ax.plot(x_inputs, inputs_compartments[:, i],
                label='Input', linestyle=':', linewidth=2, color='gray')
        # Plot true labels as solid line
        ax.plot(x_pred, labels_compartments[:, i],
                label='Labels', linewidth=2, color='blue')
        # Plot predictions as dashed line
        ax.plot(x_pred, pred_compartments[:, i],
                label='Predictions', linestyle='--', linewidth=2, color='red')
        ax.set_title(compartments[i])
        ax.set_xlabel('Day')
        ax.set_ylabel('Value')
        ax.set_yscale('log')
        ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(CWD, 'predictions_vs_labels.png'))
    plt.show()


def plot_germany_initial_conditions(inputs, compartment, num_samples=50, log_scale=False, relative=False, size_tikz=30, size_labels=22, size_title=22, size_legend=22):
    """
    For each of num_samples (default: 50) samples, plots an individual Germany map 
    that shows the state of the specified compartment at day 0.

    The process for each sample:
      - Reverse the log-transformation using np.expm1.
      - Transpose and reshape the data from an initial shape of 
          (NUM_Compartments, INPUT_WIDTH, NUMBER_OF_NODES)
        into a shape of (NUMBER_OF_NODES, INPUT_WIDTH, 6, 8), where the 48 compartments 
        are reshaped into 6 age groups and 8 compartments.
      - Extract day 0 (the first of the INPUT_WIDTH historical days).
      - For the chosen compartment (provided as a string or index), sum over the 6 age groups,
        resulting in a scalar per region.
      - Optionally, if relative is True, the values are normalized by the region's population 
        (obtained from an external JSON file) and scaled per 100,000 inhabitants.
      - Plot the resulting regional values on a Germany map (using a provided shapefile).

    The shapefile is expected to be located at:
      os.path.join(os.getcwd(), 'tools', 'vg2500_12-31.utm32s.shape', 'vg2500', 'VG2500_KRS.shp')

    Parameters:
      inputs: Input data with initial shape (NUM_Compartments, INPUT_WIDTH, NUMBER_OF_NODES).
      compartment: The compartment to be plotted (either as a string or an index).
      num_samples: The number of samples to plot (default 50).
      log_scale: If True, use a logarithmic color scale for the plots.
      relative: If True, normalize the compartment values by the regional population 
                (per 100,000 inhabitants) using external population data.
    """

    # Define the list of compartments
    compartments = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                    'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical',
                    'Recovered', 'Dead']

    # Load population data for Germany to normalize the values (if relative scaling is used)
    df_population = pd.read_json(
        os.path.join(CWD, 'data', 'pydata', 'Germany', "county_current_population.json"))
    df_population = df_population[["ID_County", "Population"]]

    # Determine the compartment index if a string is provided
    if isinstance(compartment, str):
        try:
            compartment_index = compartments.index(compartment)
        except ValueError:
            raise ValueError(
                f"Compartment '{compartment}' not found. Available compartments: {compartments}")
    elif isinstance(compartment, int):
        compartment_index = compartment
        if compartment_index < 0 or compartment_index >= len(compartments):
            raise ValueError(
                f"Compartment index {compartment_index} is invalid. It must be between 0 and {len(compartments)-1}.")
    else:
        raise ValueError(
            "The 'compartment' parameter must be either a string or an integer.")

    compartment_name = compartments[compartment_index]

    # Load the shapefile for Germany
    shapefile_path = os.path.join(
        os.getcwd(), 'tools', 'vg2500_12-31.utm32s.shape', 'vg2500', 'VG2500_KRS.shp')
    if not os.path.exists(shapefile_path):
        raise FileNotFoundError(f"Shapefile not found: {shapefile_path}")
    map_data = gpd.read_file(shapefile_path)
    # Use only the first NUMBER_OF_NODES regions if the shapefile contains more regions
    if len(map_data) > NUMBER_OF_NODES:
        map_data = map_data.iloc[:NUMBER_OF_NODES].copy()
    else:
        map_data = map_data.copy()

    # Select the first num_samples from the dataset
    if inputs.shape[0] < num_samples:
        num_samples = inputs.shape[0]
        print(f"Only {num_samples} samples available in the dataset.")
    selected_samples = inputs[:num_samples]

    # Determine global min and max values for a unified color scale across all samples,
    # applying relative scaling if needed.
    all_values = []
    for sample in selected_samples:
        # Reverse the log transformation
        sample = np.expm1(sample)
        # Initial shape is (48, INPUT_WIDTH, NUMBER_OF_NODES)
        # Transpose to shape (NUMBER_OF_NODES, INPUT_WIDTH, 48)
        sample = sample.transpose(2, 1, 0)
        # Reshape: (NUMBER_OF_NODES, INPUT_WIDTH, 48) -> (NUMBER_OF_NODES, INPUT_WIDTH, 6, 8)
        sample_reshaped = sample.reshape(NUMBER_OF_NODES, INPUT_WIDTH, 6, 8)
        # Extract day 0 (the first historical day)
        day0 = sample_reshaped[:, 0, :, :]  # shape: (NUMBER_OF_NODES, 6, 8)
        # For the chosen compartment, sum over the 6 age groups → shape: (NUMBER_OF_NODES,)
        comp_values = np.sum(day0[:, :, compartment_index], axis=1)
        if relative:
            # Apply relative scaling: normalize by the region's population (per 100,000 inhabitants)
            for idx in range(len(comp_values)):
                region_code = int(map_data.loc[idx, 'ARS'])
                population = df_population.loc[df_population['ID_County']
                                               == region_code, 'Population'].values[0]
                comp_values[idx] = comp_values[idx] / population * 100000
        all_values.append(comp_values)
    global_min = min(np.min(arr) for arr in all_values)
    global_max = max(np.max(arr) for arr in all_values)

    norm = LogNorm(vmin=global_min, vmax=global_max) if log_scale else plt.Normalize(
        vmin=global_min, vmax=global_max)

    # Create a separate plot for each sample
    save_dir = "/localdata1/gnn_paper_2024/images/with_spatial_res/initial_conditions/nodes_with_variance/"
    os.makedirs(save_dir, exist_ok=True)
    for i, sample in enumerate(selected_samples):
        # Reverse the log transformation
        sample = np.expm1(sample)
        # Initial shape is (48, INPUT_WIDTH, NUMBER_OF_NODES)
        # Transpose to shape (NUMBER_OF_NODES, INPUT_WIDTH, 48)
        sample = sample.transpose(2, 1, 0)
        # Reshape: (NUMBER_OF_NODES, INPUT_WIDTH, 48) -> (NUMBER_OF_NODES, INPUT_WIDTH, 6, 8)
        sample_reshaped = sample.reshape(NUMBER_OF_NODES, INPUT_WIDTH, 6, 8)
        # Extract day 0 (the first historical day)
        day0 = sample_reshaped[:, 0, :, :]  # shape: (NUMBER_OF_NODES, 6, 8)
        # For the chosen compartment, sum over the 6 age groups → shape: (NUMBER_OF_NODES,)
        comp_values = np.sum(day0[:, :, compartment_index], axis=1)

        # Create a copy of the GeoDataFrame and add the computed values
        map_data_copy = map_data.copy()
        map_data_copy['comp_value'] = comp_values

        # Normalize the values by regional population if relative is True
        if relative:
            for row in range(len(map_data_copy)):
                region_code = int(map_data_copy.loc[row, 'ARS'])
                population = df_population.loc[df_population['ID_County']
                                               == region_code, 'Population'].values[0]
                map_data_copy.loc[row, 'comp_value'] = map_data_copy.loc[row,
                                                                         'comp_value'] / population * 100000

        # Create the plot without colorbar
        fig, ax = plt.subplots(figsize=(10, 10))
        map_data_copy.plot(column='comp_value', ax=ax,
                           cmap='plasma', edgecolor='black', linewidth=0.5)
        ax.set_axis_off()
        plt.tight_layout()
        file_path = os.path.join(
            save_dir, f'initial_condition_sample_{i+1}_{compartment_name}.png')
        plt.savefig(file_path)
        plt.close(fig)

    # Create a vertical colorbar figure:
    sm = plt.cm.ScalarMappable(cmap='plasma', norm=norm)
    sm.set_array([])  # Dummy array for ScalarMappable
    fig_cb_vert = plt.figure(figsize=(2, 8))
    # Adjust position and size as needed
    ax_cb_vert = fig_cb_vert.add_axes([0.05, 0.05, 0.2, 0.9])
    cbar_vert = plt.colorbar(sm, cax=ax_cb_vert)
    cbar_vert.set_label(
        'Number of Infected Symptomatic per 100k', fontsize=size_labels)
    # Ensure all ticks are in the correct size
    cbar_vert.ax.tick_params(labelsize=size_tikz)
    for label in cbar_vert.ax.get_yticklabels():
        label.set_fontsize(size_tikz)
    fig_cb_vert.savefig(os.path.join(
        save_dir, f'colorbar_vertical_{compartment_name}.png'), dpi=300)
    plt.close(fig_cb_vert)

    # Create a horizontal colorbar figure:
    sm = plt.cm.ScalarMappable(cmap='plasma', norm=norm)
    sm.set_array([])  # Dummy array for ScalarMappable
    fig_cb_horiz = plt.figure(figsize=(8, 2))
    # Adjust position and size as needed
    ax_cb_horiz = fig_cb_horiz.add_axes([0.05, 0.4, 0.9, 0.2])
    cbar_horiz = plt.colorbar(sm, cax=ax_cb_horiz, orientation='horizontal')
    cbar_horiz.set_label(
        'Number of Infected Symptomatic per 100k', fontsize=size_labels)
    # increase tikz size
    cbar_horiz.ax.tick_params(labelsize=size_tikz)
    fig_cb_horiz.savefig(os.path.join(
        save_dir, f'colorbar_horizontal_{compartment_name}.png'), dpi=300)
    plt.close(fig_cb_horiz)


def main():
    # Create the adjacency matrix and compute the normalized, scaled Laplacian.
    adjacency_matrix = get_adjacency_matrix(
        NUMBER_OF_NODES, COMMUTER_FILE_PATH)
    a = rescale_laplacian(normalized_laplacian(adjacency_matrix))
    a = a.astype('float32')

    # Load the model including its weights.
    model = load_model_and_weights(a)

    # --- Sample Prediction using prepared input data ---
    # print("\n--- Sample Prediction using prepared input data ---")
    # x_sample = prepare_input_data()
    # a_sample = np.repeat(np.expand_dims(a, axis=0), x_sample.shape[0], axis=0)
    # sample_prediction = model([x_sample, a_sample])
    # print("Sample prediction (Shape):", sample_prediction.shape)
    # print("Sample prediction (Raw values):")
    # print(sample_prediction.numpy())

    # # --- Load dataset and compute metrics ---
    # print("\n--- Loading dataset and computing metrics ---")
    test_inputs, test_labels = load_dataset(DAY)
    mean_mape_per_region, mean_mape_reversed_per_region = compute_metrics(
        model, test_inputs, test_labels, a)

    # # --- Plot the mean MAPE per region on a map ---
    plot_map(mean_mape_per_region, mean_mape_reversed_per_region)

    # --- Predict and plot for a single dataset ---
    # real_dataset_path = "/localdata1/zunk_he/memilio/saves/eextrapolated_GNN_data_30days_nodeswithvariance_single.pickle"
    # predict_and_plot(model, real_dataset_path, a)

    # ---Plot initial conditions for some samples (day 0) using the prepared inputs ---
    # plot_germany_initial_conditions(
    #     test_inputs, compartment="InfectedSymptoms", num_samples=2, log_scale=True, relative=True)


if __name__ == '__main__':
    main()
