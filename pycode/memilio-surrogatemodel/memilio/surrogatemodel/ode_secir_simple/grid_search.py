# from memilio.surrogatemodel.ode_secir_simple import network_architectures
# from memilio.surrogatemodel.ode_secir_simple.model import split_data, get_test_statistic
import os
import tensorflow as tf
import pickle
import pandas as pd
import time
from sklearn.model_selection import KFold
import numpy as np
import memilio.surrogatemodel.ode_secir_simple.network_architectures as architectures


label_width = 30
num_outputs = 8

# Set up the data path 
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
    os.path.dirname(os.path.realpath(path)))), 'data')

filename = "data_secir_simple_30days_10k.pickle"
filename_df = "dataframe_30days_10k_nodamp.csv"

# General parameters for the training process:
early_stop = [100]
max_epochs = [2]
losses=[tf.keras.losses.MeanAbsolutePercentageError()]
optimizers= ["AdamW"]
metrics=[[tf.keras.metrics.MeanAbsoluteError(), tf.keras.metrics.MeanAbsolutePercentageError()]]

# Define grid search parameters
hidden_layers = [1,2]
neurons_in_hidden_layer = [32]
activation_function = ["relu"]
models = ["CNN" ]

# Collecting parameters
training_parameters = [(early, epochs, loss, optimizer, metric) 
                       for early in early_stop for epochs in max_epochs for loss in losses 
                       for optimizer in optimizers for metric in metrics]

model_parameters = [(label_width, num_outputs, layer, neuron_number, activation, modelname)
              for layer in hidden_layers for neuron_number in neurons_in_hidden_layer 
              for activation in activation_function for modelname in models]

### Preparing data 
# Loading data 
if not os.path.isfile(os.path.join(path_data, filename)):
    raise ValueError(f"No dataset found in path: {path_data}")

with open(os.path.join(path_data, filename), 'rb') as file:
    data = pickle.load(file)

# Split the data: 80% will be used for grid search with cross-validation, and the remaining 20% is withheld for testing
inputs_grid_search = data['inputs'][:int((0.8 * len(data['inputs'])))]
labels_grid_search = data['labels'][:int((0.8 * len(data['labels'])))]
inputs_withhold = data['inputs'][int((0.8 * len(data['inputs']))):]
labels_withhold = data['labels'][int((0.8 * len(data['labels']))):]

# Function to train and evaluate the model using cross-validation
def train_and_evaluate_model(param, inputs, labels, training_parameter, Print = False):
    """ Training and evaluating a model with given architecture using 5-fold cross validation, returning a dictionary with the main training statistics. 

    :param param: tuple of parameters describing the model architecture, it should be of the form 
            (num_days_per_output, num_outputs, num_hidden_layers, neurons_per_layer, name_activation, name_architecture)
    :param inputs: training inputs 
    :param labels: training output labels 
    :param training_parameter: tuple of parameters used for the training process, it should be of the form
        (early_stop, max_epochs, loss, optimizer, metrics), where loss is a loss-function implemented in keras, optimizer is the name of the used optimizer, 
        metrics is a list of used training metrics, e.g. [tf.keras.metrics.MeanAbsoluteError(), tf.keras.metrics.MeanAbsolutePercentageError()]
    :param Print:  Boolean, whether or not the evaluation results are printed. 
    :returns: a dictionary of training statistics of the form 
        {"model", "activation","optimizer","mean_train_loss_kfold","mean_val_loss_kfold","training_time", "train_losses", "val_losses"}
    
    """
    # Unpacking parameters 
    _, _, _, _, activation, modelname = param
    early_stop, max_epochs, loss, optimizer, metrics = training_parameter
    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                      patience=early_stop,
                                                      mode='min')
    
    # Preparing K-Fold Cross-Validation
    kf = KFold(n_splits=5)
    train_losses = []
    val_losses = []

    losses_history_all = []
    val_losses_history_all = []

    start = time.perf_counter()
    for train_idx, val_idx in kf.split(inputs):
        # Clearing any information about the previous model 
        tf.keras.backend.clear_session()

        # Gather training and validation data based on the fold
        train_inputs = tf.gather(inputs, indices=train_idx)
        train_labels = tf.gather(labels, indices=train_idx)
        valid_inputs = tf.gather(inputs, indices=val_idx)
        valid_labels = tf.gather(labels, indices=val_idx)

        # Initializing model
        model = initialize_model(param)
        # Compile the model
        model.compile(loss=loss,
                      optimizer=optimizer,
                      metrics=metrics)
        
        # Train the model
        history = model.fit(train_inputs, train_labels, epochs=max_epochs,
                            validation_data=(valid_inputs, valid_labels),
                            callbacks=[early_stopping])

        train_losses.append(np.min(history.history['loss']))
        val_losses.append(np.min(history.history['val_loss']))
        losses_history_all.append(history.history['loss'])
        val_losses_history_all.append(history.history['val_loss'])


    elapsed = time.perf_counter() - start

    # Print out the results
    if Print: 
        elapsed = time.perf_counter() - start
        print(f"Best train losses: {train_losses}")
        print(f"Best validation losses: {val_losses}")
        print("--------------------------------------------")
        print(f"K-Fold Train Score: {np.mean(train_losses)}")
        print(f"K-Fold Validation Score: {np.mean(val_losses)}")
        print(f"Time for training: {elapsed:.4f} seconds")
        print(f"Time for training: {elapsed / 60:.4f} minutes")

    # After cross-validation, we can test on the withhold dataset (outside of the loop) ? 
    return {
        "model": modelname, 
        "activation": activation,
        "optimizer": optimizer,
        "mean_train_loss_kfold": np.mean(train_losses), 
        "mean_val_loss_kfold": np.mean(val_losses), 
        "training_time": elapsed/60,
        "train_losses": [losses_history_all],
        "val_losses": [val_losses_history_all]
    }



def initialize_model(parameters): 
    """ Initialize model from given list of parameters

    :param parameters: tuple of parameters describing the model architecture, it should be of the form 
            (num_days_per_output, num_outputs, num_hidden_layers, neurons_per_layer, name_activation, name_architecture)
    :returns: tensor flow keras model with the given architecture 
     """
    label_width, num_outputs, layer, neuron_number, activation, modelname = parameters
    if modelname == "Dense":
        return architectures.mlp_multi_input_multi_output(
            label_width = label_width, 
            num_outputs = num_outputs,
            num_hidden_layers = layer, 
            num_neurons_per_layer = neuron_number,
            activation = activation
        )
    elif modelname == "LSTM":
        return architectures.lstm_multi_input_multi_output(
            label_width = label_width, 
            num_outputs = num_outputs,
            num_hidden_layers = layer, 
            num_neurons_per_layer = neuron_number,
            activation = activation
        )
    elif modelname == "CNN":
        return architectures.cnn_multi_input_multi_output(
            label_width = label_width, 
            num_outputs = num_outputs,
            num_hidden_layers = layer, 
            num_neurons_per_layer = neuron_number,
            activation = activation
        )
    
def perform_grid_search(model_parameters, inputs, labels, training_parameters, filename_df):
    """ Performing grid search for a given set of model parameters

    The results are stored in directory 'secir_simple_grid_search', each row has the form 
    ['model', 'optimizer', 'number_of_hidden_layers', 'number_of_neurons', 'activation',
                                    'mean_test_MAPE', 'kfold_train', 'kfold_val',
                                    'kfold_test', 'training_time', 'train_losses', 'val_losses']

    :param model_parameters: List of tuples of model parameters, each entry should be a tuple of the form 
        (num_days_per_output, num_outputs, num_hidden_layers, neurons_per_layer, name_activation, name_architecture)
    :param inputs: training input data 
    :param labels: training label data 
    :param training_parameters: List of tuples of parameters used for the training process, each should be of the form
        (early_stop, max_epochs, loss, optimizer, metrics), where loss is a loss-function implemented in keras, optimizer is the name of the used optimizer, 
        metrics is a list of used training metrics, e.g. [tf.keras.metrics.MeanAbsoluteError(), tf.keras.metrics.MeanAbsolutePercentageError()]
    :param filename_df: String, giving name of the file, where the data is stored, actual filename is given by filename_df + ".csv"
    
    """
    # Create a DataFrame to store the results
    df_results = pd.DataFrame(columns=['model','optimizer', 'number_of_hidden_layers', 'number_of_neurons', 'activation',
                                    'mean_test_MAPE', 'kfold_train', 'kfold_val',
                                    'kfold_test', 'training_time', 'train_losses', 'val_losses'])
    
    # Iterate the different model architectures and save the training results
    for param in model_parameters:
        for training_parameter in training_parameters:
            _,_, layer, neuron_number, activation, modelname = param
            results = train_and_evaluate_model(param, inputs, labels, training_parameter)
            df_results.loc[len(df_results.index)] = [
                modelname, results["optimizer"], layer, neuron_number, activation, np.nan,  # Placeholder for test score
                results["mean_train_loss_kfold"],
                results["mean_val_loss_kfold"],
                np.nan,  # Placeholder for final test score
                results["training_time"],
                results["train_losses"],
                results["val_losses"]
            ]

    # Save the results in .csv file 
    path = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(os.path.dirname(
        os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'secir_simple_grid_search')
    if not os.path.isdir(file_path):
        os.mkdir(file_path)
    file_path = os.path.join(file_path, filename_df)
    df_results.to_csv(file_path)



# Performing grid search 
start_hyper = time.perf_counter()

perform_grid_search(model_parameters, inputs_grid_search, labels_grid_search, training_parameters, filename_df)

elapsed_hyper = time.perf_counter() - start_hyper
print(
    "Time for hyperparameter testing: {:.4f} minutes".format(
        elapsed_hyper / 60))
print(
    "Time for hyperparameter testing: {:.4f} hours".format(
        elapsed_hyper / 60 / 60))