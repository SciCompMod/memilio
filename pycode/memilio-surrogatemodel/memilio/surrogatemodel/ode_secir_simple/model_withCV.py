#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
#
# Authors: Agatha Schmidt, Henrik Zunker
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
import os
import pickle
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.model_selection import KFold

from memilio.simulation.secir import InfectionState
from memilio.surrogatemodel.ode_secir_simple import network_architectures


def network_fit(path,filename, modelname,  max_epochs=30, early_stop=100, plot=False):
    """! Training and evaluation of a given model with mean squared error loss and Adam optimizer using the mean absolute error as a metric.

    @param path path of the dataset. 
    @param model Keras sequential model.
    @param max_epochs int maximum number of epochs in training. 
    @param early_stop Integer that forces an early stop of training if the given number of epochs does not give a significant reduction of validation loss. 

    """

    if not os.path.isfile(os.path.join(path, 'data_secir_simple.pickle')):
        ValueError("no dataset found in path: " + path)

    file = open(os.path.join(path,filename), 'rb')

    data = pickle.load(file)

    df_results  = pd.DataFrame(
    columns=['model', 'mean_test_MAPE', 'kfold_train',
             'kfold_val', 'kfold_test', 'training_time'])     

    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                      patience=early_stop,
                                                      mode='min')
    

    kf = KFold(n_splits=5)
    train_idxs = []
    test_idxs = []
    for i, (train_index, test_index) in enumerate(kf.split(data['inputs'])):

        train_idxs.append(train_index)
        test_idxs.append(test_index)

    test_scores = []
    train_losses = []
    val_losses = []

    losses_history_all = []
    val_losses_history_all = []

    start = time.perf_counter()
    for train_idx, test_idx in zip(train_idxs, test_idxs):

        model = network_architectures.lstm_multi_input_multi_output(label_width)
        
        train_inputs = tf.gather(data['inputs'], indices = train_idx[:(int(0.8*len(train_idx)))])
        train_labels = tf.gather(data['labels'], indices = train_idx[:(int(0.8*len(train_idx)))])
        valid_inputs = tf.gather(data['inputs'], indices = train_idx[(int(0.8*len(train_idx))):])
        valid_labels = tf.gather(data['labels'], indices = train_idx[(int(0.8*len(train_idx))):])
        test_inputs = tf.gather(data['inputs'], indices = test_idx)
        test_labels = tf.gather(data['labels'], test_idx)


        model.compile(
            loss=tf.keras.losses.MeanAbsolutePercentageError(),
            optimizer=tf.keras.optimizers.Adam(),
            metrics=[tf.keras.metrics.MeanAbsoluteError()])

        history = model.fit(train_inputs, train_labels, epochs=max_epochs,
                            validation_data=(valid_inputs, valid_labels),
                            callbacks=[early_stopping])
        
        df = get_test_statistic(test_inputs, test_labels, model)
        print(df)
        print('mean: ',  df.mean())

        test_scores.append(df.mean()[0])
        train_losses.append(np.asarray(history.history['loss']).min())
        val_losses.append(np.asarray(history.history['val_loss']).min())
        losses_history_all.append(history.history['loss'])
        val_losses_history_all.append(history.history['val_loss'])

    # print out stats
    elapsed = time.perf_counter() - start
    print("Best train losses: {} ".format(train_losses))
    print("Best validation losses: {}".format(val_losses))
    print("Test values: {}".format(test_scores))
    print("--------------------------------------------")
    print("K-Fold Train Score:{}".format(np.mean(train_losses)))
    print("K-Fold Validation Score:{}".format(np.mean(val_losses)))
    print("K-Fold Test Score: {}".format(np.mean(test_scores)))

    print("Time for training: {:.4f} seconds".format(elapsed))
    print("Time for training: {:.4f} minutes".format(elapsed/60))


        

    
    # save the model
    # path = os.path.dirname(os.path.realpath(__file__))
    # path_models = os.path.join(
    #     os.path.dirname(
    #         os.path.realpath(os.path.dirname(os.path.realpath(path)))),
    #     'saved_models_secir_simple_30days_baseline')
    # if not os.path.isdir(path_models):
    #     os.mkdir(path_models)

    # model.save(path_models, 'LSTM.h5')


    
    df_results.loc[len(df_results.index)] = [modelname, df.mean()[0] , np.mean(train_losses),
                             np.mean(val_losses),
                             np.mean(test_scores),
                             (elapsed / 60)]

    
    path = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(
            os.path.dirname(
                os.path.realpath(os.path.dirname(os.path.realpath(path)))),
            'secir_simple_dataframes')
    if not os.path.isdir(file_path):
            os.mkdir(file_path)
    file_path = file_path+'secir_simple_baseline'+modelname
    df_results.to_csv(file_path)



    return history


def plot_losses(history):
    """! Plots the losses of the model training.  

    @param history model training history. 

    """
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('model loss')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper left')
    if os.path.isdir("plots") == False:
        os.mkdir("plots")
    plt.savefig('plots/losses_plot.png')
    plt.show()


def get_test_statistic(test_inputs, test_labels, model):
    """! Calculates the mean absolute percentage error based on the test dataset.   

    @param test_inputs inputs from test data.
    @param test_labels labels (output) from test data.
    @param model trained model. 

    """

    pred = model(test_inputs)
    pred = pred.numpy()
    test_labels = np.array(test_labels)

    diff = pred - test_labels
    relative_err = (abs(diff))/abs(test_labels)
    # reshape [batch, time, features] -> [features, time * batch]
    relative_err_transformed = relative_err.transpose(2, 0, 1).reshape(8, -1)
    relative_err_means_percentage = relative_err_transformed.mean(axis=1) * 100
        
    # delete the two confirmed compartments from InfectionStates
    compartment_array = []
    for compartment in InfectionState.values():
        compartment_array.append(compartment) 
    index = [3,5]
    compartments_cleaned= np.delete(compartment_array, index)

    mean_percentage = pd.DataFrame(
        data=relative_err_means_percentage,
        index=[str(compartment).split('.')[1]
               for compartment in compartment_array],
        columns=['Percentage Error'])

    return mean_percentage


def split_data(inputs, labels, split_train=0.7,
               split_valid=0.2, split_test=0.1):
    """! Split data set in training, validation and testing data sets.

   @param inputs input dataset
   @param labels label dataset
   @param split_train Share of training data sets.
   @param split_valid Share of validation data sets.
   @param split_test Share of testing data sets.
   """

    if split_train + split_valid + split_test > 1 + 1e-10:
        raise ValueError(
            "Summed data set shares are greater than 1. Please adjust the values.")
    elif inputs.shape[0] != labels.shape[0] or inputs.shape[2] != labels.shape[2]:
        raise ValueError(
            "Number of batches or features different for input and labels")

    n = inputs.shape[0]
    n_train = int(n * split_train)
    n_valid = int(n * split_valid)
    n_test = n - n_train - n_valid

    inputs_train, inputs_valid, inputs_test = tf.split(
        inputs, [n_train, n_valid, n_test], 0)
    labels_train, labels_valid, labels_test = tf.split(
        labels, [n_train, n_valid, n_test], 0)

    data = {
        'train_inputs': inputs_train,
        'train_labels': labels_train,
        'valid_inputs': inputs_valid,
        'valid_labels': labels_valid,
        'test_inputs': inputs_test,
        'test_labels': labels_test
    }

    return data


if __name__ == "__main__":
    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')
    
    filename = "data_secir_simple.pickle"
    max_epochs = 1500
    label_width = 30 

    model_output = network_fit(
        path_data, filename,  modelname  = 'LSTM_30',
        max_epochs=max_epochs)
