import numpy as np
import pandas as pd
from datetime import date
from math import ceil
import random
import os
from progress.bar import Bar # pip install progess
from sklearn.model_selection import train_test_split
from sklearn.compose import make_column_transformer
from sklearn.preprocessing import MinMaxScaler
import tensorflow as tf
import matplotlib.pyplot as plt
from tensorflow import keras
from keras import Sequential
from keras.layers import Dense
from keras.optimizers import Adam
from keras import backend as K

import seaborn as sns # plot after normalization
from data_secir_simple import generate_data






def splitdata(data, split_train=0.7, 
              split_valid=0.2, split_test=0.1,
              normalize=False):
    if split_train + split_valid + split_test != 1:
        ValueError("Split ratio not equal 1! Please adjust the values")
    n = len(data)
    train_data = data[0:int(n*split_train)]
    val_data = data[int(n*split_train):int(n*(1-split_test))]
    test_data = data[int(n*(1-split_test)):]

    return train_data, val_data, test_data






### normalization is done now when generating the data 


# for input_width = 5 and label_width = 20
def plotCol_Dense(inputs, labels, model=None, plot_col='Susceptible', max_subplots=3):

    input_width = 5
    label_width = 20
    plt.figure(figsize=(12, 8))
    cols = np.array(['Susceptible', 'Exposed', 'Carrier',
                        'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead'])  
    plot_col_index = np.where(cols == plot_col)[0][0]
    max_n = min(max_subplots, len(inputs)/100)
    for n in range(max_n):
        plt.subplot(max_n, 1, n+1)
        plt.ylabel(f'{plot_col}')
        # plt.plot(np.arange(1,6), inputs.iloc[0:5, plot_col_index],
        #         label='Inputs', marker='.', zorder=-10)

        plt.plot(np.arange(0,19), inputs.iloc[0:19, plot_col_index],
                label='Inputs', marker='.', zorder=-10)

        # plt.scatter(np.arange(6,26), labels.iloc[0:20, plot_col_index],
        #             edgecolors='k', label='Labels', c='#2ca02c', s=64)
            
        plt.scatter(np.arange(1,20), labels.iloc[(n*100):(n*100)+19, plot_col_index],
                    edgecolors='k', label='Labels', c='#2ca02c', s=64)
        
        if model is not None:
            predictions = model(inputs.to_numpy())
            # plt.scatter(np.arange(6,26), predictions[0:20],
            #         marker='X', edgecolors='k', label='Predictions',
            #         c='#ff7f0e', s=64)

            plt.scatter(np.arange(1,20), predictions[(n*100):(n*100) + 19, plot_col_index],
                    marker='X', edgecolors='k', label='Predictions',
                    c='#ff7f0e', s=64)

    plt.xlabel('days')
    plt.show()

def plotCol(inputs, inputs_scaled, labels, model=None, plot_col='Susceptible', max_subplots=3):

    input_width = 4
    label_width = 1
    
    plt.figure(figsize=(12, 8))
    cols = np.array(['Susceptible', 'Exposed', 'Carrier',
                        'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead'])  
    plot_col_index = np.where(cols == plot_col)[0][0]
    max_n = min(max_subplots, len(inputs)/100)

    # predictions = model(inputs) # for just one input: input_series = tf.expand_dims(inputs[n], axis=0) -> model(input_series)
    
    for n in range(max_n):
        plt.subplot(max_n, 1, n+1)
        plt.ylabel(f'{plot_col}')
        
        # plt.plot(np.arange(1,6), inputs.iloc[0:5, plot_col_index],
        #         label='Inputs', marker='.', zorder=-10)
        input_array = inputs[n * 20].numpy()
        label_array = labels[n * 20].numpy()
        plt.plot(np.arange(0,input_width), input_array[0:input_width, plot_col_index],
                label='Inputs', marker='.', zorder=-10)

        # plt.scatter(np.arange(6,26), labels.iloc[0:20, plot_col_index],
        #             edgecolors='k', label='Labels', c='#2ca02c', s=64)
            
        plt.scatter(np.arange(input_width,input_width+label_width), label_array[0:label_width,plot_col_index],
                    edgecolors='k', label='Labels', c='#2ca02c', s=64)
        
        if model is not None:
            input_series = tf.expand_dims(inputs_scaled[n * 20], axis=0)
            pred = model(input_series)
            pred = pred.numpy()

            input_series


            if (np.arange(input_width,input_width+label_width)).shape != pred[0,:,plot_col_index].shape[0]:
                plt.scatter(np.arange(input_width,input_width+label_width), 
                                    pred[0,-1,plot_col_index], # :
                                    marker='X', edgecolors='k', label='Predictions',
                                    c='#ff7f0e', s=64)
            else:
                plt.scatter(np.arange(2,input_width+1), 
                                    pred[0,1:,plot_col_index],
                                    marker='X', edgecolors='k', label='Predictions',
                                    c='#ff7f0e', s=64)

    plt.xlabel('days')
    plt.show()
    plt.savefig('evaluation_secir_simple_' + plot_col + '.pdf')

def visualize_predictions(
               testlabels,
               predictions,
               save_evaluation_pdf,
               path_data):
    """! Visualize the predictions of the trained model to evaluate the performance.

    @param testdata Testdata thats different than the traindata and not known for the model.
    @param testlabels Same as testdata but with labels.
    @param predictions Output of the model with testdata as input
    """
    labels = testlabels.to_numpy()
    titel = np.array(['Susceptible', 'Exposed', 'Carrier',
                    'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead'])              
    x = list(range(1, 101))

    for i in range(0,titel.size):
            plt.plot(x,labels[:,i] , "b" , label="Secir simple")
            plt.plot(x,predictions[:,i], "-r", label="Model prediction")
            plt.xlabel('Days')
            plt.ylabel('population')
            plt.title(titel[i])
            plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=2, fancybox=True, shadow=True)

            # max error in %
            err = abs(labels[:,i] - predictions[:,i])
            indx_max_error = np.where(err == np.amax(err))[0][0]
            if predictions[indx_max_error,i] > 0:
                e = err.max() / predictions[indx_max_error,i] * 100
            else:
                raise ValueError("Population values of zero or lower are outside the defined area!")

            # Display the maximal error , rounded to the 4st decimal, under the actual plot
            plt.figtext(0.5, 0.01, 'Maximal Error in "%%" = (%s)'%(format(e, '.4f')), 
                ha="center", fontsize=15, bbox={"facecolor":"orange", "alpha":0.5, "pad":5})

            if save_evaluation_pdf:
                path_save = os.path.join(path_data, 'evaluation_secir_simple')
                if not os.path.isdir(path_save):
                    os.mkdir(path_save)
                plt.savefig(os.path.join(path_save,'evaluation_secir_simple_' + titel[i] + '.pdf')) 
            # only show plot for Susceptible compartment. Others are saved
            if not i :
                plt.show()   

            # clean currect figure
            plt.clf() 
    return None

# create and train an neural network based on the secir simple example. 
# if no dataset is already create, we build one with default 500 runs
def network_secir_simple(path, epochs=30, num_runs_traindata=500, save_evaluation_pdf=False, path_save=""):
    """! Generate the model and train with the created dataset.

    If the dataset is not created yet, we create a new one with default value
    'nums_runs_traindata = 50'. Since the number of training runs is freely selectable with
    'epochs', unit tests are easier to implement afterwards.

   @param path Path to dataset
   @param epochs Number of epochs in the train process. 
   @param number of runs to create dataset. Only used, if no dataset is found.
   """

    if not os.path.isfile(os.path.join(path, 'traindata_secir_simple.txt')) or \
        not os.path.isfile(os.path.join(path, 'traindata_secir_simple.txt')):

        generate_data(num_runs_traindata, path)

    input = pd.read_csv(
            os.path.join(os.path.dirname(os.path.realpath(path)), 'data35', 'traindata_secir_simple.txt'),
            sep=' ')
    
    labels = pd.read_csv(
            os.path.join(os.path.dirname(os.path.realpath(path)), 'data35', 'labels_secir_simple.txt'),
            sep=' ')

    # split in train and test/valid data (ratio 80/20).
    # shuffling is also done by this function
    # TODO: Also valid data needed
    X_train, X_test, Y_train, Y_test = train_test_split(input, labels, test_size=0.2, random_state=42)

    # define transformation to prepare data for training. Mainly min-max scaling, OneHotEncoder for classification params
    transformer = make_column_transformer(
    (MinMaxScaler(), 
        ['Susceptible', 'Exposed', 'Carrier',
                    'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead', 'Damping Factor']))

    transformer.fit(X_train)
    X_train = transformer.transform(X_train)
    X_test = transformer.transform(X_test)

    # TODO: Plot afterwards

    # Regarding the network, we use a fully connected one with several hidden layers.
    # Since we have a regression problem, we use the mse loss function.
    # use the leakyRelu as activation function, except in the output layer.
    
    # set seed for pseudo-randomization
    tf.random.set_seed(42)

    # we want to predict the population in the diferent compartments. 
    # Todo: Maybe add drop out or other regularization against overfitting.
    model = Sequential([
        Dense(128, activation=tf.keras.layers.LeakyReLU(alpha=0.01)),
        Dense(128, activation=tf.keras.layers.LeakyReLU(alpha=0.01)),
        Dense(128, activation=tf.keras.layers.LeakyReLU(alpha=0.01)),
        Dense(64, activation=tf.keras.layers.LeakyReLU(alpha=0.01)),
        Dense(8)
    ])

    model.compile(
        loss=tf.keras.losses.mse,
        optimizer=Adam(),
        metrics=['mae']
    )

    model.fit(X_train, Y_train, epochs=epochs)

    model.evaluate(X_test, Y_test)

    # Use the orderer dataset to evaluate the performance of the trained model
    X_test = transformer.transform(input.tail(100))
    preds = model.predict(X_test)
    visualize_predictions(labels.tail(100), preds, save_evaluation_pdf, path_data=path)


# create and train an neural network based on the secir simple example. 
# if no dataset is already create, we build one with default 500 runs
def ml_network(path,  max_epochs=30, early_stop=4, num_runs_traindata=500, save_evaluation_pdf=False):

    input_width = 1
    label_width = 1

    if not os.path.isfile(os.path.join(path, 'data_secir_simple.txt')):
        generate_data(num_runs_traindata, path)

    data = pd.read_csv(
        os.path.join(os.path.dirname(os.path.realpath(path)), 'data35', 'data_secir_simple.txt'),
        sep=' ')

    train_data, valid_data, test_data = splitdata(data, normalize=False)

    train_inputs, train_labels = getPartitions(input_width = input_width, label_width = label_width, data=train_data)
    valid_inputs, valid_labels = getPartitions(input_width = input_width, label_width = label_width, data=valid_data)
    test_inputs, test_labels = getPartitions(input_width = input_width, label_width = label_width, data=test_data)


    train_inputs, valid_inputs, test_inputs = normalizeData(train_inputs, 
                                                    valid_inputs,
                                                    test_inputs,
                                                    data)

    model = tf.keras.Sequential([
                # tf.keras.layers.Flatten(),
                tf.keras.layers.Dense(units=32, activation='relu'),
                tf.keras.layers.Dense(units=32, activation='relu'),
                tf.keras.layers.Dense(units=train_labels.shape[1])])


    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                    patience=early_stop,
                                                    mode='min')

    model.compile(loss=tf.keras.losses.MeanSquaredError(),
                optimizer=tf.keras.optimizers.Adam(),
                metrics=[tf.keras.metrics.MeanAbsoluteError()])

    history = model.fit(train_inputs, train_labels, epochs=max_epochs,
                      validation_data=(valid_inputs, valid_labels),
                      callbacks=[early_stopping])
    
    plotCol_Dense(test_inputs, test_labels, model=model, plot_col='Susceptible', max_subplots=3)

    return history

def ml_network_multi_input(path,  max_epochs=30, early_stop=4, num_runs_traindata=500, save_evaluation_pdf=False):

    input_width = 5
    label_width = 1

    if not os.path.isfile(os.path.join(path, 'data_secir_simple.txt')):
        generate_data(num_runs_traindata, path)

    data = pd.read_csv(
        os.path.join(os.path.dirname(os.path.realpath(path)), 'data35', 'data_secir_simple.txt'),
        sep=' ')

    # normalization is kind of difficult for data nested in lists. So start with normalize all data
    # and then repeat progress again w/o normalization to get the correct labels
    data_mean = data.mean()
    data_std = data.std()
    data = (data - data_mean) / data_std

    train_data, valid_data, test_data = splitdata(data, normalize=False)

    train_inputs, train_labels = getPartitions(input_width = input_width, label_width = label_width, data=train_data)
    valid_inputs, valid_labels = getPartitions(input_width = input_width, label_width = label_width, data=valid_data)
    test_inputs, test_labels = getPartitions(input_width = input_width, label_width = label_width, data=test_data)

    # data = pd.read_csv(
    #     os.path.join(os.path.dirname(os.path.realpath(path)), 'data', 'data_secir_simple.txt'),
    #     sep=' ')

    # againt w/o normalization for the labels

    # train_data, valid_data, test_data = splitdata(data, normalize=False)

    # _, train_labels = getPartitions(input_width = input_width, label_width = label_width, data=train_data)
    # _, valid_labels = getPartitions(input_width = input_width, label_width = label_width, data=valid_data)
    # test_inputs_no_scale, test_labels = getPartitions(input_width = input_width, label_width = label_width, data=test_data)

    # convert list to tensors
    train_inputs = tf.stack(train_inputs)
    train_labels = tf.stack(train_labels)
    valid_inputs = tf.stack(valid_inputs)
    valid_labels = tf.stack(valid_labels)
    test_inputs = tf.stack(test_inputs)
    test_labels = tf.stack(test_labels)

    model = tf.keras.Sequential([
                tf.keras.layers.Flatten(),
                tf.keras.layers.Dense(units=32, activation='relu'),
                tf.keras.layers.Dense(units=32, activation='relu'),
                tf.keras.layers.Dense(units=train_labels[0].shape[1]),
                tf.keras.layers.Reshape([1, -1]),])


    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                    patience=early_stop,
                                                    mode='min')

    model.compile(loss=tf.keras.losses.MeanSquaredError(),
                optimizer=tf.keras.optimizers.Adam(),
                metrics=[tf.keras.metrics.MeanAbsoluteError()])

    history = model.fit(train_inputs, train_labels, epochs=max_epochs,
                      validation_data=(valid_inputs, valid_labels),
                      callbacks=[early_stopping])
    
    
    # plotCol(test_inputs, test_inputs, test_labels, model=model, plot_col='Susceptible', max_subplots=3)

    return history

def lstm_network_multi_input(path,  max_epochs=30, early_stop=4, num_runs_traindata=500, save_evaluation_pdf=False):

    input_width = 4
    label_width = 1

    if not os.path.isfile(os.path.join(path, 'data_secir_simple.txt')):
        generate_data(num_runs_traindata, path)

    data = pd.read_csv(
        os.path.join(os.path.dirname(os.path.realpath(path)), 'data35', 'data_secir_simple.txt'),
        sep=' ')

    # normalization is kind of difficult for data nested in lists. So start with normalize all data
    # and then repeat progress again w/o normalization to get the correct labels
    data_mean = data.mean()
    data_std = data.std()
    data = (data - data_mean) / data_std

    train_data, valid_data, test_data = splitdata(data, normalize=False)

    train_inputs, train_labels = getPartitions(input_width = input_width, label_width = label_width, data=train_data)
    valid_inputs, valid_labels = getPartitions(input_width = input_width, label_width = label_width, data=valid_data)
    test_inputs, test_labels = getPartitions(input_width = input_width, label_width = label_width, data=test_data)

    # data = pd.read_csv(
    #     os.path.join(os.path.dirname(os.path.realpath(path)), 'data', 'data_secir_simple.txt'),
    #     sep=' ')

    # againt w/o normalization for the labels

    # train_data, valid_data, test_data = splitdata(data, normalize=False)

    # _, train_labels = getPartitions(input_width = input_width, label_width = label_width, data=train_data)
    # _, valid_labels = getPartitions(input_width = input_width, label_width = label_width, data=valid_data)
    # test_inputs_no_scale, test_labels = getPartitions(input_width = input_width, label_width = label_width, data=test_data)

    # convert list to tensors
    train_inputs = tf.stack(train_inputs)
    train_labels = tf.stack(train_labels)
    valid_inputs = tf.stack(valid_inputs)
    valid_labels = tf.stack(valid_labels)
    test_inputs = tf.stack(test_inputs)
    test_labels = tf.stack(test_labels)
    # test_inputs_no_scale = tf.stack(test_inputs_no_scale)

    if input_width == 1:
        train_inputs = tf.expand_dims(train_inputs, axis=2)
        train_labels = tf.expand_dims(train_labels, axis=2)
        valid_inputs = tf.expand_dims(valid_inputs, axis=2)
        valid_labels = tf.expand_dims(valid_labels, axis=2)
        test_inputs = tf.expand_dims(test_inputs, axis=2)
        test_labels = tf.expand_dims(test_labels, axis=2)


    model = tf.keras.models.Sequential([
            # Shape [batch, time, features] => [batch, time, lstm_units]
            tf.keras.layers.LSTM(32, return_sequences=True),
            # Shape => [batch, time, features]
            tf.keras.layers.Dense(units=8)
        ])



    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                    patience=early_stop,
                                                    mode='min')

    model.compile(loss=tf.keras.losses.MeanSquaredError(),
                optimizer=tf.keras.optimizers.Adam(),
                metrics=[tf.keras.metrics.MeanAbsoluteError()])

    history = model.fit(train_inputs, train_labels, epochs=max_epochs,
                      validation_data=(valid_inputs, valid_labels),
                      callbacks=[early_stopping])
    
    plotCol(test_inputs, test_inputs, test_labels, model=model, plot_col='Susceptible', max_subplots=3)

    return history

def cnn_network_multi_output(path,  max_epochs=30, early_stop=4, num_runs_traindata=500, save_evaluation_pdf=False):

    input_width = 5
    label_width = 20
    num_outputs = 8

    if not os.path.isfile(os.path.join(path, 'data_secir_simple.txt')):
        generate_data(num_runs_traindata, path)

    data = pd.read_csv(
        os.path.join(os.path.dirname(os.path.realpath(path)), 'data35', 'data_secir_simple.txt'),
        sep=' ')

    # normalization is kind of difficult for data nested in lists. So start with normalize all data
    # and then repeat progress again w/o normalization to get the correct labels
    data_mean = data.mean()
    data_std = data.std()
    data = (data - data_mean) / data_std

    train_data, valid_data, test_data = splitdata(data, normalize=False)

    train_inputs, train_labels = getPartitions(input_width = input_width, label_width = label_width, data=train_data)
    valid_inputs, valid_labels = getPartitions(input_width = input_width, label_width = label_width, data=valid_data)
    test_inputs, test_labels = getPartitions(input_width = input_width, label_width = label_width, data=test_data)

    # convert list to tensors
    train_inputs = tf.stack(train_inputs)
    train_labels = tf.stack(train_labels)
    valid_inputs = tf.stack(valid_inputs)
    valid_labels = tf.stack(valid_labels)
    test_inputs = tf.stack(test_inputs)
    test_labels = tf.stack(test_labels)

    CONV_WIDTH = 3
    model = tf.keras.Sequential([
        # Shape [batch, time, features] => [batch, CONV_WIDTH, features]
        tf.keras.layers.Lambda(lambda x: x[:, -CONV_WIDTH:, :]),
        # Shape => [batch, 1, conv_units]
        tf.keras.layers.Conv1D(256, activation='relu', kernel_size=(CONV_WIDTH)),
        # Shape => [batch, 1,  out_steps*features]
        tf.keras.layers.Dense(label_width*num_outputs,
                            kernel_initializer=tf.initializers.zeros()),
        # Shape => [batch, out_steps, features]
        tf.keras.layers.Reshape([label_width, num_outputs])
    ])

    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                    patience=early_stop,
                                                    mode='min')

    model.compile(loss=tf.keras.losses.MeanSquaredError(),
                optimizer=tf.keras.optimizers.Adam(),
                metrics=[tf.keras.metrics.MeanAbsoluteError()])

    history = model.fit(train_inputs, train_labels, epochs=max_epochs,
                      validation_data=(valid_inputs, valid_labels),
                      callbacks=[early_stopping])
    
    plotCol(test_inputs, test_inputs, test_labels, model=model, plot_col='Susceptible', max_subplots=3)

    return history

def lstm_network_multi_output(path,  max_epochs=30, early_stop=4, num_runs_traindata=500, save_evaluation_pdf=False):

    input_width = 5
    label_width = 20
    num_outputs = 8

    if not os.path.isfile(os.path.join(path, 'data_secir_simple.txt')):
        generate_data(num_runs_traindata, path)

    data = pd.read_csv(
        os.path.join(os.path.dirname(os.path.realpath(path)), 'data35', 'data_secir_simple.txt'),
        sep=' ')

    # normalization is kind of difficult for data nested in lists. So start with normalize all data
    # and then repeat progress again w/o normalization to get the correct labels
    data_mean = data.mean()
    data_std = data.std()
    data = (data - data_mean) / data_std

    train_data, valid_data, test_data = splitdata(data, normalize=False)

    train_inputs, train_labels = getPartitions(input_width = input_width, label_width = label_width, data=train_data)
    valid_inputs, valid_labels = getPartitions(input_width = input_width, label_width = label_width, data=valid_data)
    test_inputs, test_labels = getPartitions(input_width = input_width, label_width = label_width, data=test_data)

    # convert list to tensors
    train_inputs = tf.stack(train_inputs)
    train_labels = tf.stack(train_labels)
    valid_inputs = tf.stack(valid_inputs)
    valid_labels = tf.stack(valid_labels)
    test_inputs = tf.stack(test_inputs)
    test_labels = tf.stack(test_labels)

    model = tf.keras.Sequential([
    # Shape [batch, time, features] => [batch, lstm_units].
    # Adding more `lstm_units` just overfits more quickly.
    tf.keras.layers.LSTM(32, return_sequences=False),
    # Shape => [batch, out_steps*features].
    tf.keras.layers.Dense(label_width*num_outputs,
                          kernel_initializer=tf.initializers.zeros()),
    # Shape => [batch, out_steps, features].
    tf.keras.layers.Reshape([label_width, num_outputs])])


    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                    patience=early_stop,
                                                    mode='min')

    model.compile(loss=tf.keras.losses.MeanSquaredError(),
                optimizer=tf.keras.optimizers.Adam(),
                metrics=[tf.keras.metrics.MeanAbsoluteError()])

    history = model.fit(train_inputs, train_labels, epochs=max_epochs,
                      validation_data=(valid_inputs, valid_labels),
                      callbacks=[early_stopping])
    
    plotCol(test_inputs, test_inputs, test_labels, model=model, plot_col='Susceptible', max_subplots=3)

    return history
    
# simple benchmarking
def timereps(reps, model, input):
    from time import time
    start = time()
    for i in range(0, reps):
        _ = model(input)
    end = time()
    time_passed = end - start
    print(time_passed)
    return (end - start) / reps

# Plot Performance
def plot_histories(histories):
    model_names = ["LSTM", "Dense"]
    count_names = 0
    interval = np.arange(len(histories))
    for x in interval:
        history = histories[x]
        width = 0.3

        train_mae = history.history['mean_absolute_error'][-1]
        valid_mae = history.history['val_mean_absolute_error'][-1]

        plt.bar(x - 0.17, train_mae, width, label='Train')
        plt.bar(x + 0.17, valid_mae, width, label='Validation')
        name = model_names[x]
        # plt.xticks(ticks=x, labels=[name,name],
        #         rotation=45)
        plt.ylabel(f'MAE')
        _ = plt.legend()
        count_names = count_names + 1
    
    plt.show()
    plt.savefig('evaluation_single_shot.pdf')

x = 5

print('x')

if __name__ == "__main__":
    # TODO: Save contact matrix depending on the damping.
    # In the actual state it might be enough to save the regular one and the damping

    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(os.path.dirname(os.path.realpath(path)))), 'data')

    max_epochs = 1

    ### Data Generation ###
    num_runs = 100
    generate_data(num_runs, path_data)

    ### Plot Data ###
    # data = pd.read_csv(
    #     os.path.join(os.path.dirname(os.path.realpath(path)), 'data35', 'data_secir_simple.txt'),
    #     sep=' ')

    # train_data, val_data, test_data = splitdata(data)

    # inputs, labels = getPartitions(input_width = 1, label_width = 1, data=train_data)
    # plotCol(inputs, labels, plot_col='Susceptible', max_subplots=3)

    ### Models ###
    # single input
    # ml_network(path_data, max_epochs=max_epochs)

    # # multi input
    # lstm_hist = lstm_network_multi_input(path_data, max_epochs=max_epochs)
    # # ml_hist = ml_network_multi_input(path_data, max_epochs=max_epochs)

    # # Multi output
    # # cnn_output = cnn_network_multi_output(path_data, max_epochs=max_epochs)
    # # lstm_hist_multi = lstm_network_multi_output(path_data, max_epochs=max_epochs)
    
    # histories = [ lstm_hist, ml_hist]
    # plot_histories(histories)
