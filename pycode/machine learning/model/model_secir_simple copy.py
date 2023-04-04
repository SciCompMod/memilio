import numpy as np
import pandas as pd
from datetime import date
from math import ceil
import random
import os
import pickle
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
from data_secir_simple import generate_data, splitdata
from different_models import *




def plotCol(inputs, labels, model=None, plot_col='Susceptible', max_subplots=5):

    input_width = inputs.shape[1]
    label_width = labels.shape[1]
    
    plt.figure(figsize=(12, 8))
    cols = np.array(['Susceptible', 'Exposed', 'Carrier',
                        'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead'])  
    plot_col_index = np.where(cols == plot_col)[0][0]
    max_n = min(max_subplots, inputs.shape[0])

    # predictions = model(inputs) # for just one input: input_series = tf.expand_dims(inputs[n], axis=0) -> model(input_series)
    for n in range(max_n):
        plt.subplot(max_n, 1, n+1)
        plt.ylabel(f'{plot_col}')

        input_array = inputs[n].numpy()
        label_array = labels[n].numpy()
        plt.plot(np.arange(0,input_width), input_array[:, plot_col_index],
                label='Inputs', marker='.', zorder=-10)
        plt.scatter(np.arange(input_width,input_width+label_width), label_array[:,plot_col_index],
                    edgecolors='k', label='Labels', c='#2ca02c', s=64)
        
        if model is not None:
            input_series = tf.expand_dims(inputs[n], axis=0)
            pred = model(input_series)
            pred = pred.numpy()
            plt.scatter(np.arange(input_width,input_width+pred.shape[-2]), 
                                pred[0,:,plot_col_index],
                                marker='X', edgecolors='k', label='Predictions',
                                c='#ff7f0e', s=64)

    plt.xlabel('days')
    plt.show()
    plt.savefig('evaluation_secir_simple_' + plot_col + '.pdf')

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



def network_fit(path, model,  MAX_EPOCHS=30, early_stop=4, save_evaluation_pdf=False):

    if not os.path.isfile(os.path.join(path, 'data_secir_simple.pickle')):
        ValueError("no dataset found in path: " + path)

   
    file = open(os.path.join(path, 'data_secir_simple.pickle'),'rb')


    data = pickle.load(file)

    train_data, valid_data, test_data = splitdata(data, normalize=False)

    train = list(map(list, zip(*train_data)))
    valid = list(map(list, zip(*valid_data)))
    test = list(map(list, zip(*test_data)))


    train_inputs = tf.constant(train[0])    
    train_labels = tf.constant(train[1])
    valid_inputs = tf.constant(valid[0])
    valid_labels = tf.constant(valid[1])
    test_inputs  = tf.constant(test[0])
    test_labels  = tf.constant(test[1])

    


    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                    patience=early_stop,
                                                    mode='min')

    model.compile(loss=tf.keras.losses.MeanSquaredError(),
                optimizer=tf.keras.optimizers.Adam(),
                metrics=[tf.keras.metrics.MeanAbsoluteError()])

    history = model.fit(train_inputs, train_labels, epochs=MAX_EPOCHS,
                      validation_data=(valid_inputs, valid_labels),
                      callbacks=[early_stopping])

    plotCol(test_inputs, test_labels, model=model, plot_col='Infected', max_subplots=3)
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




print('x')

if __name__ == "__main__":
    # TODO: Save contact matrix depending on the damping.
    # In the actual state it might be enough to save the regular one and the damping

    

    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(path)), 'data_simple')

    MAX_EPOCHS = 50

    ### Models ###
    # single input
    # hist = network_fit(path_data, model=single_output(), MAX_EPOCHS=MAX_EPOCHS)

    # # multi input
    # lstm_hist = network_fit(path_data, model=lstm_network_multi_input(), MAX_EPOCHS=MAX_EPOCHS)
    # ml_hist = network_fit(path_data, model=multilayer_multi_input(), MAX_EPOCHS=MAX_EPOCHS)

    # # Multi output 
    cnn_output = network_fit(path_data, model=cnn_multi_output(), MAX_EPOCHS=MAX_EPOCHS)
    # lstm_hist_multi = network_fit(path_data, model=lstm_multi_output(), MAX_EPOCHS=MAX_EPOCHS)
    
    # histories = [ lstm_hist, ml_hist]
    # plot_histories(histories)
