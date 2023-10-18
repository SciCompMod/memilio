from memilio.surrogatemodel.ode_secir_simple import network_architectures
from memilio.surrogatemodel.ode_secir_simple.model import split_data, get_test_statistic
import os 
import tensorflow as tf
import pickle
import pandas as pd
import time
from sklearn.model_selection import KFold
import numpy as np



path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')
    
filename = "data_secir_simple.pickle"
filename_df = "dataframe"
max_epochs = 2
label_width = 30 
early_stop = 100

hidden_layers = [0,1,2,3,4]
neurons_in_hidden_layer = [32, 62, 128, 512, 1024]
models = ["Dense", "CNN", "LSTM"]

label_width = 30


parameters = []
for layer in hidden_layers:
    for neuron_number in neurons_in_hidden_layer:
        for modelname in models:
            parameters.append((layer, neuron_number, modelname))



df_results  = pd.DataFrame(
    columns=['model', 'number_of_hidden_layers', 'number_of_neurons',
             'mean_test_MAPE', 'kfold_train',
             'kfold_val', 'kfold_test', 'training_time',
             'train_losses', 'val_losses'])            


#for param in parameters: 

def train_and_evaluate_model(param):
    layer =param[0]
    neuron_number = param[1]
    modelname = param[2]
    

    
    if not os.path.isfile(os.path.join(path, 'data_secir_simple.pickle')):
        ValueError("no dataset found in path: " + path)

    file = open(os.path.join(path_data,filename), 'rb')

    data = pickle.load(file)


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

        if modelname == 'Dense':
                
            model = tf.keras.Sequential([
                tf.keras.layers.Flatten(),
                tf.keras.layers.Dense(units=neuron_number, activation='relu')])
                
            for i in range(layer):
                        model.add(tf.keras.layers.Dense(units=neuron_number, activation='relu'))
                    
            model.add(tf.keras.layers.Dense(units=label_width*8))
            model.add(tf.keras.layers.Reshape([label_width,8]))

        elif modelname == 'CNN':
            conv_size=3
            num_outputs = 8
            model = tf.keras.Sequential([
            tf.keras.layers.Lambda(lambda x: x[:, -conv_size:, :]),
            tf.keras.layers.Conv1D(neuron_number, activation='relu',
                                    kernel_size=(conv_size))])
            for i in range(layer):
                model.add(tf.keras.layers.Dense(units=neuron_number, activation='relu'))
                
            model.add(tf.keras.layers.Dense(label_width*num_outputs,
                                    kernel_initializer=tf.initializers.zeros()))
            model.add(tf.keras.layers.Reshape([label_width, num_outputs]))

        elif modelname == "LSTM":
                
            num_outputs = 8
            model = tf.keras.Sequential([
            tf.keras.layers.LSTM(neuron_number, return_sequences=False)])
            
            for i in range(layer):
                    model.add(tf.keras.layers.Dense(units=neuron_number, activation='relu'))
            
            model.add(tf.keras.layers.Dense(label_width*num_outputs,
                                    kernel_initializer=tf.initializers.zeros()))
            model.add(tf.keras.layers.Reshape([label_width, num_outputs]))
        #import torch
        #x = torch.take(data['inputs'], torch.tensor(train_idx[:(int(0.8*len(train_idx)))]))
         
        train_inputs = tf.gather(data['inputs'], indices = train_idx[:(int(0.8*len(train_idx)))])
        train_labels = tf.gather(data['labels'], indices = train_idx[:(int(0.8*len(train_idx)))])
        valid_inputs = tf.gather(data['inputs'], indices = train_idx[(int(0.8*len(train_idx))):])
        valid_labels = tf.gather(data['labels'], indices = train_idx[(int(0.8*len(train_idx))):])
        test_inputs = tf.gather(data['inputs'], indices = test_idx)
        test_labels = tf.gather(data['labels'], test_idx)
         
        #start = time.perf_counter()

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

    df_results.loc[len(df_results.index)] = [modelname, layer, neuron_number,  df.mean()[0] , np.mean(train_losses),
                             np.mean(val_losses),
                             np.mean(test_scores),
                             (elapsed / 60),
                             [losses_history_all],
                             [val_losses_history_all]]
    
    path = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(
        os.path.dirname(
            os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'secir_simple_grid_search')
    if not os.path.isdir(file_path):
        os.mkdir(file_path)
    file_path = file_path+filename_df
    df_results.to_csv(file_path)











    
   
 
      
        
     


