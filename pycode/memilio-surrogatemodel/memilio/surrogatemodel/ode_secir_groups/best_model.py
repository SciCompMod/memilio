from memilio.surrogatemodel.ode_secir_groups import network_architectures
from memilio.surrogatemodel.ode_secir_groups.model import split_data, get_test_statistic
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
    
filename = "data_secir_groups_30days_nodamp.pickle"
filename_df = "dataframe_0_512_LSTM_opt"

label_width = 30 
early_stop = 100

hidden_layers = 0
neurons_in_hidden_layer = 512
optimizer = ['Adam', 'Nadam','SGD', 'RMSProp', 'Adagrad']
#activations = ['relu', 'elu', 'sigmoid', 'leaky_relu', 'tanh', 'linear', 'softmax']
model= "LSTM"

label_width = 30

df_results  = pd.DataFrame(
    columns=['model', 'number_of_hidden_layers', 'number_of_neurons', 'optimizer', 
             'mean_test_MAPE', 'kfold_train',
             'kfold_val', 'kfold_test', 'training_time',
             'train_losses', 'val_losses'])       





parameters = []
for o in optimizer:
    
        parameters.append(o)


#parameters = []
#for o in optimizer:
#    for a in activations:
#        parameters.append((o,a))

#param = [hidden_layers, neurons_in_hidden_layer, model, optimizer, activations]
#for param in parameters: 

def train_and_evaluate_model(param, max_epochs):
        layer =hidden_layers
        neuron_number = neurons_in_hidden_layer
        modelname = "LSTM"
        #optimizer = param[0]
        optimizer = param
        #activation = param[1]
        

        
        if not os.path.isfile(os.path.join(path_data, 'data_secir_simple.pickle')):
            ValueError("no dataset found in path: " + path_data)

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

            if modelname == "LSTM":
                    
                num_outputs = 8*6
                model = tf.keras.Sequential([
                tf.keras.layers.LSTM(neuron_number, return_sequences=False)])
                
                for i in range(layer):
                        #model.add(tf.keras.layers.Dense(units=neuron_number, activation=activation))
                        model.add(tf.keras.layers.Dense(units=neuron_number))

                
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
                #optimizer=tf.keras.optimizers.Adam(),
                optimizer = optimizer ,
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

        df_results.loc[len(df_results.index)] = [modelname, layer, neuron_number, optimizer,  df.mean()[0] , np.mean(train_losses),
                                np.mean(val_losses),
                                np.mean(test_scores),
                                (elapsed / 60),
                                [losses_history_all],
                                [val_losses_history_all]]
        
        path = os.path.dirname(os.path.realpath(__file__))
        file_path = os.path.join(
            os.path.dirname(
                os.path.realpath(os.path.dirname(os.path.realpath(path)))),
            'secir_groups_LSTM_hyperparamter')
        if not os.path.isdir(file_path):
            os.mkdir(file_path)
        file_path = os.path.join(file_path,filename_df)
        df_results.to_csv(file_path)




start_hyper = time.perf_counter()
max_epochs = 1500

#train_and_evaluate_model(param, max_epochs)

for param in parameters:
     train_and_evaluate_model(param, max_epochs)

elapsed_hyper = time.perf_counter() - start_hyper
print(
    "Time for hyperparameter testing: {:.4f} minutes".format(
        elapsed_hyper / 60))
print(
    "Time for hyperparameter testing: {:.4f} hours".format(
        elapsed_hyper / 60 / 60))



    
   
 
      
        
     


