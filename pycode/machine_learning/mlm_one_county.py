# Machine Learning model ==> classification problem

import os
import datetime
import itertools
import numpy as np
import pandas as pd
import tensorflow as tf

from matplotlib import pyplot
from collections import Counter

from keras.utils import to_categorical
from keras.models import Model
from keras.models import Sequential
from keras.models import load_model
from keras.layers import Dense
from keras.layers import Conv1D
from keras.layers import Flatten
from keras.layers import concatenate
from keras.callbacks import ModelCheckpoint
from keras.optimizers import Adam

from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split



class ML_Model:

    __prepared_data = {}

    def __init__(self, file_name):
        # load and prepocess data
        t1 = datetime.datetime.now()
        self.__prepared_data = self.load_data(file_name)
        __loading_status = True
        t2 = datetime.datetime.now()
        print("Loading duration: " + str(t2-t1))


    def deactivateGPU(self):
        # deactivate gpu because its slower
        os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

    def activateGPU(self):
        os.environ["CUDA_VISIBLE_DEVICES"] = "0"


    # data loading and preprocessing
    # returns:  trainX, trainY, testX, testY    type==numpy.array
    def load_data(self, file_name):
        df = pd.read_csv(file_name)
        # CSV format: damping_matrix, dampingday, rki_acute_infected[age_groups], rki_recoverd[age_groups], rki_deaths[age_groups], divi_hospital, divi_icu, max_infected[age_groups], max_infected_time[age_groups], last_day_infected[age_groups], infected_sum[age_groups], last_day_deaths[age_groups]

        data = df.values
        # number of dampings
        n_dampings = 3
        # calculating the count of X values
        divider = 48 + (n_dampings * 2)
        # seperate i/o
        X, Y =  data[:,:divider], data[:,divider:]     # :2 => Indizies bis (exklusive 2), also 0&1   | 2 => Index 2


        # split data into test and train set
        trainX, testX, trainY, testY = train_test_split(X, Y, test_size=0.2, random_state=42)

        trainX_cnn, testX_cnn = [trainX[:,i] for i in range(0,2*n_dampings,2)], [testX[:,i] for i in range(0,2*n_dampings,2)]
        trainX_mlp, testX_mlp = [list([item[i] for i in range(1,2*n_dampings,2)]) + list(item[n_dampings*2:]) for item in trainX], [list([item[i] for i in range(1,2*n_dampings,2)]) + list(item[n_dampings*2:]) for item in testX] ###trainX[:,n_dampings*2-1:] + [trainX[:,i] for i in range(1,2*n_dampings,2)], testX[:,n_dampings*2-1:] + [testX[:,i] for i in range(1,2*n_dampings,2)]

        tmp_trainX_cnn = []
        for item in trainX_cnn:
            damping_matrix_train = []
            for damping_matrix_string in item:
                # conversion from string format to 2d list
                tmp = [[float(num) for num in row.replace('[','').replace(']','').replace(',','').split(' ') if len(num) > 0] for row in damping_matrix_string.split('], ')]
                damping_matrix_train.append(np.array(tmp))
            tmp_trainX_cnn.append(np.array(damping_matrix_train))
        trainX_cnn = tmp_trainX_cnn


        tmp_testX_cnn = []
        for item in testX_cnn:
            damping_matrix_test = []
            for damping_matrix_string in item:
                # conversion from string format to 2d list
                tmp = [[float(num) for num in row.replace('[','').replace(']','').replace(',','').split(' ') if len(num) > 0] for row in damping_matrix_string.split('], ')]
                damping_matrix_test.append(np.array(tmp))
            tmp_testX_cnn.append(np.array(damping_matrix_test))
        testX_cnn = tmp_testX_cnn


        # normalize data
        sc = StandardScaler()
        sc.fit(trainX_mlp)
        trainX_mlp = sc.transform(trainX_mlp)
        testX_mlp = sc.transform(testX_mlp)
        sc.fit(trainY)
        scaleTrainY = sc.transform(trainY)
        scaleTestY = sc.transform(testY)


        return {"trainX_mlp":trainX_mlp, "trainX_cnn":trainX_cnn, "trainY":trainY, "testX_mlp":testX_mlp, "testX_cnn":testX_cnn, "testY":testY, "scaleTrainY":scaleTrainY, "scaleTestY":scaleTestY}


    # creation of neural network
    # returns:  NN-model    type==Sequential
    def define_mlp_model(self, n_input):
        # using Sequential, cause every layer got exactly on i/o-tensor
        model = Sequential()
        #input layer
        model.add(Dense(51, input_dim=n_input, kernel_initializer='he_uniform', activation='relu'))
        #hidden layer
        model.add(Dense(512, activation='relu'))
        return model


    def define_cnn_model(self):
        model = Sequential()

        model.add(Conv1D(72, kernel_size=(1), activation='relu', input_shape=(6,6)))

        model.add(Flatten())

        model.add(Dense(512, activation='relu'))

        return model



    def define_combined_model(self, n_output, n_dampings, eta=0.001, loss_func='mse', metrics_funcs=['mse', 'mae']):
        # mlp => dampingdays, rki
        mlp = self.define_mlp_model(n_input=48+n_dampings)
        # cnn => damping matrices
        cnn = [self.define_cnn_model() for i in range(n_dampings)]

        combinedInput = concatenate([mlp.output]+[item.output for item in cnn])

        x = Dense(1024, activation='relu')(combinedInput)
        x = Dense(1024, activation='relu')(x)
        x = Dense(n_output, activation='linear')(x)

        model = Model(inputs=[mlp.input]+[item.input for item in cnn], outputs=x)

        opt = Adam(lr=eta)
        model.compile(loss=loss_func, optimizer=opt, metrics=metrics_funcs)
        return model

    # training and evaluation process
    # returns: scores, histories   type==list 
    def evaluate_model(self, data, n_epochs=40, n_batch_size=512, name="best_model",verbose=0, metrics_monitor='mse', eta=0.001, loss_func='mae', metrics_funcs=['mae', 'mse']):

        # save time to calculate run-time
        t1=datetime.datetime.now()

        # define model
        model = self.define_combined_model(n_output=len(data["trainY"][0]), n_dampings=len(data["trainX_cnn"]), eta=eta, loss_func=loss_func, metrics_funcs=metrics_funcs)

        mc = ModelCheckpoint(name + '.h5', monitor=metrics_monitor, mode='min', verbose=verbose, save_best_only=True)

        
        # train model
        history = model.fit(x=[data["trainX_mlp"]]+ data["trainX_cnn"], y=data["scaleTrainY"], epochs=n_epochs, batch_size=n_batch_size, verbose=verbose, callbacks=[mc])
        # evaluate model
        #_, acc = model.evaluate(testdataX, testdataY, verbose=0)
        # load the saved model
        saved_model = load_model(name + '.h5')
        error_train = saved_model.evaluate([data["trainX_mlp"]]+ data["trainX_cnn"], data["scaleTrainY"], verbose=0)
        error_test= saved_model.evaluate([data["testX_mlp"]]+ data["testX_cnn"], data["scaleTestY"], verbose=0)
        # save time to calculate run-time
        t2=datetime.datetime.now()
        print("LM: " + str(len(data["testX_mlp"][0])))
        print("LC: " + str(len(data["testX_cnn"][0])))
        x_input = [data["testX_mlp"]]+ data["testX_cnn"]
        print(len(x_input))
        predtest = saved_model.predict(x_input)
        #print(len(x_input))
        #print(np.shape(x_input))
        #print(0/0)
        #return
        sc = StandardScaler()
        sc.fit(data["trainY"])
        predtest = sc.inverse_transform(predtest)
        print('Calculating time: ', str(t2-t1))
        print('Train [MAE, MSE]: ', str(error_train[1:]))
        print('Test [MAE, MSE]: ', str(error_test[1:]))
        print("* shown infos are just correct if default isn't change")

        return [error_train[1:], error_test[1:]], history, predtest, t2-t1


    # plot diagnostic learning curves
    def summarize_diagnostics(self, history, scores, predtest, data):
        # plot accuracy to epochs
        pyplot.subplot(3,2,1)
        pyplot.title('Mean Squared Error')
        pyplot.plot(history.history['mse'], color='blue', label='train')
        pyplot.subplot(3,2,2)
        pyplot.title('Mean Absolute Error')
        pyplot.plot(history.history['mae'], color='blue', label='train')

        max_infected = [sum([day[i] for i in range(0,6,1)]) for day in data["testY"]]
        last_day_infected = [sum([day[i] for i in range(12,18,1)]) for day in data["testY"]]
        sum_infected = [sum([day[i] for i in range(18,24,1)]) for day in data["testY"]]

        pred_max_infected = [sum([day[i] for i in range(0,6,1)]) for day in predtest]
        pred_last_day_infected = [sum([day[i] for i in range(12,18,1)]) for day in predtest]
        pred_sum_infected = [sum([day[i] for i in range(18,24,1)]) for day in predtest]

        x_ax = range(len(data["testX_mlp"]))
        pyplot.subplot(3,2,3)
        pyplot.title('Maximum Infected (Groups 2-4)')
        pyplot.plot(x_ax, max_infected,  label="real",color="c")
        pyplot.plot(x_ax, pred_max_infected, label="pred")
        pyplot.legend()
        pyplot.subplot(3,2,4)
        pyplot.title('Day when Maximum Infected (Groups 2-4)')
        pyplot.plot(x_ax, data["testY"][:,6],  label="real",color="m")
        pyplot.plot(x_ax, predtest[:,6], label="pred")
        pyplot.legend()
        pyplot.subplot(3,2,5)
        pyplot.title('Infected on last day of simulation')
        pyplot.plot(x_ax, last_day_infected,  label="real",color="r")
        pyplot.plot(x_ax, pred_last_day_infected, label="pred")
        pyplot.legend()
        pyplot.subplot(3,2,6)
        pyplot.title('Infected sum over simulation')
        pyplot.plot(x_ax, sum_infected,  label="real",color="g")
        pyplot.plot(x_ax, pred_sum_infected, label="pred")
        pyplot.legend()

        pyplot.show()


    def run_Model(self, n_epochs, n_batch_size, name, eta, loss_func, metrics_funcs):
        data = self.__prepared_data
        # train and evaluate model
        scores, history, predtest, calculation_time = self.evaluate_model(data, n_epochs=n_epochs, n_batch_size=n_batch_size, name=name, eta=eta, loss_func=loss_func, metrics_funcs=metrics_funcs)
        # show results
        #self.summarize_diagnostics(history, scores, predtest, data)
        histories = []
        for item in metrics_funcs:
            histories.append(history.history[item])
        return scores, histories, predtest, calculation_time

    def getTestY(self):
        return self.__prepared_data["testY"]

    def getYScaler(self):
        return StandardScaler().fit(self.__prepared_data["trainY"])

    def getXMLPScaler(self):
        return StandardScaler().fit(self.__prepared_data["trainX_mlp"])
