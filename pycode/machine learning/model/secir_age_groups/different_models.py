import tensorflow as tf
from tensorflow import keras
from keras.utils import to_categorical
from keras.models import Model
from keras.models import Sequential
from keras.models import load_model
from keras.layers import Dense, BatchNormalization
from keras.layers import Dropout
from keras.layers import Conv1D
from keras.layers import Dense, SimpleRNN
from keras.layers import Flatten
from keras.layers import MaxPooling1D
from keras.layers import GaussianNoise
from keras.layers import concatenate
from keras.callbacks import ModelCheckpoint
from keras.optimizers import Adam


def mlp_model():

    # using Sequential, cause every layer got exactly on i/o-tensor
    model = Sequential()
    # input layer
    model.add(Dense(312, input_dim=312,
              kernel_initializer='he_uniform', activation='relu'))
    # hidden layer
    model.add(Dense(512,
                    kernel_initializer='he_uniform', activation='relu'))

    model.add(Dense(512,
                    kernel_initializer='he_uniform', activation='relu'))

    model.add(Dense(960,  # 1440
              kernel_initializer='he_uniform', activation='linear'))
    return model


def cnn_model():

    model = Sequential()
    model.add(
        Conv1D(
            filters=64, kernel_size=3, activation='relu',
            input_shape=(312, 1)))  # 312
    model.add(Conv1D(filters=64, kernel_size=3, activation='relu'))
    # model.add(Dropout(0.5))
    model.add(MaxPooling1D(pool_size=2))
    model.add(Flatten())
    model.add(GaussianNoise(0.35))

    model.add(Dense(512, activation='relu'))

    model.add(Dense(512, activation='relu'))

    model.add(Dense(960, activation='linear'))  # 1440
    return model


def lstm_multi_output():
    output_dim = 960
    model = tf.keras.Sequential(
        [tf.keras.layers.LSTM(
            64, kernel_initializer=tf.keras.initializers.GlorotUniform(
                seed=42),
            input_shape=(5, 120),
            return_sequences=False, recurrent_activation='relu',
            go_backwards=True),
         tf.keras.layers.Dense(
             output_dim, kernel_initializer=tf.initializers.zeros())])

    return model
