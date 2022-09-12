import tensorflow as tf
from tensorflow import keras
from keras.utils import to_categorical
from keras.models import Model
from keras.models import Sequential
from keras.models import load_model
from keras.layers import Dense
from keras.layers import Dropout
from keras.layers import Conv1D
from keras.layers import Flatten
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

    model.add(Dense(1080,  # 1440
              kernel_initializer='he_uniform', activation='linear'))
    return model
