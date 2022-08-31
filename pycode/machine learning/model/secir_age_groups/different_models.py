import tensorflow as tf
from tensorflow import keras
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


def define_mlp_model(n_input):

    # using Sequential, cause every layer got exactly on i/o-tensor
    model = Sequential()
    # input layer
    model.add(Flatten())
    model.add(Dense(49, input_dim=n_input * 5,
              kernel_initializer='he_uniform', activation='relu'))
    # hidden layer
    model.add(Dense(512, activation='relu'))
    return model


def multilayer_multi_input():
    model = tf.keras.Sequential([
        tf.keras.layers.Flatten(),
        tf.keras.layers.Dense(
            units=128,  kernel_initializer='he_uniform', activation='relu'),
        tf.keras.layers.Dense(units=512, activation='relu')
        # tf.keras.layers.Reshape([1, -1]),
    ])
    return model


def define_cnn_model():
    model = Sequential()
    model.add(Conv1D(72, kernel_size=(1), activation='relu', input_shape=(6, 6)))
    model.add(Flatten())
    model.add(Dense(512, activation='relu'))
    return model


def define_combined_model(n_output, n_dampings=1):
    # mlp => dampingdays, rki
    mlp = multilayer_multi_input()  # define_mlp_model(n_input=48+n_dampings)
    # cnn => damping matrices
    cnn = [define_cnn_model() for _ in range(n_dampings)]
    combinedInput = concatenate([mlp.output]+[item.output for item in cnn])
    x = Dense(1024, activation='relu')(combinedInput)
    x = Dense(1024, activation='relu')(x)
    x = Dense(n_output, activation='linear')(x)
    model = Model(inputs=[mlp.input]+[item.input for item in cnn], outputs=x)
    return model
