
import tensorflow as tf
from keras.models import Sequential
from keras.layers import Dense, Conv1D, MaxPooling1D, Flatten, GaussianNoise


def mlp_model(label_width):

    model = Sequential()
    # input layer
    model.add(Dense(276, input_dim=276,
              kernel_initializer='he_uniform', activation='relu'))
    # hidden layer
    model.add(Dense(512,
                    kernel_initializer='he_uniform', activation='relu'))

    model.add(Dense(512,
                    kernel_initializer='he_uniform', activation='relu'))

    model.add(Dense(48*label_width,
              kernel_initializer='he_uniform', activation='linear'))
    return model


def cnn_model():

    model = Sequential()
    model.add(
        Conv1D(
            filters=64, kernel_size=3, activation='relu',
            input_shape=(276, 1)))
    model.add(Conv1D(filters=64, kernel_size=3, activation='relu'))
    # model.add(Dropout(0.5))
    model.add(MaxPooling1D(pool_size=2))
    model.add(Flatten())
    model.add(GaussianNoise(0.35))

    model.add(Dense(512, activation='relu'))

    model.add(Dense(512, activation='relu'))

    model.add(Dense(1440, activation='linear'))
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
