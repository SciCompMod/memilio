import tensorflow as tf
from tensorflow import keras
from keras.utils import to_categorical
from keras.models import Model
from keras.models import Sequential
from keras.models import load_model
from keras.layers import Dense, BatchNormalization
from keras.layers import Dropout
from keras.layers import Conv1D
from keras.layers import LSTM
from keras.layers import Dense, SimpleRNN
from keras.layers import Flatten
from keras.layers import MaxPooling1D
from keras.layers import GaussianNoise
from keras.layers import concatenate
from keras.callbacks import ModelCheckpoint
from keras.optimizers import Adam


def mlp_model(input_dim, output_dim):

    # using Sequential, cause every layer got exactly on i/o-tensor
    model = Sequential()
    # input layer
    model.add(Dense(input_dim, input_dim=input_dim,
              kernel_initializer='he_uniform', activation='relu'))
    # hidden layer
    model.add(Dense(512,
                    kernel_initializer='he_uniform', activation='relu'))

    model.add(Dense(512,
                    kernel_initializer='he_uniform', activation='relu'))

    model.add(Dense(output_dim,  # 1440
              kernel_initializer='he_uniform', activation='linear'))
    return model


def cnn_model(input_dim, output_dim):

    model = Sequential()
    model.add(
        Conv1D(
            filters=64, kernel_size=3, activation='relu',
            input_shape=(input_dim, 1)))  # 312
    model.add(Conv1D(filters=64, kernel_size=3, activation='relu'))
    # model.add(Dropout(0.5))
    model.add(MaxPooling1D(pool_size=2))
    model.add(Flatten())
    model.add(GaussianNoise(0.35))

    model.add(Dense(512, activation='relu'))

    model.add(Dense(512, activation='relu'))

    model.add(Dense(output_dim, activation='linear'))  # 1440
    return model


def lstm_multi_output(input_dim, output_dim):
    
    model = tf.keras.Sequential(
        [tf.keras.layers.LSTM(
            64, kernel_initializer=tf.keras.initializers.GlorotUniform(
                seed=42),
            input_shape=(input_dim),
            return_sequences=False, recurrent_activation='relu',
            go_backwards=True),
         tf.keras.layers.Dense(
             output_dim, kernel_initializer=tf.initializers.zeros())])

    return model




# def cnn_lstm_hybrid(input_dim, output_dim):
    
#     model = Sequential()
#     model.add(
#         (Conv1D(
#             filters=64, kernel_size=3, activation='relu',
#             input_shape=(input_dim, 1))))  # 312
#     #model.add(Conv1D(filters=64, kernel_size=3, activation='relu'))
#     # model.add(Dropout(0.5))
#     model.add(((MaxPooling1D(pool_size=2))))
#     # model.add((Flatten()))
#     model.add(GaussianNoise(0.35))

#     model.add(LSTM(512, activation='relu'))
#     model.add(BatchNormalization()),
#     model.add(Dense(512, activation='relu')),
#     model.add(BatchNormalization()),
#     model.add(Dense(512, activation='relu'))
#     model.add(Dense(512, activation='relu'))

#     model.add(Dense(output_dim, activation='linear'))  # 1440
#     return model

# for lstm data


def cnn_lstm_hybrid_2(input_dim, output_dim):

    model = Sequential()
    model.add(
        (Conv1D(
            filters=64, kernel_size=3, activation='relu',
            input_shape=(input_dim))))  # 312
    #model.add(Conv1D(filters=64, kernel_size=3, activation='relu'))
    model.add(Dropout(0.5))
    model.add(((MaxPooling1D(pool_size=2))))
    # model.add((Flatten()))
    model.add(GaussianNoise(0.35))

    model.add(LSTM(512, activation='relu'))
    model.add(BatchNormalization()),
    model.add(Dense(512, activation='relu'))
    # model.add(BatchNormalization()),
    # model.add(Dense(512, activation='relu'))
    # model.add(BatchNormalization()),
    # model.add(Dense(512, activation='relu'))
    model.add(BatchNormalization()),
    model.add(Dense(output_dim, activation='linear'))  # 1440
    return model


def cnn_lstm_hybrid_3damp(input_dim, output_dim):

    model = Sequential()
    model.add((
        Conv1D(
            filters=64, kernel_size=3, activation='relu',
            input_shape=(5, input_dim),
            kernel_initializer=tf.keras.initializers.GlorotUniform(
                seed=42))))

    model.add(BatchNormalization()),
    model.add(((MaxPooling1D(pool_size=3))))
    # model.add((Flatten()))
    model.add(GaussianNoise(0.2))

    model.add(LSTM(512, activation='relu'))
    model.add(BatchNormalization()),
    model.add(Dropout(0.35))
    model.add(Dense(512, activation='relu'))
    model.add(BatchNormalization()),
    model.add(Dense(output_dim, activation='linear'))  # 1440
    return model


def cnn_lstm_hybrid_3damp_new(input_dim, output_dim):

    model = Sequential()
    model.add((
        Conv1D(
            filters=64, kernel_size=3, activation='relu',
            input_shape=(5, input_dim),
            kernel_initializer=tf.keras.initializers.GlorotUniform(
                seed=42))))

    model.add(BatchNormalization()),
    model.add(((MaxPooling1D(pool_size=3))))
    # model.add((Flatten()))
    model.add(GaussianNoise(0.2))

    model.add(LSTM(512, activation='relu'))
    model.add(BatchNormalization()),
    model.add(Dropout(0.35))
    model.add(Dense(512, activation='relu'))
    model.add(BatchNormalization()),
    model.add(Dropout(0.35))
    model.add(Dense(512, activation='relu'))
    model.add(BatchNormalization()),
    model.add(Dense(output_dim, activation='linear'))  # 1440
    return model


def cnn_lstm_hybrid_4damp(input_dim, output_dim):

    model = Sequential()
    model.add(
        (Conv1D(
            filters=64, kernel_size=3, activation='relu',
            input_shape=(5, input_dim))))

    model.add(BatchNormalization()),
    model.add(((MaxPooling1D(pool_size=2))))
    # model.add((Flatten()))
    model.add(GaussianNoise(0.35))

    model.add(LSTM(512, activation='relu'))
    model.add(BatchNormalization()),
    model.add(Dropout(0.35))
    model.add(Dense(512, activation='relu'))
    model.add(BatchNormalization()),
    model.add(Dense(output_dim, activation='linear'))  # 1440
    return model


def cnn_lstm_hybrid_5damp(input_dim, output_dim):

    model = Sequential()
    model.add((
        Conv1D(
            filters=64, kernel_size=3, activation='relu',
            input_shape=(5, input_dim),
            kernel_initializer=tf.keras.initializers.GlorotUniform(
                seed=42))))

    model.add(BatchNormalization()),
    model.add(((MaxPooling1D(pool_size=2))))
    # model.add((Flatten()))
    model.add(GaussianNoise(0.35))

    model.add(LSTM(512, activation='relu'))
    model.add(BatchNormalization()),
    model.add(Dropout(0.35))
    model.add(Dense(512, activation='relu'))
    # model.add(BatchNormalization()),
    #model.add(Dense(512, activation='relu'))
    model.add(BatchNormalization()),
    model.add(Dense(output_dim, activation='linear'))  # 1440
    return model

# def lstm_multi_output():
#     #output_dim = 960
#     output_dim = 7200
#     model = tf.keras.Sequential(
#         [tf.keras.layers.LSTM(
#             1024, kernel_initializer=tf.keras.initializers.GlorotUniform(
#                 seed=42),
#             input_shape=(5, 233),
#             return_sequences=False, recurrent_activation='relu',
#             go_backwards=True),
#          # tf.keras.layers.GaussianNoise(0.25),
#          tf.keras.layers.Dense(
#              output_dim, kernel_initializer=tf.initializers.zeros()),
#          ])

#     return model


# def cnn_model(input_dim, output_dim):

#     model = Sequential()
#     model.add(
#         Conv1D(
#             filters=64, kernel_size=3, activation='relu',
#             input_shape=(input_dim, 1)))  # 312
#     model.add(Conv1D(filters=64, kernel_size=3, activation='relu'))
#     # model.add(Dropout(0.5))
#     model.add(MaxPooling1D(pool_size=2))
#     model.add(Flatten())
#     model.add(GaussianNoise(0.35))

#     model.add(Dense(512, activation='relu'))

#     model.add(Dense(512, activation='relu'))

#     model.add(Dense(output_dim, activation='linear'))  # 1440
#     return model
