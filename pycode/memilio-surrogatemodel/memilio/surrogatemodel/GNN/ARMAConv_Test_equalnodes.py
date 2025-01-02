import os
import pickle
import spektral
import time
import json
import pandas as pd
import numpy as np
import tensorflow as tf

from keras.layers import Dense
from keras.losses import MeanAbsolutePercentageError
from keras.metrics import mean_absolute_percentage_error
from keras.models import Model
from keras.optimizers import Adam
from spektral.data import MixedLoader
from spektral.layers import ARMAConv
from spektral.transforms.normalize_adj import NormalizeAdj
from spektral.utils.convolution import normalized_laplacian, rescale_laplacian

# load and prepare data
# path = os.path.dirname(os.path.realpath(__file__))
# path_data = os.path.join(
#    os.path.dirname(
#        os.path.realpath(os.path.dirname(os.path.realpath(path)))),
#    'data_GNN_paper')

# file = open(os.path.join(path_data, 'GNN_data_30days_10k.pickle'), 'rb')
file = open(
    '/hpc_data/schm_a45/data_paper/GNN_data_60days_samenodes_1k.pickle', 'rb')
data_secir = pickle.load(file)

len_dataset = data_secir['inputs'].shape[0]
number_of_nodes = data_secir['inputs'].shape[3]

# transpose to num_runs, num_nodes, input_days, (num_compartments*num_age_groups) e.g. 1000,400,5,48
data_secir['inputs'] = data_secir['inputs'].transpose(0, 3, 2, 1)
data_secir['labels'] = data_secir['labels'].transpose(0, 3, 2, 1)

n_input_days = data_secir['inputs'].shape[2]
n_label_days = data_secir['labels'].shape[2]
n_age_groups = 6
n_compartments = 8

train_inputs = data_secir['inputs'][:int(0.8*len(data_secir['inputs']))]
train_labels = data_secir['labels'][:int(0.8*len(data_secir['labels']))]
test_inputs = data_secir['inputs'][int(0.8*len(data_secir['inputs'])):]
test_labels = data_secir['labels'][int(0.8*len(data_secir['labels'])):]

# # our dataset only contain 1000 samples. The variance in these datasets is relatively high
# # variance. We want to augment the trainig dataset by copying the data (first without any changes)
# # Repeat along the first axis 5 times
# train_inputs = np.tile(data_secir['inputs'][:int(
#     0.8*len(data_secir['inputs']))], (5, 1, 1, 1))
# train_labels = np.tile(data_secir['labels'][:int(
#     0.8*len(data_secir['labels']))], (5, 1, 1, 1))

new_inputs_train = train_inputs.reshape(
    len(train_inputs), number_of_nodes, n_input_days*n_age_groups*n_compartments)
new_labels_train = train_labels.reshape(
    len(train_labels), number_of_nodes, n_label_days*n_age_groups*n_compartments)

new_inputs_test = test_inputs.reshape(
    len(test_inputs), number_of_nodes, n_input_days*n_age_groups*n_compartments)
new_labels_test = test_labels.reshape(
    len(test_labels), number_of_nodes, n_label_days*n_age_groups*n_compartments)

# new_inputs = np.asarray(
#     data_secir['inputs']).reshape(
#     len_dataset, number_of_nodes, n_input_days*n_age_groups*n_compartments)
# new_labels = np.asarray(data_secir['labels']).reshape(
#     len_dataset, number_of_nodes, n_label_days*n_age_groups*n_compartments)


######## open commuter data #########
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(path)))))),
                         'data/mobility')
commuter_file = open(os.path.join(
    path_data, 'commuter_mobility.txt'), 'rb')
commuter_data = pd.read_csv(commuter_file, sep=" ", header=None)
sub_matrix = commuter_data.iloc[:number_of_nodes, 0:number_of_nodes]

adjacency_matrix = np.asarray(sub_matrix)
adjacency_matrix[adjacency_matrix > 0] = 1

node_features_train = new_inputs_train
node_labels_train = new_labels_train

node_features_test = new_inputs_test
node_labels_test = new_labels_test

layer = ARMAConv
number_of_layers = 3
number_of_channels = 512
parameters = [layer, number_of_layers, number_of_channels]

df = pd.DataFrame(
    columns=['layer', 'number_of_layers', 'channels', 'train_MAPE',
             'val_MAPE', 'test_MAPE', 'notlog_test_MAPE' ,'training_time',
             'train_losses', 'val_losses'])


def get_test_statistic(loader, model):
    output = []
    step = 0
    while step < loader.steps_per_epoch:
        step += 1
        inputs, target = loader.__next__()
        pred = model(inputs, training=False)

    pred = pred.numpy()
    #test_labels = np.array(target)
    print('Target',target)
    test_labels = np.asarray(target)
    diff = pred - test_labels
    relative_err = (abs(diff))/abs(test_labels)
    # reshape [batch, time, features] -> [features, time * batch]
    relative_err_transformed = relative_err.transpose(2, 0, 1).reshape(8, -1)
    relative_err_means_percentage = relative_err_transformed.mean(axis=1) * 100

    # same evaluation for rescaled data:
    pred = np.expm1(pred)
    test_labels = np.expm1(np.array(test_labels))

    diff_rescaled = pred - test_labels
    relative_err_rescaled = (abs(diff_rescaled))/abs(test_labels)
    # reshape [batch, time, features] -> [features, time * batch]
    relative_err_transformed_rescaled = relative_err_rescaled.transpose(
        2, 0, 1).reshape(8, -1)
    relative_err_means_percentage_rescaled = relative_err_transformed_rescaled.mean(
        axis=1) * 100

    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']
    mean_percentage = pd.DataFrame(
        data=relative_err_means_percentage,
        index=infectionstates,
        columns=['MAPE_Scaled'])
    mean_percentage['MAPE_rescaled'] = relative_err_means_percentage_rescaled

    print('MAPE scaled data: ', mean_percentage['MAPE_Scaled'].mean(), "%")
    print('MAPE rescaled data: ', mean_percentage['MAPE_rescaled'].mean(), "%")

    # save the results as csv

    path = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'testing_paper')
    # file_path = '/localdata1/gnn_paper_2024/data_Iteration2/results/testing'
    path_models = '/hpc_data/schm_a45/data_paper/'

    if not os.path.isdir(file_path):
        os.mkdir(file_path)
    file_path = os.path.join(file_path, filename_df)
    mean_percentage.to_csv(file_path)
    return mean_percentage


def train_and_evaluate_model(
        epochs, learning_rate, param, save_name, filename, filename_df):

    layer_name, number_of_layer, channels = param

    class MyDataset(spektral.data.dataset.Dataset):
        def read(self):

            self.a = rescale_laplacian(normalized_laplacian(adjacency_matrix))

            return [spektral.data.Graph(x=x, y=y) for x, y in zip(node_features_train, node_labels_train)]

            super().__init__(**kwargs)

    # data = MyDataset()
    data_train = MyDataset(transforms=NormalizeAdj())

    class MyDataset(spektral.data.dataset.Dataset):
        def read(self):

            self.a = rescale_laplacian(normalized_laplacian(adjacency_matrix))

            return [spektral.data.Graph(x=x, y=y) for x, y in zip(node_features_test, node_labels_test)]

            super().__init__(**kwargs)

    data_test = MyDataset(transforms=NormalizeAdj())

    batch_size = 32
    epochs = epochs
    es_patience = 200  # Patience for early stopping

    class Net(Model):
        def __init__(self):
            super().__init__()
            self.conv1 = layer(channels, order = 2, activation='elu')
            self.conv2 = layer(channels, order = 2, activation='elu')
            self.conv3 = layer(channels, order = 2, activation='elu')
            # self.conv4 = layer(channels, order=2, iterations=2,
            #                   share_weights=True,  activation='elu')
            # not important if data_train or data_test, as we only need label size
            self.dense = Dense(data_train.n_labels, activation="linear")

        def call(self, inputs):
            x, a = inputs
            a = np.asarray(a)
            x = self.conv1([x, a])
            x = self.conv2([x, a])
            x = self.conv3([x, a])
            # x = self.conv4([x, a])
            output = self.dense(x)

            return output

    optimizer = Adam(learning_rate=learning_rate)
    loss_fn = MeanAbsolutePercentageError()

    def train_step(inputs, target):
        with tf.GradientTape() as tape:
            predictions = model(inputs, training=True)
            loss = loss_fn(target, predictions) + sum(model.losses)
        gradients = tape.gradient(loss, model.trainable_variables)
        optimizer.apply_gradients(zip(gradients, model.trainable_variables))
        acc = tf.reduce_mean(
            mean_absolute_percentage_error(
                target, predictions))
        return loss, acc

    def evaluate(loader):
        output = []
        step = 0
        while step < loader.steps_per_epoch:
            step += 1
            inputs, target = loader.__next__()
            pred = model(inputs, training=False)
            outs = (
                loss_fn(target, pred),
                tf.reduce_mean(mean_absolute_percentage_error(target, pred)),
                len(target),  # Keep track of batch size
            )
            output.append(outs)
            if step == loader.steps_per_epoch:
                output = np.array(output)
                return np.average(output[:, :-1], 0, weights=output[:, -1])

    def evaluate_r(loader):
        output = []
        step = 0
        while step < loader.steps_per_epoch:
            step += 1
            inputs, target = loader.__next__()
            pred = model(inputs, training=False)
            outs = (
                loss_fn(np.expm1(target), np.expm1(pred)),
                tf.reduce_mean(mean_absolute_percentage_error(
                    np.expm1(target), np.expm1(pred))),
                len(target),  # Keep track of batch size
            )
            output.append(outs)
            if step == loader.steps_per_epoch:
                output = np.array(output)
                return np.average(output[:, :-1], 0, weights=output[:, -1])

    test_scores = []
    test_scores_r = []
    train_losses = []
    val_losses = []

    losses_history_all = []
    val_losses_history_all = []

    start = time.perf_counter()

    model = Net()

    # 80% for training --> 20% if that for validation
    idxs = np.arange(len(data_train))
    split_va = int(0.8 * len(data_train))
    idx_tr, idx_va = np.split(idxs, [split_va])

    data_tr = data_train[idx_tr]
    data_va = data_train[idx_va]
    data_te = data_test

    # Data loaders
    loader_tr = MixedLoader(
        data_tr, batch_size=batch_size, epochs=epochs)
    loader_va = MixedLoader(data_va,  batch_size=batch_size)
    loader_te = MixedLoader(data_te, batch_size=data_te.n_graphs)

    epoch = step = 0
    best_val_loss = np.inf
    best_weights = None
    patience = es_patience
    results = []
    losses_history = []
    val_losses_history = []

    start = time.perf_counter()

    for batch in loader_tr:
        step += 1
        loss, acc = train_step(*batch)
        results.append((loss, acc))
        if step == loader_tr.steps_per_epoch:
            step = 0
            epoch += 1

            # Compute validation loss and accuracy
            val_loss, val_acc = evaluate(loader_va)
            print(
                "Ep. {} - Loss: {:.3f} - Acc: {:.3f} - Val loss: {:.3f} - Val acc: {:.3f}".format(
                    epoch, *np.mean(results, 0), val_loss, val_acc
                )
            )

            # Check if loss improved for early stopping
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                patience = es_patience
                print("New best val_loss {:.3f}".format(val_loss))
                best_weights = model.get_weights()
            else:
                patience -= 1
                if patience == 0:
                    print(
                        "Early stopping (best val_loss: {})".format(
                            best_val_loss))
                    break
            results = []
            losses_history.append(loss)
            val_losses_history.append(val_loss)

    ################################################################################
    # Evaluate model
    ################################################################################
    model.set_weights(best_weights)  # Load best model
    test_loss, test_acc = evaluate(loader_te)
    test_loss_r, test_acc_r = evaluate_r(loader_te)
    #mape = get_test_statistic(loader_te, model)

    print(
        "Done. Test loss: {:.4f}. Test acc: {:.2f}".format(
            test_loss, test_acc))
    test_scores.append(test_loss)
    test_scores_r.append(test_loss_r)
    train_losses.append(np.asarray(losses_history).min())
    val_losses.append(np.asarray(val_losses_history).min())
    losses_history_all.append(np.asarray(losses_history))
    val_losses_history_all.append(val_losses_history)

    elapsed = time.perf_counter() - start

    # print out stats
    print("Best train losses: {} ".format(train_losses))
    print("Best validation losses: {}".format(val_losses))
    print("Test values: {}".format(test_scores))
    print("--------------------------------------------")
    print("Train Score:{}".format(np.mean(train_losses)))
    print("Validation Score:{}".format(np.mean(val_losses)))
    print("Test Score (log): {}".format(np.mean(test_scores)))
    print("Test Score (orig.): {}".format(np.mean(test_scores_r)))

    print("Time for training: {:.4f} seconds".format(elapsed))
    print("Time for training: {:.4f} minutes".format(elapsed/60))


    # Ensure that save_name has the .pickle extension
    if not save_name.endswith('.pickle'):
        save_name += '.pickle'

    # Save best weights as pickle
    path = os.path.dirname(os.path.realpath(__file__))
    file_path_w = os.path.join(path, 'saved_weights_paper')

    # Ensure the directory exists
    if not os.path.isdir(file_path_w):
        os.mkdir(file_path_w)

    # Construct the full file path by joining the directory with save_name
    file_path_w = os.path.join(file_path_w, save_name)

    # Save the weights to the file
    with open(file_path_w, 'wb') as f:
        pickle.dump(best_weights, f)

    # safe df
    df.loc[len(df.index)] = [layer_name, number_of_layer, channels,
                             np.mean(train_losses),
                             np.mean(val_losses),
                             np.mean(test_scores),
                             np.mean(test_scores_r), 
                             (elapsed / 60),
                             [losses_history_all],
                             [val_losses_history_all]]
    print(df)

    path = os.path.dirname(os.path.realpath(__file__))
    file_path_df = os.path.join(
        os.path.dirname(
            os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'model_evaluations_paper')
    if not os.path.isdir(file_path_df):
        os.mkdir(file_path_df)
    file_path_df = file_path_df+filename
    df.to_csv(file_path_df)


start_hyper = time.perf_counter()
epochs = 1500
filename = '/GNN_60days_samenodes_1k_test.csv'  # name for df
filename_df = 'GNN_60days_samenodes_1k_test.csv'  # name for df
save_name = 'GNN_60days_samenodes_1k_test'  # name for model
# for param in parameters:
train_and_evaluate_model(epochs, 0.001, parameters,
                         save_name, filename, filename_df)

elapsed_hyper = time.perf_counter() - start_hyper
print(
    "Time for hyperparameter testing: {:.4f} minutes".format(
        elapsed_hyper / 60))
print(
    "Time for hyperparameter testing: {:.4f} hours".format(
        elapsed_hyper / 60 / 60))
