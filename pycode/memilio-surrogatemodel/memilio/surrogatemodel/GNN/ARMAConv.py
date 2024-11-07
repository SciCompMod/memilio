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
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(
    os.path.dirname(
        os.path.realpath(os.path.dirname(os.path.realpath(path)))),
    'data_GNN_paper')

file = open(os.path.join(path_data, 'GNN_data_30days_1k.pickle'), 'rb')
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

# our dataset only contain 1000 samples. The variance in these datasets is relatively high
# variance. We want to augment the trainig dataset by copying the data (first without any changes)
test_inputs = data_secir['inputs'][int(0.8*len(data_secir['inputs'])):]
test_labels = data_secir['labels'][int(0.8*len(data_secir['labels'])):]

# Repeat along the first axis 5 times
train_inputs = np.tile(data_secir['inputs'][:int(
    0.8*len(data_secir['inputs']))], (5, 1, 1, 1))
train_labels = np.tile(data_secir['labels'][:int(
    0.8*len(data_secir['labels']))], (5, 1, 1, 1))

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
number_of_layers = 2
number_of_channels = 1024
parameters = [layer, number_of_layers, number_of_channels]

df = pd.DataFrame(
    columns=['layer', 'number_of_layers', 'channels', 'kfold_train',
             'kfold_val', 'kfold_test', 'training_time',
             'train_losses', 'val_losses'])


def train_and_evaluate_model(
        epochs, learning_rate, param, save_name, filename):

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
    es_patience = 100  # Patience for early stopping

    class Net(Model):
        def __init__(self):
            super().__init__()
            self.conv1 = layer(channels, order=2, iterations=2,
                               share_weights=True,  activation='elu')
            self.conv2 = layer(channels, order=2, iterations=2,
                               share_weights=True,  activation='elu')
            # not important if data_train or data_test, as we only need label size
            self.dense = Dense(data_train.n_labels, activation="linear")

        def call(self, inputs):
            x, a = inputs
            a = np.asarray(a)
            x = self.conv1([x, a])
            x = self.conv2([x, a])
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

    test_scores = []
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

    print(
        "Done. Test loss: {:.4f}. Test acc: {:.2f}".format(
            test_loss, test_acc))
    test_scores.append(test_loss)
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
    print("K-Fold Train Score:{}".format(np.mean(train_losses)))
    print("K-Fold Validation Score:{}".format(np.mean(val_losses)))
    print("K-Fold Test Score: {}".format(np.mean(test_scores)))

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
filename = '/GNN_data_30days_1k_dataaug.csv'  # name for df
save_name = 'GNN_data_30days_1k_dataaug'  # name for model
# for param in parameters:
train_and_evaluate_model(epochs, 0.001, parameters, save_name, filename)

elapsed_hyper = time.perf_counter() - start_hyper
print(
    "Time for hyperparameter testing: {:.4f} minutes".format(
        elapsed_hyper / 60))
print(
    "Time for hyperparameter testing: {:.4f} hours".format(
        elapsed_hyper / 60 / 60))
