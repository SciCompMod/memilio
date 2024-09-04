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
    'data_GNN_with_1_dampings')

file = open(os.path.join(path_data, 'data_secir_age_groups.pickle'), 'rb')
data_secir = pickle.load(file)

len_dataset = data_secir['inputs'].shape[0]
number_of_nodes = data_secir['inputs'].shape[1]

shape_input_flat = np.asarray(
    data_secir['inputs']).shape[2]*np.asarray(data_secir['inputs']).shape[3]
shape_labels_flat = np.asarray(
    data_secir['labels']).shape[2]*np.asarray(data_secir['labels']).shape[3]


damping_factors = data_secir['damping_coeff']  # one for every run
damping_days = data_secir['damping_day']  # one for every rum
contact_matrices = np.concatenate(data_secir['damped_matrix'])

n_dampings = len(data_secir['damping_day'][0])


new_inputs = np.asarray(
    data_secir['inputs']).reshape(
    len_dataset, number_of_nodes, shape_input_flat)
new_labels = np.asarray(data_secir['labels']).reshape(
    len_dataset, number_of_nodes, shape_labels_flat)

inputs_with_damp = np.dstack((new_inputs, np.repeat(np.asarray(damping_factors)[:, np.newaxis, :], 400, axis=1),
                              (np.repeat(np.asarray(damping_days)
                               [:, np.newaxis, :], 400, axis=1)),
                              (np.repeat(np.asarray(contact_matrices)[:, np.newaxis, :], 400, axis=1).reshape([len_dataset, number_of_nodes, 36*n_dampings]))))


######## open commuter data #########
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(path,
                         'data')
commuter_file = open(os.path.join(
    path_data, 'commuter_migration_scaled.txt'), 'rb')
commuter_data = pd.read_csv(commuter_file, sep=" ", header=None)
sub_matrix = commuter_data.iloc[:number_of_nodes, 0:number_of_nodes]

adjacency_matrix = np.asarray(sub_matrix)
adjacency_matrix[adjacency_matrix > 0] = 1
node_features = inputs_with_damp

node_labels = new_labels

layer = ARMAConv
number_of_layers = 1
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
            self.a = rescale_laplacian(
                normalized_laplacian(adjacency_matrix))

            return [spektral.data.Graph(x=x, y=y) for x, y in zip(node_features, node_labels)]

            super().__init__(**kwargs)

    # data = MyDataset()
    data = MyDataset(transforms=NormalizeAdj())
    batch_size = 32
    epochs = epochs
    es_patience = 100  # Patience for early stopping

    class Net(Model):
        def __init__(self):
            super().__init__()
            self.conv1 = layer(channels, activation='elu')
            self.dense = Dense(data.n_labels, activation="linear")

        def call(self, inputs):
            x, a = inputs

            a = np.asarray(a)
            x = self.conv1([x, a])
            output = self.dense(x)

            return output

    optimizer = Adam(learning_rate=learning_rate)
    # optimizer = optimizer(learning_rate=learning_rate)
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

    n_days = int(new_labels.shape[2]/48)

    test_scores = []
    train_losses = []
    val_losses = []

    losses_history_all = []
    val_losses_history_all = []

    start = time.perf_counter()

    model = Net()
    # idxs = np.random.permutation(len(data)) #random
    idxs = np.arange(len(data))  # same for each run
    split_va, split_te = int(0.8 * len(data)), int(0.9 * len(data))
    idx_tr, idx_va, idx_te = np.split(idxs, [split_va, split_te])

    data_tr = data[idx_tr]
    data_va = data[idx_va]
    data_te = data[idx_te]

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

    ##########################################################################
    ##########################################################################
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

    # save df
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
filename = '/GNN_paper_onedamp.csv'  # name for df
save_name = 'GNN_paper_onedamp'  # name for model
# for param in parameters:
train_and_evaluate_model(epochs, 0.001, parameters, save_name, filename)

elapsed_hyper = time.perf_counter() - start_hyper
print(
    "Time for hyperparameter testing: {:.4f} minutes".format(
        elapsed_hyper / 60))
print(
    "Time for hyperparameter testing: {:.4f} hours".format(
        elapsed_hyper / 60 / 60))
