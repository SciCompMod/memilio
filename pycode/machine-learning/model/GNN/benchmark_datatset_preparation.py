import os
import numpy as np
import pickle
import numpy as np
# import bz2file as bz2

path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(
    os.path.dirname(
        os.path.realpath(os.path.dirname(os.path.realpath(path)))),
    'data/data_GNN_400pop_one_var_damp_30days_1k')

file = open(os.path.join(path_data, 'data_secir_age_groups.pickle'), 'rb')
data_secir = pickle.load(file)

inputs = np.asarray(data_secir['inputs'])
labels = np.asarray(data_secir['labels'])
damping_days = np.asarray(data_secir['damping_day'])
damping_coeff = np.asarray(data_secir['damping_coeff'])

len_dataset = inputs.shape[1]

input_reshaped = inputs.reshape(1000, 400, 5, 48)
# labels_reshaped = labels.reshape(31, 400, 1000, 48)
# labels_30 = labels_reshaped[:30]
# labels_right_shape = labels_30.reshape(1000, 400, 30, 48)

new_labels = []
for nodes in labels:
    for dataset in nodes:
        new_labels.append(dataset[:30])


labels_right_shape = np.asarray(new_labels).reshape(1000, 400, 30, 48)

idxs = np.random.permutation(len_dataset)

split_te = int(0.9 * len_dataset)
idx_tr,  idx_te = np.split(idxs, [split_te])
inputs_train = input_reshaped[idx_tr]
labels_train = labels_right_shape[idx_tr]
damping_days_train = damping_days[0][idx_tr]
damping_coeff_train = damping_coeff[0][idx_tr]


inputs_test = input_reshaped[idx_te]
labels_test = labels_right_shape[idx_te]
damping_days_test = damping_days[0][idx_te]
damping_coeff_test = damping_coeff[0][idx_te]


data = {
    "train_inputs": inputs_train,
    "train_labels": labels_train,
    "test_inputs": inputs_test,
    "test_labels": labels_test,
    "train_damping_day": damping_days_train,
    "test_dampings_day": damping_days_test,
    "train_damping_coeff": damping_coeff_train,
    "test_damping_coeff": damping_coeff_test

}


path = os.path.dirname(os.path.realpath(__file__))

if not os.path.isdir(path):
    os.mkdir(path)

# save dict to json file
with open(os.path.join(path, 'benchmark_data_splitted_30days_onedamp.pickle'), 'wb') as f:

    pickle.dump(data, f)
