from data_generator import *


def generate_data(num_runs, input_width,  days, number_of_dampings,
                  save_data=True):
    """! Generate dataset by calling run_secir_simulation (num_runs)-often

   @param num_runs Number of times, the function run_secir_simulation is called.
   @param path Path, where the datasets are stored.
   @param input_width number of time steps used for model input.
   @param label_width number of time steps (days) used as model output/label.  
   @param save_data Option to deactivate the save of the dataset. Per default true.
   """
    data = {
        "inputs": [],
        "labels": [],
        "damping_days": [],
        "contact_matrix": []
    }

    #days = damping_days[-1]+20

    population = get_population()

    # show progess in terminal for longer runs
    # Due to the random structure, theres currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)

    for _ in range(num_runs):

        damping_days = [0]*number_of_dampings
        damping_days[0] = randrange(5, 20)
        for i in range(1, number_of_dampings):
            damping_days[i] = damping_days[(i-1)]+randrange(10, 20)

        data_run, damping_matrix = run_secir_groups_simulation(
            days, damping_days, population[random.randint(0, len(population) - 1)])

        # for i in damping_matrix:
        #    data["contact_matrix"].append(i)

        inputs = data_run[:input_width]
        data["inputs"].append(inputs)

        data["labels"].append(data_run[input_width:])
        data['contact_matrix'].append(damping_matrix)
        data['damping_days'].append(damping_days)

        bar.next()

    bar.finish()

    if save_data:

        transformer = FunctionTransformer(np.log1p, validate=True)
        inputs = np.asarray(data['inputs']).transpose(2, 0, 1).reshape(48, -1)
        scaled_inputs = transformer.transform(inputs)
        scaled_inputs = scaled_inputs.transpose().reshape(num_runs, input_width, 48)
        scaled_inputs_list = scaled_inputs.tolist()

        labels = np.asarray(data['labels']).transpose(2, 0, 1).reshape(48, -1)
        scaled_labels = transformer.transform(labels)
        scaled_labels = scaled_labels.transpose().reshape(
            num_runs, np.asarray(data['labels']).shape[1], 48)
        scaled_labels_list = scaled_labels.tolist()

        # cast dfs to tensors
        data['inputs'] = (scaled_inputs_list)
        data['labels'] = (scaled_labels_list)

    return data


input_width = 5
num_runs = 10
days = 200
data_set_names = ['data1', 'data2', 'data3', 'data4',
                  'data5', 'data6', 'data7', 'data8', 'data9', 'data10']
damping_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

datasets = []
for damp, name in zip(damping_list, data_set_names):
    name = generate_data(num_runs, input_width,
                         days, damp)
    datasets.append(name)


inputs = []
for i in range(len(datasets)):
    inputs.append(datasets[i]['inputs'])
new_inputs = np.asarray(inputs).reshape(num_runs * 10, input_width, 48)


labels = []
for i in range(len(datasets)):
    labels.append(datasets[i]['labels'])
new_labels = np.asarray(labels).reshape(num_runs * 10, days-input_width, 48)

matrices = []
for i in range(len(datasets)):
    matrices.append(datasets[i]['contact_matrix'])
new_matrices = np.asarray(matrices, dtype=object).reshape(num_runs*10).tolist()

damping_days = []
for i in range(len(datasets)):
    damping_days.append(datasets[i]['damping_days'])
new_damping_days = np.asarray(damping_days, dtype=object).reshape(num_runs*10)


full_dataset = {
    "inputs": new_inputs,
    "labels": new_labels,
    "damping_days": new_damping_days,
    "contact_matrix": new_matrices
}

path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
    os.path.dirname(os.path.realpath(path)))), 'data_multiple_dampings')

if not os.path.isdir(path_data):
    os.mkdir(path_data)

# save dict to json file
with open(os.path.join(path_data, 'data_multiple_dampings.pickle'), 'wb') as f:
    pickle.dump(full_dataset, f)
