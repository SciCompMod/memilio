import time 
import os
import numpy as np
import pickle
import random 
import pandas as pd 
from progress.bar import Bar
import tensorflow as tf

############# secir simple: ODE performance #############
from memilio.surrogatemodel.ode_secir_simple.data_generation import run_secir_simulation
from memilio.surrogatemodel.ode_secir_simple.model import split_data

label_widths = [25,55,85,115,145]
df_secir_simple_ODE = pd.DataFrame(columns = ['days', 'time_per_run'])
for lw in label_widths:


        input_width = 5 
        label_width = lw
        num_runs = 1000 
        data = {
                "inputs": [],
                "labels": []
        }
        days = input_width + label_width - 1

        times = []
        for i in range(100):
                start = time.time()

                # show progess in terminal for longer runs
                # Due to the random structure, theres currently no need to shuffle the data
                bar = Bar('Number of Runs done', max=num_runs)
                for _ in range(0, num_runs):
                        data_run = run_secir_simulation(days)
                        data['inputs'].append(data_run[:input_width])
                        data['labels'].append(data_run[input_width:])
                        bar.next()
                bar.finish()

                end = time.time()
                times.append(end - start)
        print('SECIR SIMPLE : Mean time per rum for ' + str(input_width + label_width) + 'days: ' + str(np.asarray(times).mean()))
        df_secir_simple_ODE.loc[len(df_secir_simple_ODE)]=[days+1, np.asarray(times).mean()]
print(df_secir_simple_ODE)



###### secir simple LSTM performance #############
model_names = ['saved_models_secir_simple_bestLSTM_2024', 'saved_models_secir_simple_bestLSTM_2024_60days', 'saved_models_secir_simple_bestLSTM_2024_90days',
               'saved_models_secir_simple_bestLSTM_2024_120days', 'saved_models_secir_simple_bestLSTM_2024_150days']
filenames = ['data_secir_simple.pickle', 'data_secir_simple_60days.pickle', 'data_secir_simple_90days.pickle', 'data_secir_simple_120days.pickle',
             'data_secir_simple_150days.pickle']

df_secir_simple_LSTM = pd.DataFrame(columns = ['days', 'time_per_run'])
for modelname, filename in zip(model_names, filenames):
       
        secirsimple_model = tf.keras.models.load_model(os.path.join('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/', modelname))

        #path = os.path.dirname(os.path.realpath(__file__))
        path_data = ('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/data')
        filename = filename

        file = open(os.path.join(path_data,filename), 'rb')

        data = pickle.load(file)
        data_splitted = split_data(data['inputs'], data['labels'])

        test_inputs = data_splitted['test_inputs']
        test_labels = data_splitted['test_labels']
        days = test_labels.shape[1]

        times = []
        for i in range(100):

                start = time.time()
                pred = secirsimple_model.predict(test_inputs)
                end = time.time()
                times.append(end - start)
        print('SECIR SIMPLE : Mean time per rum for ' + str(days) + 'days: ' + str(np.asarray(times).mean()))
        df_secir_simple_LSTM.loc[len(df_secir_simple_LSTM)]=[days, np.asarray(times).mean()]
print(df_secir_simple_LSTM)


############# secir groups ODE performance #############

from memilio.surrogatemodel.ode_secir_groups.data_generation_nodamp import run_secir_groups_simulation, get_population

label_widths = [25,55,85,115,145]
df_secir_groups_ODE = pd.DataFrame(columns = ['days', 'time_per_run'])
for lw in label_widths:
        input_width = 5 
        label_width = lw
        num_runs = 1000 
        data = {
                "inputs": [],
                "labels": []
        }
        days = input_width + label_width - 1

        path_population = os.path.abspath(
                r"data//pydata//Germany//county_population.json")

        # Load population data
        population = get_population(path_population)

        times = []
        for i in range(100):
                start = time.time()

                # show progess in terminal for longer runs
                # Due to the random structure, theres currently no need to shuffle the data
                bar = Bar('Number of Runs done', max=num_runs)
                for _ in range(0, num_runs):
                        data_run = run_secir_groups_simulation(
                                        days, population[random.randint(0, len(population) - 1)])
                        data['inputs'].append(data_run[:input_width])
                        data['labels'].append(data_run[input_width:])
                        bar.next()
                bar.finish()

                end = time.time()
                times.append(end - start)
        #print(np.asarray(times).mean())
        print('SECIR GROUPS : Mean time per rum for ' + str(input_width + label_width) + 'days: ' + str(np.asarray(times).mean()))
        df_secir_groups_ODE.loc[len(df_secir_groups_ODE)]=[days+1, np.asarray(times).mean()]
print(df_secir_groups_ODE)


############# secir groups LSTM performance #############
model_names = ['saved_models_secir_groups_best_LSTM', 'saved_models_secir_groups_best_LSTM_60days', 'saved_models_secir_groups_best_LSTM_90days',
               'saved_models_secir_groups_best_LSTM_120days', 'saved_models_secir_groups_best_LSTM_150days']
filenames = ['data_secir_groups_30days_nodamp.pickle', 'data_secir_groups_60days_nodamp.pickle', 'data_secir_groups_90days_nodamp.pickle', 
             'data_secir_groups_120days_nodamp.pickle','data_secir_groups_150days_nodamp.pickle']

df_secir_groups_LSTM = pd.DataFrame(columns = ['days', 'time_per_run'])
for modelname, filename in zip(model_names, filenames):

        secirgroups_model =tf.keras.models.load_model( os.path.join('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/', modelname))

        path = os.path.dirname(os.path.realpath(__file__))
        path_data = os.path.join(os.path.dirname(os.path.realpath(
                os.path.dirname(os.path.realpath(path)))), 'data')
        
        filename = filename

        file = open(os.path.join(path_data,filename), 'rb')

        data = pickle.load(file)
        data_splitted = split_data(data['inputs'], data['labels'])

        test_inputs = data_splitted['test_inputs']
        test_labels = data_splitted['test_labels']
        days = test_labels.shape[1]

        times = []
        for i in range(100):
                start = time.time()
                pred = secirgroups_model.predict(test_inputs)
                end = time.time()
                times.append(end - start)
        print('SECIR GROUPS : Mean time per rum for ' + str(days) + 'days: ' + str(np.asarray(times).mean()))
        df_secir_groups_LSTM.loc[len(df_secir_groups_LSTM)]=[days, np.asarray(times).mean()]
print(df_secir_groups_LSTM)



######## ODE secir groups one damp #################

from memilio.surrogatemodel.ode_secir_groups.data_generation import run_secir_groups_simulation, get_population

input_width = 5 
label_width = 25
num_runs = 1000 
data = {
        "inputs": [],
        "labels": [], 
        "contact_matrix":[],
        "damping_day":[]

    }
days = input_width + label_width - 1

path_population = os.path.abspath(
        r"data//pydata//Germany//county_population.json")

# Load population data
population = get_population(path_population)

times = []
for i in range(100):
    start = time.time()


    bar = Bar('Number of Runs done', max=num_runs)

    
    for i in range(0, num_runs):

        # Generate a random damping day
        damping_day = random.randrange(
            input_width, input_width+label_width)
        

        data_run, damped_contact_matrix = run_secir_groups_simulation(
            days, damping_day, population[random.randint(0, len(population) - 1)])
        #data_run = run_secir_groups_simulation(
        #    days, population[random.randint(0, len(population) - 1)])
        data['inputs'].append(data_run[:input_width])
        data['labels'].append(data_run[input_width:])
        data['contact_matrix'].append(np.array(damped_contact_matrix))
        data['damping_day'].append(damping_day)
        bar.next()
    bar.finish()
    end = time.time()
    times.append(end - start)
print('SECIR GROUPS ONE DAMP: Mean time per rum for ' + str(days+1) + 'days: ' + str(np.asarray(times).mean()))



############### one damp secir groups LSTM #########################
from memilio.surrogatemodel.ode_secir_groups.model import split_contact_matrices, flat_input

secirgroups_model = tf.keras.models.load_model('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/saved_models_secir_groups_onedamp_best_LSTM__30days')

path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')
    
filename = "data_secir_groups_30days_onevardamp_baseline.pickle"

file = open(os.path.join(path_data,filename), 'rb')

data = pickle.load(file)
data_splitted = split_data(data['inputs'], data['labels'])

test_inputs_compartments = (data_splitted["test_inputs"])
test_labels = (data_splitted["test_labels"])
contact_matrices = split_contact_matrices(
            tf.stack(data["contact_matrix"]))
contact_matrices_test = flat_input(contact_matrices['test'])
n = np.array(data['damping_day']).shape[0]
test_days = data['damping_day'][int(n*0.9):]

new_contact_test = []
for i in contact_matrices_test:
            new_contact_test.extend([i for j in range(5)])

new_contact_test = tf.reshape(tf.stack(new_contact_test), [
            contact_matrices_test.shape[0], 5, contact_matrices_test.shape[1]])

new_damping_days_test = []
for i in test_days:
            new_damping_days_test.extend([i for j in range(5)])
new_damping_days_test = tf.reshape(
            tf.stack(new_damping_days_test),
            [test_inputs_compartments.shape[0],
             5, 1])

test_inputs = tf.concat(
            (tf.cast(test_inputs_compartments, tf.float16),
             tf.cast(new_contact_test, tf.float16),
             tf.cast(new_damping_days_test, tf.float16)),
            axis=2)
times = []
for i in range(100):
        start = time.time()
        pred = secirgroups_model.predict(test_inputs)
        end = time.time()
        times.append(end - start)
       
print('SECIR GROUPS ONE  DAMP  Mean time per rum for 30 days ' + str(np.asarray(times).mean()))





############## ODE multiple dampings ####################
from memilio.surrogatemodel.ode_secir_many_dampings.data_generation_with_adjusted_dampdrawing import run_secir_groups_simulation, generate_dampings_withshadowdamp, get_population

df_secir_groups_damps = pd.DataFrame(columns = ['dampings', 'time_per_run']) 
number_of_dampings = [2,3,4,5]
for nd in number_of_dampings:
        data = {
                "inputs": [],
                "labels": [],
                "contact_matrix": [],
                "damping_day": []
                }

        input_width = 5
        label_width = 95
        num_runs = 1000
        # The number of days is the same as the sum of input and label width.
        # Since the first day of the input is day 0, we still need to subtract 1.
        days = input_width + label_width - 1


        path_population = os.path.abspath(
                r"data//pydata//Germany//county_population.json")
        # Load population data
        population = get_population(path_population)

        damping_days = generate_dampings_withshadowdamp(number_of_dampings = nd, days= days, min_distance=7, min_damping_day=input_width, n_runs = num_runs)

        # show progess in terminal for longer runs
        # Due to the random structure, theres currently no need to shuffle the data
        times = []
        for i in range(100):
                start = time.time()
                bar = Bar('Number of Runs done', max=num_runs)
                for i in range(0, num_runs):
                        data_run, damped_contact_matrix, damping_days_s = run_secir_groups_simulation(
                        days, number_of_dampings,  damping_days[i], population[random.randint(0, len(population) - 1)])
                        #data_run = run_secir_groups_simulation(
                        #    days, population[random.randint(0, len(population) - 1)])
                        data['inputs'].append(data_run[:input_width])
                        data['labels'].append(data_run[input_width:])
                        data['contact_matrix'].append(np.array(damped_contact_matrix))
                        data['damping_day'].append(damping_days_s)
                        bar.next()
                bar.finish()
                end = time.time()
                times.append(end - start)
        print('SECIR GROUPS DAMPINGS : Mean time per rum for ' + str(input_width + label_width+1) + 'days: ' + str(np.asarray(times).mean()))
        df_secir_groups_damps.loc[len(df_secir_groups_damps)]=[int(nd), np.asarray(times).mean()]
print(df_secir_groups_damps)

########### secir groups multiple dampings LSTM #######
from memilio.surrogatemodel.ode_secir_many_dampings.model_2024 import split_contact_matrices, split_data, flat_input

model_names = ['saved_models_secir_grpups_2damp_LSTM_100days', 'saved_models_secir_grpups_3damp_LSTM_100days', 'saved_models_secir_grpups_4damp_LSTM_100days',
               'saved_models_secir_grpups_5damp_LSTM_100days']
filenames = ['data_secir_groups_100days_2damp.pickle', 'data_secir_groups_100days_3damp.pickle', 'data_secir_groups_100days_4damp.pickle', 
             'data_secir_groups_100days_5damp.pickle']

df_secir_groups_dampings_LSTM = pd.DataFrame(columns = ['dampings', 'time_per_run'])
for modelname, filename in zip(model_names, filenames):

        secirgroups_model =tf.keras.models.load_model( os.path.join('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/', modelname))

        filename = filename
        path_data = os.path.join('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/data', filename)
        
        file = open(os.path.join(path_data), 'rb')

        data = pickle.load(file)
        data_splitted = split_data(data['inputs'], data['labels'])

        test_inputs_compartments = (data_splitted["test_inputs"])
        test_labels = (data_splitted["test_labels"])
        test_labels = tf.reshape(test_labels, [test_labels.shape[0], test_labels.shape[1]* test_labels.shape[2]])

        contact_matrices = split_contact_matrices(
            tf.stack(data["contact_matrix"]))
        contact_matrices_test = flat_input(contact_matrices['test'])

        n = np.array(data['damping_day']).shape[0]
        test_days = data['damping_day'][int(n*0.9):]

        new_contact_test = []
        for i in contact_matrices_test:
            new_contact_test.extend([i for j in range(5)])

        new_contact_test = tf.reshape(tf.stack(new_contact_test), [
            contact_matrices_test.shape[0], 5, contact_matrices_test.shape[1]])

        new_damping_days_test = []
        for i in test_days:
            new_damping_days_test.extend([i for j in range(5)])

        number_of_dampings = new_damping_days_test[0].shape[0]

        new_damping_days_test = tf.reshape(
            tf.stack(new_damping_days_test),
            [test_inputs_compartments.shape[0],
             5, number_of_dampings])

        test_inputs = tf.concat(
            (tf.cast(test_inputs_compartments, tf.float16),
             tf.cast(new_contact_test, tf.float16),
             tf.cast(new_damping_days_test, tf.float16)),
            axis=2)

        times = []
        for i in range(100):
                start = time.time()
                pred = secirgroups_model.predict(test_inputs)
                end = time.time()
                times.append(end - start)
        print('SECIR GROUPS MULTIPLE DAMP LSTM : Mean time per rum for ' + str(number_of_dampings) + 'dampings: ' + str(np.asarray(times).mean()))
        df_secir_groups_dampings_LSTM.loc[len(df_secir_groups_dampings_LSTM)]=[number_of_dampings, np.asarray(times).mean()]
print(df_secir_groups_dampings_LSTM)





############### anlyizing results ###############
