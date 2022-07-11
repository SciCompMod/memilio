from memilio.simulation import UncertainContactMatrix, ContactMatrix, Damping
from memilio.simulation.secir import SecirModel, simulate, AgeGroup, Index_InfectionState, SecirSimulation
from memilio.simulation.secir import InfectionState as State
import numpy as np
import pandas as pd
from datetime import date
from math import ceil
import random
import os
from progress.bar import Bar # pip install progess
from skimpy import skim # show overview over dataset at the end. Sometimes really helpful
from sklearn.model_selection import train_test_split
from sklearn.compose import make_column_transformer
from sklearn.preprocessing import MinMaxScaler
import tensorflow as tf
import matplotlib.pyplot as plt
from tensorflow.keras import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.optimizers import Adam
from tensorflow.keras import backend as K



def run_secir_simulation():
    """! Slightly modified example of Secir Simple.

   Here, the initial values are choosen randomly, so the model is no longer deterministic.
   Necessary to create the dataset for later training process.
   This method is called within a loop in the function 'generate_data', which also sets the number of runs.

   """

    # Define Comartment names
    compartments = ['Susceptible', 'Exposed', 'Carrier',
                    'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead']

    # Define random population number in intervall (0,500000)
    populations = [500000 * random.random()]

    days = 100  # number of days to simulate

    # since we defined the number of days, we already know the dimension of our dataset.
    # We aim to save the compartment population for each day.
    # i.e. for days = 100 (+1 for t0) and len(compartments = 8) -> 101 x  8 array
    # we also add one column for the damping factor
    data = np.zeros((days+1,len(compartments) + 1))

    start_day = 1
    start_month = 1
    start_year = 2019
    dt = 0.1
    num_groups = 1
    num_compartments = len(compartments)

    # Initialize Parameters
    model = SecirModel(1)

    A0 = AgeGroup(0)


    # Define random values for damping coef and start date
    damping_factor = random.random()
    damping_date = ceil(random.random() * days) # scale with days since random is in Interval (0,1)

    # save damping factor in vector. -1 for index correction
    data[damping_date-1:,-1] = damping_factor

    # Set parameters

    # Compartment transition duration
    model.parameters.IncubationTime[A0] = 5.2  # R_2^(-1)+R_3^(-1)
    model.parameters.InfectiousTimeMild[A0] = 6.  # 4-14  (=R4^(-1))
    # 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
    model.parameters.SerialInterval[A0] = 4.2
    model.parameters.HospitalizedToHomeTime[A0] = 12.  # 7-16 (=R5^(-1))
    model.parameters.HomeToHospitalizedTime[A0] = 5.  # 2.5-7 (=R6^(-1))
    model.parameters.HospitalizedToICUTime[A0] = 2.  # 1-3.5 (=R7^(-1))
    model.parameters.ICUToHomeTime[A0] = 8.  # 5-16 (=R8^(-1))
    model.parameters.ICUToDeathTime[A0] = 5.  # 3.5-7 (=R5^(-1))

    # Initial number of people in each compartment with random numbers 
    model.populations[A0, Index_InfectionState(State.Exposed)] = 250 * random.random()
    model.populations[A0, Index_InfectionState(State.Carrier)] = 120 * random.random()
    model.populations[A0, Index_InfectionState(State.Infected)] = 120 * random.random()
    model.populations[A0, Index_InfectionState(State.Hospitalized)] = 50 * random.random()
    model.populations[A0, Index_InfectionState(State.ICU)] = 30 * random.random()
    model.populations[A0, Index_InfectionState(State.Recovered)] = 30 * random.random()
    model.populations[A0, Index_InfectionState(State.Dead)] = 0
    model.populations.set_difference_from_total(
        (A0, Index_InfectionState(State.Susceptible)), populations[0])

    # Compartment transition propabilities
    model.parameters.RelativeCarrierInfectability[A0] = 0.67
    model.parameters.InfectionProbabilityFromContact[A0] = 1.0
    model.parameters.AsymptomaticCasesPerInfectious[A0] = 0.09  # 0.01-0.16
    model.parameters.RiskOfInfectionFromSymptomatic[A0] = 0.25  # 0.05-0.5
    model.parameters.HospitalizedCasesPerInfectious[A0] = 0.2  # 0.1-0.35
    model.parameters.ICUCasesPerHospitalized[A0] = 0.25  # 0.15-0.4
    model.parameters.DeathsPerICU[A0] = 0.3  # 0.15-0.77
    # twice the value of RiskOfInfectionFromSymptomatic
    model.parameters.MaxRiskOfInfectionFromSymptomatic[A0] = 0.5

    model.parameters.StartDay = (
        date(start_year, start_month, start_day) - date(start_year, 1, 1)).days

    # model.parameters.ContactPatterns.cont_freq_mat[0] = ContactMatrix(np.r_[0.5])
    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
        (num_groups, num_groups)) * 1
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.ones(
        (num_groups, num_groups)) * 0
    model.parameters.ContactPatterns.cont_freq_mat.add_damping(
        Damping(coeffs=np.r_[damping_factor], t=damping_date, level=0, type=0))

    # Apply mathematical constraints to parameters
    model.apply_constraints()

    # Run Simulation
    # store data for t = 0
    data[0,:-1] = model.populations.get_compartments()
    # TODO: Do we need to build the whole Model for new, if we only change the number of days?
    for day in range(1, days+1):
        result = simulate(0, day, dt, model)

        # Maybe round data
        # data[day,:] = np.ceil(result_day)
        data[day,:-1] = result.get_last_value()

    # create labels
    # the labels are the population inside the copartmens on the next time step.
    labels = data[1:,:-1]

    return data[:-1,:], labels



def generate_data(num_runs, save_data=True):
    """! Generate dataset by calling run_secir_simulation (num_runs)-often

   @param num_runs Number of times, the function run_secir_simulation is called.
   @param save_data Option to deactivate the save of the dataset. Per default true.
   """
    data = []
    labels = []

    # show progess in terminal for longer runs
    bar = Bar('Number of Runs done', max=num_runs)
    for _ in range(0, num_runs):
        data_run, labels_run = run_secir_simulation()
        data.append(data_run)
        labels.append(labels_run)
        bar.next()
    
    bar.finish()

    if save_data:
    
        columns = ['Susceptible', 'Exposed', 'Carrier',
                        'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead', 'Damping Factor']

        dataset = pd.DataFrame(np.concatenate(data), columns=columns)

        labels_set = pd.DataFrame(np.concatenate(labels), columns=columns[:-1])

        # short overview over dataset
        skim(dataset)

        # save train dataset as csv file
        path = os.path.dirname(os.path.realpath(__file__))
        dataset.to_csv(
            os.path.join(os.path.dirname(os.path.realpath(path)), 'data', 'traindata_secir_simple.txt'),
            sep=' ', index=False,
            header=True)

        labels_set.to_csv(
            os.path.join(os.path.dirname(os.path.realpath(path)), 'data', 'labels_secir_simple.txt'),
            sep=' ', index=False,
            header=True)


def visualize_predictions(
               testdata ,
               testlabels,
               predictions):
    """! Visualize the predictions of the trained model to evaluate the performance.

    @param testdata Testdata thats different than the traindata and not known for the model.
    @param testlabels Same as testdata but with labels.
    @param predictions Output of the model with testdata as input
    """
    # Todo: Visualization. We need the data in correct order.
    return None

# create and train an neural network based on the secir simple example. 
# if no dataset is already create, we build one with default 50 runs
def network_secir_simple(path, epochs=30, num_runs_traindata=50):
    """! Generate the model and train with the created dataset.

    If the dataset is not created yet, we create a new one with default value
    'nums_runs_traindata = 50'. Since the number of training runs is freely selectable with
    'epochs', unit tests are easier to implement afterwards.

   @param path Path to dataset
   @param epochs Number of epochs in the train process. 
   @param number of runs to create dataset. Only used, if no dataset is found.
   """

    if not os.path.isfile(os.path.join(path, 'traindata_secir_simple.txt')) or \
        not os.path.isfile(os.path.join(path, 'traindata_secir_simple.txt')):

        generate_data(num_runs)

    input = pd.read_csv(
            os.path.join(os.path.dirname(os.path.realpath(path)), 'data', 'traindata_secir_simple.txt'),
            sep=' ')
    
    labels = pd.read_csv(
            os.path.join(os.path.dirname(os.path.realpath(path)), 'data', 'labels_secir_simple.txt'),
            sep=' ')

    # split in train and test/valid data (ratio 80/20).
    # shuffling is also done by this function
    X_train, X_test, Y_train, Y_test = train_test_split(input, labels, test_size=0.2, random_state=42)

    # define transformation to prepare data for training. Mainly min-max scaling, OneHotEncoder for classification params
    transformer = make_column_transformer(
    (MinMaxScaler(), 
        ['Susceptible', 'Exposed', 'Carrier',
                    'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead', 'Damping Factor']))

    transformer.fit(X_train)
    X_train = transformer.transform(X_train)
    X_test = transformer.transform(X_test)

    # Regarding the network, we use a fully connected one with several hidden layers.
    # Since we have a regression problem, we use the mse loss function.
    # use the leakyRelu as activation function, except in the output layer.
    
    # set seed for pseudo-randomization
    tf.random.set_seed(42)

    # we want to predict the population in the diferent compartments. 
    # Todo: Maybe add drop out or other mechanics against overfitting.
    model = Sequential([
        Dense(128, activation=tf.keras.layers.LeakyReLU(alpha=0.01)),
        Dense(128, activation=tf.keras.layers.LeakyReLU(alpha=0.01)),
        Dense(128, activation=tf.keras.layers.LeakyReLU(alpha=0.01)),
        Dense(64, activation=tf.keras.layers.LeakyReLU(alpha=0.01)),
        Dense(8)
    ])

    model.compile(
        loss=tf.keras.losses.mse,
        optimizer=Adam(),
        metrics=['mae']
    )

    model.fit(X_train, Y_train, epochs=epochs)

    preds = model.predict(X_test)
    visualize_predictions(X_test,Y_test,preds)


if __name__ == "__main__":
    # TODO: Save contact matrix depending on the damping.
    # In the actual state it might be enough to save the regular one and the damping
    num_runs = 50
    generate_data(num_runs)

    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(path)), 'data')
    network_secir_simple(path_data)
    
    
