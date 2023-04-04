from cProfile import label
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
from sklearn.model_selection import train_test_split
from sklearn.compose import make_column_transformer
from sklearn.preprocessing import MinMaxScaler
import tensorflow as tf
import matplotlib.pyplot as plt
import memilio.simulation as mio
import pickle
import seaborn as sns # plot after normalization

# very simple data generator 
# no damping considered 
# sclaing instead of normalizing 


def run_secir_simulation():
    """! Slightly modified example of Secir Simple.

   Here, the initial values are choosen randomly, so the model is no longer deterministic.
   Necessary to create the dataset for later training process.
   This method is called within a loop in the function 'generate_data', which also sets the number of runs.

   """
    mio.set_log_level(mio.LogLevel.Off)

    # Define Comartment names
    compartments = ['Susceptible', 'Exposed', 'Carrier',
                    'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead']

    # Define random population number in intervall (0,500000)
    populations = [500000 * random.random()]

    days = 35  # number of days to simulate

    # since we defined the number of days, we already know the dimension of our dataset.
    # We aim to save the compartment population for each day.
    # i.e. for days = 100 (+1 for t0) and len(compartments = 8) -> 101 x  8 array
    # we also add one column for the damping factor
    data = np.zeros((days,len(compartments) + 1))

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
    #data[damping_date-1:,-1] = damping_factor

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
    for day in range(1, days):
        result = simulate(0, day, dt, model)

        # Maybe round data
        # data[day,:] = np.ceil(result_day)
        data[day,:-1] = result.get_last_value()

    
    df_data = pd.DataFrame(data = data)
    df_data.drop(8, axis = 1, inplace = True)
    damping_factor_array =  [damping_factor]* days
    damping_date_array = [damping_date]* days

    df_data['damping_factor'] = damping_factor_array
    df_data['damping_date'] = damping_date_array

    data = df_data.values 





    return data


def get_scaler(num_runs):
    """! Generates a big dataset ordered by features in order to define the values needed for normalization"""


    # get the DataFrame 
    data = pd.DataFrame(columns =['G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'factor', 'date'])

    for _ in range(0, num_runs):
        data_run = run_secir_simulation()
        for i in range (0, len(data_run)):
            data.loc[len(data)] = data_run[i].tolist()

    # build scaler based on the data 
    scaler = MinMaxScaler()
    scaler.fit(data.values)



    return scaler




def generate_data(num_runs, path, save_data=True):
    """! Generate dataset by calling run_secir_simulation (num_runs)-often

   @param num_runs Number of times, the function run_secir_simulation is called.
   @param path Path, where the datasets are stored.
   @param save_data Option to deactivate the save of the dataset. Per default true.
   """

    scaler = get_scaler(1000)
    data =[]


    # how many timesteps should be input for the ML model? 
    n_inputs  = 5  



    # show progess in terminal for longer runs
    bar = Bar('Number of Runs done', max=num_runs)
    for _ in range(0, num_runs):
        data_run = run_secir_simulation()


        to_scale = data_run[:n_inputs]
        inputs_run = scaler.transform(to_scale)


        df_data = pd.DataFrame(data = data_run )
        df_data.drop([8,9], axis = 1, inplace = True )
        data_nodamp = df_data.values
        labels_run = data_nodamp[n_inputs:]

        data.append(inputs_run.reshape(1,n_inputs, 10).tolist())
        data[_].append(labels_run.tolist())


        bar.next()
    
    bar.finish()

    if save_data:
    
        columns = ['Inputs', 'Labels']

        dataset = pd.DataFrame(data, columns=columns)


        # short overview over dataset
        #print(dataset.describe().transpose())



        # check if data directory exists. If necessary create it.
        if not os.path.isdir(path): 
            os.mkdir(path)

        # save dict to json file
        with open(os.path.join(path, 'data_secir_simple.pickle'), 'wb') as f:
            pickle.dump(data, f)




        print("Data saved in " + path)


print('x')
if __name__ == "__main__":
    # TODO: Save contact matrix depending on the damping.
    # In the actual state it might be enough to save the regular one and the damping

    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(path)), 'data_withdamping')



    num_runs = 100000
    generate_data(num_runs, path_data, save_data=True)
    
