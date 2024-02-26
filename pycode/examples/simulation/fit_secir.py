import argparse
from datetime import date, datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import json

from memilio.simulation import AgeGroup, ContactMatrix, Damping, UncertainContactMatrix
from memilio.simulation.secir import Index_InfectionState
from memilio.simulation.secir import InfectionState as State
from memilio.simulation.secir import (Model, Simulation,
                                      interpolate_simulation_result, simulate)
import lmfit
import os
import math


def linearize(x, xleft, xright, yleft, yright):
    if xleft == xright:
         return yleft
    return yleft+(x-xleft)*(yright-yleft)/(xright-xleft)
     

def fit_func_secir(x,RelativeTransmissionNoSymptoms,
                    TransmissionProbabilityOnContact,
                    RecoveredPerInfectedNoSymptoms,
                    RiskOfInfectionFromSymptomatic,
                    SeverePerInfectedSymptoms,
                    CriticalPerSevere,
                    DeathsPerCritical,
                    S0,E0,C0,I0,H0,U0,R0,D0):
    print("Function called")
    print("Parameters: ")
    print(RelativeTransmissionNoSymptoms)
    print(TransmissionProbabilityOnContact)
    print(RecoveredPerInfectedNoSymptoms)
    print(RiskOfInfectionFromSymptomatic)
    print(SeverePerInfectedSymptoms)
    print(CriticalPerSevere)
    print(DeathsPerCritical)
    print(S0,E0,C0,I0,H0,U0,R0,D0)
    # Define Comartment names
    compartments = [
        'Susceptible', 'Exposed', 'InfectedNoSymptoms', 'InfectedSymptoms',
        'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

    days = 100  # number of days to simulate
    start_day = 1
    start_month = 1
    start_year = 2019
    dt = 0.1
    num_groups = 1
    num_compartments = len(compartments)

    A0 = AgeGroup(0)

    # Initialize Parameters
    model = Model(1)

    model.parameters.IncubationTime[A0] = 5.2
    model.parameters.TimeInfectedSymptoms[A0] = 6.
    # 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
    model.parameters.SerialInterval[A0] = 4.2
    model.parameters.TimeInfectedSevere[A0] = 12.  # 7-16 (=R5^(-1))
    model.parameters.TimeInfectedCritical[A0] = 8.

    model.parameters.StartDay = (
        date(start_year, start_month, start_day) - date(start_year, 1, 1)).days

     # Compartment transition propabilities
    model.parameters.RelativeTransmissionNoSymptoms[A0] = RelativeTransmissionNoSymptoms
    model.parameters.TransmissionProbabilityOnContact[A0] = TransmissionProbabilityOnContact
    model.parameters.RecoveredPerInfectedNoSymptoms[A0] = RecoveredPerInfectedNoSymptoms
    model.parameters.RiskOfInfectionFromSymptomatic[A0] = RiskOfInfectionFromSymptomatic
    model.parameters.SeverePerInfectedSymptoms[A0] = SeverePerInfectedSymptoms
    model.parameters.CriticalPerSevere[A0] = CriticalPerSevere
    model.parameters.DeathsPerCritical[A0] = DeathsPerCritical
    # twice the value of RiskOfInfectionFromSymptomatic
    model.parameters.MaxRiskOfInfectionFromSymptomatic[A0] = 2*RiskOfInfectionFromSymptomatic

    # Initial number of people in each compartment
    model.populations[A0, State.Susceptible] = S0
    model.populations[A0, State.Exposed] = E0
    model.populations[A0, State.InfectedNoSymptoms] = C0
    model.populations[A0, State.InfectedNoSymptomsConfirmed] = 0
    model.populations[A0, State.InfectedSymptoms] = I0
    model.populations[A0, State.InfectedSymptomsConfirmed] = 0
    model.populations[A0, State.InfectedSevere] = H0
    model.populations[A0, State.InfectedCritical] = U0
    model.populations[A0, State.Recovered] = R0
    model.populations[A0, State.Dead] = D0

    # model.parameters.ContactPatterns.cont_freq_mat[0] = ContactMatrix(np.r_[0.5])
    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
        (num_groups, num_groups)) * 1
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.ones(
        (num_groups, num_groups)) * 0
    model.parameters.ContactPatterns.cont_freq_mat.add_damping(
        Damping(coeffs=np.r_[0.9], t=30.0, level=0, type=0))

    # Apply mathematical constraints to parameters
    model.apply_constraints()

    # Run Simulation
    result_sim = simulate(0, days, dt, model)

    # The resulting time series has a number of timepoints different than 100 -> format of our data
    num_time_points = result_sim.get_num_time_points()

    # Make x iterable if it is only a float/int
    if not hasattr(x, '__iter__'):
        x = np.array([x])

    result_arr = np.zeros(len(x))

    # rescale x
    x *= (num_time_points-1)/99

    counter = 0
    for x0 in x:
        if 0 <= x0 < num_time_points: # Number of Deaths is queried
                result_arr[counter] = linearize(x0,math.floor(x0),math.ceil(x0),result_sim.get_value(math.floor(x0))[9], result_sim.get_value(math.ceil(x0))[9])
        if num_time_points <= x0 < 2*num_time_points: # Number of Recovered is queried
                x0 = x0-num_time_points 
                result_arr[counter] = linearize(x0,math.floor(x0),math.ceil(x0),result_sim.get_value(math.floor(x0))[8], result_sim.get_value(math.ceil(x0))[8])
        counter +=1

    return result_arr



# Read in data from cases_all_county_ma7.json for County 1001.0 for the first hundred days: 24.04.2020-01.08.2020
def read_data():
        # First, extract the subset of data that we need
        df = pd.read_json('pycode/examples/simulation/cases_all_county_ma7.json',)
        df_subset = df[df['Date'] <= '2020-08-01']
        df_sub = df_subset[df_subset['ID_County'] == 1001.0]
        # Next, turn dataframe into two-dim numpy array
        array = np.empty((2,100))
        array[0] = df_sub['Recovered'] # Recovered
        array[1] = df_sub['Deaths'] # Deaths
        return np.concatenate([array[0],array[1]])

def fit_model():
    lmfit_params = lmfit.Parameters()

    # Initial populations as fixed parameters
    lmfit_params.add('S0',value=82760) # To be done
    lmfit_params.add('E0',value=100)
    lmfit_params.add('C0',value=50)
    lmfit_params.add('I0',value=50)
    lmfit_params.add('H0',value=20)
    lmfit_params.add('U0',value=10)
    lmfit_params.add('R0',value=10)
    lmfit_params.add('D0',value=0)
    lmfit_params['S0'].vary = False
    lmfit_params['E0'].vary = False
    lmfit_params['C0'].vary = False
    lmfit_params['I0'].vary = False
    lmfit_params['H0'].vary = False
    lmfit_params['U0'].vary = False
    lmfit_params['R0'].vary = False
    lmfit_params['D0'].vary = False

    # Variable parameters varied in fitting process
    lmfit_params.add('RelativeTransmissionNoSymptoms',value=0.67)
    lmfit_params.add('TransmissionProbabilityOnContact',value=0.8, min = 0.0, max = 1.0)
    lmfit_params.add('RecoveredPerInfectedNoSymptoms',value=0.09,min=0.01,max=0.16)
    lmfit_params.add('RiskOfInfectionFromSymptomatic',value=0.25,min=0.05,max=0.5)
    lmfit_params.add('SeverePerInfectedSymptoms',value=0.2)
    lmfit_params.add('CriticalPerSevere',value=0.25)
    lmfit_params.add('DeathsPerCritical',value=0.1)
    '''
    lmfit_params['RecoveredPerInfectedNoSymptoms'].vary = False
    lmfit_params['RiskOfInfectionFromSymptomatic'].vary = False
    lmfit_params['SeverePerInfectedSymptoms'].vary = False
    lmfit_params['CriticalPerSevere'].vary = False
    lmfit_params['DeathsPerCritical'].vary = False
    '''
    # Don't need MaxRiskOfInfectionFromSymptomatic, just take 2*RiskOfInfectionFromSymptomatic
    # lmfit_params.add('MaxRiskOfInfectionFromSymptomatic',value=0.5)

    total_data_array = read_data()

    secir_model = lmfit.Model(fit_func_secir)

    results = secir_model.fit(total_data_array,params=lmfit_params,x = np.arange(0,200),method='leastsq')

    best_params = results.best_values

    print("THE BEST VALUES ARE: ")
    print(fit_func_secir(np.arange(0.0,100.0),best_params['RelativeTransmissionNoSymptoms'],
                    best_params['TransmissionProbabilityOnContact'],
                    best_params['RecoveredPerInfectedNoSymptoms'],
                    best_params['RiskOfInfectionFromSymptomatic'],
                    best_params['SeverePerInfectedSymptoms'],
                    best_params['CriticalPerSevere'],
                    best_params['DeathsPerCritical'],
                    1000,500,0,0,0,0,0,0))
    
    print("THIS IS THE TOTAL DATA ARRAY: ")
    print(total_data_array)

if __name__ == '__main__':
    fit_model()