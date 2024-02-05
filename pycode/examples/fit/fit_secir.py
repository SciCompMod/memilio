import argparse
from datetime import date, datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from memilio.simulation import AgeGroup, ContactMatrix, Damping, UncertainContactMatrix
from memilio.simulation.secir import Index_InfectionState
from memilio.simulation.secir import InfectionState as State
from memilio.simulation.secir import (Model, Simulation,
                                      interpolate_simulation_result, simulate)
import lmfit

# First read in the existing data for one county and try to fit that

def fit_func_secir(params, initial_populations):
    # Define Comartment names
    compartments = [
        'Susceptible', 'Exposed', 'InfectedNoSymptoms', 'InfectedSymptoms',
        'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']
    # Define population of age groups
    populations = [83000]

    days = 100  # number of days to simulate
    start_day = 1
    start_month = 1
    start_year = 2019
    dt = 0.1
    num_groups = 1
    num_compartments = len(compartments)

    # Initialize Parameters
    model = Model(1)
    model.parameters = params

    # Initial number of people in each compartment
    model.populations[A0, State.Exposed] = 100
    model.populations[A0, State.InfectedNoSymptoms] = 50
    model.populations[A0, State.InfectedNoSymptomsConfirmed] = 0
    model.populations[A0, State.InfectedSymptoms] = 50
    model.populations[A0, State.InfectedSymptomsConfirmed] = 0
    model.populations[A0, State.InfectedSevere] = 20
    model.populations[A0, State.InfectedCritical] = 10
    model.populations[A0, State.Recovered] = 10
    model.populations[A0, State.Dead] = 0
    model.populations.set_difference_from_total(
        (A0, State.Susceptible), populations[0])

    # Apply mathematical constraints to parameters
    model.apply_constraints()

    # Run Simulation
    result = simulate(0, days, dt, model)
    # interpolate results
    result = interpolate_simulation_result(result)
    return result

def fit_model():
    lmfit_params = lmfit.Parameters()

    # Initial populations as fixed parameters
    lmfit_params.add('S0',value=...) # To be done
    lmfit_params.add('E0',value=...)
    lmfit_params.add('C0',value=...)
    lmfit_params.add('I0',value=...)
    lmfit_params.add('H0',value=...)
    lmfit_params.add('U0',value=...)
    lmfit_params.add('R0',value=...)
    lmfit_params.add('D0',value=...)

    # Variable parameters varied in fitting process
    lmfit_params.add('RelativeTransmissionNoSymptoms',value=0.67)
    lmfit_params.add('TransmissionProbabilityOnContact',value=1.0)
    lmfit_params.add('RecoveredPerInfectedNoSymptoms',value=0.09,min=0.01,max=0.16)
    lmfit_params.add('RiskOfInfectionFromSymptomatic',value=0.25,min=0.05,max=0.5)
    lmfit_params.add('SeverePerInfectedSymptoms',value=0.2)
    lmfit_params.add('CriticalPerSevere',value=0.25)
    lmfit_params.add('DeathsPerCritical',value=0.3)
    lmfit_params.add('MaxRiskOfInfectionFromSymptomatic',value=0.5)