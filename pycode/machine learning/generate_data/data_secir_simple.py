from memilio.simulation import UncertainContactMatrix, ContactMatrix, Damping
from memilio.simulation.secir import SecirModel, simulate, AgeGroup, Index_InfectionState, SecirSimulation
from memilio.simulation.secir import InfectionState as State
import numpy as np
from datetime import date
from math import ceil
import random


def run_secir_simulation():

    # Define Comartment names
    compartments = ['Susceptible', 'Exposed', 'Carrier',
                    'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead']

    # Define random population number in intervall (0,500000)
    populations = [500000 * random.random()]

    days = 100  # number of days to simulate
    # since we defined the number of days, we already know the dimension of our dataset.
    # We aim to save the compartment population for each day.
    # i.e. for days = 100 (+1 for t0) and len(compartments = 8) -> 101 x  8 array
    data = np.zeros((days+1,len(compartments)))
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
    data[0,:] = model.populations.get_compartments()
    # TODO: Do we need to build the whole Model for new, if we only change the number of days?
    for day in range(1, days+1):
        print("Simulate with day =" , day)
        result = simulate(0, day, dt, model)

        # Maybe round data
        # data[day,:] = np.ceil(result_day)
        data[day,:] = result.get_last_value()

    return data

def getDataAndLabels(data):
    rows = np.shape(data)[0]
    cols = np.shape(data)[1]
    X = data[:100,:]
    Y = data[1:,:]
    return X,Y



if __name__ == "__main__":
    # TODO: Save contact matrix depending on the damping.
    # In the actual state it might be enough to save the normal one and the damping-

    # data_ = 
    data = run_secir_simulation()
    X,Y = getDataAndLabels(run_secir_simulation())
