#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Martin J. Kuehn, Wadim Koslow
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################
from epidemiology.secir import (UncertainContactMatrix, ContactMatrix, Damping, SecirModel,
                                simulate, AgeGroup, Index_InfectionState, SecirSimulation)
from epidemiology.secir import InfectionState as State
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

def run_secir_simulation():
    """
    Runs the c++ secir model using one age group 
    and plots the results
    """

    # Define Comartment names
    compartments = ['Susceptible', 'Exposed', 'Carrier', 'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead']
    # Define population of age groups
    populations = [83000]

    days = 100 # number of days to simulate
    start_day = 1
    start_month = 1
    start_year = 2019
    dt = 0.1
    num_groups = 1
    num_compartments = len(compartments)

    # Initialize Parameters
    model = SecirModel(1)

    A0 = AgeGroup(0)

    # Set parameters

    # Compartment transition duration
    model.parameters.IncubationTime[A0] = 5.2  # R_2^(-1)+R_3^(-1)
    model.parameters.InfectiousTimeMild[A0] =  6.  # 4-14  (=R4^(-1))
    model.parameters.SerialInterval[A0] = 4.2   # 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
    model.parameters.HospitalizedToHomeTime[A0] = 12.  # 7-16 (=R5^(-1))
    model.parameters.HomeToHospitalizedTime[A0] = 5.  # 2.5-7 (=R6^(-1))
    model.parameters.HospitalizedToICUTime[A0] = 2.  # 1-3.5 (=R7^(-1))
    model.parameters.ICUToHomeTime[A0] = 8.  # 5-16 (=R8^(-1))
    model.parameters.ICUToDeathTime[A0] = 5.  # 3.5-7 (=R5^(-1))

    # Initial number of people in each compartment
    model.populations[A0, Index_InfectionState(State.Exposed)] = 100
    model.populations[A0, Index_InfectionState(State.Carrier)] = 50
    model.populations[A0, Index_InfectionState(State.Infected)] = 50
    model.populations[A0, Index_InfectionState(State.Hospitalized)] = 20
    model.populations[A0, Index_InfectionState(State.ICU)] = 10
    model.populations[A0, Index_InfectionState(State.Recovered)] = 10
    model.populations[A0, Index_InfectionState(State.Dead)] = 0
    model.populations.set_difference_from_total((A0, Index_InfectionState(State.Susceptible)), populations[0])
    
    # Compartment transition propabilities
    model.parameters.RelativeCarrierInfectability[A0] = 0.67 #TODO: must this parameter be set?
    model.parameters.InfectionProbabilityFromContact[A0] = 1.0
    model.parameters.AsymptoticCasesPerInfectious[A0] = 0.09  # 0.01-0.16
    model.parameters.RiskOfInfectionFromSympomatic[A0] = 0.25  # 0.05-0.5
    model.parameters.HospitalizedCasesPerInfectious[A0] = 0.2  # 0.1-0.35
    model.parameters.ICUCasesPerHospitalized[A0] = 0.25  # 0.15-0.4
    model.parameters.DeathsPerHospitalized[A0] = 0.3  # 0.15-0.77
    #TODO what about MaxRiskOfInfectionFromSympomatic?

    #TODO how to do this: model.parameters.set_start_day(start_day + start_month * 30) # TODO: start day has to adapted more precisely!

    #TODO use "contacts = ContactMatrix(np.r_[0.5])" or differ between baseline and minimum?
    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones((num_groups, num_groups)) * 1
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.ones((num_groups, num_groups)) * 0
    model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(coeffs = np.r_[0.9], t = 30.0, level = 0, type = 0))

    # Apply mathematical constraints to parameters
    model.apply_constraints()

    # Run Simulation
    result = simulate(0, days, dt, model)
    print(result.get_last_value())

    num_time_points = result.get_num_time_points()
    result_array = result.as_ndarray()
    t = result_array[0, :]
    group_data = np.transpose(result_array[1:, :])

    #sum over all groups
    data = np.zeros((num_time_points,num_compartments))
    for i in range(num_groups):
        data += group_data[:, i * num_compartments : (i + 1) * num_compartments]

    # Plot Results
    datelist = np.array(pd.date_range(datetime(start_year, start_month, start_day), periods=days, freq='D').strftime('%m-%d').tolist())

    tick_range = (np.arange(int(days / 10) + 1 ) * 10)
    tick_range[-1] -=1
    fig, ax = plt.subplots()
    ax.plot(t, data[:,0], label='#Susceptible')
    ax.plot(t, data[:,1], label='#Exposed')
    ax.plot(t, data[:,2], label='#Carrying')
    ax.plot(t, data[:,3], label='#Infected')
    ax.plot(t, data[:,4], label='#Hospitalzed')
    ax.plot(t, data[:,5], label='#ICU')
    ax.plot(t, data[:,6], label='#Recovered')
    ax.plot(t, data[:,7], label='#Died')
    ax.set_title("SECIR model simulation")
    ax.set_xticks(tick_range)
    ax.set_xticklabels(datelist[tick_range],rotation=45)
    ax.legend()
    fig.tight_layout
    fig.savefig('Secir_simple.pdf')

    plt.show()
    plt.close()

if __name__ == "__main__":
    run_secir_simulation()
