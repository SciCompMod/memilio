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
import epidemiology.secir as secir
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime


def run_secir_groups_simulation():
    """
    Runs the c++ secir model using mulitple age groups 
    and plots the results
    """

    # Define Comartment names
    compartments = ['Susceptible', 'Exposed', 'Carrier', 'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead']
    # Define age Groups
    groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80+']
    # Define population of age groups
    populations = [30000, 60000, 150000, 200000, 120000, 100000] #secir simple uses 820000 as pop size, this adds up to 660000

    days = 100 # number of days to simulate
    start_day = 1
    start_month = 1
    start_year = 2019
    dt = 0.1
    num_groups = len(groups)
    num_compartments = len(compartments)

    # set contact frequency matrix
    baseline_contact_matrix = "../../data/contacts/baseline_school_pf_eig.txt"

    # Initialize Parameters
    model = secir.SecirModel(len(populations))

    # set parameters
    for i in range(num_groups):
        # Compartment transition duration
        model.parameters.times[i].set_incubation(5.2)
        model.parameters.times[i].set_infectious_mild(6)
        model.parameters.times[i].set_serialinterval(4.2)
        model.parameters.times[i].set_hospitalized_to_home(12)
        model.parameters.times[i].set_home_to_hospitalized(5)
        model.parameters.times[i].set_hospitalized_to_icu(2)
        model.parameters.times[i].set_icu_to_home(8)
        model.parameters.times[i].set_icu_to_death(5)

        # Initial number of peaople in each compartment
        model.populations[secir.AgeGroup(i), secir.Index_InfectionState(secir.InfectionState.Exposed)] = 100
        model.populations[secir.AgeGroup(i), secir.Index_InfectionState(secir.InfectionState.Carrier)] = 40
        model.populations[secir.AgeGroup(i), secir.Index_InfectionState(secir.InfectionState.Infected)] = 80
        model.populations[secir.AgeGroup(i), secir.Index_InfectionState(secir.InfectionState.Hospitalized)] = 40
        model.populations[secir.AgeGroup(i), secir.Index_InfectionState(secir.InfectionState.ICU)] = 20
        model.populations[secir.AgeGroup(i), secir.Index_InfectionState(secir.InfectionState.Recovered)] = 7
        model.populations[secir.AgeGroup(i), secir.Index_InfectionState(secir.InfectionState.Dead)] = 3
        model.populations.set_difference_from_group_total_AgeGroup((secir.AgeGroup(i), 
            secir.Index_InfectionState(secir.InfectionState.Susceptible)), populations[i])

        # Compartment transition propabilities
        model.parameters.probabilities[0].set_infection_from_contact(1.0)
        model.parameters.probabilities[0].set_carrier_infectability(0.67)
        model.parameters.probabilities[0].set_asymp_per_infectious(0.09)
        model.parameters.probabilities[0].set_risk_from_symptomatic(0.25)
        model.parameters.probabilities[0].set_hospitalized_per_infectious(0.2)
        model.parameters.probabilities[0].set_icu_per_hospitalized(0.25)
        model.parameters.probabilities[0].set_dead_per_icu(0.3)
        
    
    # set contact rates and emulate some mitigations
    # set contact frequency matrix
    model.parameters.get_contact_patterns().cont_freq_mat[0].baseline = np.loadtxt(baseline_contact_matrix)
    model.parameters.get_contact_patterns().cont_freq_mat[0].minimum = np.ones((num_groups, num_groups)) * 0

    # Define Damping on Contacts
    model.parameters.get_contact_patterns().cont_freq_mat.add_damping(secir.Damping(np.ones((num_groups, num_groups)) * 0.9, 30, 0, 0))

     # Apply mathematical constraints to parameters
    model.apply_constraints()

    # Run Simulation
    result = secir.simulate(0, days, dt, model)
    #print(result.get_last_value())

    num_time_points = result.get_num_time_points()
    result_array = result.as_ndarray()
    t = result_array[0, :]
    group_data = np.transpose(result_array[1:, :])

    # sum over all groups
    data = np.zeros((num_time_points, num_compartments))
    for i in range(num_groups):
        data += group_data[:, i * num_compartments : (i + 1) * num_compartments]

    # Plot Results
    datelist = np.array(pd.date_range(datetime(start_year, start_month, start_day), periods=days, freq='D').strftime('%m-%d').tolist())

    tick_range = (np.arange(int(days / 10) + 1 ) * 10)
    tick_range[-1] -=1
    fig, ax = plt.subplots()
    ax.plot(t, data[:,0], label='#Susceptible')
    ax.plot(t, data[:,1], label='#Exposed')
    ax.plot(t, data[:,2], label='#Carrier')
    ax.plot(t, data[:,3], label='#Infected')
    ax.plot(t, data[:,4], label='#Hospitalzed')
    ax.plot(t, data[:,5], label='#ICU')
    ax.plot(t, data[:,6], label='#Recovered')
    ax.plot(t, data[:,7], label='#Dead')
    ax.set_title("SECIR simulation results (entire population)")
    ax.set_xticks(tick_range)
    ax.set_xticklabels(datelist[tick_range],rotation=45)
    ax.legend()
    fig.tight_layout
    fig.savefig('Secir_by_compartments.pdf')

    
    # plot dynamics in each comparment by age group
    fig, ax = plt.subplots(4, 2, figsize=(12, 15))

    for i, title in zip(range(num_compartments), compartments):
        
        for j, group in enumerate(groups):
            ax[int(np.floor(i / 2)), int(i % 2)].plot(t, group_data[:, j*num_compartments+i], label = group)

        ax[int(np.floor(i / 2)), int(i % 2)].set_title(title, fontsize=10)
        ax[int(np.floor(i / 2)), int(i % 2)].legend()

        ax[int(np.floor(i / 2)), int(i % 2)].set_xticks(tick_range)
        ax[int(np.floor(i / 2)), int(i % 2)].set_xticklabels(datelist[tick_range], rotation=45)
    plt.subplots_adjust(hspace=0.5, bottom=0.1, top=0.9)
    fig.suptitle('SECIR simulation results by age group in each compartment')
    fig.savefig('Secir_age_groups_in_compartments.pdf')


    fig, ax = plt.subplots(4, 2, figsize=(12, 15))
    for i, title in zip(range(num_compartments), compartments):
        ax[int(np.floor(i / 2)), int(i % 2)].plot(t, data[:, i])
        ax[int(np.floor(i / 2)), int(i % 2)].set_title(title, fontsize=10)

        ax[int(np.floor(i / 2)), int(i % 2)].set_xticks(tick_range)
        ax[int(np.floor(i / 2)), int(i % 2)].set_xticklabels(datelist[tick_range], rotation=45)
    plt.subplots_adjust(hspace=0.5, bottom=0.1, top=0.9)
    fig.suptitle('SECIR simulation results by compartment (entire population)')
    fig.savefig('Secir_all_parts.pdf')

    plt.show()
    plt.close()

    # return data


if __name__ == "__main__":
    run_secir_groups_simulation()
