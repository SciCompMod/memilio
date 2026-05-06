#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Martin J. Kuehn, Wadim Koslow, Annalena Lange, Khoa Nguyen
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
import argparse
import os
from datetime import date, datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from memilio.simulation import AgeGroup, ContactMatrix, Damping, UncertainContactMatrix
from memilio.simulation.osecir import Index_InfectionState
from memilio.simulation.osecir import InfectionState as State
from memilio.simulation.osecir import (Model, Simulation,
                                       interpolate_simulation_result, simulate)


def run_ode_secir_groups_simulation(show_plot=True):
    """Runs the c++ ODE SECIHURD model using mulitple age groups
    and plots the results

    :param show_plot:  (Default value = True)

    """

    # Define Comartment names
    compartments = [
        'Susceptible', 'Exposed', 'InfectedNoSymptoms', 'InfectedSymptoms',
        'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']
    # Define age Groups
    groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80+']
    # Define population of age groups
    populations = [40000, 70000, 190000, 290000, 180000, 60000]

    days = 100  # number of days to simulate
    start_day = 1
    start_month = 1
    start_year = 2019
    dt = 0.1
    num_groups = len(groups)
    num_compartments = len(compartments)

    # set contact frequency matrix
    data_dir = os.path.join(os.path.dirname(
        __file__), "..", "..", "..", "data")
    baseline_contact_matrix0 = os.path.join(
        data_dir, "contacts/baseline_home.txt")
    baseline_contact_matrix1 = os.path.join(
        data_dir, "contacts/baseline_school_pf_eig.txt")
    baseline_contact_matrix2 = os.path.join(
        data_dir, "contacts/baseline_work.txt")
    baseline_contact_matrix3 = os.path.join(
        data_dir, "contacts/baseline_other.txt")

    # Initialize Parameters
    model = Model(len(populations))

    # set parameters
    for i in range(num_groups):
        # Compartment transition duration
        model.parameters.TimeExposed[AgeGroup(i)] = 3.2
        model.parameters.TimeInfectedNoSymptoms[AgeGroup(i)] = 2.
        model.parameters.TimeInfectedSymptoms[AgeGroup(i)] = 6.
        model.parameters.TimeInfectedSevere[AgeGroup(i)] = 12.
        model.parameters.TimeInfectedCritical[AgeGroup(i)] = 8.

        # Initial number of peaople in each compartment
        model.populations[AgeGroup(i), State.Exposed] = 100
        model.populations[AgeGroup(i), State.InfectedNoSymptoms] = 50
        model.populations[AgeGroup(i), State.InfectedNoSymptomsConfirmed] = 0
        model.populations[AgeGroup(i), State.InfectedSymptoms] = 50
        model.populations[AgeGroup(i), State.InfectedSymptomsConfirmed] = 0
        model.populations[AgeGroup(i), State.InfectedSevere] = 20
        model.populations[AgeGroup(i), State.InfectedCritical] = 10
        model.populations[AgeGroup(i), State.Recovered] = 10
        model.populations[AgeGroup(i), State.Dead] = 0
        model.populations.set_difference_from_group_total_AgeGroup(
            (AgeGroup(i), State.Susceptible), populations[i])

        # Compartment transition propabilities

        model.parameters.RelativeTransmissionNoSymptoms[AgeGroup(i)] = 0.67
        model.parameters.TransmissionProbabilityOnContact[AgeGroup(i)] = 1.0
        model.parameters.RecoveredPerInfectedNoSymptoms[AgeGroup(
            i)] = 0.09  # 0.01-0.16
        model.parameters.RiskOfInfectionFromSymptomatic[AgeGroup(
            i)] = 0.25  # 0.05-0.5
        model.parameters.SeverePerInfectedSymptoms[AgeGroup(
            i)] = 0.2  # 0.1-0.35
        model.parameters.CriticalPerSevere[AgeGroup(
            i)] = 0.25  # 0.15-0.4
        model.parameters.DeathsPerCritical[AgeGroup(i)] = 0.3  # 0.15-0.77
        # twice the value of RiskOfInfectionFromSymptomatic
        model.parameters.MaxRiskOfInfectionFromSymptomatic[AgeGroup(i)] = 0.5

    model.parameters.StartDay = (
        date(start_year, start_month, start_day) - date(start_year, 1, 1)).days

    # set contact rates and emulate some mitigations
    # set contact frequency matrix
    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.loadtxt(baseline_contact_matrix0) \
        + np.loadtxt(baseline_contact_matrix1) + \
        np.loadtxt(baseline_contact_matrix2) + \
        np.loadtxt(baseline_contact_matrix3)
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.ones(
        (num_groups, num_groups)) * 0
    model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=np.ones((num_groups, num_groups)) * 0.9, t=30.0, level=0, type=0))

    # Apply mathematical constraints to parameters
    model.apply_constraints()

    # Run Simulation
    result = simulate(0, days, dt, model)

    # interpolate results
    result = interpolate_simulation_result(result)

    # print(result.get_last_value())
    num_time_points = result.get_num_time_points()
    result_array = result.as_ndarray()
    t = result_array[0, :]
    group_data = np.transpose(result_array[1:, :])

    # sum over all groups
    data = np.zeros((num_time_points, num_compartments))
    for i in range(num_groups):
        data += group_data[:, i * num_compartments: (i + 1) * num_compartments]

    # Plot Results
    datelist = np.array(
        pd.date_range(
            datetime(start_year, start_month, start_day),
            periods=days, freq='D').strftime('%m-%d').tolist())

    tick_range = (np.arange(int(days / 10) + 1) * 10)
    tick_range[-1] -= 1
    fig, ax = plt.subplots()
    ax.plot(t, data[:, 0], label='#Susceptible')
    ax.plot(t, data[:, 1], label='#Exposed')
    ax.plot(t, data[:, 2], label='#InfectedNoSymptoms')
    ax.plot(t, data[:, 3], label='#InfectedSymptoms')
    ax.plot(t, data[:, 4], label='#Hospitalzed')
    ax.plot(t, data[:, 5], label='#InfectedCritical')
    ax.plot(t, data[:, 6], label='#Recovered')
    ax.plot(t, data[:, 7], label='#Dead')
    ax.set_title("ODE SECIR simulation results (entire population)")
    ax.set_xticks(tick_range)
    ax.set_xticklabels(datelist[tick_range], rotation=45)
    ax.legend()
    fig.tight_layout
    fig.savefig('osecir_by_compartments.pdf')

    # plot dynamics in each comparment by age group
    fig, ax = plt.subplots(4, 2, figsize=(12, 15))

    for i, title in zip(range(num_compartments), compartments):

        for j, group in enumerate(groups):
            ax[int(np.floor(i / 2)), int(i % 2)].plot(t,
                                                      group_data[:, j*num_compartments+i], label=group)

        ax[int(np.floor(i / 2)), int(i % 2)].set_title(title, fontsize=10)
        ax[int(np.floor(i / 2)), int(i % 2)].legend()

        ax[int(np.floor(i / 2)), int(i % 2)].set_xticks(tick_range)
        ax[int(np.floor(i / 2)), int(i % 2)
           ].set_xticklabels(datelist[tick_range], rotation=45)
    plt.subplots_adjust(hspace=0.5, bottom=0.1, top=0.9)
    fig.suptitle(
        'ODE SECIR simulation results by age group in each compartment')
    fig.savefig('osecir_age_groups_in_compartments.pdf')

    fig, ax = plt.subplots(4, 2, figsize=(12, 15))
    for i, title in zip(range(num_compartments), compartments):
        ax[int(np.floor(i / 2)), int(i % 2)].plot(t, data[:, i])
        ax[int(np.floor(i / 2)), int(i % 2)].set_title(title, fontsize=10)

        ax[int(np.floor(i / 2)), int(i % 2)].set_xticks(tick_range)
        ax[int(np.floor(i / 2)), int(i % 2)
           ].set_xticklabels(datelist[tick_range], rotation=45)
    plt.subplots_adjust(hspace=0.5, bottom=0.1, top=0.9)
    fig.suptitle(
        'ODE SECIR simulation results by compartment (entire population)')
    fig.savefig('osecir_all_parts.pdf')

    if show_plot:
        plt.show()
        plt.close()

    # return data


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'ode_secir_groups',
        description='Simple example demonstrating the setup and simulation of the ODE SECIHURD model with multiple age groups.')
    arg_parser.add_argument('-p', '--show_plot',
                            action='store_const', const=True, default=False)
    args = arg_parser.parse_args()
    run_ode_secir_groups_simulation(**args.__dict__)
