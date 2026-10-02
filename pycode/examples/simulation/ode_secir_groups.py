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

# This example requires memilio.simulation and memilio.plot to be installed!

import argparse
import os
from datetime import date

import matplotlib.pyplot as plt
import numpy as np

from memilio.plot.plotTimeSeries import plot_time_series
from memilio.simulation import AgeGroup, Damping
from memilio.simulation.osecir import InfectionState as State
from memilio.simulation.osecir import (Model, interpolate_simulation_result,
                                       simulate)


def run_ode_secir_groups_simulation(show_plot=True):
    """Runs the c++ ODE SECIHURD model using mulitple age groups
    and plots the results

    :param show_plot: Whether to show the figures interactively.
        (Default value = True)

    """

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

    # set contact frequency matrix
    data_dir = os.path.join(os.path.dirname(
        __file__), "..", "..", "..", "data", "Germany")
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

    start_date = date(start_year, start_month, start_day)

    # Plot results summed over all age groups, one line per compartment. The
    # compartment names are taken from the InfectionState enum of the model.
    ax = plot_time_series(
        result, labels=State.values(), groups=groups, start_date=start_date,
        title='ODE SECIR simulation results (entire population)')
    ax.figure.savefig('osecir_by_compartments.pdf')

    # One panel per compartment with one line per age group.
    fig, axes = plt.subplots(5, 2, figsize=(16, 18), layout='constrained')
    for state, ax in zip(State.values(), axes.flat):
        plot_time_series(
            result, labels=State.values(), groups=groups, sum_groups=False,
            select=state, start_date=start_date, ax=ax, title=state.name)
    fig.suptitle(
        'ODE SECIR simulation results by age group in each compartment')
    fig.savefig('osecir_age_groups_in_compartments.pdf')

    # One panel per compartment, summed over all age groups.
    fig, axes = plt.subplots(5, 2, figsize=(16, 18), layout='constrained')
    for state, ax in zip(State.values(), axes.flat):
        plot_time_series(
            result, labels=State.values(), groups=groups, select=state,
            start_date=start_date, ax=ax, title=state.name)
    fig.suptitle(
        'ODE SECIR simulation results by compartment (entire population)')
    fig.savefig('osecir_all_parts.pdf')

    if show_plot:
        plt.show()
    plt.close('all')


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'ode_secir_groups',
        description='Simple example demonstrating the setup and simulation of the ODE SECIHURD model with multiple age groups.')
    arg_parser.add_argument('-p', '--show_plot',
                            action='store_const', const=True, default=False)
    args = arg_parser.parse_args()
    run_ode_secir_groups_simulation(**args.__dict__)
