#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Maximilian Betz
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

from memilio.simulation import ContactMatrix, Damping, UncertainContactMatrix, AgeGroup
from memilio.simulation.osecir import Index_InfectionState
from memilio.simulation.osecir import InfectionState as State
from memilio.simulation.osecir import Model, Simulation, simulate


def run_ode_secir_ageres_simulation(show_plot=True):
    """Runs the c++ ODE SECIHURD model using multiple age groups
    and plots the results

    :param show_plot:  (Default value = True)

    """

    model = Model(3)

    t0 = 0
    tmax = 50
    dt = 0.1

    cont_freq = 10
    nb_total_t0, nb_exp_t0, nb_inf_t0, nb_car_t0, nb_hosp_t0, nb_icu_t0, nb_rec_t0, nb_dead_t0 = 10000, 100, 50, 50, 20, 10, 10, 0

    nb_comp = 10
    nb_groups = 3
    fact = 1.0/nb_groups

    model.parameters.StartDay = 60
    model.parameters.Seasonality.value = 0.2
    model.parameters.TestAndTraceCapacity.value = 35

    for i in range(0, 3):
        Ai = AgeGroup(i)

        model.parameters.TimeExposed[Ai] = 3.2
        model.parameters.TimeInfectedNoSymptoms[Ai] = 2.0
        model.parameters.TimeInfectedSymptoms[Ai] = 5.8
        model.parameters.TimeInfectedSevere[Ai] = 9.5
        model.parameters.TimeInfectedCritical[Ai] = 7.1

        model.populations[Ai, State.Exposed] = fact * nb_exp_t0
        model.populations[Ai, State.InfectedNoSymptoms] = fact * nb_car_t0
        model.populations[Ai, State.InfectedSymptoms] = fact * nb_inf_t0
        model.populations[Ai, State.InfectedSevere] = fact * nb_hosp_t0
        model.populations[Ai, State.InfectedCritical] = fact * nb_icu_t0
        model.populations[Ai, State.Recovered] = fact * nb_rec_t0
        model.populations[Ai, State.Dead] = fact * nb_dead_t0
        model.populations.set_difference_from_group_total_AgeGroup(
            (Ai, State.Susceptible), fact * nb_total_t0)

        model.parameters.TransmissionProbabilityOnContact[Ai] = 0.05
        model.parameters.RelativeTransmissionNoSymptoms[Ai] = 0.7
        model.parameters.RecoveredPerInfectedNoSymptoms[Ai] = 0.09
        model.parameters.RiskOfInfectionFromSymptomatic[Ai] = 0.25
        model.parameters.MaxRiskOfInfectionFromSymptomatic[Ai] = 0.45
        model.parameters.SeverePerInfectedSymptoms[Ai] = 0.2
        model.parameters.CriticalPerSevere[Ai] = 0.25
        model.parameters.DeathsPerCritical[Ai] = 0.3

    model.apply_constraints()

    contacts = ContactMatrix(
        np.full((nb_groups, nb_groups), fact * cont_freq))
    contacts.add_damping(
        Damping(coeffs=np.full((nb_groups, nb_groups), 0.7), t=30.0, level=0, type=0))
    model.parameters.ContactPatterns.cont_freq_mat[0] = contacts

    model = model

    result = simulate(t0=t0, tmax=tmax,
                      dt=dt, model=model)

    print(result.get_last_value())

    num_time_points = result.get_num_time_points()
    result_array = result.as_ndarray()
    t = result_array[0, :]
    group_data = np.transpose(result_array[1:, :])

    # sum over all groups
    data = np.zeros((num_time_points, nb_comp))
    for i in range(nb_groups):
        data += group_data[:, i * nb_comp: (i + 1) * nb_comp]

    start_day = 1
    start_month = 1
    start_year = 2019

    if (show_plot):
        # Plot Results
        datelist = np.array(
            pd.date_range(
                datetime(start_year, start_month, start_day),
                periods=tmax, freq='D').strftime('%m-%d').tolist())

        tick_range = (np.arange(int(tmax / 10) + 1) * 10)
        tick_range[-1] -= 1
        fig, ax = plt.subplots()
        ax.plot(t, data[:, 0], label='#Susceptible')
        ax.plot(t, data[:, 1], label='#Exposed')
        ax.plot(t, data[:, 2], label='#Carrying')
        ax.plot(t, data[:, 3], label='#InfectedSymptoms')
        ax.plot(t, data[:, 4], label='#Hospitalzed')
        ax.plot(t, data[:, 5], label='#InfectedCritical')
        ax.plot(t, data[:, 6], label='#Recovered')
        ax.plot(t, data[:, 7], label='#Died')
        ax.set_title("ODE SECIR model simulation")
        ax.set_xticks(tick_range)
        ax.set_xticklabels(datelist[tick_range], rotation=45)
        ax.legend()
        fig.tight_layout
        fig.savefig('osecir_simple.pdf')

        plt.show()
        plt.close()

    for timestep in range(num_time_points):
        print(timestep, ": ", result.get_time(timestep), "\n")


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'ode_secir_ageres',
        description='Example demonstrating the setup and simulation of the ODE SECIHURD model with multiple age groups.')
    arg_parser.add_argument('-p', '--show_plot',
                            action='store_const', const=True, default=False)
    args = arg_parser.parse_args()
    run_ode_secir_ageres_simulation(**args.__dict__)
