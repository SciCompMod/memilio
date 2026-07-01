#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Martin J. Kuehn, Maximilian Betz
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

from memilio.simulation import AgeGroup, Damping, SimulationDay
from memilio.simulation.osecirvvs import InfectionState, interpolate_simulation_result
from memilio.simulation.osecirvvs import Model, simulate


def run_ode_secirvvs_simulation(show_plot=True):
    """Runs the c++ ODE SECIRVVS model using a single age group
    and plots the results

    :param show_plot:  (Default value = True)

    """

    t0 = 0
    tmax = 30  # number of days to simulate
    dt = 0.1
    num_groups = 1

    # Initialize Parameters
    model = Model(num_groups)

    # set parameters
    for i in range(num_groups):
        # Initial number of peaople in each compartment
        model.populations[AgeGroup(i), InfectionState.ExposedNaive] = 10
        model.populations[AgeGroup(
            i), InfectionState.ExposedImprovedImmunity] = 11
        model.populations[AgeGroup(
            i), InfectionState.ExposedPartialImmunity] = 12
        model.populations[AgeGroup(
            i), InfectionState.InfectedNoSymptomsNaive] = 13
        model.populations[AgeGroup(
            i), InfectionState.InfectedNoSymptomsNaiveConfirmed] = 13
        model.populations[AgeGroup(
            i), InfectionState.InfectedNoSymptomsPartialImmunity] = 14
        model.populations[AgeGroup(
            i), InfectionState.InfectedNoSymptomsPartialImmunityConfirmed] = 14
        model.populations[AgeGroup(
            i), InfectionState.InfectedNoSymptomsImprovedImmunity] = 15
        model.populations[AgeGroup(
            i), InfectionState.InfectedNoSymptomsImprovedImmunityConfirmed] = 15
        model.populations[AgeGroup(
            i), InfectionState.InfectedSymptomsNaive] = 5
        model.populations[AgeGroup(
            i), InfectionState.InfectedSymptomsNaiveConfirmed] = 5
        model.populations[AgeGroup(
            i), InfectionState.InfectedSymptomsPartialImmunity] = 6
        model.populations[AgeGroup(
            i), InfectionState.InfectedSymptomsPartialImmunityConfirmed] = 6
        model.populations[AgeGroup(
            i), InfectionState.InfectedSymptomsImprovedImmunity] = 7
        model.populations[AgeGroup(
            i), InfectionState.InfectedSymptomsImprovedImmunityConfirmed] = 7
        model.populations[AgeGroup(i), InfectionState.InfectedSevereNaive] = 8
        model.populations[AgeGroup(
            i), InfectionState.InfectedSevereImprovedImmunity] = 1
        model.populations[AgeGroup(
            i), InfectionState.InfectedSeverePartialImmunity] = 2
        model.populations[AgeGroup(
            i), InfectionState.InfectedCriticalNaive] = 3
        model.populations[AgeGroup(
            i), InfectionState.InfectedCriticalPartialImmunity] = 4
        model.populations[AgeGroup(
            i), InfectionState.InfectedCriticalImprovedImmunity] = 5
        model.populations[AgeGroup(
            i), InfectionState.SusceptibleImprovedImmunity] = 6
        model.populations[AgeGroup(
            i), InfectionState.SusceptiblePartialImmunity] = 7
        model.populations[AgeGroup(i), InfectionState.DeadNaive] = 0
        model.populations[AgeGroup(i), InfectionState.DeadPartialImmunity] = 0
        model.populations[AgeGroup(i), InfectionState.DeadImprovedImmunity] = 0
        model.populations.set_difference_from_group_total_AgeGroup(
            (AgeGroup(i), InfectionState.SusceptibleNaive), 1000)

    model.parameters.ICUCapacity.value = 100
    model.parameters.TestAndTraceCapacity.value = 0.0143
    model.parameters.DailyPartialVaccinations.resize_SimulationDay(
        SimulationDay(tmax + 1))
    model.parameters.DailyFullVaccinations.resize_SimulationDay(
        SimulationDay(tmax + 1))
    daily_vaccinations = 10
    for i, num_vaccinations in enumerate(range(0, daily_vaccinations * (tmax + 1), daily_vaccinations)):
        model.parameters.DailyPartialVaccinations[AgeGroup(
            0), SimulationDay(i)] = num_vaccinations
        model.parameters.DailyFullVaccinations[AgeGroup(
            0), SimulationDay(i)] = num_vaccinations

    # contact patterns
    baseline = np.ones((num_groups, num_groups)) * 0.5
    np.fill_diagonal(baseline, 5.0)
    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = baseline
    model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=np.ones((num_groups, num_groups)) * 0.3, t=5.0, level=0, type=0))

    # times
    model.parameters.TimeInfectedSymptoms[AgeGroup(0)] = 7
    model.parameters.TimeInfectedSevere[AgeGroup(0)] = 6
    model.parameters.TimeInfectedCritical[AgeGroup(0)] = 7

    # probabilities
    model.parameters.TransmissionProbabilityOnContact[AgeGroup(0)] = 0.15
    model.parameters.RelativeTransmissionNoSymptoms[AgeGroup(0)] = 0.5
    # The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
    # depends on incidence and test and trace capacity
    model.parameters.RiskOfInfectionFromSymptomatic[AgeGroup(0)] = 0.0
    model.parameters.MaxRiskOfInfectionFromSymptomatic[AgeGroup(0)] = 0.4
    model.parameters.RecoveredPerInfectedNoSymptoms[AgeGroup(0)] = 0.2
    model.parameters.SeverePerInfectedSymptoms[AgeGroup(0)] = 0.1
    model.parameters.CriticalPerSevere[AgeGroup(0)] = 0.1
    model.parameters.DeathsPerCritical[AgeGroup(0)] = 0.1

    model.parameters.ReducExposedPartialImmunity[AgeGroup(0)] = 0.8
    model.parameters.ReducExposedImprovedImmunity[AgeGroup(0)] = 0.331
    model.parameters.ReducInfectedSymptomsPartialImmunity[AgeGroup(0)] = 0.65
    model.parameters.ReducInfectedSymptomsImprovedImmunity[AgeGroup(0)] = 0.243
    model.parameters.ReducInfectedSevereCriticalDeadPartialImmunity[AgeGroup(
        0)] = 0.1
    model.parameters.ReducInfectedSevereCriticalDeadImprovedImmunity[AgeGroup(
        0)] = 0.091
    model.parameters.ReducTimeInfectedMild[AgeGroup(0)] = 0.9

    model.parameters.Seasonality.value = 0.2

    model.apply_constraints()

    # Run Simulation
    result = simulate(t0, tmax, dt, model)

    # # interpolate results
    result = interpolate_simulation_result(result)

    print(result.get_last_value())


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'ode_secirvvs_simple',
        description='Simple example demonstrating the setup and simulation of the ODE SECIRVVS model with a single age group.')
    args = arg_parser.parse_args()
    run_ode_secirvvs_simulation(**args.__dict__)
