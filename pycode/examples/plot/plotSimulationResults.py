#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Kilian Volmer
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
"""
Example demonstrating the standard TimeSeries plot of the MEmilio plot
package on the result of an ODE SEIR simulation with two age groups.
"""
import argparse
import os
from datetime import date

import matplotlib.pyplot as plt
import numpy as np

from memilio.plot.plotTimeSeries import (plot_time_series,
                                         time_series_to_dataframe)
from memilio.simulation import AgeGroup, Damping
from memilio.simulation.oseir import InfectionState, Model, simulate

AGE_GROUPS = ['0-19', '20+']
START_DATE = date(2020, 3, 1)


def run_ode_seir_simulation(days=100, dt=0.1):
    """ Runs the ODE SEIR model with two age groups.

    :param days: Number of days to simulate. (Default value = 100)
    :param dt: Initial time step. (Default value = 0.1)
    :returns: Simulation result as TimeSeries.
    """
    group_populations = [15000, 68000]
    num_groups = len(group_populations)
    model = Model(num_groups)

    for i, population in enumerate(group_populations):
        group = AgeGroup(i)
        model.parameters.TimeExposed[group] = 5.2
        model.parameters.TimeInfected[group] = 6.
        model.parameters.TransmissionProbabilityOnContact[group] = 1. * i
        model.populations[group, InfectionState.Exposed] = 100
        model.populations[group, InfectionState.Infected] = 50
        model.populations[group, InfectionState.Recovered] = 10
        model.populations.set_difference_from_group_total_AgeGroup(
            (group, InfectionState.Susceptible), population)

    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
        (num_groups, num_groups))
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.zeros(
        (num_groups, num_groups))
    model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=np.ones((num_groups, num_groups)) * 0.9, t=30.0, level=0,
        type=0))

    model.check_constraints()
    return simulate(0, days, dt, model)


def plot_results(result, output_path='.', show_plot=False):
    """ Creates the standard plots of the simulation result.

    :param result: TimeSeries returned by the simulation.
    :param output_path: Directory the figures are written to.
        (Default value = '.')
    :param show_plot: Whether to show the figures interactively.
        (Default value = False)
    """
    # All compartments, summed over the age groups. The compartment names
    # are taken from the InfectionState enum of the model.
    ax = plot_time_series(
        result, labels=InfectionState.values(), groups=AGE_GROUPS,
        title='ODE SEIR simulation')
    ax.figure.savefig(
        os.path.join(output_path, 'seir_compartments.png'), dpi=150)

    # Two panels in one figure: a selection of compartments per age group
    # on a date axis, and a selection on a logarithmic axis.
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout='constrained')
    plot_time_series(
        result, labels=InfectionState.values(), groups=AGE_GROUPS,
        sum_groups=False, select=['Exposed', 'Infected'],
        start_date=START_DATE, ax=axes[0],
        title='Exposed and infected per age group')
    plot_time_series(
        result, labels=InfectionState.values(), groups=AGE_GROUPS,
        select=['Exposed', 'Infected', 'Recovered'], log_scale=True,
        start_date=START_DATE, ax=axes[1], title='Logarithmic scale')
    fig.savefig(os.path.join(output_path, 'seir_selection.pdf'))

    # The same data as a tidy data frame, e.g. for other plotting libraries.
    df = time_series_to_dataframe(
        result, labels=InfectionState.values(), groups=AGE_GROUPS,
        start_date=START_DATE)
    print(df.head())
    print(df.pivot_table(index='Date', columns='Compartments', values='Values',
                         aggfunc='sum', observed=False).tail())

    if show_plot:
        plt.show()
    plt.close('all')


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(
        'plotSimulationResults',
        description='Plots the result of an ODE SEIR simulation with the '
        'standard TimeSeries plot of the MEmilio plot package.')
    arg_parser.add_argument('-p', '--show_plot', action='store_true',
                            help='Show the figures interactively.')
    arg_parser.add_argument('-o', '--output_path', default='.',
                            help='Directory the figures are written to.')
    args = arg_parser.parse_args()
    plot_results(run_ode_seir_simulation(), args.output_path, args.show_plot)
