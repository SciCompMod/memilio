#############################################################################
# Copyright (C) 2020-2026 MEmilio
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

import matplotlib.pyplot as plt
import numpy as np

import memilio.simulation as mio
import memilio.simulation.osecir as osecir
from memilio.plot.plotTimeSeries import plot_time_series


def run_ode_secir_mobility_simulation(plot_results=True):
    """

    :param plot_results:  (Default value = True)

    """
    mio.set_log_level(mio.LogLevel.Warning)

    t0 = 0
    tmax = 50

    # setup basic parameters
    model = osecir.Model(1)

    model.parameters.TimeExposed[mio.AgeGroup(0)] = 3.2
    model.parameters.TimeInfectedNoSymptoms[mio.AgeGroup(0)] = 2.
    model.parameters.TimeInfectedSymptoms[mio.AgeGroup(0)] = 6
    model.parameters.TimeInfectedSevere[mio.AgeGroup(0)] = 12
    model.parameters.TimeInfectedCritical[mio.AgeGroup(0)] = 8

    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.r_[0.5]
    model.parameters.ContactPatterns.cont_freq_mat[0].add_damping(
        mio.Damping(np.r_[0.3], t=0.3))

    model.parameters.TransmissionProbabilityOnContact[mio.AgeGroup(0)] = 1.0
    model.parameters.RelativeTransmissionNoSymptoms[mio.AgeGroup(0)] = 0.67
    model.parameters.RecoveredPerInfectedNoSymptoms[mio.AgeGroup(0)] = 0.09
    model.parameters.RiskOfInfectionFromSymptomatic[mio.AgeGroup(0)] = 0.25
    model.parameters.SeverePerInfectedSymptoms[mio.AgeGroup(0)] = 0.2
    model.parameters.CriticalPerSevere[mio.AgeGroup(0)] = 0.25
    model.parameters.DeathsPerCritical[mio.AgeGroup(0)] = 0.3

    # two regions with different populations and with some mobility between them
    graph = osecir.MobilityGraph()
    model.populations[mio.AgeGroup(0), osecir.InfectionState.Exposed] = 100
    model.populations[mio.AgeGroup(
        0), osecir.InfectionState.InfectedNoSymptoms] = 50
    model.populations[mio.AgeGroup(
        0), osecir.InfectionState.InfectedSymptoms] = 50
    model.populations[mio.AgeGroup(
        0), osecir.InfectionState.InfectedSevere] = 20
    model.populations[mio.AgeGroup(
        0), osecir.InfectionState.InfectedCritical] = 10
    model.populations[mio.AgeGroup(0), osecir.InfectionState.Recovered] = 10
    model.populations[mio.AgeGroup(0), osecir.InfectionState.Dead] = 0
    model.populations.set_difference_from_group_total_AgeGroup((
        mio.AgeGroup(0),
        osecir.InfectionState.Susceptible),
        10000)
    model.apply_constraints()
    graph.add_node(id=0, model=model, t0=t0)  # copies the model into the graph
    model.populations[mio.AgeGroup(0), osecir.InfectionState.Exposed] = 0
    model.populations[mio.AgeGroup(
        0), osecir.InfectionState.InfectedNoSymptoms] = 0
    model.populations[mio.AgeGroup(
        0), osecir.InfectionState.InfectedNoSymptomsConfirmed] = 0
    model.populations[mio.AgeGroup(
        0), osecir.InfectionState.InfectedSymptoms] = 0
    model.populations[mio.AgeGroup(
        0), osecir.InfectionState.InfectedSymptomsConfirmed] = 0
    model.populations[mio.AgeGroup(
        0), osecir.InfectionState.InfectedSevere] = 0
    model.populations[mio.AgeGroup(
        0), osecir.InfectionState.InfectedCritical] = 0
    model.populations[mio.AgeGroup(0), osecir.InfectionState.Recovered] = 0
    model.populations[mio.AgeGroup(0), osecir.InfectionState.Dead] = 0
    model.populations.set_difference_from_group_total_AgeGroup((
        mio.AgeGroup(0),
        osecir.InfectionState.Susceptible),
        2000)
    model.apply_constraints()
    graph.add_node(id=1, model=model, t0=t0)
    mobility_coefficients = 0.1 * np.ones(model.populations.numel())
    mobility_coefficients[osecir.InfectionState.Dead] = 0
    mobility_params = mio.MobilityParameters(mobility_coefficients)
    # one coefficient per (age group x compartment)
    graph.add_edge(0, 1, mobility_params)
    graph.add_edge(1, 0, mobility_params)

    # run simulation
    sim = osecir.MobilitySimulation(graph, t0, dt=0.5)
    sim.advance(tmax)

    # process results
    region0_result = osecir.interpolate_simulation_result(
        sim.graph.get_node(0).property.result)
    region1_result = osecir.interpolate_simulation_result(
        sim.graph.get_node(1).property.result)

    if plot_results:
        region_results = [region0_result, region1_result]
        region_labels = ['Region 0', 'Region 1']

        # All compartments of each region on a logarithmic axis.
        fig, axes = plt.subplots(1, 2, figsize=(16, 5), layout='constrained')
        for region_result, region_label, ax in zip(
                region_results, region_labels, axes):
            plot_time_series(
                region_result, labels=osecir.InfectionState.values(),
                ax=ax, title=region_label)
        fig.suptitle('ODE SECIR simulation results for both regions')
        fig.savefig('osecir_mobility_by_compartments.pdf')

        # Stack the results of the regions into one array with the elements
        # ordered region by region (the interpolated results share the same
        # time points), so that the regions can be passed as groups.
        results = np.vstack([region0_result.as_ndarray(),
                             region1_result.as_ndarray()[1:]])

        # One panel per compartment with one line per region.
        fig, axes = plt.subplots(5, 2, figsize=(14, 16), layout='constrained')
        for state, ax in zip(osecir.InfectionState.values(), axes.flat):
            plot_time_series(
                results, labels=osecir.InfectionState.values(),
                groups=region_labels, sum_groups=False, select=state, ax=ax,
                title=state.name)
        fig.suptitle('Simulation results for each region in each compartment')
        fig.savefig('osecir_region_results_compartments.pdf')
        plt.close('all')


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'ode_secir_mobility',
        description='Example demonstrating the setup and simulation of a geographically resolved ODE SECIHURD model with mobility.')
    args = arg_parser.parse_args()
    run_ode_secir_mobility_simulation()
