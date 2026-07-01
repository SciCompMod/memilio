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

import numpy as np
import matplotlib.pyplot as plt

import memilio.simulation as mio
import memilio.simulation.osecir as osecir


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

    if (plot_results):
        results = [region0_result.as_ndarray(), region1_result.as_ndarray()]
        t = results[0][0, :]
        tick_range = (np.arange(int((len(t) - 1) / 10) + 1) * 10)
        tick_range[-1] -= 1

        fig, ax = plt.subplots(figsize=(10, 6))
        for idx, result_region in enumerate(results):
            region_label = f'Region {idx}'
            ax.plot(t, result_region[1, :],
                    label=f'{region_label} - #Susceptible')
            ax.plot(t, result_region[2, :], label=f'{region_label} - #Exposed')
            ax.plot(t, result_region[3, :] + result_region[4, :],
                    label=f'{region_label} - #InfectedNoSymptoms')
            ax.plot(t, result_region[5, :] + result_region[6, :],
                    label=f'{region_label} - #InfectedSymptoms')
            ax.plot(t, result_region[7, :],
                    label=f'{region_label} - #Hospitalzed')
            ax.plot(t, result_region[8, :],
                    label=f'{region_label} - #InfectedCritical')
            ax.plot(t, result_region[9, :],
                    label=f'{region_label} - #Recovered')
            ax.plot(t, result_region[10, :], label=f'{region_label} - #Dead')

        ax.set_title(
            "ODE SECIR simulation results for both regions (entire population)")
        ax.set_xticks(tick_range)
        ax.legend(loc='upper right', bbox_to_anchor=(1, 0.6))
        plt.yscale('log')
        fig.tight_layout
        fig.savefig('osecir_mobility_by_compartments.pdf')

        fig, ax = plt.subplots(5, 2, figsize=(12, 15))
        compartments = [
            'Susceptible', 'Exposed', 'InfectedNoSymptoms',
            'InfectedNoSymptomsConfirmed', 'InfectedSymptoms',
            'InfectedSymptomsConfirmed', 'InfectedSevere', 'InfectedCritical',
            'Recovered', 'Dead']
        num_compartments = len(compartments)

        for i, title in zip(range(num_compartments), compartments):
            ax[int(np.floor(i / 2)), int(i % 2)].plot(t,
                                                      results[0][i+1, :], label="Region 0")
            ax[int(np.floor(i / 2)), int(i % 2)].plot(t,
                                                      results[1][i+1, :], label="Region 1")
            ax[int(np.floor(i / 2)), int(i % 2)].set_title(title, fontsize=10)
            ax[int(np.floor(i / 2)), int(i % 2)].legend()
            ax[int(np.floor(i / 2)), int(i % 2)].set_xticks(tick_range)

        plt.subplots_adjust(hspace=0.5, bottom=0.1, top=0.9)
        fig.suptitle('Simulation results for each region in each compartment')
        fig.savefig('osecir_region_results_compartments.pdf')


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'ode_secir_mobility',
        description='Example demonstrating the setup and simulation of a geographically resolved ODE SECIHURD model with mobility.')
    args = arg_parser.parse_args()
    run_ode_secir_mobility_simulation()
