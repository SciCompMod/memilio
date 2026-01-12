#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors:
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

import memilio.simulation as mio
import memilio.simulation.osecir as osecir


def run_mobility_parameter_study():
    """ """
    mio.set_log_level(mio.LogLevel.Warning)

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
    graph = osecir.ModelGraph()
    model.populations[mio.AgeGroup(0), osecir.InfectionState.Exposed] = 100
    model.populations[mio.AgeGroup(
        0), osecir.InfectionState.InfectedNoSymptoms] = 50
    model.populations[mio.AgeGroup(
        0), osecir.InfectionState.InfectedNoSymptomsConfirmed] = 0
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
    graph.add_node(id=0, model=model)  # copies the model into the graph
    model.populations[mio.AgeGroup(0), osecir.InfectionState.Exposed] = 0
    model.populations[mio.AgeGroup(
        0), osecir.InfectionState.InfectedNoSymptoms] = 0
    model.populations[mio.AgeGroup(
        0), osecir.InfectionState.InfectedSymptoms] = 0
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
    graph.add_node(id=1, model=model)
    mobility_coefficients = 0.1 * np.ones(model.populations.numel())
    mobility_coefficients[osecir.InfectionState.Dead] = 0
    mobility_params = mio.MobilityParameters(mobility_coefficients)
    # one coefficient per (age group x compartment)
    graph.add_edge(0, 1, mobility_params)
    graph.add_edge(1, 0, mobility_params)

    # process the result of one run
    def handle_result(graph, run_idx):
        """

        :param graph: 
        :param run_idx: 

        """
        print(f'run {handle_result.c}')
        handle_result.c = handle_result.c + 1
        for node_idx in range(graph.num_nodes):
            node = graph.get_node(node_idx)
            result = node.property.result
            model = node.property.model
            print(f"  node {node_idx}")
            print(
                "  initial InfectedNoSymptoms count {}.".format(
                    model.populations
                    [mio.AgeGroup(0),
                     osecir.InfectionState.InfectedNoSymptoms].value))
            print(f"  compartments at t = {result.get_time(0)}:")
            print("  ", result.get_value(0))
            print(f"  compartments at t = {result.get_last_time()}:")
            print("  ", result.get_last_value())
    handle_result.c = 0

    # study with unknown number of undetected InfectedNoSymptoms
    carrier_distribution = mio.ParameterDistributionNormal(
        50, 2000, 200, 100, 2.5758)
    graph.get_node(0).property.populations[mio.AgeGroup(
        0), osecir.InfectionState.InfectedNoSymptoms].set_distribution(carrier_distribution)

    t0 = 0
    tmax = 50
    study = osecir.GraphParameterStudy(graph, t0, tmax, dt=1.0, num_runs=3)
    study.run(handle_result)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'mobility_parameter_study',
        description='Example demonstrating setup of ensemble runs of a geographically resolved ODE SECIHURD model with mobility.')
    args = arg_parser.parse_args()
    run_mobility_parameter_study()
