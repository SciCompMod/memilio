#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
import memilio.simulation.secir as secir


def parameter_study():
    mio.set_log_level(mio.LogLevel.Warning)

    # setup basic parameters
    model = secir.Model(1)

    model.parameters.IncubationTime[secir.AgeGroup(0)] = 5.2
    model.parameters.SerialInterval[secir.AgeGroup(0)] = 4.2
    model.parameters.TimeInfectedSymptoms[secir.AgeGroup(0)] = 6
    model.parameters.TimeInfectedSevere[secir.AgeGroup(0)] = 12
    model.parameters.TimeInfectedCritical[secir.AgeGroup(0)] = 8

    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.r_[0.5]
    model.parameters.ContactPatterns.cont_freq_mat[0].add_damping(
        mio.Damping(np.r_[0.3], t=0.3))

    model.parameters.TransmissionProbabilityOnContact[secir.AgeGroup(0)] = 1.0
    model.parameters.RelativeTransmissionNoSymptoms[secir.AgeGroup(0)] = 0.67
    model.parameters.RecoveredPerInfectedNoSymptoms[secir.AgeGroup(0)] = 0.09
    model.parameters.RiskOfInfectionFromSymptomatic[secir.AgeGroup(0)] = 0.25
    model.parameters.SeverePerInfectedSymptoms[secir.AgeGroup(0)] = 0.2
    model.parameters.CriticalPerSevere[secir.AgeGroup(0)] = 0.25
    model.parameters.DeathsPerCritical[secir.AgeGroup(0)] = 0.3

    # two regions with different populations and with some migration between them
    graph = secir.ModelGraph()
    model.populations[secir.AgeGroup(0), secir.InfectionState.Exposed] = 100
    model.populations[secir.AgeGroup(
        0), secir.InfectionState.InfectedNoSymptoms] = 50
    model.populations[secir.AgeGroup(
        0), secir.InfectionState.InfectedSymptoms] = 50
    model.populations[secir.AgeGroup(
        0), secir.InfectionState.InfectedSevere] = 20
    model.populations[secir.AgeGroup(
        0), secir.InfectionState.InfectedCritical] = 10
    model.populations[secir.AgeGroup(0), secir.InfectionState.Recovered] = 10
    model.populations[secir.AgeGroup(0), secir.InfectionState.Dead] = 0
    model.populations.set_difference_from_group_total_AgeGroup((
        secir.AgeGroup(0),
        secir.InfectionState.Susceptible),
        10000)
    model.apply_constraints()
    graph.add_node(id=0, model=model)  # copies the model into the graph
    model.populations[secir.AgeGroup(0), secir.InfectionState.Exposed] = 0
    model.populations[secir.AgeGroup(
        0), secir.InfectionState.InfectedNoSymptoms] = 0
    model.populations[secir.AgeGroup(
        0), secir.InfectionState.InfectedSymptoms] = 0
    model.populations[secir.AgeGroup(
        0), secir.InfectionState.InfectedSevere] = 0
    model.populations[secir.AgeGroup(
        0), secir.InfectionState.InfectedCritical] = 0
    model.populations[secir.AgeGroup(0), secir.InfectionState.Recovered] = 0
    model.populations[secir.AgeGroup(0), secir.InfectionState.Dead] = 0
    model.populations.set_difference_from_group_total_AgeGroup((
        secir.AgeGroup(0),
        secir.InfectionState.Susceptible),
        2000)
    model.apply_constraints()
    graph.add_node(id=1, model=model)
    migration_coefficients = 0.1 * np.ones(model.populations.numel())
    migration_params = mio.MigrationParameters(migration_coefficients)
    # one coefficient per (age group x compartment)
    graph.add_edge(0, 1, migration_params)
    # directed graph -> add both directions so coefficients can be different
    graph.add_edge(1, 0, migration_params)

    # process the result of one run
    def handle_result(graph, run_idx):
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
                    [secir.AgeGroup(0),
                     secir.InfectionState.InfectedNoSymptoms].value))
            print(f"  compartments at t = {result.get_time(0)}:")
            print("  ", result.get_value(0))
            print(f"  compartments at t = {result.get_last_time()}:")
            print("  ", result.get_last_value())
    handle_result.c = 0

    # study with unknown number of undetected InfectedNoSymptoms
    carrier_distribution = mio.ParameterDistributionNormal(50, 2000, 200, 100)
    graph.get_node(0).property.populations[secir.AgeGroup(
        0), secir.InfectionState.InfectedNoSymptoms].set_distribution(carrier_distribution)

    t0 = 0
    tmax = 50
    study = secir.ParameterStudy(graph, t0, tmax, dt=1.0, num_runs=3)
    study.run(handle_result)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'migration_parameter_study',
        description='Example demonstrating setup of ensemble runs of a geographically resolved SECIR model with travel.')
    args = arg_parser.parse_args()
    parameter_study()
