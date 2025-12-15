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


def run_ode_secir_parameter_study():
    """ """
    # setup basic parameters
    num_groups = 6
    model = osecir.Model(num_groups)

    for i in range(num_groups):
        group = mio.AgeGroup(i)

        model.parameters.TimeExposed[group] = 3.2
        model.parameters.TimeInfectedNoSymptoms[group] = 2.
        model.parameters.TimeInfectedSymptoms[group] = 6
        model.parameters.TimeInfectedSevere[group] = 12
        model.parameters.TimeInfectedCritical[group] = 8

        model.populations[group, osecir.InfectionState.Exposed] = 100
        model.populations[group, osecir.InfectionState.InfectedNoSymptoms] = 50
        model.populations[group,
                          osecir.InfectionState.InfectedNoSymptomsConfirmed] = 0
        model.populations[group, osecir.InfectionState.InfectedSymptoms] = 50
        model.populations[group,
                          osecir.InfectionState.InfectedSymptomsConfirmed] = 0
        model.populations[group, osecir.InfectionState.InfectedSevere] = 20
        model.populations[group, osecir.InfectionState.InfectedCritical] = 10
        model.populations[group, osecir.InfectionState.Recovered] = 10
        model.populations[group, osecir.InfectionState.Dead] = 0
        model.populations.set_difference_from_group_total_AgeGroup(
            (group, osecir.InfectionState.Susceptible), 10000)

        model.parameters.TransmissionProbabilityOnContact[group].set_distribution(
            mio.ParameterDistributionUniform(0.1, 0.2))
        model.parameters.RelativeTransmissionNoSymptoms[group] = 0.67
        model.parameters.RecoveredPerInfectedNoSymptoms[group] = 0.09
        model.parameters.RiskOfInfectionFromSymptomatic[group] = 0.25
        model.parameters.SeverePerInfectedSymptoms[group] = 0.2
        model.parameters.CriticalPerSevere[group] = 0.25
        model.parameters.DeathsPerCritical[group] = 0.3

    model.parameters.ContactPatterns.cont_freq_mat = mio.ContactMatrixGroup(
        4, num_groups)
    model.parameters.ContactPatterns.cont_freq_mat[0] = mio.ContactMatrix(
        np.ones((num_groups, num_groups))*0.5)
    model.parameters.ContactPatterns.cont_freq_mat[1] = mio.ContactMatrix(
        np.ones((num_groups, num_groups))*0.5)
    model.parameters.ContactPatterns.cont_freq_mat[2] = mio.ContactMatrix(
        np.ones((num_groups, num_groups))*0.5)
    model.parameters.ContactPatterns.cont_freq_mat[3] = mio.ContactMatrix(
        np.ones((num_groups, num_groups))*0.5)
    model.parameters.ContactPatterns.cont_freq_mat.add_damping(
        mio.Damping(np.ones((num_groups, num_groups))*0.7, 7.0))
    print(model.parameters.ContactPatterns.cont_freq_mat[1].baseline)
    # process the result of one run

    def handle_result(graph, run_idx):
        """

        :param graph: 
        :param run_idx: 

        """
        group = mio.AgeGroup(0)
        print("run {} with infection rate {:.2G}".format(handle_result.c, graph.get_node(
            0).property.model.parameters.TransmissionProbabilityOnContact[group].value))
        print("compartments at t = {}:".format(
            graph.get_node(0).property.result.get_time(0)))
        print(graph.get_node(0).property.result.get_value(0))
        print("compartments at t = {}:".format(
            graph.get_node(0).property.result.get_last_time()))
        print(graph.get_node(0).property.result.get_last_value())
        handle_result.c += 1
    handle_result.c = 0

    # study the effect of different infection rates

    model.apply_constraints()

    graph = osecir.ModelGraph()
    graph.add_node(0, model)
    graph.add_node(1, model)
    mobility_coefficients = 0.01 * np.ones(model.populations.numel())
    for i in range(num_groups):
        flat_index = model.populations.get_flat_index(
            osecir.MultiIndex_PopulationsArray(mio.AgeGroup(i), osecir.InfectionState.Dead))
        mobility_coefficients[flat_index] = 0
    graph.add_edge(0, 1, mobility_coefficients)
    graph.add_edge(1, 0, mobility_coefficients)

    study = osecir.GraphParameterStudy(
        graph, t0=0, tmax=10, dt=0.5, num_runs=3)
    study.run(handle_result)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'ode_secir_parameter_study',
        description='Example demonstrating ensemble runs of a ODE SECIHURD model.')
    args = arg_parser.parse_args()
    run_ode_secir_parameter_study()
