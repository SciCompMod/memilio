#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Martin J. Kuehn, Wadim Koslow, Daniel Abele, Khoa Nguyen
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
from datetime import date

import matplotlib.pyplot as plt
import numpy as np

from memilio.plot.plotTimeSeries import plot_time_series
from memilio.simulation import AgeGroup, Damping
from memilio.simulation.osecir import InfectionState as State
from memilio.simulation.osecir import (Model, interpolate_simulation_result,
                                       simulate)

import openturns as ot
import openturns.viewer as otv


def run_ode_secir_simulation(vars):
    """Runs the c++ ODE SECIHURD model using one age group
    """

    (TimeExposed, TimeInfectedNoSymptoms, TimeInfectedSymptoms,
     TimeInfectedSevere, TimeInfectedCritical, TransmissionProbabilityOnContact,
     RelativeTransmissionNoSymptoms, RiskOfInfectionFromSymptomatic,
     RecoveredPerInfectedNoSymptoms, SeverePerInfectedSymptoms, CriticalPerSevere, DeathsPerCritical) = (
        float(vars[i]) for i in range(len(vars))
    )

    # Define population of age groups
    populations = [83000]

    days = 30  # number of days to simulate
    start_day = 1
    start_month = 1
    start_year = 2019
    dt = 0.1
    num_groups = 1

    # Initialize Parameters
    model = Model(1)

    A0 = AgeGroup(0)

    # Set parameters
    # Compartment transition duration
    model.parameters.TimeExposed[A0] = TimeExposed
    model.parameters.TimeInfectedNoSymptoms[A0] = TimeInfectedNoSymptoms
    model.parameters.TimeInfectedSymptoms[A0] = TimeInfectedSymptoms
    model.parameters.TimeInfectedSevere[A0] = TimeInfectedSevere
    model.parameters.TimeInfectedCritical[A0] = TimeInfectedCritical

    # Initial number of people in each compartment
    model.populations[A0, State.Exposed] = 10
    model.populations[A0, State.InfectedNoSymptoms] = 0
    model.populations[A0, State.InfectedNoSymptomsConfirmed] = 0
    model.populations[A0, State.InfectedSymptoms] = 0
    model.populations[A0, State.InfectedSymptomsConfirmed] = 0
    model.populations[A0, State.InfectedSevere] = 0
    model.populations[A0, State.InfectedCritical] = 0
    model.populations[A0, State.Recovered] = 0
    model.populations[A0, State.Dead] = 0
    model.populations.set_difference_from_total(
        (A0, State.Susceptible), populations[0])

    # Compartment transition propabilities
    model.parameters.RelativeTransmissionNoSymptoms[A0] = RelativeTransmissionNoSymptoms
    model.parameters.TransmissionProbabilityOnContact[A0] = TransmissionProbabilityOnContact
    model.parameters.RecoveredPerInfectedNoSymptoms[A0] = RecoveredPerInfectedNoSymptoms
    model.parameters.RiskOfInfectionFromSymptomatic[A0] = RiskOfInfectionFromSymptomatic
    model.parameters.SeverePerInfectedSymptoms[A0] = SeverePerInfectedSymptoms
    model.parameters.CriticalPerSevere[A0] = CriticalPerSevere
    model.parameters.DeathsPerCritical[A0] = DeathsPerCritical

    model.parameters.StartDay = (
        date(start_year, start_month, start_day) - date(start_year, 1, 1)).days

    # model.parameters.ContactPatterns.cont_freq_mat[0] = ContactMatrix(np.r_[0.5])
    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
        (num_groups, num_groups)) * 10

    # Run Simulation
    result = simulate(0, days, dt, model)
    # interpolate results
    result = interpolate_simulation_result(result)

    # print(result.get_last_value())

    deaths = float(result.get_last_value()[int(State.Dead)])
    total_infections = populations[0] - \
        float(result.get_last_value()[int(State.Susceptible)])

    # return [total_infections]
    return [deaths]


def main():

    # Define parameters to be drawn
    variable_list = [
        ("TimeExposed", 3.335),
        ("TimeInfectedNoSymptoms", 2.58916),
        ("TimeInfectedSymptoms", 6.94547),
        ("TimeInfectedSevere", 7.28196),
        ("TimeInfectedCritical", 13.066),
        ("TransmissionProbabilityOnContact", 0.07333),
        ("RelativeTransmissionNoSymptoms", 1),
        ("RiskOfInfectionFromSymptomatic", 0.25),
        ("RecoveredPerInfectedNoSymptoms", 1-0.79310),
        ("SeverePerInfectedSymptoms", 0.07864),
        ("CriticalPerSevere", 0.17318),
        ("DeathsPerCritical", 0.21718),
    ]

    # convert names and define marginals
    names = [name for name, _ in variable_list]
    marginals = [ot.Uniform(0.5*mean, 1.5*mean) for _, mean in variable_list]

    # define joint distribution
    distribution = ot.JointDistribution(marginals)
    distribution.setDescription(names)

    # define seed and sample size
    ot.RandomGenerator.SetSeed(10)
    size = 1000

    # compute parameter samples
    X = distribution.getSample(size)
    X.setDescription([
        r"$T_E$", r"$T_{I_{NS}}$", r"$T_{I_{Sy}}$", r"$T_{I_{Sev}}$", r"$T_{I_{Cr}}$", r"$\rho$", r"$\xi_{NS}$", r"$\xi_{Sy}$",
        r"$\mu_{I_{NS}}^R$", r"$\mu_{I_{Sy}}^{I_{Sev}}$", r"$\mu_{I_{Sev}}^{I_{Cr}}$", r"$\mu_{I_{Cr}}^{D}$"])

    # define and evaluate MEmilio model
    model = ot.PythonFunction(
        len(variable_list), 1, run_ode_secir_simulation
    )
    Y = model(X)
    # Y.setDescription(["Total infections"])
    Y.setDescription(["Deaths"])

    sobol = ot.RankSobolSensitivityAlgorithm(X, Y)
    sobol.setBootstrapSize(200)
    sobol.setConfidenceLevel(0.95)

    indices = sobol.getFirstOrderIndices()
    intervals = sobol.getFirstOrderIndicesInterval()

    lower = intervals.getLowerBound()
    upper = intervals.getUpperBound()

    # Print the sensitivity ranking
    print("First-order Sobol indices:")
    for i in sorted(
        range(len(names)), key=lambda j: indices[j], reverse=True
    ):
        print(
            f"  {names[i]:36s} S1={indices[i]:7.3f} "
            f"CI=[{lower[i]:7.3f}, {upper[i]:7.3f}]"
        )

    # Draw the Sobol indices and confidence intervals
    graph = sobol.draw()
    graph.setTitle("")
    view = otv.View(graph, figure_kw={"figsize": (11, 5)})

    ax = view.getAxes()[0]
    # ax.set_title("Total infections after 30 days", fontsize=24)
    ax.set_title("Deaths after 30 days", fontsize=24)
    ax.legend(fontsize=20)
    ax.tick_params(axis="x", labelsize=20)
    ax.tick_params(axis="y", labelsize=20)
    ax.set_ylabel("Index value", fontsize=20)
    ax.set_xlabel("Input", fontsize=20)

    # Remove parameter names displayed inside the plot
    for txt in ax.texts:
        txt.set_visible(False)

    # Add parameter names as real x-axis ticks
    names = list(X.getDescription())
    positions = range(1, len(names) + 1)

    ax.set_xticks(positions)
    ax.set_xticklabels(names)

    # Set maximum y value
    ax.set_ylim(-0.1, 0.5)

    view.save("memilio_sobol_rank.png", dpi=300, bbox_inches="tight")
    otv.View.ShowAll()


if __name__ == "__main__":
    main()
