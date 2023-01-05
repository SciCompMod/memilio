#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Martin J. Kuehn, Wadim Koslow, Annalena Lange
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

from enum import Enum, auto

from memilio.simulation import ContactMatrix, Damping, UncertainContactMatrix
from memilio.simulation.secir import AgeGroup, Index_InfectionState
from memilio.simulation.secir import InfectionState as State
from memilio.simulation.secir import (Model, Simulation,
                                      interpolate_simulation_result, simulate)


class ContactLocation(Enum):
    Home = auto()
    School = auto()
    Work = auto()
    Other = auto()


class Intervention(Enum):
    Home = auto()
    SchoolClosure = auto()
    HomeOffice = auto()
    GatheringBanFacilitiesClosure = auto()
    PhysicalDistanceAndMasks = auto()
    SeniorAwareness = auto()


class InterventionLevel(Enum):
    Main = auto()
    PhysicalDistanceAndMasks = auto()
    SeniorAwareness = auto()
    Holidays = auto()


class SecirSimulation:

    def __init__(self, data_dir):
        t0 = 0.0
        self.num_age_groups = 6
        self.model = Model(self.num_age_groups)
        self.parameters = self.model.parameters

    def set_covid_parameters(self):

        # times
        incubationTime = 5.2
        serialIntervalMin = 0.5 * 2.67 + 0.5 * 5.2
        serialIntervalMax = 0.5 * 4.0 + 0.5 * 5.2
        timeInfectedSymptomsMin = [5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465]
        timeInfectedSymptomsMax = [8.427, 8.427, 8.4684, 8.3139, 8.169, 8.085]
        timeInfectedSevereMin = [3.925, 3.925, 4.85, 6.4, 7.2, 9.0]
        timeInfectedSevereMax = [6.075, 6.075, 7.0, 8.7, 9.8, 13.0]
        timeInfectedCriticalMin = [4.95, 4.95, 4.86, 14.14, 14.4, 10.0]
        timeInfectedCriticalMax = [8.95, 8.95, 8.86, 20.58, 19.8, 13.2]

        # probabilities
        transmissionProbabilityOnContactMin = [0.02, 0.05, 0.05, 0.05, 0.08, 0.15]
        transmissionProbabilityOnContactMax = [0.04, 0.07, 0.07, 0.07, 0.10, 0.20]
        relativeTransmissionNoSymptomsMin = 1.0
        relativeTransmissionNoSymptomsMax = 1.0
        riskOfInfectionFromSymptomaticMin = 0.1
        riskOfInfectionFromSymptomaticMax = 0.3
        maxRiskOfInfectionFromSymptomaticMin = 0.3
        maxRiskOfInfectionFromSymptomaticMax = 0.5
        recoveredPerInfectedNoSymptomsMin = [0.2, 0.2, 0.15, 0.15, 0.15, 0.15]
        recoveredPerInfectedNoSymptomsMax = [0.3, 0.3, 0.25, 0.25, 0.25, 0.25]
        severePerInfectedSymptomsMin = [0.006, 0.006, 0.015, 0.049, 0.15, 0.20]
        severePerInfectedSymptomsMax = [0.009, 0.009, 0.023, 0.074, 0.18, 0.25]
        criticalPerSevereMin = [0.05, 0.05, 0.05, 0.10, 0.25, 0.35]
        criticalPerSevereMax = [0.10, 0.10, 0.10, 0.20, 0.35, 0.45]
        deathsPerCriticalMin = [0.0, 0.0, 0.10, 0.10, 0.3, 0.5]
        deathsPerCriticalMax = [0.1, 0.1, 0.18, 0.18, 0.5, 0.7]

        # seasonality
        seasonalityMin = 0.1
        seasonalityMax = 0.3
