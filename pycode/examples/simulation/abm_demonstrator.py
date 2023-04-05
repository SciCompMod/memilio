#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Martin J. Kuehn, Wadim Koslow, Daniel Abele
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
import pandas as pd
import numpy as np

from memilio.simulation import abm
from memilio.simulation.abm import AgeGroup
from memilio.simulation.abm import VaccinationState

class LocationMapping:

    def __init__(self):
        self.inputId = None
        self.modelId = None

def set_infection_parameters():
    infection_params = abm.GlobalInfectionParameters()

    #0-4
    infection_params.IncubationPeriod[AgeGroup.Age0to4, VaccinationState.Unvaccinated] = 4
    infection_params.SusceptibleToExposedByCarrier[AgeGroup.Age0to4, VaccinationState.Unvaccinated] = 0.05
    infection_params.SusceptibleToExposedByInfected[AgeGroup.Age0to4, VaccinationState.Unvaccinated] = 0.05
    infection_params.CarrierToInfected[AgeGroup.Age0to4, VaccinationState.Unvaccinated] = 0.276
    infection_params.CarrierToRecovered[AgeGroup.Age0to4, VaccinationState.Unvaccinated] = 0.092
    infection_params.InfectedToRecovered[AgeGroup.Age0to4, VaccinationState.Unvaccinated] = 0.142
    infection_params.InfectedToSevere[AgeGroup.Age0to4, VaccinationState.Unvaccinated] = 0.001
    infection_params.SevereToRecovered[AgeGroup.Age0to4, VaccinationState.Unvaccinated] = 0.186
    infection_params.SevereToCritical[AgeGroup.Age0to4, VaccinationState.Unvaccinated] = 0.015
    infection_params.CriticalToRecovered[AgeGroup.Age0to4, VaccinationState.Unvaccinated] = 0.143
    infection_params.CriticalToDead[AgeGroup.Age0to4, VaccinationState.Unvaccinated] = 0.001
    infection_params.RecoveredToSusceptible[AgeGroup.Age0to4, VaccinationState.Unvaccinated] = 0

    #5-14
    infection_params.IncubationPeriod[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 4
    infection_params.SusceptibleToExposedByCarrier[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.1
    infection_params.SusceptibleToExposedByInfected[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.1
    infection_params.CarrierToInfected[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.276
    infection_params.CarrierToRecovered[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.092
    infection_params.InfectedToRecovered[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.142
    infection_params.InfectedToSevere[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.001
    infection_params.SevereToRecovered[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.186
    infection_params.SevereToCritical[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.015
    infection_params.CriticalToRecovered[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.143
    infection_params.CriticalToDead[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.001
    infection_params.RecoveredToSusceptible[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.

    #15-34
    infection_params.IncubationPeriod[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 4
    infection_params.SusceptibleToExposedByCarrier[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.13
    infection_params.SusceptibleToExposedByInfected[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.13
    infection_params.CarrierToInfected[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.315
    infection_params.CarrierToRecovered[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.079
    infection_params.InfectedToRecovered[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.139
    infection_params.InfectedToSevere[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.003
    infection_params.SevereToRecovered[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.157
    infection_params.SevereToCritical[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.013
    infection_params.CriticalToRecovered[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.126
    infection_params.CriticalToDead[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.021
    infection_params.RecoveredToSusceptible[AgeGroup.Age5to14, VaccinationState.Unvaccinated] = 0.

    #35-59
    infection_params.IncubationPeriod[AgeGroup.Age35to59, VaccinationState.Unvaccinated] = 4
    infection_params.SusceptibleToExposedByCarrier[AgeGroup.Age35to59, VaccinationState.Unvaccinated] = 0.11
    infection_params.SusceptibleToExposedByInfected[AgeGroup.Age35to59, VaccinationState.Unvaccinated] = 0.11
    infection_params.CarrierToInfected[AgeGroup.Age35to59, VaccinationState.Unvaccinated] = 0.315
    infection_params.CarrierToRecovered[AgeGroup.Age35to59, VaccinationState.Unvaccinated] = 0.079
    infection_params.InfectedToRecovered[AgeGroup.Age35to59, VaccinationState.Unvaccinated] = 0.136
    infection_params.InfectedToSevere[AgeGroup.Age35to59, VaccinationState.Unvaccinated] = 0.009
    infection_params.SevereToRecovered[AgeGroup.Age35to59, VaccinationState.Unvaccinated] = 0.113
    infection_params.SevereToCritical[AgeGroup.Age35to59, VaccinationState.Unvaccinated] = 0.02
    infection_params.CriticalToRecovered[AgeGroup.Age35to59, VaccinationState.Unvaccinated] = 0.05
    infection_params.CriticalToDead[AgeGroup.Age35to59, VaccinationState.Unvaccinated] = 0.008
    infection_params.RecoveredToSusceptible[AgeGroup.Age35to59, VaccinationState.Unvaccinated] = 0.

    #60-79
    infection_params.IncubationPeriod[AgeGroup.Age60to79, VaccinationState.Unvaccinated] = 4
    infection_params.SusceptibleToExposedByCarrier[AgeGroup.Age60to79, VaccinationState.Unvaccinated] = 0.04
    infection_params.SusceptibleToExposedByInfected[AgeGroup.Age60to79, VaccinationState.Unvaccinated] = 0.04
    infection_params.CarrierToInfected[AgeGroup.Age60to79, VaccinationState.Unvaccinated] = 0.315
    infection_params.CarrierToRecovered[AgeGroup.Age60to79, VaccinationState.Unvaccinated] = 0.079
    infection_params.InfectedToRecovered[AgeGroup.Age60to79, VaccinationState.Unvaccinated] = 0.123
    infection_params.InfectedToSevere[AgeGroup.Age60to79, VaccinationState.Unvaccinated] = 0.024
    infection_params.SevereToRecovered[AgeGroup.Age60to79, VaccinationState.Unvaccinated] = 0.083
    infection_params.SevereToCritical[AgeGroup.Age60to79, VaccinationState.Unvaccinated] = 0.035
    infection_params.CriticalToRecovered[AgeGroup.Age60to79, VaccinationState.Unvaccinated] = 0.035
    infection_params.CriticalToDead[AgeGroup.Age60to79, VaccinationState.Unvaccinated] = 0.023
    infection_params.RecoveredToSusceptible[AgeGroup.Age60to79, VaccinationState.Unvaccinated] = 0.

    

    return infection_params

def read_txt(path):
    input = pd.read_csv(path, sep='\t')
    return input

def add_households(world, distribution, num_inhabitants):
    household_sizes = np.zeros(len(distribution))
    while num_inhabitants > 0:
        size = np.random.choice(np.arange(0, len(distribution)), p = distribution)
        household_sizes[size] += 1
        num_inhabitants -= (size + 1)

    #one-person household
    one_person_household_member = abm.HouseholdMember()

def create_locations_from_input(world, input_areas):
    has_school = False
    has_hospital = False
    for index, area in input_areas.iterrows():
        if('residential' in area.type):
            add_households(world, [0., 0.1, 0., 0., 0.9], area.inhabitants)
        elif(area.type == 'recreational'):
            continue
        elif(area.type == 'shopping_business'):
            continue
        elif(area.type == 'university'):
            continue
        elif(area.type == 'mixed'):
            continue

def run_abm_simulation():

    LocationIds = []

    infection_params = set_infection_parameters()
    world = abm.World(infection_params)

    t0 = abm.TimePoint(0)
    tmax = t0 + abm.days(14)

    areas = read_txt('C:/Users/bick_ju/Documents/INSIDe/Demonstrator/INSIDeDemonstrator/INSIDe_Demonstrator_AreaList.txt')
    create_locations_from_input(world, areas)
    #print(areas)

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'abm demonstrator',
        description='Example demonstrating the agent-based model for a synthetic population.')
    args = arg_parser.parse_args()
    run_abm_simulation(**args.__dict__)