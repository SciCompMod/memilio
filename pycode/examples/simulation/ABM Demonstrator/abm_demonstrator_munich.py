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
import numpy as np
import os
import sys
import h5py
import time
import random

import memilio.simulation as mio
from memilio.simulation import abm
from memilio.simulation import AgeGroup
from memilio.simulation.abm import VirusVariant
from memilio.simulation.abm import History
from memilio.simulation.abm import Infection

import pandas as pd

# class used to map input areas to abm locations


class LocationMapping:

    def __init__(self):
        self.inputId = None
        self.modelId = []


# number of age groups
num_age_groups = 6
age_group_0_to_4 = AgeGroup(0)
age_group_5_to_14 = AgeGroup(1)
age_group_15_to_34 = AgeGroup(2)
age_group_35_to_59 = AgeGroup(3)
age_group_60_to_79 = AgeGroup(4)
age_group_80_plus = AgeGroup(5)


def age_group_to_string(age_group):
    if (age_group == age_group_0_to_4):
        return "0"
    elif (age_group == age_group_5_to_14):
        return "1"
    elif (age_group == age_group_15_to_34):
        return "2"
    elif (age_group == age_group_35_to_59):
        return "3"
    elif (age_group == age_group_60_to_79):
        return "4"
    elif (age_group == age_group_80_plus):
        return "5"
    else:
        return "AgeGroup not found"


def set_infection_parameters(parameters):

    infection_params = abm.Parameters(num_age_groups)

    # AgeGroup 0-4
    abm.set_incubationPeriod(
        infection_params, VirusVariant.Wildtype, age_group_0_to_4, parameters.loc["Age0to4_IncubationPeriod"].value, parameters.loc["Age0to4_IncubationPeriod"].dev)
    abm.set_TimeInfectedNoSymptomsToSymptoms(infection_params, VirusVariant.Wildtype, age_group_0_to_4,
                                             parameters.loc["Age0to4_InfectedNoSymptomsToSymptoms"].value, parameters.loc["Age0to4_InfectedNoSymptomsToSymptoms"].dev)
    abm.set_TimeInfectedNoSymptomsToRecovered(infection_params, VirusVariant.Wildtype, age_group_0_to_4,
                                              parameters.loc["Age0to4_InfectedNoSymptomsToRecovered"].value, parameters.loc["Age0to4_InfectedNoSymptomsToRecovered"].dev)
    abm.set_TimeInfectedSymptomsToRecovered(infection_params, VirusVariant.Wildtype, age_group_0_to_4,
                                            parameters.loc["Age0to4_InfectedSymptomsToRecovered"].value, parameters.loc["Age0to4_InfectedSymptomsToRecovered"].dev)
    abm.set_TimeInfectedSymptomsToSevere(infection_params, VirusVariant.Wildtype, age_group_0_to_4,
                                         parameters.loc["Age0to4_InfectedSymptomsToSevere"].value, parameters.loc["Age0to4_InfectedSymptomsToSevere"].dev)
    abm.set_TimeInfectedSevereToRecovered(infection_params, VirusVariant.Wildtype, age_group_0_to_4,
                                          parameters.loc["Age0to4_SevereToRecovered"].value, parameters.loc["Age0to4_SevereToRecovered"].dev)
    abm.set_TimeInfectedSevereToCritical(infection_params, VirusVariant.Wildtype, age_group_0_to_4,
                                         parameters.loc["Age0to4_SevereToCritical"].value, parameters.loc["Age0to4_SevereToCritical"].dev)
    abm.set_TimeInfectedCriticalToRecovered(infection_params, VirusVariant.Wildtype, age_group_0_to_4,
                                            parameters.loc["Age0to4_CriticalToRecovered"].value, parameters.loc["Age0to4_CriticalToRecovered"].dev)
    abm.set_TimeInfectedCriticalToDead(infection_params, VirusVariant.Wildtype, age_group_0_to_4,
                                       parameters.loc["Age0to4_CriticalToDead"].value, parameters.loc["Age0to4_CriticalToDead"].dev)
    abm.set_viral_load_parameters(infection_params, VirusVariant.Wildtype, age_group_0_to_4,
                                  parameters.loc["peak_max"].value, parameters.loc["peak_max"].value,
                                  parameters.loc["incline"].value, parameters.loc["incline"].value,
                                  parameters.loc["decline"].value, parameters.loc["decline"].value)
    abm.set_infectivity_parameters(
        infection_params, VirusVariant.Wildtype, age_group_0_to_4,
        parameters.loc["alpha"].value, parameters.loc["alpha"].value,
        parameters.loc["beta"].value, parameters.loc["beta"].value)
    infection_params.SymptomaticPerInfectedNoSymptoms[VirusVariant.Wildtype,
                                                      age_group_0_to_4] = parameters.loc["Age0to4_SymptomsPerInfectedNoSymptoms"].value
    infection_params.SeverePerInfectedSymptoms[VirusVariant.Wildtype,
                                               age_group_0_to_4] = parameters.loc["Age0to4_SeverePerInfectedSymptoms"].value
    infection_params.CriticalPerInfectedSevere[VirusVariant.Wildtype,
                                               age_group_0_to_4] = parameters.loc["Age0to4_CriticalPerInfectedSevere"].value
    infection_params.DeathsPerInfectedCritical[VirusVariant.Wildtype,
                                               age_group_0_to_4] = parameters.loc["Age0to4_DeathsPerInfectedCritical"].value

    # AgeGroup 5-14
    abm.set_incubationPeriod(
        infection_params, VirusVariant.Wildtype, age_group_5_to_14, parameters.loc["Age5to14_IncubationPeriod"].value, parameters.loc["Age5to14_IncubationPeriod"].dev)
    abm.set_TimeInfectedNoSymptomsToSymptoms(infection_params, VirusVariant.Wildtype, age_group_5_to_14,
                                             parameters.loc["Age5to14_InfectedNoSymptomsToSymptoms"].value, parameters.loc["Age5to14_InfectedNoSymptomsToSymptoms"].dev)
    abm.set_TimeInfectedNoSymptomsToRecovered(infection_params, VirusVariant.Wildtype, age_group_5_to_14,
                                              parameters.loc["Age5to14_InfectedNoSymptomsToRecovered"].value, parameters.loc["Age5to14_InfectedNoSymptomsToRecovered"].dev)
    abm.set_TimeInfectedSymptomsToRecovered(infection_params, VirusVariant.Wildtype, age_group_5_to_14,
                                            parameters.loc["Age5to14_InfectedSymptomsToRecovered"].value, parameters.loc["Age5to14_InfectedSymptomsToRecovered"].dev)
    abm.set_TimeInfectedSymptomsToSevere(infection_params, VirusVariant.Wildtype, age_group_5_to_14,
                                         parameters.loc["Age5to14_InfectedSymptomsToSevere"].value, parameters.loc["Age5to14_InfectedSymptomsToSevere"].dev)
    abm.set_TimeInfectedSevereToRecovered(infection_params, VirusVariant.Wildtype, age_group_5_to_14,
                                          parameters.loc["Age5to14_SevereToRecovered"].value, parameters.loc["Age5to14_SevereToRecovered"].dev)
    abm.set_TimeInfectedSevereToCritical(infection_params, VirusVariant.Wildtype, age_group_5_to_14,
                                         parameters.loc["Age5to14_SevereToCritical"].value, parameters.loc["Age5to14_SevereToCritical"].dev)
    abm.set_TimeInfectedCriticalToRecovered(infection_params, VirusVariant.Wildtype, age_group_5_to_14,
                                            parameters.loc["Age5to14_CriticalToRecovered"].value, parameters.loc["Age5to14_CriticalToRecovered"].dev)
    abm.set_TimeInfectedCriticalToDead(infection_params, VirusVariant.Wildtype, age_group_5_to_14,
                                       parameters.loc["Age5to14_CriticalToDead"].value, parameters.loc["Age5to14_CriticalToDead"].dev)
    abm.set_viral_load_parameters(infection_params, VirusVariant.Wildtype, age_group_5_to_14,
                                  parameters.loc["peak_max"].value, parameters.loc["peak_max"].value,
                                  parameters.loc["incline"].value, parameters.loc["incline"].value,
                                  parameters.loc["decline"].value, parameters.loc["decline"].value)
    abm.set_infectivity_parameters(
        infection_params, VirusVariant.Wildtype, age_group_5_to_14,
        parameters.loc["alpha"].value, parameters.loc["alpha"].value,
        parameters.loc["beta"].value, parameters.loc["beta"].value)
    infection_params.SymptomaticPerInfectedNoSymptoms[VirusVariant.Wildtype,
                                                      age_group_5_to_14] = parameters.loc["Age5to14_SymptomsPerInfectedNoSymptoms"].value
    infection_params.SeverePerInfectedSymptoms[VirusVariant.Wildtype,
                                               age_group_5_to_14] = parameters.loc["Age5to14_SeverePerInfectedSymptoms"].value
    infection_params.CriticalPerInfectedSevere[VirusVariant.Wildtype,
                                               age_group_5_to_14] = parameters.loc["Age5to14_CriticalPerInfectedSevere"].value
    infection_params.DeathsPerInfectedCritical[VirusVariant.Wildtype,
                                               age_group_5_to_14] = parameters.loc["Age5to14_DeathsPerInfectedCritical"].value

    # AgeGroup 15-34
    abm.set_incubationPeriod(
        infection_params, VirusVariant.Wildtype, age_group_15_to_34, parameters.loc["Age15to34_IncubationPeriod"].value, parameters.loc["Age15to34_IncubationPeriod"].dev)
    abm.set_TimeInfectedNoSymptomsToSymptoms(infection_params, VirusVariant.Wildtype, age_group_15_to_34,
                                             parameters.loc["Age15to34_InfectedNoSymptomsToSymptoms"].value, parameters.loc["Age15to34_InfectedNoSymptomsToSymptoms"].dev)
    abm.set_TimeInfectedNoSymptomsToRecovered(infection_params, VirusVariant.Wildtype, age_group_15_to_34,
                                              parameters.loc["Age15to34_InfectedNoSymptomsToRecovered"].value, parameters.loc["Age15to34_InfectedNoSymptomsToRecovered"].dev)
    abm.set_TimeInfectedSymptomsToRecovered(infection_params, VirusVariant.Wildtype, age_group_15_to_34,
                                            parameters.loc["Age15to34_InfectedSymptomsToRecovered"].value, parameters.loc["Age15to34_InfectedSymptomsToRecovered"].dev)
    abm.set_TimeInfectedSymptomsToSevere(infection_params, VirusVariant.Wildtype, age_group_15_to_34,
                                         parameters.loc["Age15to34_InfectedSymptomsToSevere"].value, parameters.loc["Age15to34_InfectedSymptomsToSevere"].dev)
    abm.set_TimeInfectedSevereToRecovered(infection_params, VirusVariant.Wildtype, age_group_15_to_34,
                                          parameters.loc["Age15to34_SevereToRecovered"].value, parameters.loc["Age15to34_SevereToRecovered"].dev)
    abm.set_TimeInfectedSevereToCritical(infection_params, VirusVariant.Wildtype, age_group_15_to_34,
                                         parameters.loc["Age15to34_SevereToCritical"].value, parameters.loc["Age15to34_SevereToCritical"].dev)
    abm.set_TimeInfectedCriticalToRecovered(infection_params, VirusVariant.Wildtype, age_group_15_to_34,
                                            parameters.loc["Age15to34_CriticalToRecovered"].value, parameters.loc["Age15to34_CriticalToRecovered"].dev)
    abm.set_TimeInfectedCriticalToDead(infection_params, VirusVariant.Wildtype, age_group_15_to_34,
                                       parameters.loc["Age15to34_CriticalToDead"].value, parameters.loc["Age15to34_CriticalToDead"].dev)
    abm.set_viral_load_parameters(infection_params, VirusVariant.Wildtype, age_group_15_to_34,
                                  parameters.loc["peak_max"].value, parameters.loc["peak_max"].value,
                                  parameters.loc["incline"].value, parameters.loc["incline"].value,
                                  parameters.loc["decline"].value, parameters.loc["decline"].value)
    abm.set_infectivity_parameters(
        infection_params, VirusVariant.Wildtype, age_group_15_to_34,
        parameters.loc["alpha"].value, parameters.loc["alpha"].value,
        parameters.loc["beta"].value, parameters.loc["beta"].value)
    infection_params.SymptomaticPerInfectedNoSymptoms[VirusVariant.Wildtype,
                                                      age_group_15_to_34] = parameters.loc["Age15to34_SymptomsPerInfectedNoSymptoms"].value
    infection_params.SeverePerInfectedSymptoms[VirusVariant.Wildtype,
                                               age_group_15_to_34] = parameters.loc["Age15to34_SeverePerInfectedSymptoms"].value
    infection_params.CriticalPerInfectedSevere[VirusVariant.Wildtype,
                                               age_group_15_to_34] = parameters.loc["Age15to34_CriticalPerInfectedSevere"].value
    infection_params.DeathsPerInfectedCritical[VirusVariant.Wildtype,
                                               age_group_15_to_34] = parameters.loc["Age15to34_DeathsPerInfectedCritical"].value

    # AgeGroup 35-59
    abm.set_incubationPeriod(
        infection_params, VirusVariant.Wildtype, age_group_35_to_59, parameters.loc["Age35to59_IncubationPeriod"].value, parameters.loc["Age35to59_IncubationPeriod"].dev)
    abm.set_TimeInfectedNoSymptomsToSymptoms(infection_params, VirusVariant.Wildtype, age_group_35_to_59,
                                             parameters.loc["Age35to59_InfectedNoSymptomsToSymptoms"].value, parameters.loc["Age35to59_InfectedNoSymptomsToSymptoms"].dev)
    abm.set_TimeInfectedNoSymptomsToRecovered(infection_params, VirusVariant.Wildtype, age_group_35_to_59,
                                              parameters.loc["Age35to59_InfectedNoSymptomsToRecovered"].value, parameters.loc["Age35to59_InfectedNoSymptomsToRecovered"].dev)
    abm.set_TimeInfectedSymptomsToRecovered(infection_params, VirusVariant.Wildtype, age_group_35_to_59,
                                            parameters.loc["Age35to59_InfectedSymptomsToRecovered"].value, parameters.loc["Age35to59_InfectedSymptomsToRecovered"].dev)
    abm.set_TimeInfectedSymptomsToSevere(infection_params, VirusVariant.Wildtype, age_group_35_to_59,
                                         parameters.loc["Age35to59_InfectedSymptomsToSevere"].value, parameters.loc["Age35to59_InfectedSymptomsToSevere"].dev)
    abm.set_TimeInfectedSevereToRecovered(infection_params, VirusVariant.Wildtype, age_group_35_to_59,
                                          parameters.loc["Age35to59_SevereToRecovered"].value, parameters.loc["Age35to59_SevereToRecovered"].dev)
    abm.set_TimeInfectedSevereToCritical(infection_params, VirusVariant.Wildtype, age_group_35_to_59,
                                         parameters.loc["Age35to59_SevereToCritical"].value, parameters.loc["Age35to59_SevereToCritical"].dev)
    abm.set_TimeInfectedCriticalToRecovered(infection_params, VirusVariant.Wildtype, age_group_35_to_59,
                                            parameters.loc["Age35to59_CriticalToRecovered"].value, parameters.loc["Age35to59_CriticalToRecovered"].dev)
    abm.set_TimeInfectedCriticalToDead(infection_params, VirusVariant.Wildtype, age_group_35_to_59,
                                       parameters.loc["Age35to59_CriticalToDead"].value, parameters.loc["Age35to59_CriticalToDead"].dev)
    abm.set_viral_load_parameters(infection_params, VirusVariant.Wildtype, age_group_35_to_59,
                                  parameters.loc["peak_max"].value, parameters.loc["peak_max"].value,
                                  parameters.loc["incline"].value, parameters.loc["incline"].value,
                                  parameters.loc["decline"].value, parameters.loc["decline"].value)
    abm.set_infectivity_parameters(
        infection_params, VirusVariant.Wildtype, age_group_35_to_59,
        parameters.loc["alpha"].value, parameters.loc["alpha"].value,
        parameters.loc["beta"].value, parameters.loc["beta"].value)
    infection_params.SymptomaticPerInfectedNoSymptoms[VirusVariant.Wildtype,
                                                      age_group_35_to_59] = parameters.loc["Age35to59_SymptomsPerInfectedNoSymptoms"].value
    infection_params.SeverePerInfectedSymptoms[VirusVariant.Wildtype,
                                               age_group_35_to_59] = parameters.loc["Age35to59_SeverePerInfectedSymptoms"].value
    infection_params.CriticalPerInfectedSevere[VirusVariant.Wildtype,
                                               age_group_35_to_59] = parameters.loc["Age35to59_CriticalPerInfectedSevere"].value
    infection_params.DeathsPerInfectedCritical[VirusVariant.Wildtype,
                                               age_group_35_to_59] = parameters.loc["Age35to59_DeathsPerInfectedCritical"].value

    # AgeGroup 60-79
    abm.set_incubationPeriod(
        infection_params, VirusVariant.Wildtype, age_group_60_to_79, parameters.loc["Age60to79_IncubationPeriod"].value, parameters.loc["Age60to79_IncubationPeriod"].dev)
    abm.set_TimeInfectedNoSymptomsToSymptoms(infection_params, VirusVariant.Wildtype, age_group_60_to_79,
                                             parameters.loc["Age60to79_InfectedNoSymptomsToSymptoms"].value, parameters.loc["Age60to79_InfectedNoSymptomsToSymptoms"].dev)
    abm.set_TimeInfectedNoSymptomsToRecovered(infection_params, VirusVariant.Wildtype, age_group_60_to_79,
                                              parameters.loc["Age60to79_InfectedNoSymptomsToRecovered"].value, parameters.loc["Age60to79_InfectedNoSymptomsToRecovered"].dev)
    abm.set_TimeInfectedSymptomsToRecovered(infection_params, VirusVariant.Wildtype, age_group_60_to_79,
                                            parameters.loc["Age60to79_InfectedSymptomsToRecovered"].value, parameters.loc["Age60to79_InfectedSymptomsToRecovered"].dev)
    abm.set_TimeInfectedSymptomsToSevere(infection_params, VirusVariant.Wildtype, age_group_60_to_79,
                                         parameters.loc["Age60to79_InfectedSymptomsToSevere"].value, parameters.loc["Age60to79_InfectedSymptomsToSevere"].dev)
    abm.set_TimeInfectedSevereToRecovered(infection_params, VirusVariant.Wildtype, age_group_60_to_79,
                                          parameters.loc["Age60to79_SevereToRecovered"].value, parameters.loc["Age60to79_SevereToRecovered"].dev)
    abm.set_TimeInfectedSevereToCritical(infection_params, VirusVariant.Wildtype, age_group_60_to_79,
                                         parameters.loc["Age60to79_SevereToCritical"].value, parameters.loc["Age60to79_SevereToCritical"].dev)
    abm.set_TimeInfectedCriticalToRecovered(infection_params, VirusVariant.Wildtype, age_group_60_to_79,
                                            parameters.loc["Age60to79_CriticalToRecovered"].value, parameters.loc["Age60to79_CriticalToRecovered"].dev)
    abm.set_TimeInfectedCriticalToDead(infection_params, VirusVariant.Wildtype, age_group_60_to_79,
                                       parameters.loc["Age60to79_CriticalToDead"].value, parameters.loc["Age60to79_CriticalToDead"].dev)
    abm.set_viral_load_parameters(infection_params, VirusVariant.Wildtype, age_group_60_to_79,
                                  parameters.loc["peak_max"].value, parameters.loc["peak_max"].value,
                                  parameters.loc["incline"].value, parameters.loc["incline"].value,
                                  parameters.loc["decline"].value, parameters.loc["decline"].value)
    abm.set_infectivity_parameters(
        infection_params, VirusVariant.Wildtype, age_group_60_to_79,
        parameters.loc["alpha"].value, parameters.loc["alpha"].value,
        parameters.loc["beta"].value, parameters.loc["beta"].value)
    infection_params.SymptomaticPerInfectedNoSymptoms[VirusVariant.Wildtype,
                                                      age_group_60_to_79] = parameters.loc["Age60to79_SymptomsPerInfectedNoSymptoms"].value
    infection_params.SeverePerInfectedSymptoms[VirusVariant.Wildtype,
                                               age_group_60_to_79] = parameters.loc["Age60to79_SeverePerInfectedSymptoms"].value
    infection_params.CriticalPerInfectedSevere[VirusVariant.Wildtype,
                                               age_group_60_to_79] = parameters.loc["Age60to79_CriticalPerInfectedSevere"].value
    infection_params.DeathsPerInfectedCritical[VirusVariant.Wildtype,
                                               age_group_60_to_79] = parameters.loc["Age60to79_DeathsPerInfectedCritical"].value

    # AgeGroup 80+
    abm.set_incubationPeriod(
        infection_params, VirusVariant.Wildtype, age_group_80_plus, parameters.loc["Age80plus_IncubationPeriod"].value, parameters.loc["Age80plus_IncubationPeriod"].dev)
    abm.set_TimeInfectedNoSymptomsToSymptoms(infection_params, VirusVariant.Wildtype, age_group_80_plus,
                                             parameters.loc["Age80plus_InfectedNoSymptomsToSymptoms"].value, parameters.loc["Age80plus_InfectedNoSymptomsToSymptoms"].dev)
    abm.set_TimeInfectedNoSymptomsToRecovered(infection_params, VirusVariant.Wildtype, age_group_80_plus,
                                              parameters.loc["Age80plus_InfectedNoSymptomsToRecovered"].value, parameters.loc["Age80plus_InfectedNoSymptomsToRecovered"].dev)
    abm.set_TimeInfectedSymptomsToRecovered(infection_params, VirusVariant.Wildtype, age_group_80_plus,
                                            parameters.loc["Age80plus_InfectedSymptomsToRecovered"].value, parameters.loc["Age80plus_InfectedSymptomsToRecovered"].dev)
    abm.set_TimeInfectedSymptomsToSevere(infection_params, VirusVariant.Wildtype, age_group_80_plus,
                                         parameters.loc["Age80plus_InfectedSymptomsToSevere"].value, parameters.loc["Age80plus_InfectedSymptomsToSevere"].dev)
    abm.set_TimeInfectedSevereToRecovered(infection_params, VirusVariant.Wildtype, age_group_80_plus,
                                          parameters.loc["Age80plus_SevereToRecovered"].value, parameters.loc["Age80plus_SevereToRecovered"].dev)
    abm.set_TimeInfectedSevereToCritical(infection_params, VirusVariant.Wildtype, age_group_80_plus,
                                         parameters.loc["Age80plus_SevereToCritical"].value, parameters.loc["Age80plus_SevereToCritical"].dev)
    abm.set_TimeInfectedCriticalToRecovered(infection_params, VirusVariant.Wildtype, age_group_80_plus,
                                            parameters.loc["Age80plus_CriticalToRecovered"].value, parameters.loc["Age80plus_CriticalToRecovered"].dev)
    abm.set_TimeInfectedCriticalToDead(infection_params, VirusVariant.Wildtype, age_group_80_plus,
                                       parameters.loc["Age80plus_CriticalToDead"].value, parameters.loc["Age80plus_CriticalToDead"].dev)
    abm.set_viral_load_parameters(infection_params, VirusVariant.Wildtype, age_group_80_plus,
                                  parameters.loc["peak_max"].value, parameters.loc["peak_max"].value,
                                  parameters.loc["incline"].value, parameters.loc["incline"].value,
                                  parameters.loc["decline"].value, parameters.loc["decline"].value)
    abm.set_infectivity_parameters(
        infection_params, VirusVariant.Wildtype, age_group_80_plus,
        parameters.loc["alpha"].value, parameters.loc["alpha"].value,
        parameters.loc["beta"].value, parameters.loc["beta"].value)
    infection_params.SymptomaticPerInfectedNoSymptoms[VirusVariant.Wildtype,
                                                      age_group_80_plus] = parameters.loc["Age80plus_SymptomsPerInfectedNoSymptoms"].value
    infection_params.SeverePerInfectedSymptoms[VirusVariant.Wildtype,
                                               age_group_80_plus] = parameters.loc["Age80plus_SeverePerInfectedSymptoms"].value
    infection_params.CriticalPerInfectedSevere[VirusVariant.Wildtype,
                                               age_group_80_plus] = parameters.loc["Age80plus_CriticalPerInfectedSevere"].value
    infection_params.DeathsPerInfectedCritical[VirusVariant.Wildtype,
                                               age_group_80_plus] = parameters.loc["Age80plus_DeathsPerInfectedCritical"].value

    return infection_params


def read_txt(path):
    # read input file and save it in a pd.Dataframe
    return pd.read_csv(path, sep='\t')


def insert_locations_to_map(mapping, inputId, locations):
    map = LocationMapping()
    map.inputId = inputId
    for loc in locations:
        type = str(int(loc[1]))
        index = str(loc[0].index())
        if int(loc[1]) < 10:
            type = "0" + type
        if loc[0].index() < 10:
            index = "0" + index
        map.modelId.append(type+index)
    mapping.append(map)
    return mapping


def assign_infection_states(model, t0, exposed_pct, infected_no_symptoms_pct, infected_symptoms_pct,
                            infected_severe_pct, infected_critical_pct, recovered_pct):
    susceptible_pct = 1 - exposed_pct - infected_no_symptoms_pct - \
        infected_symptoms_pct - infected_severe_pct - \
        infected_critical_pct - recovered_pct
    for person in model.persons:
        # draw infection state from distribution for every agent
        infection_state = np.random.choice(np.arange(0, int(abm.InfectionState.Count)),
                                           p=[susceptible_pct, exposed_pct, infected_no_symptoms_pct,
                                               infected_symptoms_pct, infected_severe_pct, infected_critical_pct, recovered_pct, 0.0])
        if (abm.InfectionState(infection_state) != abm.InfectionState.Susceptible):
            person.add_new_infection(Infection(
                model, person, VirusVariant.Wildtype, t0, abm.InfectionState(infection_state), False), t0)


def find_all_locations_of_type(model, type):
    locations = []
    for loc in model.locations:
        if (loc.type == type):
            locations.append(loc.id)
    return locations


def convert_loc_id_to_string(loc):
    type = str(int(loc[1]))
    index = str(loc[0].index())
    if int(loc[1]) < 10:
        type = "0" + type
    if int(loc[0].index()) < 10:
        index = "0" + index

    return type + index


def get_agents_per_location(loc_id, agents):
    agents_per_loc = []
    for a in agents:
        if (int(a[1]) == int(loc_id[1]) and a[0].index() == loc_id[0].index()):
            agents_per_loc.append(a)
    return agents_per_loc


def write_results_to_file(path, log):
    location_ids = log[1][0]
    time_points = log[0]
    agents = log[2]
    with open(path, 'w') as f:
        for location_id_index in range(len(location_ids)):
            line = convert_loc_id_to_string(
                location_ids[location_id_index]) + " " + str(len(time_points))
            for t in range(len(time_points)):
                agents_per_loc = get_agents_per_location(
                    location_ids[location_id_index], agents[t])
                line += " " + str(time_points[t]) + \
                    " " + str(len(agents_per_loc))
                for a in agents_per_loc:
                    if (a[3].days > 1000):
                        time_since_transmission = -1.0
                    else:
                        time_since_transmission = a[3].hours
                    line += " " + str(a[2].index()) + " " + \
                        str(time_since_transmission)
            f.write(line)
            f.write('\n')
    f.close()


def convert_infection_state_to_string(infection_state):
    if (infection_state == abm.InfectionState.Susceptible):
        return "S"
    elif (infection_state == abm.InfectionState.Exposed):
        return "E"
    elif (infection_state == abm.InfectionState.InfectedNoSymptoms):
        return "I_ns"
    elif (infection_state == abm.InfectionState.InfectedSymptoms):
        return "I_sy"
    elif (infection_state == abm.InfectionState.InfectedSevere):
        return "I_sev"
    elif (infection_state == abm.InfectionState.InfectedCritical):
        return "I_cri"
    elif (infection_state == abm.InfectionState.Recovered):
        return "R"
    elif (infection_state == abm.InfectionState.Dead):
        return "D"
    else:
        raise Exception("Infection state not found")


def write_infection_paths_to_file_states(path, log):
    agent_ids = [log[2][0][i][1] for i in range(len(log[2][0]))]
    with open(path, 'w') as f:
        for id in agent_ids:
            line = str(id) + " "
            for t in range(len(log[2])):
                line += convert_infection_state_to_string(
                    log[2][t][id][3]) + " "
            f.write(line)
            f.write('\n')
    f.close()


def write_infection_paths_to_file(path, model, tmax):
    with open(path, 'w') as f:
        f.write("Agent_id S E I_ns I_sy I_sev I_cri R D\n")
        for person in model.persons:
            line = str(int(person.id.index())) + " "
            if person.infection_state(tmax) == abm.InfectionState.Susceptible:
                line += str(int(tmax.hours)) + " "
                for i in range(int(abm.InfectionState.Count)-1):
                    line += "0 "
            else:
                time_S = max(
                    person.infection.get_infection_start() - abm.TimePoint(0), abm.TimeSpan(0))
                time_E = person.infection.get_time_in_state(
                    abm.InfectionState.Exposed)
                time_INS = person.infection.get_time_in_state(
                    abm.InfectionState.InfectedNoSymptoms)
                time_ISy = person.infection.get_time_in_state(
                    abm.InfectionState.InfectedSymptoms)
                time_ISev = person.infection.get_time_in_state(
                    abm.InfectionState.InfectedSevere)
                time_ICri = person.infection.get_time_in_state(
                    abm.InfectionState.InfectedCritical)
                time_R = abm.TimePoint(0)
                time_D = abm.TimePoint(0)
                if (person.infection_state(tmax) == abm.InfectionState.Recovered):
                    time_infected = time_E + time_INS + time_ISy + time_ISev + time_ICri
                    if (time_S.hours == 0):
                        time_R = tmax - \
                            (time_infected +
                             (person.infection.get_infection_start() - abm.TimePoint(0)))
                    else:
                        time_R = tmax - time_S - time_infected
                if (person.infection_state(tmax) == abm.InfectionState.Dead):
                    time_infected = time_E + time_INS + time_ISy + time_ISev + time_ICri
                    if (time_S.hours == 0):
                        time_R = tmax - \
                            (time_infected +
                             (person.infection.get_infection_start() - abm.TimePoint(0)))
                    else:
                        time_R = tmax - time_S - time_infected
                line += str(time_S.hours) + " " + str(time_E.hours) + " " + str(time_INS.hours) + " " + str(time_ISy.hours) + " " \
                    + str(time_ISev.hours) + " " + str(time_ICri.hours) + \
                    " " + str(time_R.hours) + " " + \
                    str(time_D.hours) + " "
            f.write(line)
            f.write('\n')
    f.close()


def write_location_mapping_to_file(path, mapping):
    with open(path, 'w') as f:
        for id in mapping:
            line = id.inputId + " "
            for modelId in id.modelId:
                line += modelId + " "
            f.write(line)
            f.write('\n')
        f.close()


def set_sim_result_at_start(sim):
    for location in sim.model.locations:
        result = sim.result.get_last_value()
        result += location.population.get_last_value()


def write_age_and_hh(model, path):
    with open(path, 'w') as f:
        for person in model.persons:
            line = str(person.id) + " " + age_group_to_string(person.age) + " " + \
                str(person.assigned_location(abm.LocationType.Home))
            f.write(line)
            f.write('\n')
        f.close()


def write_compartments_to_file(model, path, timepoints):
    with open(path, 'w') as f:
        f.write("t S E Ins Isy Isev Icri R D\n")
        for t in range(len(timepoints)):
            tp = abm.TimePoint(0) + abm.hours(t)
            line = str(timepoints[t]) + " "
            comps = np.zeros(int(abm.InfectionState.Count))
            for person in model.persons:
                state = person.infection_state(tp)
                comps[int(state)] += 1
            for c in comps:
                line += str(c) + " "
            f.write(line)
            f.write('\n')
        f.close()


def convert_time_since_transmission(time):
    if (time.days > 1000):
        return -1.0
    else:
        return time.hours


def write_results_to_h5(path, log):
    file = h5py.File(path, 'w')
    result_list = []
    for agent in log[3][0]:
        gr = file.create_group(str(agent.index()))
        loc_ids = np.array([convert_loc_id_to_string((log[2][t][agent.index()][0], log[2][t][agent.index()][1]))
                           for t in range(len(log[2]))], dtype=np.str_)
        time_since_transm = np.array([convert_time_since_transmission(
            log[2][t][agent.index()][3]) for t in range(len(log[2]))], dtype=np.float64)
        ds = gr.create_dataset('loc_ids', shape=len(
            loc_ids), dtype=h5py.string_dtype())
        ds[:] = loc_ids
        gr.create_dataset("time_since_transm",
                          data=time_since_transm, dtype=np.float64)
    file.close()


def age_to_age_group(age):
    if (age < 5):
        return age_group_0_to_4
    elif (age < 15):
        return age_group_5_to_14
    elif (age < 35):
        return age_group_15_to_34
    elif (age < 60):
        return age_group_35_to_59
    elif (age < 80):
        return age_group_60_to_79
    else:
        return age_group_80_plus


def initialize_model2(model, trip_file):
    start = time.time()
    # This map has the puid from the trip df as key and the ABM PersonIds as values
    person_map = {}
    # Map for the locations. The keys are the location types. For every location type we have:
    # - a dictionary with traffic zone ids as keys and as value
    #   - a dictionary with huid/location_id from trip data as key and the ABM LocationId as value
    location_map = {0: {}, 1: {}, 2: {}, 3: {}, 4: {}, 5: {}}

    # Load the data
    trip_df = pd.read_csv(trip_file)

    # Preprocess data to avoid repeated filtering
    # To create a person, the puid, the home location and a person's age is needed
    pids_df = trip_df[['puid', 'huid', 'age']].drop_duplicates()
    # Get the number of trips per location for each person sorted by traffic zone and location type
    locs = trip_df.groupby(['puid', 'location_type', 'end_zone'])[
        'loc_id_end'].value_counts().reset_index(name='count')
    # Get all shops in the df
    all_shops = trip_df[trip_df['location_type']
                        == 4][['end_zone', 'loc_id_end']].drop_duplicates()
    # Get all events in the df
    all_events = trip_df[trip_df['location_type']
                         == 3][['end_zone', 'loc_id_end']].drop_duplicates()
    # Get all work places in the df
    all_works = trip_df[trip_df['location_type']
                        == 2][['end_zone', 'loc_id_end']].drop_duplicates()
    # Get all schools in the df
    all_schools = trip_df[trip_df['location_type']
                          == 1][['end_zone', 'loc_id_end']].drop_duplicates()

    # Create fast lookup dictionaries for home zones
    # i.e. creates dictionary with loc_id as key and traffic zone id as value
    start_zone_lookup = trip_df.set_index(
        'loc_id_start')['start_zone'].to_dict()
    end_zone_lookup = trip_df.set_index('loc_id_end')['end_zone'].to_dict()
    for index, row in pids_df.iterrows():
        if row['home_in_munich'] == 0:
            continue
        # Find home zone
        home_zone = start_zone_lookup.get(
            row['huid']) or end_zone_lookup.get(row['huid'])
        if home_zone is None:
            # print("Error: Home zone was not found")
            continue

        # Ensure home_zone is in integer format
        # If it's a series take the first entry
        if isinstance(home_zone, pd.Series):
            home_zone = home_zone.iloc[0]

        # Add home location to the model
        if home_zone in location_map[0]:
            # Home zone is in the location map already
            home = location_map[0][home_zone].get(row['huid'])
            if home is None:
                # Home id is not in the location map/model
                home = model.add_location(abm.LocationType.Home)
                location_map[0][home_zone][row['huid']] = home
        else:
            # Home zone (and therefore also home id) is not in the location map/model
            home = model.add_location(abm.LocationType.Home)
            location_map[0][home_zone] = {row['huid']: home}

        # Add the person to the model
        p = model.add_person(home, age_to_age_group(row['age']))
        person_map[row['puid']] = p

        ### Assign shop ###
        # Get all shops the person visits sorted by number of trips to that shop
        shops = locs[(locs['puid'] == row['puid']) & (
            locs['location_type'] == 4)].sort_values(by='count', ascending=False)

        if shops.empty:
            # No trips to shops
            # Assign a shop in home zone
            if home_zone in location_map[4]:
                # Choose shop uniquely distributed from all shops in home zone
                shop = random.choice(list(location_map[4][home_zone].values()))
            else:
                # Get all shops in home zone
                shops_in_hz = all_shops[all_shops['end_zone'] == home_zone]
                if not shops_in_hz.empty:
                    # Shop in home zone exists
                    shop_id = random.choice(list(shops_in_hz['loc_id_end']))
                    shop = model.add_location(abm.LocationType.BasicsShop)
                    location_map[4][home_zone] = {shop_id: shop}
                else:
                    print('Info: No shop in home zone was found in the data.')
                    shop = model.add_location(abm.LocationType.BasicsShop)
                    location_map[4][home_zone] = {'xx': shop}
        else:
            # Use the person's most frequently visited shop
            shop_zone = shops.iloc[0]['end_zone']
            shop_id = shops.iloc[0]['loc_id_end']
            if shop_zone in location_map[4]:
                # Traffic zone of most frequently visited shop already exists in location map
                shop = location_map[4][shop_zone].get(shop_id)
                if shop is None:
                    # Shop id is not in the location map/model
                    shop = model.add_location(abm.LocationType.BasicsShop)
                    location_map[4][shop_zone][shop_id] = shop
            else:
                # Shop zone (and therefore also shop id) is not in the location map/model
                shop = model.add_location(abm.LocationType.BasicsShop)
                location_map[4][shop_zone] = {shop_id: shop}

        # Assign the shop to the person
        model.persons[p.index()].set_assigned_location(
            abm.LocationType.BasicsShop, shop)

        ### Assign event ###
        # Get all events the person visits sorted by number of trips to that event
        events = locs[(locs['puid'] == row['puid']) & (
            locs['location_type'] == 3)].sort_values(by='count', ascending=False)

        if events.empty:
            # No trips to events
            # Assign an event in home zone
            if home_zone in location_map[3]:
                # Choose event uniquely distributed from all events in home zone
                event = random.choice(
                    list(location_map[3][home_zone].values()))
            else:
                # Get all events in home zone
                events_in_hz = all_events[all_events['end_zone'] == home_zone]
                if not events_in_hz.empty:
                    # Event in home zone exists
                    event_id = random.choice(list(events_in_hz['loc_id_end']))
                    event = model.add_location(abm.LocationType.SocialEvent)
                    location_map[3][home_zone] = {event_id: event}
                else:
                    print('Info: No event in home zone was found in the data.')
                    event = model.add_location(abm.LocationType.SocialEvent)
                    location_map[3][home_zone] = {'xx': event}
        else:
            # Use the person's most frequently visited event
            event_zone = events.iloc[0]['end_zone']
            event_id = events.iloc[0]['loc_id_end']
            if event_zone in location_map[3]:
                # Traffic zone of most frequently visited event already exists in location map
                event = location_map[3][event_zone].get(event_id)
                if event is None:
                    # Event id is not in the location map/model
                    event = model.add_location(abm.LocationType.SocialEvent)
                    location_map[3][event_zone][event_id] = event
            else:
                # Event zone (and therefore also event id) is not in the location map/model
                event = model.add_location(abm.LocationType.SocialEvent)
                location_map[3][event_zone] = {event_id: event}

        # Assign the event to the person
        model.persons[p.index()].set_assigned_location(
            abm.LocationType.SocialEvent, event)

        ### Assign work to ages 15 to 59###
        if (model.persons[p.index()].age == age_group_15_to_34 or model.persons[p.index()].age == age_group_35_to_59):
            # Get all work places the person visits sorted by number of trips to that work place
            works = locs[(locs['puid'] == row['puid']) & (
                locs['location_type'] == 2)].sort_values(by='count', ascending=False)

            if works.empty:
                # No trips to work places
                # Assign a work place in home zone
                if home_zone in location_map[2]:
                    # Choose work place uniquely distributed from all work places in home zone
                    work = random.choice(
                        list(location_map[2][home_zone].values()))
                else:
                    # Get all work places in home zone
                    works_in_hz = all_works[all_works['end_zone'] == home_zone]
                    if not works_in_hz.empty:
                        # Work place in home zone exists
                        work_id = random.chioce(
                            list(works_in_hz['loc_id_end']))
                        work = model.add_location(abm.LocationType.Work)
                        location_map[2][home_zone] = {work_id: work}
                    else:
                        print(
                            'Info: No work place in home zone was found in the data.')
                        work = model.add_location(abm.LocationType.Work)
                        location_map[2][home_zone] = {'xx': work}
            else:
                # Use the person's most frequently visited work place
                work_zone = works.iloc[0]['end_zone']
                work_id = works.iloc[0]['loc_id_end']
                if work_zone in location_map[2]:
                    # Traffic zone of most frequently visited work place already exists in location map
                    work = location_map[2][work_zone].get(work_id)
                    if work is None:
                        # Work id is not in the location map/model
                        work = model.add_location(abm.LocationType.Work)
                        location_map[2][work_zone][work_id] = work
                else:
                    # Work zone (and therefore also work id) is not in the location map/model
                    work = model.add_location(abm.LocationType.Work)
                    location_map[2][work_zone] = {work_id: work}

            # Assign the work place to the person
            model.persons[p.index()].set_assigned_location(
                abm.LocationType.Work, work)

        ### Assign school to ages 5 to 14###
        if (model.persons[p.index()].age == age_group_5_to_14):
            # Get all schools the person visits sorted by number of trips to that school
            schools = locs[(locs['puid'] == row['puid']) & (
                locs['location_type'] == 1)].sort_values(by='count', ascending=False)

            if schools.empty:
                # No trips to schools
                # Assign a school in home zone
                if home_zone in location_map[1]:
                    # Choose school uniquely distributed from all schools in home zone
                    school = random.choice(
                        list(location_map[1][home_zone].values()))
                else:
                    # Get all schools in home zone
                    schools_in_hz = all_schools[all_schools['end_zone']
                                                == home_zone]
                    if not schools_in_hz.empty:
                        # School in home zone exists
                        school_id = random.choice(
                            list(schools_in_hz['loc_id_end']))
                        school = model.add_location(abm.LocationType.School)
                        location_map[1][home_zone] = {school_id: school}
                    else:
                        print(
                            'Info: No school in home zone was found in the data.')
                        school = model.add_location(abm.LocationType.School)
                        location_map[1][home_zone] = {'xx': school}
            else:
                # Use the person's most frequently visited school
                school_zone = schools.iloc[0]['end_zone']
                school_id = schools.iloc[0]['loc_id_end']
                if school_zone in location_map[1]:
                    # Traffic zone of most frequently visited school already exists in location map
                    school = location_map[1][school_zone].get(school_id)
                    if school is None:
                        # School id is not in the location map/model
                        school = model.add_location(abm.LocationType.School)
                        location_map[1][school_zone][school_id] = school
                else:
                    # School zone (and therefore also school id) is not in the location map/model
                    school = model.add_location(abm.LocationType.School)
                    location_map[1][school_zone] = {school_id: school}

            # Assign the work place to the person
            model.persons[p.index()].set_assigned_location(
                abm.LocationType.School, school)
    end = time.time()
    print(f'Time to initialize model: {end - start} seconds')
    print("Number of persons in model: ", len(model.persons))
    print('Initialization complete.')


def save_persons(model, trip_file):
    start = time.time()
    # Map for the locations. The keys are the location types. For every location type we have:
    # - a dictionary with traffic zone ids as keys and as value
    #   - a dictionary with huid/location_id from trip data as key and the ABM LocationId as value
    location_map = {0: {}, 1: {}, 2: {}, 3: {}, 4: {}, 5: {}}

    # Load the data
    trip_df = pd.read_csv(trip_file)

    # Preprocess data to avoid repeated filtering
    # To create a person, the puid, the home location and a person's age is needed
    pids_df = trip_df[['puid', 'huid', 'age',
                       'home_in_munich']].drop_duplicates()
    # Get the number of trips per location for each person sorted by traffic zone and location type
    locs = trip_df.groupby(['puid', 'location_type', 'end_zone'])[
        'loc_id_end'].value_counts().reset_index(name='count')
    # Get all shops in the df
    all_shops = trip_df[trip_df['location_type']
                        == 4][['end_zone', 'loc_id_end']].drop_duplicates()
    # Get all events in the df
    all_events = trip_df[trip_df['location_type']
                         == 3][['end_zone', 'loc_id_end']].drop_duplicates()
    # Get all work places in the df
    all_works = trip_df[trip_df['location_type']
                        == 2][['end_zone', 'loc_id_end']].drop_duplicates()
    # Get all schools in the df
    all_schools = trip_df[trip_df['location_type']
                          == 1][['end_zone', 'loc_id_end']].drop_duplicates()

    # Create fast lookup dictionaries for home zones
    # i.e. creates dictionary with loc_id as key and traffic zone id as value
    start_zone_lookup = trip_df.set_index(
        'loc_id_start')['start_zone'].to_dict()
    end_zone_lookup = trip_df.set_index('loc_id_end')['end_zone'].to_dict()
    persons = []
    i = 1
    for index, row in pids_df.iterrows():
        print(i, '/', len(pids_df))
        i += 1
        if row['home_in_munich'] == 0:
            continue
        home_id = row['huid']
        age = row['age']
        # Find home zone
        home_zone = start_zone_lookup.get(
            row['huid']) or end_zone_lookup.get(row['huid'])
        if home_zone is None:
            print("Error: Home zone was not found")
            continue

        # Ensure home_zone is in integer format
        # If it's a series take the first entry
        if isinstance(home_zone, pd.Series):
            home_zone = home_zone.iloc[0]

        ### Assign shop ###
        # Get all shops the person visits sorted by number of trips to that shop
        shops = locs[(locs['puid'] == row['puid']) & (
            locs['location_type'] == 4)].sort_values(by='count', ascending=False)

        if shops.empty:
            shop_zone = home_zone
            # No trips to shops
            # Assign a shop in home zone
            if home_zone in location_map[4]:
                # Choose shop uniquely distributed from all shops in home zone
                shop_id = random.choice(
                    list(location_map[4][home_zone].keys()))
            else:
                # Get all shops in home zone
                shops_in_hz = all_shops[all_shops['end_zone'] == home_zone]
                if not shops_in_hz.empty:
                    # Shop in home zone exists
                    shop_id = random.choice(list(shops_in_hz['loc_id_end']))
                    location_map[4][home_zone] = {shop_id: shop_id}
                else:
                    print('Info: No shop in home zone was found in the data.')
                    shop_id = -1
        else:
            # Use the person's most frequently visited shop
            shop_zone = shops.iloc[0]['end_zone']
            shop_id = shops.iloc[0]['loc_id_end']
            if shop_zone in location_map[4]:
                # Traffic zone of most frequently visited shop already exists in location map
                shop = location_map[4][shop_zone].get(shop_id)
                if shop is None:
                    # Shop id is not in the location map/model
                    location_map[4][shop_zone][shop_id] = shop_id
            else:
                # Shop zone (and therefore also shop id) is not in the location map/model
                location_map[4][shop_zone] = {shop_id: shop_id}

        ### Assign event ###
        # Get all events the person visits sorted by number of trips to that event
        events = locs[(locs['puid'] == row['puid']) & (
            locs['location_type'] == 3)].sort_values(by='count', ascending=False)

        if events.empty:
            event_zone = home_zone
            # No trips to events
            # Assign an event in home zone
            if home_zone in location_map[3]:
                # Choose event uniquely distributed from all events in home zone
                event_id = random.choice(
                    list(location_map[3][home_zone].keys()))
            else:
                # Get all events in home zone
                events_in_hz = all_events[all_events['end_zone'] == home_zone]
                if not events_in_hz.empty:
                    # Event in home zone exists
                    event_id = random.choice(list(events_in_hz['loc_id_end']))
                    location_map[3][home_zone] = {event_id: event_id}
                else:
                    print('Info: No event in home zone was found in the data.')
                    event_id = -1
        else:
            # Use the person's most frequently visited event
            event_zone = events.iloc[0]['end_zone']
            event_id = events.iloc[0]['loc_id_end']
            if event_zone in location_map[3]:
                # Traffic zone of most frequently visited event already exists in location map
                event = location_map[3][event_zone].get(event_id)
                if event is None:
                    # Event id is not in the location map/model
                    location_map[3][event_zone][event_id] = event_id
            else:
                # Event zone (and therefore also event id) is not in the location map/model
                location_map[3][event_zone] = {event_id: event_id}

        work_zone = -2
        work_id = -2
        ### Assign work to ages 15 to 59###
        if (row['age'] > 15 and row['age'] < 60):
            # Get all work places the person visits sorted by number of trips to that work place
            works = locs[(locs['puid'] == row['puid']) & (
                locs['location_type'] == 2)].sort_values(by='count', ascending=False)

            if works.empty:
                work_zone = home_zone
                # No trips to work places
                # Assign a work place in home zone
                if home_zone in location_map[2]:
                    # Choose work place uniquely distributed from all work places in home zone
                    work_id = random.choice(
                        list(location_map[2][home_zone].keys()))
                else:
                    # Get all work places in home zone
                    works_in_hz = all_works[all_works['end_zone'] == home_zone]
                    if not works_in_hz.empty:
                        # Work place in home zone exists
                        work_id = random.choice(
                            list(works_in_hz['loc_id_end']))
                        location_map[2][home_zone] = {work_id: work_id}
                    else:
                        print(
                            'Info: No work place in home zone was found in the data.')
                        work_id = -1
            else:
                # Use the person's most frequently visited work place
                work_zone = works.iloc[0]['end_zone']
                work_id = works.iloc[0]['loc_id_end']
                if work_zone in location_map[2]:
                    # Traffic zone of most frequently visited work place already exists in location map
                    work = location_map[2][work_zone].get(work_id)
                    if work is None:
                        # Work id is not in the location map/model
                        location_map[2][work_zone][work_id] = work_id
                else:
                    # Work zone (and therefore also work id) is not in the location map/model
                    location_map[2][work_zone] = {work_id: work_id}

        school_zone = -2
        school_id = -2
        ### Assign school to ages 5 to 14###
        if (row['age'] > 4 and row['age'] < 16):
            # Get all schools the person visits sorted by number of trips to that school
            schools = locs[(locs['puid'] == row['puid']) & (
                locs['location_type'] == 1)].sort_values(by='count', ascending=False)

            if schools.empty:
                school_zone = home_zone
                # No trips to schools
                # Assign a school in home zone
                if home_zone in location_map[1]:
                    # Choose school uniquely distributed from all schools in home zone
                    school_id = random.choice(
                        list(location_map[1][home_zone].keys()))
                else:
                    # Get all schools in home zone
                    schools_in_hz = all_schools[all_schools['end_zone']
                                                == home_zone]
                    if not schools_in_hz.empty:
                        # School in home zone exists
                        school_id = random.choice(
                            list(schools_in_hz['loc_id_end']))
                        location_map[1][home_zone] = {school_id: school_id}
                    else:
                        print(
                            'Info: No school in home zone was found in the data.')
                        school_id = -1
            else:
                # Use the person's most frequently visited school
                school_zone = schools.iloc[0]['end_zone']
                school_id = schools.iloc[0]['loc_id_end']
                if school_zone in location_map[1]:
                    # Traffic zone of most frequently visited school already exists in location map
                    school = location_map[1][school_zone].get(school_id)
                    if school is None:
                        # School id is not in the location map/model
                        location_map[1][school_zone][school_id] = school_id
                else:
                    # School zone (and therefore also school id) is not in the location map/model
                    location_map[1][school_zone] = {school_id: school_id}
        persons.append({'puid': row['puid'], 'age': age, 'home_zone': home_zone, 'home_id': home_id, 'shop_zone': shop_zone, 'shop_id': shop_id,
                        'event_zone': event_zone, 'event_id': event_id, 'work_zone': work_zone, 'work_id': work_id,
                        'school_zone': school_zone, 'school_id': school_id})
    end = time.time()
    print(f'Time to create person list: {end - start} seconds')
    print("Number of persons in model: ", len(persons))
    start = time.time()
    df = pd.DataFrame(persons)
    df.to_csv('persons.csv')
    end = time.time()
    print(f'Time to create and save person df: {end - start} seconds')


def run_abm_simulation(sim_num):
    mio.set_log_level(mio.LogLevel.Off)
    input_path = sys.path[0] + '/input/'
    output_path = sys.path[0] + '/output/'
    # set seed for fixed model initialization (locations and initial infection states)
    np.random.seed(sim_num)
    # starting time point
    t0 = abm.TimePoint(0)
    # end time point: simulation will run 14 days
    tmax = t0 + abm.days(14)
    # create simulation with starting timepoint and number of age groups
    sim = abm.Simulation(t0, num_age_groups)
    # initialize model
    save_persons(sim.model, os.path.join(
        input_path, 'muenchen_microscopic_result_modified.csv'))

    # # read infection parameters
    # parameters = pd.read_csv(os.path.join(
    #     input_path, 'parameter_table.csv'), index_col=0)
    # # set seeds for simulation
    # abm.set_seeds(sim.model, sim_num)
    # # set infection parameters
    # sim.model.parameters = set_infection_parameters(parameters)
    # # set age groups that go to school and work
    # abm.set_AgeGroupGoToSchool(sim.model.parameters, age_group_5_to_14)
    # abm.set_AgeGroupGoToWork(sim.model.parameters, age_group_15_to_34)
    # abm.set_AgeGroupGoToWork(sim.model.parameters, age_group_35_to_59)
    # # assign initial infection states according to distribution
    # assign_infection_states(sim.model, t0, 0.002, 0.005,
    #                         0.0029, 0.0001, 0.0, 0.0)
    # # output object
    # history = History()
    # # just used for debugging
    # # write_age_and_hh(sim.model, os.path.join(output_path, 'age_hh.txt'))
    # # advance simulation until tmax
    # sim.advance(tmax, history)
    # # results collected during the simulation
    # log = history.log
    # # write infection paths per agent to file
    # write_infection_paths_to_file(os.path.join(
    #     output_path, str(sim_num) + '_infection_paths.txt'), sim.model, tmax)
    # # write compartment size per time step to file
    # write_compartments_to_file(sim.model, os.path.join(
    #     output_path, str(sim_num) + '_comps.csv'), log[0])
    # start = time.time()
    # write_results_to_h5(os.path.join(
    #     output_path, str(sim_num) + '_output.h5'), log)
    # end = time.time()
    # print(f'Time to write output h5: {end - start} seconds')

    print('done')


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'abm demonstrator',
        description='Example demonstrating the agent-based model for a synthetic population.')
    args = arg_parser.parse_args()
    # set LogLevel
    mio.set_log_level(mio.LogLevel.Off)
    mio.abm.set_log_level_warn()
    for i in range(1, 2):
        run_abm_simulation(i, **args.__dict__)
