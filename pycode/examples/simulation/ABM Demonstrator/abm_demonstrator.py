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


def make_one_person_households(number_of_households):
    # one-person household member
    one_person_household_member = abm.HouseholdMember(num_age_groups)
    # set weights for household member
    one_person_household_member.set_age_weight(age_group_15_to_34, 4364)
    one_person_household_member.set_age_weight(age_group_35_to_59, 7283)
    one_person_household_member.set_age_weight(age_group_60_to_79, 4100)
    one_person_household_member.set_age_weight(age_group_80_plus, 1800)

    # create one-person household group
    household_group = abm.HouseholdGroup()
    for hh in range(int(number_of_households)):
        household = abm.Household()
        household.add_members(one_person_household_member, 1)
        household_group.add_households(household, 1)

    return household_group


def make_multiple_person_households(household_size,
                                    number_of_two_parent_households, number_of_one_parent_households,
                                    number_of_other_households):

    # members for multiple person households
    child = abm.HouseholdMember(num_age_groups)
    child.set_age_weight(age_group_0_to_4, 1)
    child.set_age_weight(age_group_5_to_14, 1)

    parent = abm.HouseholdMember(num_age_groups)
    parent.set_age_weight(age_group_15_to_34, 2)
    parent.set_age_weight(age_group_35_to_59, 2)
    parent.set_age_weight(age_group_60_to_79, 1)

    other = abm.HouseholdMember(num_age_groups)
    other.set_age_weight(age_group_0_to_4, 5000)
    other.set_age_weight(age_group_5_to_14, 6000)
    other.set_age_weight(age_group_15_to_34, 14943)
    other.set_age_weight(age_group_35_to_59, 22259)
    other.set_age_weight(age_group_60_to_79, 11998)
    other.set_age_weight(age_group_80_plus, 5038)

    household_group = abm.HouseholdGroup()
    # add two parent households
    hh_two_parents = abm.Household()
    hh_two_parents.add_members(child, household_size - 2)
    hh_two_parents.add_members(parent, 2)
    household_group.add_households(
        hh_two_parents, number_of_two_parent_households)
    # add one parent households
    hh_one_parent = abm.Household()
    hh_one_parent.add_members(child, household_size - 1)
    hh_one_parent.add_members(parent, 1)
    household_group.add_households(
        hh_one_parent, number_of_one_parent_households)
    # add other households
    hh_other = abm.Household()
    hh_other.add_members(parent, 1)
    hh_other.add_members(other, household_size - 1)
    household_group.add_households(hh_other, number_of_other_households)

    return household_group


def add_households(world, distribution, num_inhabitants):
    household_sizes = np.zeros(len(distribution))
    locationIds = []
    new_index = len(world.locations)
    # distribute inhabintants to households
    while num_inhabitants > 0:
        size = np.random.choice(
            np.arange(0, len(distribution)), p=distribution)
        household_sizes[size] += 1
        num_inhabitants -= (size + 1)

    # one-person household
    one_person_household_group = make_one_person_households(household_sizes[0])
    abm.add_household_group_to_world(world, one_person_household_group)

    # two-person households
    two_person_two_parents = int(0.85 * household_sizes[1])
    two_person_one_parent = int(0.13 * household_sizes[1])
    two_person_other = int(
        household_sizes[1] - two_person_two_parents - two_person_one_parent)
    two_person_household_group = make_multiple_person_households(2, two_person_two_parents,
                                                                 two_person_one_parent, two_person_other)
    abm.add_household_group_to_world(world, two_person_household_group)

    # three-person households
    three_person_two_parents = int(0.83 * household_sizes[2])
    three_person_one_parent = int(0.13 * household_sizes[2])
    three_person_other = int(
        household_sizes[2] - three_person_two_parents - three_person_one_parent)
    three_person_household_group = make_multiple_person_households(3, three_person_two_parents,
                                                                   three_person_one_parent, three_person_other)
    abm.add_household_group_to_world(world, three_person_household_group)

    # four-person households
    four_person_two_parents = int(0.93 * household_sizes[3])
    four_person_one_parent = int(0.03 * household_sizes[3])
    four_person_other = int(
        household_sizes[3] - four_person_two_parents - four_person_one_parent)
    four_person_household_group = make_multiple_person_households(4, four_person_two_parents,
                                                                  four_person_one_parent, four_person_other)
    abm.add_household_group_to_world(world, four_person_household_group)

    # five-person households
    five_person_two_parents = int(0.88 * household_sizes[4])
    five_person_one_parent = int(0.05 * household_sizes[4])
    five_person_other = int(
        household_sizes[4] - five_person_two_parents - five_person_one_parent)
    five_person_household_group = make_multiple_person_households(5, five_person_two_parents,
                                                                  five_person_one_parent, five_person_other)
    abm.add_household_group_to_world(world, five_person_household_group)

    for hh in range(int(sum(household_sizes))):
        id = abm.LocationId(new_index, abm.LocationType.Home)
        locationIds.append(id)
        new_index += 1
    return locationIds


def insert_locations_to_map(mapping, inputId, locationIds):
    map = LocationMapping()
    map.inputId = inputId
    for id in locationIds:
        type = str(int(id.type))
        index = str(id.index)
        if int(id.type) < 10:
            type = "0" + type
        if int(id.index) < 10:
            index = "0" + index
        map.modelId.append(type+index)
    mapping[inputId] = map.modelId
    return mapping


def create_locations_from_input(world, input_areas, household_distribution):
    # map input area ids to corresponding abm location ids
    mapping = {}
    # bools to make sure the world has a school and a hospital
    has_school = False
    has_hospital = False
    for index, area in input_areas.iterrows():
        locationIds = []
        if ('residential' in area.type):
            # area 'residential' corresponds to location type 'Home'
            locationIds = add_households(
                world, household_distribution, area.inhabitants)
        elif (area.type == 'recreational'):
            # area 'recreational' corresponds to location type 'SocialEvent'
            location = world.add_location(abm.LocationType.SocialEvent)
            # set maximum contacts and capacity for social events
            world.get_individualized_location(
                location).infection_parameters.MaximumContacts = 30.
            world.get_individualized_location(location).set_capacity(30, 40, 0)
            locationIds.append(location)
        elif (area.type == 'shopping_business'):
            # area 'shopping_business' corresponds to location type School, Hospital, BasicsShops, Work
            if (not has_school):
                # if world does not have a school yet, a school is added
                location = world.add_location(abm.LocationType.School)
                # set maximum contacts and capacity for school
                world.get_individualized_location(
                    location).infection_parameters.MaximumContacts = 40.
                world.get_individualized_location(
                    location).set_capacity(500, 2000, 0)
                locationIds.append(location)
                has_school = True
            elif (not has_hospital):
                # if world does not have a hospital yet, a hospital and a icu is added
                locHosp = world.add_location(abm.LocationType.Hospital)
                # set maximum contacts and capacity for hospital
                world.get_individualized_location(
                    locHosp).infection_parameters.MaximumContacts = 5.
                world.get_individualized_location(
                    locHosp).set_capacity(300, 10000, 0)
                locICU = world.add_location(abm.LocationType.ICU)
                # set maximum contacts and capacity for icu
                world.get_individualized_location(
                    locICU).infection_parameters.MaximumContacts = 5.
                world.get_individualized_location(
                    locICU).set_capacity(30, 1000, 0)
                locationIds.append(locHosp)
                locationIds.append(locICU)
                has_hospital = True
            else:
                # when hospital and school has been added, the area 'shopping_business' is either
                # transformed to location type BasicsShop or location type Work with same probability
                type = np.random.choice(np.arange(0, 2), p=[0.5, 0.5])
                if (type):
                    location = world.add_location(abm.LocationType.BasicsShop)
                    # set maximum contacts and capacity for basics shops
                    world.get_individualized_location(
                        location).infection_parameters.MaximumContacts = 20.
                    world.get_individualized_location(
                        location).set_capacity(100, 1000, 0)
                    locationIds.append(location)
                else:
                    location = world.add_location(abm.LocationType.Work)
                    # set maximum contacts and capacity for work
                    world.get_individualized_location(
                        location).infection_parameters.MaximumContacts = 40.
                    world.get_individualized_location(
                        location).set_capacity(300, 2000, 0)
                    locationIds.append(location)
        elif (area.type == 'university'):
            # area 'university' corresponds to location type 'Work'
            location = world.add_location(abm.LocationType.Work)
            # set maximum contacts and capacity for work
            world.get_individualized_location(
                location).infection_parameters.MaximumContacts = 50.
            world.get_individualized_location(
                location).set_capacity(200, 4000, 0)
            locationIds.append(location)
        elif (area.type == 'mixed'):
            # area 'mixed' corresponds either to location type 'Work' or location type 'Home' with same probability
            type = np.random.choice(np.arange(0, 2), p=[0.5, 0.5])
            if (type):
                location = world.add_location(abm.LocationType.Work)
                # set maximum contacts and capacity for work
                world.get_individualized_location(
                    location).infection_parameters.MaximumContacts = 40.
                world.get_individualized_location(
                    location).set_capacity(100, 2000, 0)
                locationIds.append(location)
            else:
                locationIds = add_households(
                    world, household_distribution, area.inhabitants)
        insert_locations_to_map(mapping, str(area.id), locationIds)
    return mapping


def assign_infection_states(world, t0, exposed_pct, infected_no_symptoms_pct, infected_symptoms_pct,
                            infected_severe_pct, infected_critical_pct, recovered_pct, loc_list=[]):
    susceptible_pct = 1.0 - exposed_pct - infected_no_symptoms_pct - \
        infected_symptoms_pct - infected_severe_pct - \
        infected_critical_pct - recovered_pct
    # draw infection state from distribution for every agent
    for person in world.persons:
        if (len(loc_list) > 0):
            if (person.assigned_location(abm.LocationType.Home) in [int(x) for x in loc_list]):
                infection_state = np.random.choice(np.arange(0, int(abm.InfectionState.Count)),
                                                   p=[susceptible_pct, exposed_pct, infected_no_symptoms_pct,
                                                      infected_symptoms_pct, infected_severe_pct, infected_critical_pct, recovered_pct, 0.0])
                if (abm.InfectionState(infection_state) != abm.InfectionState.Susceptible):
                    person.add_new_infection(Infection(
                        world, person, VirusVariant.Wildtype, t0, abm.InfectionState(infection_state), False), t0)
        else:
            infection_state = np.random.choice(np.arange(0, int(abm.InfectionState.Count)),
                                               p=[susceptible_pct, exposed_pct, infected_no_symptoms_pct,
                                                  infected_symptoms_pct, infected_severe_pct, infected_critical_pct, recovered_pct, 0.0])
            if (abm.InfectionState(infection_state) != abm.InfectionState.Susceptible):
                person.add_new_infection(Infection(
                    world, person, VirusVariant.Wildtype, t0, abm.InfectionState(infection_state), False), t0)


def find_all_locations_of_type(world, type):
    locations = []
    for loc in world.locations:
        if (loc.type == type):
            locations.append(abm.LocationId(loc.index, type))
    return locations


def assign_locations(world):
    # get locations from world
    schools = find_all_locations_of_type(world, abm.LocationType.School)
    school_weights = [(1/len(schools)) for i in range(len(schools))]
    hospitals = find_all_locations_of_type(world, abm.LocationType.Hospital)
    hospital_weights = [(1/len(hospitals)) for i in range(len(hospitals))]
    icus = find_all_locations_of_type(world, abm.LocationType.ICU)
    icu_weights = [(1/len(icus)) for i in range(len(icus))]
    workplaces = find_all_locations_of_type(world, abm.LocationType.Work)
    workplace_weights = [(1/len(workplaces)) for i in range(len(workplaces))]
    basic_shops = find_all_locations_of_type(
        world, abm.LocationType.BasicsShop)
    shop_weights = [(1/len(basic_shops)) for i in range(len(basic_shops))]
    social_events = find_all_locations_of_type(
        world, abm.LocationType.SocialEvent)
    event_weights = [(1/len(social_events)) for i in range(len(social_events))]

    # assign locations to agents
    for person in world.persons:
        shop = np.random.choice(np.arange(0, len(basic_shops)), p=shop_weights)
        person.set_assigned_location(abm.LocationId(
            basic_shops[int(shop)].index, basic_shops[int(shop)].type))

        hospital = np.random.choice(
            np.arange(0, len(hospitals)), p=hospital_weights)
        person.set_assigned_location(abm.LocationId(
            hospitals[int(hospital)].index, hospitals[int(hospital)].type))

        icu = np.random.choice(np.arange(0, len(icus)), p=icu_weights)
        person.set_assigned_location(abm.LocationId(
            icus[int(icu)].index, icus[int(icu)].type))

        event = np.random.choice(
            np.arange(0, len(social_events)), p=event_weights)
        person.set_assigned_location(abm.LocationId(
            social_events[int(event)].index, social_events[int(event)].type))

        # assign school to agents between 5 and 14 years
        if (person.age == age_group_5_to_14):
            school = np.random.choice(
                np.arange(0, len(schools)), p=school_weights)
            person.set_assigned_location(abm.LocationId(
                schools[int(school)].index, schools[int(school)].type))
        # assign work to agents between 15 and 59
        if (person.age == age_group_15_to_34 or person.age == age_group_35_to_59):
            work = np.random.choice(
                np.arange(0, len(workplaces)), p=workplace_weights)
            person.set_assigned_location(abm.LocationId(
                workplaces[int(work)].index, workplaces[int(work)].type))


def convert_loc_id_to_string(loc_id):
    type = str(int(loc_id[0]))
    index = str(loc_id[1])
    if int(int(loc_id[0])) < 10:
        type = "0" + type
    if int(loc_id[1]) < 10:
        index = "0" + index

    return type + index


def get_agents_per_location(loc_id, agents):
    agents_per_loc = []
    for a in agents:
        if (int(a[0].type) == int(loc_id[0]) and a[0].index == loc_id[1]):
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
                    if (a[2].days > 1000):
                        time_since_transmission = -1.0
                    else:
                        time_since_transmission = a[2].hours
                    line += " " + str(a[1]) + " " + \
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


def write_infection_paths_to_file(path, world, tmax):
    with open(path, 'w') as f:
        f.write("Agent_id S E I_ns I_sy I_sev I_cri R D\n")
        for person in world.persons:
            line = str(int(person.id)) + " "
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
        for id in mapping.keys():
            line = id + " "
            for modelId in mapping[id]:
                line += modelId + " "
            f.write(line)
            f.write('\n')
        f.close()


def set_sim_result_at_start(sim):
    for location in sim.world.locations:
        result = sim.result.get_last_value()
        result += location.population.get_last_value()


def write_age_and_hh(world, path):
    with open(path, 'w') as f:
        for person in world.persons:
            line = str(person.id) + " " + age_group_to_string(person.age) + " " + \
                str(person.assigned_location(abm.LocationType.Home))
            f.write(line)
            f.write('\n')
        f.close()


def write_compartments_to_file(world, path, timepoints):
    with open(path, 'w') as f:
        f.write("t S E Ins Isy Isev Icri R D\n")
        for t in range(len(timepoints)):
            tp = abm.TimePoint(0) + abm.hours(t)
            line = str(timepoints[t]) + " "
            comps = np.zeros(int(abm.InfectionState.Count))
            for person in world.persons:
                state = person.infection_state(tp)
                comps[int(state)] += 1
            for c in comps:
                line += str(c) + " "
            f.write(line)
            f.write('\n')
        f.close()


def run_abm_simulation(sim_num):
    mio.set_log_level(mio.LogLevel.Warning)
    input_path = sys.path[0] + '/input/'
    output_path = sys.path[0] + '/output/'
    local_initial_outbreak = True
    # set seed for fixed model initialization (locations and initial infection states)
    np.random.seed(sim_num)
    # starting time point
    t0 = abm.TimePoint(0)
    # end time point: simulation will run 14 days
    tmax = t0 + abm.days(14)
    # distribution used to distribute the residential areas to one-, two-, three-, four- and five-person households
    household_distribution = [0.4084, 0.3375, 0.1199, 0.0965, 0.0377]
    # read txt file with inputs
    areas = read_txt(os.path.join(
        input_path, 'INSIDe_Demonstrator_AreaList.txt'))
    parameters = pd.read_csv(os.path.join(
        input_path, 'parameter_table.csv'), index_col=0)
    # create simulation with starting timepoint and number of age groups
    sim = abm.Simulation(t0, num_age_groups)
    # set seeds for simulation
    abm.set_seeds(sim.world, sim_num)
    # set infection parameters
    sim.world.parameters = set_infection_parameters(parameters)
    # as input areas do not fit one-to-one to abm location types, there has to be a mapping
    mapping = create_locations_from_input(
        sim.world, areas, household_distribution)
    # write location mapping to txt file
    write_location_mapping_to_file(
        os.path.join(output_path, str(sim_num) + '_location_mapping.txt'), mapping)

    start_locs = []
    E_pct = 0.002
    INS_pct = 0.005
    ISY_pct = 0.0029
    ISEV_pct = 0.0001
    if (local_initial_outbreak):
        E_pct = 0.03
        INS_pct = 0.03
        ISY_pct = 0.039
        ISEV_pct = 0.001
        start_locs = ["83", "84", "85"]

    # assign initial infection states according to distribution
    assign_infection_states(sim.world, t0, E_pct, INS_pct,
                            ISY_pct, ISEV_pct, 0.0, 0.0, [j for i in start_locs for j in mapping[i]])
    # assign locations to agents
    assign_locations(sim.world)
    # output object
    history = History()
    # just used for debugging
    # write_age_and_hh(sim.world, os.path.join(output_path, 'age_hh.txt'))
    # advance simulation until tmax
    sim.advance(tmax, history)
    # results collected during the simulation
    log = history.log
    # write infection paths per agent to file
    write_infection_paths_to_file(os.path.join(
        output_path, str(sim_num) + '_infection_paths.txt'), sim.world, tmax)
    # write compartment size per time step to file
    write_compartments_to_file(sim.world, os.path.join(
        output_path, str(sim_num) + '_comps.csv'), log[0])
    # write simulation results to txt file
    write_results_to_file(os.path.join(
        output_path, str(sim_num) + '_output.txt'), log)

    print('done')


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'abm demonstrator',
        description='Example demonstrating the agent-based model for a synthetic population.')
    args = arg_parser.parse_args()
    # set LogLevel
    mio.set_log_level(mio.LogLevel.Warning)
    for i in range(100):
        run_abm_simulation(i, **args.__dict__)
