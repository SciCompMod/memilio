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
import time
import random
import geopandas
import json
from multiprocessing import Pool
import copy

import memilio.simulation as mio
from memilio.simulation import abm
from memilio.simulation import AgeGroup
from memilio.simulation.abm import VirusVariant
from memilio.simulation.abm import History_sensitivity
from memilio.simulation.abm import Infection

import pandas as pd

# number of age groups
num_age_groups = 6
age_group_0_to_4 = AgeGroup(0)
age_group_5_to_15 = AgeGroup(1)
age_group_16_to_34 = AgeGroup(2)
age_group_35_to_59 = AgeGroup(3)
age_group_60_to_79 = AgeGroup(4)
age_group_80_plus = AgeGroup(5)


def set_infection_parameters(parameters, kappa):

    infection_params = abm.Parameters(num_age_groups)

    infection_params.InfectionRateFromViralShed[VirusVariant.Wildtype] = kappa

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
        infection_params, VirusVariant.Wildtype, age_group_5_to_15, parameters.loc["Age5to14_IncubationPeriod"].value, parameters.loc["Age5to14_IncubationPeriod"].dev)
    abm.set_TimeInfectedNoSymptomsToSymptoms(infection_params, VirusVariant.Wildtype, age_group_5_to_15,
                                             parameters.loc["Age5to14_InfectedNoSymptomsToSymptoms"].value, parameters.loc["Age5to14_InfectedNoSymptomsToSymptoms"].dev)
    abm.set_TimeInfectedNoSymptomsToRecovered(infection_params, VirusVariant.Wildtype, age_group_5_to_15,
                                              parameters.loc["Age5to14_InfectedNoSymptomsToRecovered"].value, parameters.loc["Age5to14_InfectedNoSymptomsToRecovered"].dev)
    abm.set_TimeInfectedSymptomsToRecovered(infection_params, VirusVariant.Wildtype, age_group_5_to_15,
                                            parameters.loc["Age5to14_InfectedSymptomsToRecovered"].value, parameters.loc["Age5to14_InfectedSymptomsToRecovered"].dev)
    abm.set_TimeInfectedSymptomsToSevere(infection_params, VirusVariant.Wildtype, age_group_5_to_15,
                                         parameters.loc["Age5to14_InfectedSymptomsToSevere"].value, parameters.loc["Age5to14_InfectedSymptomsToSevere"].dev)
    abm.set_TimeInfectedSevereToRecovered(infection_params, VirusVariant.Wildtype, age_group_5_to_15,
                                          parameters.loc["Age5to14_SevereToRecovered"].value, parameters.loc["Age5to14_SevereToRecovered"].dev)
    abm.set_TimeInfectedSevereToCritical(infection_params, VirusVariant.Wildtype, age_group_5_to_15,
                                         parameters.loc["Age5to14_SevereToCritical"].value, parameters.loc["Age5to14_SevereToCritical"].dev)
    abm.set_TimeInfectedCriticalToRecovered(infection_params, VirusVariant.Wildtype, age_group_5_to_15,
                                            parameters.loc["Age5to14_CriticalToRecovered"].value, parameters.loc["Age5to14_CriticalToRecovered"].dev)
    abm.set_TimeInfectedCriticalToDead(infection_params, VirusVariant.Wildtype, age_group_5_to_15,
                                       parameters.loc["Age5to14_CriticalToDead"].value, parameters.loc["Age5to14_CriticalToDead"].dev)
    abm.set_viral_load_parameters(infection_params, VirusVariant.Wildtype, age_group_5_to_15,
                                  parameters.loc["peak_max"].value, parameters.loc["peak_max"].value,
                                  parameters.loc["incline"].value, parameters.loc["incline"].value,
                                  parameters.loc["decline"].value, parameters.loc["decline"].value)
    abm.set_infectivity_parameters(
        infection_params, VirusVariant.Wildtype, age_group_5_to_15,
        parameters.loc["alpha"].value, parameters.loc["alpha"].value,
        parameters.loc["beta"].value, parameters.loc["beta"].value)
    infection_params.SymptomaticPerInfectedNoSymptoms[VirusVariant.Wildtype,
                                                      age_group_5_to_15] = parameters.loc["Age5to14_SymptomsPerInfectedNoSymptoms"].value
    infection_params.SeverePerInfectedSymptoms[VirusVariant.Wildtype,
                                               age_group_5_to_15] = parameters.loc["Age5to14_SeverePerInfectedSymptoms"].value
    infection_params.CriticalPerInfectedSevere[VirusVariant.Wildtype,
                                               age_group_5_to_15] = parameters.loc["Age5to14_CriticalPerInfectedSevere"].value
    infection_params.DeathsPerInfectedCritical[VirusVariant.Wildtype,
                                               age_group_5_to_15] = parameters.loc["Age5to14_DeathsPerInfectedCritical"].value

    # AgeGroup 15-34
    abm.set_incubationPeriod(
        infection_params, VirusVariant.Wildtype, age_group_16_to_34, parameters.loc["Age15to34_IncubationPeriod"].value, parameters.loc["Age15to34_IncubationPeriod"].dev)
    abm.set_TimeInfectedNoSymptomsToSymptoms(infection_params, VirusVariant.Wildtype, age_group_16_to_34,
                                             parameters.loc["Age15to34_InfectedNoSymptomsToSymptoms"].value, parameters.loc["Age15to34_InfectedNoSymptomsToSymptoms"].dev)
    abm.set_TimeInfectedNoSymptomsToRecovered(infection_params, VirusVariant.Wildtype, age_group_16_to_34,
                                              parameters.loc["Age15to34_InfectedNoSymptomsToRecovered"].value, parameters.loc["Age15to34_InfectedNoSymptomsToRecovered"].dev)
    abm.set_TimeInfectedSymptomsToRecovered(infection_params, VirusVariant.Wildtype, age_group_16_to_34,
                                            parameters.loc["Age15to34_InfectedSymptomsToRecovered"].value, parameters.loc["Age15to34_InfectedSymptomsToRecovered"].dev)
    abm.set_TimeInfectedSymptomsToSevere(infection_params, VirusVariant.Wildtype, age_group_16_to_34,
                                         parameters.loc["Age15to34_InfectedSymptomsToSevere"].value, parameters.loc["Age15to34_InfectedSymptomsToSevere"].dev)
    abm.set_TimeInfectedSevereToRecovered(infection_params, VirusVariant.Wildtype, age_group_16_to_34,
                                          parameters.loc["Age15to34_SevereToRecovered"].value, parameters.loc["Age15to34_SevereToRecovered"].dev)
    abm.set_TimeInfectedSevereToCritical(infection_params, VirusVariant.Wildtype, age_group_16_to_34,
                                         parameters.loc["Age15to34_SevereToCritical"].value, parameters.loc["Age15to34_SevereToCritical"].dev)
    abm.set_TimeInfectedCriticalToRecovered(infection_params, VirusVariant.Wildtype, age_group_16_to_34,
                                            parameters.loc["Age15to34_CriticalToRecovered"].value, parameters.loc["Age15to34_CriticalToRecovered"].dev)
    abm.set_TimeInfectedCriticalToDead(infection_params, VirusVariant.Wildtype, age_group_16_to_34,
                                       parameters.loc["Age15to34_CriticalToDead"].value, parameters.loc["Age15to34_CriticalToDead"].dev)
    abm.set_viral_load_parameters(infection_params, VirusVariant.Wildtype, age_group_16_to_34,
                                  parameters.loc["peak_max"].value, parameters.loc["peak_max"].value,
                                  parameters.loc["incline"].value, parameters.loc["incline"].value,
                                  parameters.loc["decline"].value, parameters.loc["decline"].value)
    abm.set_infectivity_parameters(
        infection_params, VirusVariant.Wildtype, age_group_16_to_34,
        parameters.loc["alpha"].value, parameters.loc["alpha"].value,
        parameters.loc["beta"].value, parameters.loc["beta"].value)
    infection_params.SymptomaticPerInfectedNoSymptoms[VirusVariant.Wildtype,
                                                      age_group_16_to_34] = parameters.loc["Age15to34_SymptomsPerInfectedNoSymptoms"].value
    infection_params.SeverePerInfectedSymptoms[VirusVariant.Wildtype,
                                               age_group_16_to_34] = parameters.loc["Age15to34_SeverePerInfectedSymptoms"].value
    infection_params.CriticalPerInfectedSevere[VirusVariant.Wildtype,
                                               age_group_16_to_34] = parameters.loc["Age15to34_CriticalPerInfectedSevere"].value
    infection_params.DeathsPerInfectedCritical[VirusVariant.Wildtype,
                                               age_group_16_to_34] = parameters.loc["Age15to34_DeathsPerInfectedCritical"].value

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


def get_person_age_string(age):
    if age == age_group_0_to_4:
        return "Age0to4"
    elif age == age_group_5_to_15:
        return "Age5to14"
    elif age == age_group_16_to_34:
        return "Age15to34"
    elif age == age_group_35_to_59:
        return "Age35to59"
    elif age == age_group_60_to_79:
        return "Age60to79"
    elif age == age_group_80_plus:
        return "Age80plus"
    else:
        print("Error: Age group cannot be found.")
        return " "


def assign_infection_states(model, t0, exposed_pct, infected_no_symptoms_pct, infected_symptoms_pct,
                            infected_severe_pct, infected_critical_pct, recovered_pct, parameters, loc_list):
    susceptible_pct = 1 - exposed_pct - infected_no_symptoms_pct - \
        infected_symptoms_pct - infected_severe_pct - \
        infected_critical_pct - recovered_pct
    for person in model.persons:
        # draw infection state from distribution for every agent
        if (len(loc_list) > 0):
            if (not (int(person.assigned_location(abm.LocationType.Home).index()) in [int(x[2:]) for x in loc_list])):
                continue
        infection_state = np.random.choice(np.arange(0, int(abm.InfectionState.Count)),
                                           p=[susceptible_pct, exposed_pct, infected_no_symptoms_pct,
                                               infected_symptoms_pct, infected_severe_pct, infected_critical_pct, recovered_pct, 0.0])
        if (abm.InfectionState(infection_state) != abm.InfectionState.Susceptible):
            shift = False
            shift_rate = 0.
            if (abm.InfectionState(infection_state) == abm.InfectionState.Exposed):
                shift = True
                param_string = get_person_age_string(
                    person.age) + "_IncubationPeriod"
                shift_rate = 1. / \
                    (np.exp(parameters.loc[param_string].value +
                     (parameters.loc[param_string].dev**2)/2.) / 4.)
            elif (abm.InfectionState(infection_state) == abm.InfectionState.InfectedNoSymptoms):
                shift = True
                param_string1 = get_person_age_string(
                    person.age) + "_InfectedNoSymptomsToSymptoms"
                shift_rate1 = 1. / \
                    (np.exp(parameters.loc[param_string1].value +
                     (parameters.loc[param_string1].dev**2)/2.) / 4.)
                param_string2 = get_person_age_string(
                    person.age) + "_InfectedNoSymptomsToRecovered"
                shift_rate2 = 1. / \
                    (np.exp(parameters.loc[param_string2].value +
                     (parameters.loc[param_string2].dev**2)/2.) / 4.)
                shift_rate = np.minimum(shift_rate1, shift_rate2)
            elif (abm.InfectionState(infection_state) == abm.InfectionState.InfectedSymptoms):
                shift = True
                param_string1 = get_person_age_string(
                    person.age) + "_InfectedSymptomsToRecovered"
                shift_rate1 = 1. / \
                    (np.exp(parameters.loc[param_string1].value +
                     (parameters.loc[param_string1].dev**2)/2.) / 4.)
                param_string2 = get_person_age_string(
                    person.age) + "_InfectedSymptomsToSevere"
                shift_rate2 = 1. / \
                    (np.exp(parameters.loc[param_string2].value +
                     (parameters.loc[param_string2].dev**2)/2.) / 4.)
                shift_rate = np.minimum(shift_rate1, shift_rate2)
            # shift = False
            person.add_new_infection(Infection(
                model, person, VirusVariant.Wildtype, t0, abm.InfectionState(infection_state), False, shift, shift_rate), t0)


def save_persons(trip_file):
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
        if (int(home_zone / 10000) != 9162):
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


def map_traffic_cell_to_wastewater_area(mapping_path, wastewater_path, new_file, new_file2):
    rng = np.random.default_rng(30)
    with open(mapping_path) as f:
        d = dict(x.rstrip().split(None, 1) for x in f)
    areas = geopandas.read_file(wastewater_path)
    new_dict = {}
    new_dict2 = {}
    for traffic_cell_id in d.keys():
        if traffic_cell_id[:4] != '9162':
            for loc in d[traffic_cell_id].split(' '):
                new_key = "0"
                if loc in new_dict2:
                    print('Problem: Location is assigned to two traffic cells')
                new_dict2[loc] = new_key
            continue
            # new_key = 'x' + traffic_cell_id
            # new_dict[new_key] = d[traffic_cell_id].split(' ')
        else:
            Id_tan = np.unique(
                areas[areas['id_n'] == int(traffic_cell_id)]['ID_TAN'])
            for loc in d[traffic_cell_id].split(' '):
                new_key = str(random.choice(Id_tan))
                if (new_key in new_dict):
                    new_dict[new_key].append(loc)
                else:
                    new_dict[new_key] = [loc]
                if loc in new_dict2:
                    print('Problem: Location is assigned to two traffic cells')
                new_dict2[loc] = new_key
    # save mapping from locs to TAN IDs
    with open(new_file, 'w') as f:
        for id in new_dict:
            line = id + " "
            for loc in new_dict[id]:
                line += loc + " "
            f.write(line)
            f.write('\n')
        f.close()
    # save mapping from TAN IDs to locs
    with open(new_file2, 'w') as f:
        for id in new_dict2:
            line = id + " "
            line += new_dict2[id] + " "
            f.write(line)
            f.write('\n')
        f.close()
    return new_dict


def num_locations(model):
    num_home = 0
    num_work = 0
    num_school = 0
    num_event = 0
    num_hosp = 0
    num_icu = 0
    num_shop = 0
    for loc in model.locations:
        if loc.type == abm.LocationType.Home:
            num_home += 1
        elif loc.type == abm.LocationType.Work:
            num_work += 1
        elif loc.type == abm.LocationType.School:
            num_school += 1
        elif loc.type == abm.LocationType.SocialEvent:
            num_event += 1
        elif loc.type == abm.LocationType.BasicsShop:
            num_shop += 1
        elif loc.type == abm.LocationType.Hospital:
            num_hosp += 1
        elif loc.type == abm.LocationType.ICU:
            num_icu += 1
        else:
            print('error')
    print('')


def create_home_mapping(map):
    for key, loc_list in map.items():
        map[key] = [loc for loc in loc_list if loc[:2] == "00"]
    return map


def write_time_to_file(sim_nums, init_times, sim_times, output_times, file):
    with open(file, 'w') as f:
        line = "Sim_num,init,simulation,output"
        f.write(line)
        f.write('\n')
        for i in range(len(sim_nums)):
            line = f"{sim_nums[i]},{init_times[i]},{sim_times[i]},{output_times[i]}"
            f.write(line)
            f.write('\n')
    f.close()

def initialize_and_run_sim(sim_num, run_settings, sim_time_frame, person_file, output_path, max_work_size, max_school_size, local_outbreak, location_closure_scheme_random, location_closure_scheme_maximum, parameters):
    print(f'Running {sim_num}')
    start_init = time.time()
    # set seed for initial infection states
    np.random.seed(run_settings["seed"][sim_num])
    # starting time point
    t0 = abm.TimePoint(0)
    # end time point of simulation
    tmax = t0 + abm.days(sim_time_frame)
    # create simulation with starting timepoint and number of age groups
    sim = abm.Simulation(t0, num_age_groups)
    # set seeds for simulation
    abm.set_seeds(sim.model, run_settings["seed"][sim_num])
    # initialize model
    abm.initialize_model(sim.model, person_file, os.path.join(
        input_path, 'hospitals.csv'), os.path.join(
        output_path, str(sim_num) + '_mapping.txt'), max_work_size, max_school_size)
    
    # set infection parameters
    sim.model.parameters = set_infection_parameters(parameters, run_settings["kappa"][sim_num])  
    # set age groups that go to school and work
    abm.set_AgeGroupGoToSchool(sim.model.parameters, age_group_5_to_15)
    abm.set_AgeGroupGoToWork(sim.model.parameters, age_group_16_to_34)
    abm.set_AgeGroupGoToWork(sim.model.parameters, age_group_35_to_59)
    # set age groups that go to shop
    abm.set_AgeGroupGoToShop(sim.model.parameters, age_group_0_to_4)
    abm.set_AgeGroupGoToShop(sim.model.parameters, age_group_5_to_15)
    abm.set_AgeGroupGoToShop(sim.model.parameters, age_group_16_to_34)
    abm.set_AgeGroupGoToShop(sim.model.parameters, age_group_35_to_59)
    abm.set_AgeGroupGoToShop(sim.model.parameters, age_group_60_to_79)
    abm.set_AgeGroupGoToShop(sim.model.parameters, age_group_80_plus)
    # add dampings
    sim.model.add_infection_rate_damping(
        abm.TimePoint(abm.days(int(run_settings["damp_time"][sim_num])).seconds), run_settings["damp_lvl"][sim_num])
    # add closure for work, event, shop and school locations
    sim.model.add_location_closure(abm.TimePoint(
        abm.days(19).seconds), abm.LocationType.Work, 0.37, location_closure_scheme_random)
    sim.model.add_location_closure(abm.TimePoint(
        abm.days(14).seconds), abm.LocationType.School, 1.0, location_closure_scheme_random)
    sim.model.add_location_closure(abm.TimePoint(
        abm.days(19).seconds), abm.LocationType.SocialEvent, 0.51, location_closure_scheme_maximum)
    sim.model.add_location_closure(abm.TimePoint(
        abm.days(19).seconds), abm.LocationType.BasicsShop, 0.13, location_closure_scheme_maximum)
    end_init = time.time()
    print(f'Time for model initialization: {end_init - start_init} seconds')
    # map locations to wastewater areas
    start_map = time.time()
    tan_map = map_traffic_cell_to_wastewater_area(os.path.join(
        output_path, str(sim_num) + '_mapping.txt'), os.path.join(input_path, 'Munich_shape250319/Verschnitt_DLR_TAN_Rep.shp'), os.path.join(
        output_path, str(sim_num) + '_mapping_tan.txt'), os.path.join(
        output_path, str(sim_num) + '_mapping_tan_locs.txt'))
    end_map = time.time()
    print(
        f'Time for mapping locations to TAN areas: {end_map - start_map} seconds')

    start_locs = []
    # specify starting locations if local outbreak should be simulated
    if (local_outbreak):
        start_home_map = time.time()
        home_map = create_home_mapping(tan_map)
        start_areas = ['58']
        start_locs = [loc for area in start_areas for loc in home_map[area]]
        end_home_map = time.time()
        print(
            f'Time for creating outbreak loc list: {end_home_map - start_home_map} seconds')

    # assign initial infection states according to distribution
    start_assign = time.time()
    assign_infection_states(sim.model, t0, run_settings["initE"][sim_num], 0.00005,
                            0.00002, 0.0, 0.0, 0.0, parameters, start_locs)
    end_assign = time.time()
    print(
        f'Time for assigning infection state: {end_assign - start_assign} seconds')

    # initialize tan wastewater ids for all locations
    start_ww_init = time.time()
    abm.set_wastewater_ids(os.path.join(
        output_path, str(sim_num) + '_mapping_tan_locs.txt'), sim.model)
    end_ww_init = time.time()
    print(
        f'Time for initializing TAN areas for all locations: {end_ww_init - start_ww_init} seconds')
    # write size per location
    abm.write_size_per_location(os.path.join(
        output_path, str(sim_num) + '_size_per_loc.txt'), sim.model)
    # output object
    history = History_sensitivity()
    start_advance = time.time()
    # advance simulation until tmax
    sim.advance(tmax, history)
    end_advance = time.time()
    print(
        f'Time for advancing simulation: {end_advance - start_advance} seconds')
    # write compartment size, new infections and summed shedding per time step to file
    start_o2 = time.time()
    abm.save_comp_output_sensitivity(os.path.join(
        output_path, str(sim_num) + '_output.csv'), sim.model, history)
    end_o2 = time.time()
    print(f'Time writing output csv to path {output_path}: {end_o2 - start_o2} seconds')

age_groups = ["Age0to4", "Age5to14", "Age15to34", "Age35to59", "Age60to79", "Age80plus"]
parameter_dict = {"T_E": [f"{i}_IncubationPeriod" for i in age_groups],
                  "T_infected": [f"{i}_InfectedNoSymptomsToSymptoms" for i in age_groups] + [f"{i}_InfectedNoSymptomsToRecovered" for i in age_groups] + [f"{i}_InfectedSymptomsToRecovered" for i in age_groups] + [f"{i}_InfectedSymptomsToSevere" for i in age_groups] + [f"{i}_SevereToRecovered" for i in age_groups] +[f"{i}_SevereToCritical" for i in age_groups] + [f"{i}_CriticalToRecovered" for i in age_groups] + [f"{i}_CriticalToDead" for i in age_groups],
                  "mu_Ins_Isy": [f"{i}_SymptomsPerInfectedNoSymptoms" for i in age_groups],
                  "vl_peak": ["peak_max"]}

"""
This function gets a percentage by which the mean of a lognormally distributed variable should be increased/decreased and calculates the new parameters that match that new mean while the variance stays the same.
"""
def get_new_mu_sigma(mu, sigma, p):
    mean = np.exp(mu + (sigma**2)/2.)
    var = np.exp(2 * mu + (sigma**2)) * (np.exp(sigma**2) - 1)
    new_mean = p * mean
    new_mu = np.log((new_mean * new_mean)/np.sqrt(new_mean * new_mean + var))
    new_sigma = np.sqrt(np.log(1 + var/(new_mean * new_mean)))
    return new_mu, new_sigma

"""
Calculates mean of lognorm distribution.
"""
def get_mean(mu, sigma):
    return np.exp(mu + (sigma**2)/2.)

def run_sim_in_parallel(sim):
    return initialize_and_run_sim(sim, global_run_settings, global_sim_time_frame, global_person_file, global_var_path, global_max_work_size, global_max_school_size, global_local_outbreak, global_location_closure_scheme_random, global_location_closure_scheme_maximum, copy.deepcopy(global_parameters))

"""
infection_parameter_to_vary: should be chosen from the parameter_dict.keys()
variations: are the factors the parameters are multiplied with for the sensitivity analysis. For +/-20%, +/-10% these would be [0.8, 0.9, 1, 1.1, 1.2]
"""
def run_sensitivity(num_sims, run_settings, input_path, output_path, infection_parameter_to_vary, p):
    global global_run_settings, global_sim_time_frame, global_person_file
    global global_var_path, global_max_work_size, global_max_school_size
    global global_local_outbreak, global_parameters
    global global_location_closure_scheme_random, global_location_closure_scheme_maximum
    
    person_file = input_path + 'persons_scaled.csv'
    local_outbreak = False
    max_work_size = 40
    max_school_size = 45
    sim_time_frame = 92
    # Possible schemes are 'random' and 'maximum'
    location_closure_scheme_random = 'random'
    location_closure_scheme_maximum = 'maximum'

    # read infection parameters
    parameters = pd.read_csv(os.path.join(
        input_path, 'parameter_table.csv'), index_col=0)
    
    global_run_settings = run_settings
    global_sim_time_frame = sim_time_frame
    global_person_file = person_file
    global_max_work_size = max_work_size
    global_max_school_size = max_school_size
    global_local_outbreak = local_outbreak
    global_location_closure_scheme_random = location_closure_scheme_random
    global_location_closure_scheme_maximum = location_closure_scheme_maximum
    
    print(f"p={p}", file=sys.stderr)
    var_path = output_path + f"{p}"
    print("Var path: ", var_path, file=sys.stderr)
    if not os.path.exists(var_path):
        os.makedirs(var_path)
        
    global_var_path = var_path
    print("Gobal var path: ", var_path, file=sys.stderr)
    # Iterate over all ABM parameters targeted by the sensitivity variation
    for ABM_param in parameter_dict[infection_parameter_to_vary]:
        if(infection_parameter_to_vary == "T_E"):
            mu_old = parameters.loc[ABM_param].value.copy()
            sigma_old = parameters.loc[ABM_param].dev.copy()
            mean_old = get_mean(mu_old, sigma_old)
            mu_new, sigma_new = get_new_mu_sigma(mu_old, sigma_old, p)
            mean_new = get_mean(mu_new, sigma_new)
            parameters.loc[ABM_param].value = mu_new
            parameters.loc[ABM_param].dev = sigma_new
            # Adapt non-symptomatic stay times such that the whole infected period stays the same
            ag = ABM_param.split('_')[0]
            # NoSymptoms to Symptoms
            mu_Ins_Sy_old = parameters.loc[f"{ag}_InfectedNoSymptomsToSymptoms"].value.copy()
            sigma_Ins_Sy_old = parameters.loc[f"{ag}_InfectedNoSymptomsToSymptoms"].dev.copy()
            p_Ins_Sy = 1 - (mean_new - mean_old)/get_mean(mu_Ins_Sy_old, sigma_Ins_Sy_old)
            mu_Ins_Sy_new, sigma_Ins_Sy_new = get_new_mu_sigma(mu_Ins_Sy_old, sigma_Ins_Sy_old, p_Ins_Sy)
            parameters.loc[f"{ag}_InfectedNoSymptomsToSymptoms"].value = mu_Ins_Sy_new
            parameters.loc[f"{ag}_InfectedNoSymptomsToSymptoms"].dev = sigma_Ins_Sy_new
            # Nosymptoms to Recovered
            mu_Ins_R_old = parameters.loc[f"{ag}_InfectedNoSymptomsToRecovered"].value.copy()
            sigma_Ins_R_old = parameters.loc[f"{ag}_InfectedNoSymptomsToRecovered"].dev.copy()
            p_Ins_R = 1 - (mean_new - mean_old)/get_mean(mu_Ins_R_old, sigma_Ins_R_old)
            mu_Ins_R_new, sigma_Ins_R_new = get_new_mu_sigma(mu_Ins_R_old, sigma_Ins_R_old, p_Ins_R)
            parameters.loc[f"{ag}_InfectedNoSymptomsToRecovered"].value = mu_Ins_R_new
            parameters.loc[f"{ag}_InfectedNoSymptomsToRecovered"].dev = sigma_Ins_R_new
        elif(infection_parameter_to_vary == "T_infected"):
            mu_old = parameters.loc[ABM_param].value.copy()
            sigma_old = parameters.loc[ABM_param].dev.copy()
            mu_new, sigma_new = get_new_mu_sigma(mu_old, sigma_old, p)
            parameters.loc[ABM_param].value = mu_new
            parameters.loc[ABM_param].dev = sigma_new
        else:
            old_value = parameters.loc[ABM_param].value.copy()
            parameters.loc[ABM_param].value = p * old_value
        
        global_parameters = parameters
    
    print(f"Starting advance for p={p}", file=sys.stderr)
    with Pool(processes=os.cpu_count()) as pool:  # Use os.cpu_count() for dynamic number of processes
        pool.map(run_sim_in_parallel, range(0, num_sims))
    print(f"End advance for p={p}", file=sys.stderr)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'abm demonstrator',
        description='Example demonstrating the agent-based model for a synthetic population of Munich.')
    args = arg_parser.parse_args()
    # set LogLevel
    mio.abm.set_log_level_warn()
    num_sims = 100
    parameter_to_vary = "vl_peak"
    input_path = '/hpc_data/bick_ju/INSIDeMunichPaper/IScienceRebuttal/Input/'
    output_path = '/hpc_data/bick_ju/INSIDeMunichPaper/IScienceRebuttal/sensitivity_analysis/' + parameter_to_vary + '/'
    run_settings = pd.read_csv(input_path + 'pop8_sel_particles.csv', sep=',')
    for p in [0.8, 0.9, 1.0, 1.1, 1.2]:
        run_sensitivity(num_sims, run_settings, input_path, output_path, parameter_to_vary, p, **args.__dict__)
