#############################################################################
# Copyright (C) 2020-2025 MEmilio
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
import numpy as np
from enum import Enum
import os

import memilio.simulation as mio
import memilio.simulation.osecirvvs as osecirvvs
import memilio.plot.createGIF as mp


class Location(Enum):
    """ """
    Home = 0
    School = 1
    Work = 2
    Other = 3


class Intervention(Enum):
    """ """
    Home = 0
    SchoolClosure = 1
    HomeOffice = 2
    GatheringBanFacilitiesClosure = 3
    PhysicalDistanceAndMasks = 4
    SeniorAwareness = 5


class InterventionLevel(Enum):
    """ """
    Main = 0
    PhysicalDistanceAndMasks = 1
    SeniorAwareness = 2
    Holidays = 3


class Simulation:
    """ """

    def __init__(self, data_dir, results_dir):
        self.num_groups = 6
        self.data_dir = data_dir
        self.results_dir = results_dir
        if not os.path.exists(self.results_dir):
            os.makedirs(self.results_dir)

        self.late = False
        self.masks = True
        self.test = True
        self.high = False
        self.long_time = False
        self.future = False

        if self.future:
            self.start_date = mio.Date(2021, 10, 15)
        else:
            self.start_date = mio.Date(2021, 6, 6)

    def set_covid_parameters(self, model):
        """

        :param model: 

        """
        def array_assign_uniform_distribution(param, min, max, num_groups=6):
            """

            :param param: 
            :param min: 
            :param max: 
            :param num_groups:  (Default value = 6)

            """
            if isinstance(
                    min, (int, float)) and isinstance(
                        max, (int, float)):
                min = min * np.ones(num_groups)
                max = max * np.ones(num_groups)
            elif not (isinstance(min, (list)) and isinstance(
                    max, (list))):
                raise TypeError(
                    "Invalid type for parameter 'min' or 'max. \
                            Expected a scalar or a list. Must be the same for both.")
            for i in range(num_groups):
                param[mio.AgeGroup(i)] = mio.UncertainValue(
                    0.5 * (max[i] + min[i]))
                param[mio.AgeGroup(i)].set_distribution(
                    mio.ParameterDistributionUniform(min[i], max[i]))

        # times
        timeExposedMin = 2.67
        timeExposedMax = 4.
        timeInfectedNoSymptomsMin = 1.2
        timeInfectedNoSymptomsMax = 2.53
        timeInfectedSymptomsMin = [
            5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465]
        timeInfectedSymptomsMax = [8.427, 8.427, 8.4684, 8.3139, 8.169, 8.085]
        timeInfectedSevereMin = [3.925, 3.925, 4.85, 6.4, 7.2, 9.]
        timeInfectedSevereMax = [6.075, 6.075, 7., 8.7, 9.8, 13.]
        timeInfectedCriticalMin = [4.95, 4.95, 4.86, 14.14, 14.4, 10.]
        timeInfectedCriticalMax = [8.95, 8.95, 8.86, 20.58, 19.8, 13.2]

        array_assign_uniform_distribution(
            model.parameters.TimeExposed, timeExposedMin, timeExposedMax)

        array_assign_uniform_distribution(
            model.parameters.TimeInfectedNoSymptoms, timeInfectedNoSymptomsMin,
            timeInfectedNoSymptomsMax)

        array_assign_uniform_distribution(
            model.parameters.TimeInfectedSymptoms, timeInfectedSymptomsMin,
            timeInfectedSymptomsMax)

        array_assign_uniform_distribution(
            model.parameters.TimeInfectedSevere, timeInfectedSevereMin,
            timeInfectedSevereMax)

        array_assign_uniform_distribution(
            model.parameters.TimeInfectedCritical, timeInfectedCriticalMin,
            timeInfectedCriticalMax)

        # probabilities
        variantFactor = 1.4
        transmissionProbabilityOnContactMin = [
            0.02 * variantFactor, 0.05 * variantFactor, 0.05 * variantFactor, 0.05 * variantFactor, 0.08 * variantFactor, 0.1 * variantFactor]
        transmissionProbabilityOnContactMax = [
            0.04 * variantFactor, 0.07 * variantFactor, 0.07 * variantFactor, 0.07 * variantFactor, 0.10 * variantFactor, 0.15 * variantFactor]
        relativeTransmissionNoSymptomsMin = 0.5
        relativeTransmissionNoSymptomsMax = 0.5

        # The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
        # depends on incidence and test and trace capacity
        riskOfInfectionFromSymptomaticMin = 0.0
        riskOfInfectionFromSymptomaticMax = 0.2
        maxRiskOfInfectionFromSymptomaticMin = 0.4
        maxRiskOfInfectionFromSymptomaticMax = 0.5
        recoveredPerInfectedNoSymptomsMin = [0.2, 0.2, 0.15, 0.15, 0.15, 0.15]
        recoveredPerInfectedNoSymptomsMax = [0.3, 0.3, 0.25, 0.25, 0.25, 0.25]
        severePerInfectedSymptomsMin = [0.006, 0.006, 0.015, 0.049, 0.15, 0.20]
        severePerInfectedSymptomsMax = [0.009, 0.009, 0.023, 0.074, 0.18, 0.25]
        criticalPerSevereMin = [0.05, 0.05, 0.05, 0.10, 0.25, 0.35]
        criticalPerSevereMax = [0.10, 0.10, 0.10, 0.20, 0.35, 0.45]
        deathsPerCriticalMin = [0.00, 0.00, 0.10, 0.10, 0.30, 0.5]
        deathsPerCriticalMax = [0.10, 0.10, 0.18, 0.18, 0.50, 0.7]

        reducExposedPartialImmunityMin = 0.75
        reducExposedPartialImmunityMax = 0.85
        reducExposedImprovedImmunityMin = 0.281
        reducExposedImprovedImmunityMax = 0.381
        reducInfectedSymptomsPartialImmunityMin = 0.6
        reducInfectedSymptomsPartialImmunityMax = 0.7
        reducInfectedSymptomsImprovedImmunityMin = 0.193
        reducInfectedSymptomsImprovedImmunityMax = 0.293
        reducInfectedSevereCriticalDeadPartialImmunityMin = 0.05
        reducInfectedSevereCriticalDeadPartialImmunityMax = 0.15
        reducInfectedSevereCriticalDeadImprovedImmunityMin = 0.041
        reducInfectedSevereCriticalDeadImprovedImmunityMax = 0.141

        if self.long_time:
            temp_reducTimeInfectedMild = 1.0
        else:
            temp_reducTimeInfectedMild = 0.5
        reducTimeInfectedMild = temp_reducTimeInfectedMild

        array_assign_uniform_distribution(
            model.parameters.TransmissionProbabilityOnContact,
            transmissionProbabilityOnContactMin,
            transmissionProbabilityOnContactMax)

        array_assign_uniform_distribution(
            model.parameters.RelativeTransmissionNoSymptoms,
            relativeTransmissionNoSymptomsMin,
            relativeTransmissionNoSymptomsMax)

        array_assign_uniform_distribution(
            model.parameters.RiskOfInfectionFromSymptomatic,
            riskOfInfectionFromSymptomaticMin,
            riskOfInfectionFromSymptomaticMax)

        array_assign_uniform_distribution(
            model.parameters.MaxRiskOfInfectionFromSymptomatic,
            maxRiskOfInfectionFromSymptomaticMin,
            maxRiskOfInfectionFromSymptomaticMax)

        array_assign_uniform_distribution(
            model.parameters.RecoveredPerInfectedNoSymptoms,
            recoveredPerInfectedNoSymptomsMin,
            recoveredPerInfectedNoSymptomsMax)

        array_assign_uniform_distribution(
            model.parameters.SeverePerInfectedSymptoms,
            severePerInfectedSymptomsMin,
            severePerInfectedSymptomsMax)

        array_assign_uniform_distribution(
            model.parameters.CriticalPerSevere,
            criticalPerSevereMin,
            criticalPerSevereMax)

        array_assign_uniform_distribution(
            model.parameters.DeathsPerCritical,
            deathsPerCriticalMin,
            deathsPerCriticalMax)

        array_assign_uniform_distribution(
            model.parameters.ReducExposedPartialImmunity,
            reducExposedPartialImmunityMin,
            reducExposedPartialImmunityMax)

        array_assign_uniform_distribution(
            model.parameters.ReducExposedImprovedImmunity,
            reducExposedImprovedImmunityMin,
            reducExposedImprovedImmunityMax)

        array_assign_uniform_distribution(
            model.parameters.ReducInfectedSymptomsPartialImmunity,
            reducInfectedSymptomsPartialImmunityMin,
            reducInfectedSymptomsPartialImmunityMax)

        array_assign_uniform_distribution(
            model.parameters.ReducInfectedSymptomsImprovedImmunity,
            reducInfectedSymptomsImprovedImmunityMin,
            reducInfectedSymptomsImprovedImmunityMax)

        array_assign_uniform_distribution(
            model.parameters.ReducInfectedSevereCriticalDeadPartialImmunity,
            reducInfectedSevereCriticalDeadPartialImmunityMin,
            reducInfectedSevereCriticalDeadPartialImmunityMax)

        array_assign_uniform_distribution(
            model.parameters.ReducInfectedSevereCriticalDeadImprovedImmunity,
            reducInfectedSevereCriticalDeadImprovedImmunityMin,
            reducInfectedSevereCriticalDeadImprovedImmunityMax)

        array_assign_uniform_distribution(
            model.parameters.ReducTimeInfectedMild,
            reducTimeInfectedMild,
            reducTimeInfectedMild)

        # start day is set to the n-th day of the year
        model.parameters.StartDay = self.start_date.day_in_year

        model.parameters.Seasonality.value = 0.2

        seasonality_min = 0.1
        seasonality_max = 0.3
        model.parameters.Seasonality = mio.UncertainValue(
            0.5 * (seasonality_max + seasonality_min))
        model.parameters.Seasonality.set_distribution(
            mio.ParameterDistributionUniform(seasonality_min, seasonality_max))

        model.parameters.StartDayNewVariant = mio.Date(2021, 6, 6).day_in_year

    def set_contact_matrices(self, model):
        """

        :param model: 

        """
        contact_matrices = mio.ContactMatrixGroup(
            len(list(Location)), self.num_groups)
        locations = ["home", "school_pf_eig", "work", "other"]

        for i, location in enumerate(locations):
            baseline_file = os.path.join(
                self.data_dir, "Germany", "contacts", "baseline_" + location + ".txt")
            contact_matrices[i] = mio.ContactMatrix(
                mio.read_mobility_plain(baseline_file),
                np.zeros((self.num_groups, self.num_groups))
            )
        model.parameters.ContactPatterns.cont_freq_mat = contact_matrices

    def set_npis(self, params, end_date):
        """

        :param params: 
        :param end_date: 

        """
        contacts = params.ContactPatterns
        dampings = contacts.dampings

        if self.test:
            params.commuter_nondetection = 0.85
        else:
            params.commuter_nondetection = 1.0

        group_weights_all = np.ones(self.num_groups)
        group_weights_seniors = np.vectorize(
            lambda i: 1.0 if i == 5 else (0.5 if i == 4 else 0.0))(
            np.arange(self.num_groups))

        loc_home = Location.Home.value
        loc_school = Location.School.value
        loc_work = Location.Work.value
        loc_other = Location.Other.value

        lvl_main = InterventionLevel.Main.value
        lvl_pd_and_masks = InterventionLevel.PhysicalDistanceAndMasks.value
        lvl_seniors = InterventionLevel.SeniorAwareness.value
        lvl_holidays = InterventionLevel.Holidays.value

        typ_home = Intervention.Home.value
        typ_school = Intervention.SchoolClosure.value
        typ_homeoffice = Intervention.HomeOffice.value
        typ_gathering = Intervention.GatheringBanFacilitiesClosure.value
        typ_distance = Intervention.PhysicalDistanceAndMasks.value
        typ_senior = Intervention.SeniorAwareness.value

        def damping_helper(
                t, min, max, damping_level, type, location,
                group_weights=group_weights_all):
            """

            :param t: 
            :param min: 
            :param max: 
            :param damping_level: 
            :param type: 
            :param location: 
            :param group_weights:  (Default value = group_weights_all)

            """
            v = mio.UncertainValue(0.5 * (max + min))
            v.set_distribution(mio.ParameterDistributionUniform(min, max))
            return mio.DampingSampling(
                value=v,
                level=damping_level,
                type=type,
                time=t,
                matrix_indices=location,
                group_weights=group_weights)

        def contacts_at_home(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_main, typ_home, [loc_home])

        def school_closure(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_main, typ_school, [loc_school])

        def home_office(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_main, typ_homeoffice, [loc_work])

        def social_events(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_main, typ_gathering, [loc_other])

        def social_events_work(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_main, typ_gathering, [loc_work])

        def physical_distancing_home(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_pd_and_masks, typ_distance, [loc_home])

        def physical_distancing_school(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_pd_and_masks, typ_distance, [loc_school])

        def physical_distancing_work(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_pd_and_masks, typ_distance, [loc_work])

        def physical_distancing_other(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_pd_and_masks, typ_distance, [loc_other])

        def senior_awareness(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_seniors, typ_senior, [loc_home, loc_other],
                group_weights_seniors)

        # OPEN SCENARIO SPRING
        start_year = mio.Date(
            2021, 1, 1)
        narrow = 0.05

        if start_year < end_date:
            static_open_scenario_spring = start_year - self.start_date
            dampings.append(contacts_at_home(
                static_open_scenario_spring, 0.0, 0.0))
            dampings.append(school_closure(
                static_open_scenario_spring, 0.0, 0.0))
            dampings.append(home_office(static_open_scenario_spring, 0.0, 0.0))
            dampings.append(social_events(
                static_open_scenario_spring, 0.0, 0.0))
            dampings.append(social_events_work(
                static_open_scenario_spring, 0.0, 0.0))
            dampings.append(physical_distancing_home(
                static_open_scenario_spring, 0.0, 0.0))
            dampings.append(physical_distancing_school(
                static_open_scenario_spring, 0.2 + narrow, 0.4 - narrow))
            dampings.append(physical_distancing_work(
                static_open_scenario_spring, 0.2 + narrow, 0.4 - narrow))
            dampings.append(physical_distancing_other(
                static_open_scenario_spring, 0.2 + narrow, 0.4 - narrow))
            dampings.append(senior_awareness(
                static_open_scenario_spring, 0.0, 0.0))

        # OPEN SCENARIO
        if self.late:
            month_open = 8
        else:
            month_open = 7

        if self.masks:
            masks_low_school = 0.2
            masks_high_school = 0.4
            masks_low = 0.2
            masks_high = 0.4
            masks_narrow = narrow
        else:
            masks_low_school = 0.0
            masks_high_school = 0.0
            masks_low = 0.0
            masks_high = 0.0
            masks_narrow = 0.0

        start_open = mio.Date(
            2021, month_open, 1)
        start_summer = start_open - self.start_date
        params.end_dynamic_npis = start_summer

        if start_open < end_date:
            dampings.append(contacts_at_home(start_summer, 0.0, 0.0))
            dampings.append(school_closure(start_summer, 0.0, 0.0))
            dampings.append(home_office(start_summer, 0.0, 0.0))
            dampings.append(social_events(start_summer, 0.0, 0.0))
            dampings.append(social_events_work(start_summer, 0.0, 0.0))
            dampings.append(physical_distancing_home(start_summer, 0.0, 0.0))
            dampings.append(physical_distancing_school(
                start_summer, masks_low_school + masks_narrow, masks_high_school - masks_narrow))
            dampings.append(physical_distancing_work(
                start_summer, masks_low + masks_narrow, masks_high - masks_narrow))
            dampings.append(physical_distancing_other(
                start_summer, masks_low + masks_narrow, masks_high - masks_narrow))
            dampings.append(senior_awareness(start_summer, 0.0, 0.0))

        start_autumn = mio.Date(2021, 10, 1) - self.start_date
        dampings.append(contacts_at_home(start_autumn, 0.0, 0.0))
        dampings.append(school_closure(
            start_autumn, 0.3 + narrow, 0.5 - narrow))
        # dampings.append(home_office(start_autumn, 0.3 + narrow, 0.5 - narrow)) #S3F only
        dampings.append(social_events(
            start_autumn, 0.3 + narrow, 0.5 - narrow))
        dampings.append(social_events_work(start_autumn, 0.0, 0.0))

        dampings.append(home_office(
            start_autumn, 0.0 + narrow, 0.2 - narrow))  # S2F

        # dampings.append(school_closure(start_autumn, 0.0 + narrow, 0.2 - narrow)) #S1F
        # dampings.append(home_office(start_autumn, 0.0 + narrow, 0.2 - narrow)) #S1F
        # dampings.append(social_events(start_autumn, 0.0 + narrow, 0.2 - narrow)) #S1F

        params.ContactPatterns.dampings = dampings

        narrow = 0.0
        # local dynamic NPIs
        dynamic_npis = params.DynamicNPIsInfectedSymptoms

        dynamic_npi_dampings = []
        dynamic_npi_dampings.append(
            contacts_at_home(0, 0.1 + narrow, 0.3 - narrow))
        dynamic_npi_dampings.append(school_closure(
            0, 0.2 + narrow, 0.4 - narrow))  # 0.25 - 0.25 in autumn
        dynamic_npi_dampings.append(home_office(0, 0.1 + narrow, 0.3 - narrow))
        dynamic_npi_dampings.append(
            social_events(0, 0.2 + narrow, 0.4 - narrow))
        dynamic_npi_dampings.append(social_events_work(0, 0.0, 0.0))
        dynamic_npi_dampings.append(physical_distancing_home(0, 0.0, 0.0))
        dynamic_npi_dampings.append(
            physical_distancing_school(0, 0.2 + narrow, 0.4 - narrow))
        dynamic_npi_dampings.append(
            physical_distancing_work(0, 0.2 + narrow, 0.4 - narrow))
        dynamic_npi_dampings.append(
            physical_distancing_other(0, 0.2 + narrow, 0.4 - narrow))
        dynamic_npi_dampings.append(senior_awareness(0, 0.0, 0.0))

        dynamic_npi_dampings2 = []
        dynamic_npi_dampings2.append(
            contacts_at_home(0, 0.5 + narrow, 0.7 - narrow))
        dynamic_npi_dampings2.append(school_closure(
            0, 0.4 + narrow, 0.6 - narrow))  # 0.25 - 0.25 in autumn
        dynamic_npi_dampings2.append(
            home_office(0, 0.2 + narrow, 0.4 - narrow))
        dynamic_npi_dampings2.append(
            social_events(0, 0.7 + narrow, 0.9 - narrow))
        dynamic_npi_dampings2.append(social_events_work(0, 0.0, 0.0))
        dynamic_npi_dampings2.append(
            physical_distancing_home(0, 0.0 + narrow, 0.2 - narrow))
        dynamic_npi_dampings2.append(
            physical_distancing_school(0, 0.2 + narrow, 0.4 - narrow))
        dynamic_npi_dampings2.append(
            physical_distancing_work(0, 0.2 + narrow, 0.4 - narrow))
        dynamic_npi_dampings2.append(
            physical_distancing_other(0, 0.2 + narrow, 0.4 - narrow))
        dynamic_npi_dampings2.append(senior_awareness(0, 0.0, 0.0))

        dynamic_npis.interval = 1.0
        dynamic_npis.duration = 14.0
        dynamic_npis.base_value = 100000
        dynamic_npis.set_threshold(35.0, dynamic_npi_dampings)
        dynamic_npis.set_threshold(100.0, dynamic_npi_dampings2)

        params.DynamicNPIsInfectedSymptoms = dynamic_npis

        # school holidays (holiday periods are set per node)
        school_holiday_value = mio.UncertainValue(0.5 * (1.0 + 1.0))
        school_holiday_value.set_distribution(
            mio.ParameterDistributionUniform(1.0, 1.0))
        params.ContactPatterns.school_holiday_damping = mio.DampingSampling(school_holiday_value, InterventionLevel.Holidays.value,
                                                                            Intervention.SchoolClosure.value, 0.0,
                                                                            [Location.School.value], group_weights_all)

    def get_graph(self, end_date):
        """

        :param end_date: 

        """
        model = osecirvvs.Model(self.num_groups)
        self.set_covid_parameters(model)
        self.set_contact_matrices(model)
        self.set_npis(model.parameters, end_date)

        graph = osecirvvs.ModelGraph()

        scaling_factor_infected = np.ones(self.num_groups)
        scaling_factor_icu = 1.0
        tnt_capacity_factor = 1.43 / 100000.

        data_dir_Germany = os.path.join(self.data_dir, "Germany")
        mobility_data_file = os.path.join(
            data_dir_Germany, "mobility", "commuter_mobility_2022.txt")
        pydata_dir = os.path.join(data_dir_Germany, "pydata")

        graph = mio.osecir.create_graph_german_county(
            model.parameters, self.start_date, end_date,
            scaling_factor_infected, scaling_factor_icu, pydata_dir, tnt_capacity_factor)

        osecirvvs.set_edges(
            mobility_data_file, graph, len(Location))

        return graph

    def run(self, num_days_sim, num_runs=10, save_graph=True, create_gif=True):
        """

        :param num_days_sim: 
        :param num_runs:  (Default value = 10)
        :param save_graph:  (Default value = True)
        :param create_gif:  (Default value = True)

        """
        mio.set_log_level(mio.LogLevel.Warning)
        end_date = self.start_date + num_days_sim

        graph = self.get_graph(end_date)

        if save_graph:
            path_graph = os.path.join(self.results_dir, "graph")
            if not os.path.exists(path_graph):
                os.makedirs(path_graph)
            osecirvvs.write_graph(graph, path_graph)

        study = osecirvvs.GraphParameterStudy(
            graph, 0., num_days_sim, 0.5, num_runs)
        ensemble = study.run(self.high)

        ensemble_results = []
        ensemble_params = []
        for run in range(num_runs):
            graph_run = ensemble[run]
            ensemble_results.append(
                osecirvvs.interpolate_simulation_result(graph_run))
            ensemble_params.append(
                [graph_run.get_node(node_indx).property.model
                 for node_indx in range(graph.num_nodes)])

        node_ids = [graph.get_node(i).id for i in range(graph.num_nodes)]

        save_percentiles = True
        save_single_runs = False

        osecirvvs.save_results(
            ensemble_results, ensemble_params, node_ids, self.results_dir,
            save_single_runs, save_percentiles)
        if create_gif:
            # any compartments in the model (see InfectionStates)
            compartments = [c for c in range(2, 23)]
            mp.create_gif_map_plot(
                self.results_dir + "/p75", self.results_dir, compartments)


if __name__ == "__main__":
    file_path = os.path.dirname(os.path.abspath(__file__))
    sim = Simulation(
        data_dir=os.path.join(file_path, "../../../data"),
        results_dir=os.path.join(file_path, "../../../results_osecirvvs"))
    sim.run(num_days_sim=30, num_runs=5)
